#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 09:21:15 2022

@author: arest, mcorrenti

This is class wrapper around doing simple photometry on a single JWST image
"""

import argparse,glob,re,sys,os,time,copy,math

import glob as glob

import numpy as np

import urllib
import scipy

from astropy.io import fits
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.table import Table,vstack
from astropy import wcs


from photutils.detection import DAOStarFinder
from photutils.background import MMMBackground, MADStdBackgroundRMS, Background2D
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from astropy.visualization import simple_norm

import jwst
from jwst.datamodels import ImageModel
from jwst import datamodels
import astropy.units as u
import pysiaf
from astroquery.gaia import Gaia
from astropy.time import Time
import pandas as pd

from jwcf import hawki, hst
from astropy.coordinates import SkyCoord


from pdastro import pdastroclass,pdastrostatsclass,makepath4file,unique,AnotB,AorB,AandB,rmfile
from astropy.coordinates import SkyCoord, match_coordinates_sky

def get_jwst_ap_corr_table():
    
    cnames = ['filter', 'pupil', 'wave', 'r10', 'r20', 'r30', 'r40', 'r50', 'r60', 'r70', 'r80',
                                        'r85', 'r90', 'sky_flux_px', 'apcorr10', 'apcorr20', 'apcorr30', 'apcorr40',
                                        'apcorr50', 'apcorr60', 'apcorr70', 'apcorr80', 'apcorr85', 'apcorr90', 'sky_in',
                                        'sky_out']
    
    if not os.path.isfile('./jwst_aperture_corrections.txt'):
        if os.path.isfile('./aperture_correction_table.txt'):
            ap_tab = './aperture_correction_table.txt'
        else:
            print("Downloading the nircam aperture correction table")

            boxlink_apcorr_table = 'https://data.science.stsci.edu/redirect/JWST/jwst-data_analysis_tools/stellar_photometry/aperture_correction_table.txt'
            boxfile_apcorr_table = './aperture_correction_table.txt'
            urllib.request.urlretrieve(boxlink_apcorr_table, boxfile_apcorr_table)
            ap_tab = './aperture_correction_table.txt'

        aper_table1 = Table.read(ap_tab,format='ascii',names=cnames)
        if os.path.isfile('./jwst_miri_apcorr_0008.fits'):
            ap_tab = Table.read('jwst_miri_apcorr_0008.fits',format='fits')
        else:
            print("Downloading the miri aperture correction table")
            urllib.request.urlretrieve('https://jwst-crds.stsci.edu/unchecked_get/references/jwst/jwst_miri_apcorr_0008.fits',
                                       'jwst_miri_apcorr_0008.fits')
        tab = Table.read('jwst_miri_apcorr_0008.fits',format='fits')
        rows = {k:[] for k in cnames}
        tab=tab[tab['subarray']=='FULL']
        tab.rename_column('subarray','pupil')
        tab.rename_column('skyin','sky_in')
        tab.rename_column('skyout','sky_out')
        for i in range(len(np.unique(tab['filter']))):
            temp_tab = tab[tab['filter']==np.unique(tab['filter'])[i]]
            for c in cnames:
                if c in tab.colnames:
                    rows[c].append(temp_tab[0][c])
                elif c=='wave':
                    rows[c].append(int(temp_tab[0]['filter'][1:-1]))
                elif c=='sky_flux_px':
                    rows[c].append(np.nan)
                elif c.startswith('r'):
                    ind = np.where(temp_tab['eefraction']==float(c[1:])/100)[0]
                    if len(ind)>0:
                        rows[c].append(temp_tab[ind]['radius'])
                        rows['apcorr'+c[1:]].append(temp_tab[ind]['apcorr'])
                    else:
                        rows[c].append(np.nan)
                        rows['apcorr'+c[1:]].append(np.nan)
        
        aper_table2 = Table(rows)
        for k in aper_table1.colnames:
            aper_table1[k] = np.array(aper_table1[k]).flatten()
            aper_table2[k] = np.array(aper_table2[k]).flatten()
        aper_table = vstack([aper_table1,aper_table2])
        
        aper_table.write('./jwst_aperture_corrections.txt',format='ascii')
    else:
        aper_table = Table.read('./jwst_aperture_corrections.txt',format='ascii')
    return(aper_table)

def get_GAIA_sources_NP(f,pm=True,radius_factor=2):
    """"Nor Pirzkal: Return a GAIA catalog/table that is location and time matched to an observation.
    It applies the GAIA proper motion by default (pm=True).
    The (x,y) pixel coordinates of the GAIA source are also returned in  'x' and 'y' column.
    radius_factor is 2 by default and hence about twice that of the FOV but can be set."""

    h = fits.open(f)[1].header
    d = fits.open(f)[1].data

    direct_with_wcs = datamodels.open(f)
    world_to_pix = direct_with_wcs.meta.wcs.get_transform('world','detector')
    pix_to_world = direct_with_wcs.meta.wcs.get_transform('detector','world')

    sy,sx = np.shape(d)
    ra0,dec0 = pix_to_world(sx/2,sy/2)
    radius = np.abs(pix_to_world(sx,sy)[1] - dec0)*radius_factor # Approx radius needed in deg
        
    job5 = Gaia.launch_job_async("SELECT * \
                    FROM gaiadr2.gaia_source WHERE CONTAINS(POINT('ICRS',\
                    gaiadr2.gaia_source.ra,gaiadr2.gaia_source.dec),\
                    CIRCLE('ICRS',{},{} ,{}))=1;".format(ra0,dec0,radius))
    tb_gaia = job5.get_results() 
    print("Number of stars:",len(tb_gaia))
    
    if pm is True:
        time_gaia = Time(tb_gaia['ref_epoch'], format = 'jyear')[0]
        time_obs = Time(h['MJD-AVG'], format ='mjd')

        dRA = ((time_obs - time_gaia).to(u.yr).value * tb_gaia['pmra'].data * u.mas / np.cos(np.deg2rad(tb_gaia['dec']))).to(u.deg).value
        dDec = ((time_obs - time_gaia).to(u.yr).value * tb_gaia['pmdec'].data * u.mas).to(u.deg).value

        tb_gaia['ra'] = tb_gaia['ra'] + dRA
        tb_gaia['dec'] = tb_gaia['dec'] + dDec
    
    tb_gaia['x'],tb_gaia['y'] = world_to_pix(tb_gaia['ra'],tb_gaia['dec'])
    
    ok = (np.isfinite(tb_gaia['x'])) & (np.isfinite(tb_gaia['y'])) 
    return tb_gaia[ok]

def get_refcat_file(refcatfilename,racol=None,deccol=None):
    cat = pdastroclass()
    cat.load(refcatfilename)
    if racol is not None:
        if not (racol in cat.t.columns):
            raise RuntimeError(f'Cannot find racol {racol} in columns {cat.t.columns}')
    if deccol is not None:
        if not (deccol in cat.t.columns):
            raise RuntimeError(f'Cannot find deccol {deccol} in columns {cat.t.columns}')
    return(cat.t)

def get_hawki_cat(ra0,dec0,radius_deg,radius_factor=1.1,
                  columns=['ID','ra','dec','ra_error_mas','dec_error_mas','J2mag','K2mag','J2_K2','q_Jmag','gdr2_source_id','sep_arcmin']):

    cat = hawki.hawki_catalog()
    cat.rename_column('ra_deg', 'ra')
    cat.rename_column('dec_deg', 'dec')
    cat.rename_column('j_2mass_extrapolated', 'j_magnitude')
    #print(cat.columns)

    pos0 = SkyCoord(ra=ra0, dec=dec0, unit=(u.deg,u.deg), frame='icrs')
    pos1 = SkyCoord(ra=cat['ra'], dec=cat['dec'], unit=(u.deg,u.deg), frame='icrs')
    sep = pos0.separation(pos1)
    print(sep.degree)
    cat['sep_deg'] = sep
    cat['sep_arcmin'] = sep/60.0
    cat['J2_K2'] = cat['J2mag'] - cat['K2mag']

    df = cat.to_pandas()
    
    # only keep the values within the radius
    ixs = df.index.values
    Ntot=len(ixs)
    (keep,) = np.where(df.loc[ixs,'sep_deg'].le(radius_deg*radius_factor))
    ixs = ixs[keep]
    print(f'Keeping {len(ixs)} out of {Ntot} of hawkI sources: all positions within {radius_deg*radius_factor:.4f} deg of RA/Dec=({ra0},{dec0})')
    df = df.loc[ixs]
     
    #print('bbb0',columns)
    if columns is not None:
        df = df[columns]

        
    return(df,'ra','dec')


def get_hst_cat(ra0,dec0,radius_deg,radius_factor=1.1,
                mjd=None,
                columns=['ID','ra','dec','ra_error_mas','dec_error_mas','pmra','epmra','pmdec','epmdec',
                        'h_mag_vega','j_mag_vega','k_mag_vega','j_k',
                        'sep_arcmin']):
    if mjd is not None:
        time_obs = Time(mjd, format ='mjd')
        decimal_year_of_observation = time_obs.decimalyear
        print(f'decimal year = {decimal_year_of_observation} for HST LMC cat proper motion correction!')
    else:
        decimal_year_of_observation = None
    
    cat = hst.hst_catalog(decimal_year_of_observation=decimal_year_of_observation)
    cat.rename_column('ra_deg', 'ra')
    cat.rename_column('dec_deg', 'dec')
    #print(cat.columns)

    pos0 = SkyCoord(ra=ra0, dec=dec0, unit=(u.deg,u.deg), frame='icrs')
    pos1 = SkyCoord(ra=cat['ra'], dec=cat['dec'], unit=(u.deg,u.deg), frame='icrs')
    sep = pos0.separation(pos1)
    #print(sep.degree)
    cat['sep_deg'] = sep
    cat['sep_arcmin'] = sep/60.0
    cat['j_k'] = cat['j_mag_vega'] - cat['k_mag_vega']

    df = cat.to_pandas()
    
    # only keep the values within the radius
    ixs = df.index.values
    Ntot=len(ixs)
    (keep,) = np.where(df.loc[ixs,'sep_deg'].le(radius_deg*radius_factor))
    ixs = ixs[keep]
    print(f'Keeping {len(ixs)} out of {Ntot} of HST LMC cat sources: all positions within {radius_deg*radius_factor:.4f} deg of RA/Dec=({ra0},{dec0})')
    df = df.loc[ixs]
     
    if columns is not None:
        df = df[columns]

    return(df,'ra','dec')


def get_GAIA_sources(ra0,dec0,radius_deg,radius_factor=1.1,
                     mjd=None, # if mjd!=None, then proper motion (pm) adjustement is done
                     pm_median = False, # if pm_median, then the median proper motion is added instead of the individual ones
                     datarelease='dr2',
                     calc_mag_errors=True,
                     rename_mag_colnames=True, # renames from f'phot_{filt}_mean_mag' to f'{filt}'
                     remove_null=True,
                     columns=['source_id','ref_epoch','ra','ra_error','dec','dec_error','pmra','pmra_error','pmdec','pmdec_error',
                              'g','g_err','bp','bp_err','rp','rp_err','bp_rp','bp_g','g_rp','bp_rp_err','bp_g_err','g_rp_err']):
#                     columns=['source_id','ref_epoch','ra','ra_error','dec','dec_error','pmra','pmra_error','pmdec','pmdec_error',
#                              'phot_g_mean_mag','phot_bp_mean_mag','phot_rp_mean_mag','bp_rp','bp_g','g_rp','phot_g_mean_flux','phot_g_mean_flux_error','phot_bp_mean_flux','phot_bp_mean_flux_error','phot_rp_mean_flux','phot_rp_mean_flux_error']):
    if datarelease is None:
        dr = 'gaiadr2'
    elif datarelease.lower() in ['dr2','gaiadr2']:
        dr = 'gaiadr2'
    elif datarelease.lower() in ['edr3','gaiaedr3']:
        dr = 'gaiaedr3'
    elif datarelease.lower() in ['dr3','gaiadr3']:
        dr = 'gaiadr3'
    else:
        raise RuntimeError(f'datarelease {datarelease} not known yet')
    
    query ="SELECT * FROM {}.gaia_source WHERE CONTAINS(POINT('ICRS',\
            {}.gaia_source.ra,{}.gaia_source.dec),\
            CIRCLE('ICRS',{},{} ,{}))=1;".format(dr,dr,dr,ra0,dec0,radius_deg)
    print(f'query:{query}')
    job5 = Gaia.launch_job_async(query)
    tb_gaia = job5.get_results() 
    print("Number of stars:",len(tb_gaia))
    
    if mjd is not None:
        print(f'### Applying proper motion correction to epoch mjd={mjd}')
        time_gaia = Time(tb_gaia['ref_epoch'], format = 'jyear')[0]
        time_obs = Time(mjd, format ='mjd')

        dRA = ((time_obs - time_gaia).to(u.yr).value * tb_gaia['pmra'].data * u.mas / np.cos(np.deg2rad(tb_gaia['dec']))).to(u.deg).value
        dDec = ((time_obs - time_gaia).to(u.yr).value * tb_gaia['pmdec'].data * u.mas).to(u.deg).value

        if pm_median:
            ok = (np.isfinite(dRA)) & (np.isfinite(dDec)) 
            dRA_median = np.median(dRA[ok])
            dDec_median = np.median(dDec[ok])
            print(f'adding median pm dRA={dRA_median} and dDec={dDec_median}')
            tb_gaia['ra1'] = tb_gaia['ra'] + dRA_median
            tb_gaia['dec1'] = tb_gaia['dec'] + dDec_median
            tb_gaia['ref_epoch1'] = time_obs.decimalyear
        else:
            tb_gaia['ra1'] = tb_gaia['ra'] + dRA
            tb_gaia['dec1'] = tb_gaia['dec'] + dDec
            tb_gaia['ref_epoch1'] = time_obs.decimalyear
        if not('ra1' in columns): columns.append('ra1')
        if not('dec1' in columns): columns.append('dec1')
        if not('ref_epoch1' in columns): columns.append('ref_epoch1')
        #tb_gaia['dRA'] = dRA
        #tb_gaia['dDec'] = dDec
        #columns.extend(['dRA','dDec'])
        racol='ra1'
        deccol='dec1'
    else:
        print(f'### NO propoer motion correction!!!')
        racol='ra'
        deccol='dec'
    
    df = tb_gaia.to_pandas()

    # renames columns from f'phot_{filt}_mean_mag' to f'{filt}'
    if rename_mag_colnames:
        for filt in ['g','bp','rp']:
            df.rename(columns={f'phot_{filt}_mean_mag':f'{filt}'},inplace=True)
                              
    if calc_mag_errors:
        for filt in ['g','bp','rp']:
            fluxcol = f'phot_{filt}_mean_flux'
            dfluxcol = f'{fluxcol}_error'
            df[f'{filt}_err'] =  2.5 / math.log(10.0) * df[dfluxcol] / df[fluxcol]

        for (filt1,filt2) in [('bp','rp'),('bp','g'),('g','rp')]:
            df[f'{filt1}_{filt2}'] = df[f'{filt1}'] - df[f'{filt2}']
            df[f'{filt1}_{filt2}_err'] = np.sqrt(np.square(df[f'{filt1}']) - np.square(df[f'{filt2}']))

    if remove_null:
        ixs = df.index.values
        for colname in [racol,deccol]:
            #print('XXX',indices)
            (notnull,) = np.where(pd.notnull(df.loc[ixs,colname]))
            ixs = ixs[notnull]
        df = df.loc[ixs]
        print("Number of stars after removing nan's:",len(df[racol]))
                 
    if columns is not None:
        df = df[columns]
    
    #print(df['source_id'])
    #sys.exit(0)

    return(df,racol,deccol)

    

def get_GAIA_sources_for_image(imagefilename,pm=True,radius_factor=1.1):
    """Return a GAIA catalog/table that is location and time matched to an observation.
    It applies the GAIA proper motion by default (pm=True).
    The (x,y) pixel coordinates of the GAIA source are also returned in  'x' and 'y' column.
    radius_factor is 1.1 by default and hence about 10% that of the FOV but can be set.
    original author: Nor Pirzkal, modified by Armin Rest
    """

    im = fits.open(imagefilename)

    hdr = im['SCI'].header
    nx = hdr['NAXIS1']
    ny = hdr['NAXIS2']

    image_model = ImageModel(im)
 
    ra0,dec0 = image_model.meta.wcs(nx/2.0-1,ny/2.0-1)
    coord0 = SkyCoord(ra0,dec0,unit=(u.deg, u.deg), frame='icrs')
    radius_deg = []
    for x in [0,nx-1]:        
        for y in [0,ny-1]:     
            ra,dec = image_model.meta.wcs(x,y)
            radius_deg.append(coord0.separation(SkyCoord(ra,dec,unit=(u.deg, u.deg), frame='icrs')).deg)
    radius_deg = np.amax(radius_deg)*radius_factor 

    df,racol,deccol = get_GAIA_sources(ra0,dec0,radius_deg,mjd=hdr['MJD-AVG'])
    return(df,racol,deccol)

def get_image_siaf_aperture(aperturename, primaryhdr, scihdr):

    instrument = primaryhdr['INSTRUME']
    apername = primaryhdr['APERNAME']
    print('GGGGGGGGGGGG',apername,aperturename)
    
    ra_ref_orig = scihdr['RA_REF']
    dec_ref_orig = scihdr['DEC_REF']
    roll_ref_orig = scihdr['ROLL_REF']
    
    siaf = pysiaf.Siaf(instrument)

    aper_orig = siaf[apername]

    v2_ref_orig = aper_orig.V2Ref     
    v3_ref_orig = aper_orig.V3Ref     

    attitude = pysiaf.utils.rotations.attitude(v2_ref_orig, v3_ref_orig, ra_ref_orig, dec_ref_orig, roll_ref_orig)

    aper_orig.set_attitude_matrix(attitude)

    image_siaf_aperture = siaf[aperturename]
    image_siaf_aperture.set_attitude_matrix(attitude)
      
    return image_siaf_aperture

def radec_to_idlX(x,y, aperturename, primaryhdr, scihdr):
    # Method from Mario
    instrument = primaryhdr['INSTRUME']
    siaf = pysiaf.Siaf(instrument)   ### this line will slow you down, better if you can pass the SIAF object directly (or read it from the main)
    apername = primaryhdr['APERNAME']
    ap= siaf.apertures[apername]
    xidl,yidl = ap.sci_to_idl(x+1,y+1)
    return xidl, yidl


def radec_to_idl(ra, dec, aperturename, primaryhdr, scihdr):

    image_siaf_aperture = get_image_siaf_aperture(aperturename, primaryhdr, scihdr)

    x_idl, y_idl      = image_siaf_aperture.convert(ra, dec, 'sky', 'idl')
  
    return x_idl, y_idl

"""
def xy_to_idl(x, y, aperturename, primaryhdr, scihdr):

    image_siaf_aperture = get_image_siaf_aperture(aperturename, primaryhdr, scihdr)

    x_idl, y_idl      = image_siaf_aperture.det(x, y, 'det', 'idl')
  
    return x_idl, y_idl
"""

def xy_to_idl(x,y, aperturename, primaryhdr, scihdr):
    # Method from Mario
    instrument = primaryhdr['INSTRUME']
    siaf = pysiaf.Siaf(instrument)   ### this line will slow you down, better if you can pass the SIAF object directly (or read it from the main)
    apername = primaryhdr['APERNAME']
    ap= siaf.apertures[apername]
    xidl,yidl = ap.sci_to_idl(x+1,y+1)
    return xidl, yidl

class jwst_photclass(pdastrostatsclass):
    def __init__(self):
        pdastrostatsclass.__init__(self)
        
        self.filters = {}
        self.filters['NIRCAM'] = ['F070W', 'F090W', 'F115W', 'F140M', 'F150W2', 'F150W', 'F162M', 'F164N', 'F182M',
                                  'F187N', 'F200W', 'F210M', 'F212N', 'F250M', 'F277W', 'F300M', 'F322W2', 'F323N',
                                  'F335M', 'F356W', 'F360M', 'F405N', 'F410M', 'F430M', 'F444W', 'F460M', 'F466N', 'F470N', 'F480M']
        self.filters['NIRISS'] = ['F090W', 'F115W', 'F140M', 'F150W', 'F158M', 'F200W', 'F277W', 'F356W', 'F380M', 'F430M', 'F444W', 'F480M']
        self.filters['FGS'] = ['NA']

        self.psf_fwhm = {}
        self.psf_fwhm['NIRCAM'] = [0.987, 1.103, 1.298, 1.553, 1.628, 1.770, 1.801, 1.494, 1.990, 2.060, 2.141, 2.304, 2.341, 1.340,
                                   1.444, 1.585, 1.547, 1.711, 1.760, 1.830, 1.901, 2.165, 2.179, 2.300, 2.302, 2.459, 2.507, 2.535, 2.574]
        self.psf_fwhm['NIRISS'] = [1.40, 1.40, 1.50, 1.50, 1.50, 1.50, 1.50, 1.60, 1.70, 1.80, 1.80, 1.80]
        self.psf_fwhm['FGS'] = [1.50]

        self.dict_utils = {}
        for instrument in self.filters:
            self.dict_utils[instrument.upper()] = {self.filters[instrument.upper()][i]: {'psf fwhm': self.psf_fwhm[instrument.upper()][i]} for i in range(len(self.filters[instrument]))}

        self.imagename = None
        self.imagetype = None
        self.im=None
        self.primaryhdr = None
        self.scihdr = None
        self.DNunits = None
        self.data=None
        self.data_bkgsub=None
        self.mask=None
        self.found_stars=None
        
        # define the radii of the aperture in units of fwhm
        self.radii_Nfwhm = [2.0]
        self.radius_Nfwhm_for_mag = self.radii_Nfwhm[0]
        self.radius_Nfwhm_sky_in = self.radii_Nfwhm[0] * 2.0
        self.radius_Nfwhm_sky_out = self.radius_Nfwhm_sky_in + 2.0
        
        # These values will be set when the photometry is run! They are
        # set depending on filter and self.radi*_Nfwhm*
        self.radii_px=None
        self.radius_sky_in_px=None
        self.radius_sky_out_px=None
        self.radius_for_mag_px=None
        
        self.instrument = None
        self.aperture = None
        
        self.refcat_short = None
        self.refcat_racol = None
        self.refcat_deccol = None
        self.refcat_xcol = None
        self.refcat_ycol = None

        
        #self.phot=pdastroclass()

    def define_options(self,parser=None,usage=None,conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)

        # default directory for output
        if 'JWST_OUTROOTDIR' in os.environ:
            outrootdir = os.environ['JWST_OUTROOTDIR']
        else:
            outrootdir = None

        parser.add_argument('image',  help='input image file')
        
        parser.add_argument('--photfilename', default='auto', help='output filename for photometry catalog. if "auto", then the fits of the image filename is substituted with phot.txt. If outrootdir and/or outsubdir is not None, then they are used for the path. (default=%(default)s)')
        parser.add_argument('--outrootdir', default=outrootdir, help='output root directory. The output directoy for the photometry file is the output root directory + the outsubdir if not None (default=%(default)s)')
        parser.add_argument('--outsubdir', default=None, help='outsubdir added to output root directory (default=%(default)s)')
        parser.add_argument('--overwrite', default=False, action='store_true', help='overwrite files if they exist.')

        parser.add_argument('-v','--verbose', default=0, action='count')

        return(parser)

        
    def load_image(self, imagename, imagetype=None, DNunits=False, use_dq=False,skip_preparing=False):
        self.imagename = imagename
        self.im = fits.open(imagename)
        self.primaryhdr = self.im['PRIMARY'].header
        self.scihdr = self.im['SCI'].header
        
        print(self.im.info())
        
        self.instrument = self.primaryhdr['INSTRUME']
        self.aperture  = self.primaryhdr['APERNAME']
        
        if self.verbose: print(f'Instrument: {self.instrument}, aperture:{self.aperture}')
        
        if imagetype is None:
            if re.search('cal\.fits$|tweakregstep\.fits$|assignwcsstep\.fits$',imagename):
                self.imagetype = 'cal'
            elif re.search('i2d\.fits$',imagename):
                self.imagetype = 'i2d'
            else:
                raise RuntimeError(f'Unknown image type for file {imagename}')
        else:
            self.imagetype = imagetype
            
        if not skip_preparing:
            # prepare the data.
            if self.imagetype == 'cal':
                #print('VVVVV',self.im.info())
                #sys.exit(0)
                dq=None
                if use_dq: 
                    dq = self.im['DQ'].data
                    print('Using DQ extension!!')
                (self.data,self.mask,self.DNunits) = self.prepare_image(self.im['SCI'].data, self.im['SCI'].header,
                                                                        area = self.im['AREA'].data,
                                                                        dq = dq,
                                                                        DNunits=DNunits)
            elif self.imagetype == 'i2d':
                (self.data,self.mask,self.DNunits) = self.prepare_image(self.im['SCI'].data, self.im['SCI'].header,
                                                                        DNunits=DNunits)
            else:
                raise RuntimeError(f'image type {self.imagetype} not yet implemented!')
        

    def prepare_image(self,data_original, imhdr, area=None, dq=None, 
                      DNunits=False, dq_ignore_bits = 2+4):
        # dq_ignore_bits contains the bits in the dq which are still ok, so they
        # should be ignored.
        # 2 = Pixel saturated during integration
        # 4 = Jump detected during integration
        
        if area is not None:
        
            if self.verbose: print('Applying Pixel Area Map')
        
            data_pam = data_original * area
            
        else:
            
            data_pam = data_original
        
        if DNunits:
            if imhdr["BUNIT"]!='MJy/sr':
                raise RuntimeError(f'imhdr["BUNIT"]={imhdr["BUNIT"]} is not MJy/sr!')
            if self.verbose: print(f'Converting units from {imhdr["BUNIT"]} to DN/s')
            data = data_pam / imhdr['PHOTMJSR']
            
        else:
            data = data_pam
        
        
        if dq is not None:
            # dq_ignore_bits are removed from the mask!
            #fits.writeto('TEST_dq_delme.fits',dq,overwrite=True,output_verify='ignore')
            mask = np.bitwise_and(dq,np.full(data_original.shape, ~dq_ignore_bits, dtype='int'))
        else:
            mask = np.zeros(data_original.shape, dtype='int')

        # hijack a few bits for our purposes...
        mask[np.isnan(data)==True] = 8
        mask[np.isfinite(data)==False] = 8
        mask[np.where(data==0)] = 16
        #fits.writeto('TEST_mask_delme.fits',mask,overwrite=True,output_verify='ignore')
               
        data[mask>0] = np.nan

        return data,mask,DNunits
    
    def get_bool_mask(self, data,  mask=None):
        if mask is None:  
            boolmask = np.full(np.shape(data), False, dtype=bool)
        else:
            boolmask = np.where(mask>0, True, False)
            
       
        boolmask[np.isnan(data)==True] = True
        boolmask[np.isfinite(data)==False] = True         
        
        return(boolmask)
        
            
    def calc_bkg(self, data, mask=None, var_bkg=False):
    
        bkgrms = MADStdBackgroundRMS()
        mmm_bkg = MMMBackground()
       
        
        if var_bkg:
            print('Using 2D Background')
            sigma_clip = SigmaClip(sigma=3.)
            
            boolmask = self.get_bool_mask(data,mask=mask)
    
            bkg = Background2D(data, (500, 500), filter_size=(3, 3), sigma_clip=sigma_clip, bkg_estimator=mmm_bkg,
                               coverage_mask=boolmask, fill_value=0.0)
    
            data_bkgsub = data.copy()
            data_bkgsub = data_bkgsub - bkg.background
    
            median = bkg.background_median
            std = bkg.background_rms_median
            if self.verbose: print('Background median and rms using Background 2D:', median, std)
            
        else:
            std = bkgrms(data)
            bkg = mmm_bkg(data)
            data_bkgsub = data.copy()
            data_bkgsub -= bkg
            if self.verbose: print('Background and rms using MMMBackground and MADStdBackgroundRMS:', bkg, std)
    
        return data_bkgsub, std
    
    def get_fwhm_psf(self,filt,pupil, instrument=None):
        # in the future, this can be changed to get the values directly from CRDS
        if instrument is None:
            instrument = self.instrument
        if instrument is None:
            raise RuntimeError('Can\'t get FWHM, instrument is not known')
            
        # NIRISS is special: it has some filters in the pupil wheel
        if instrument.upper() == 'NIRISS':
            if re.search('^F',filt) is None:
                if re.search('^F',pupil) is None:
                    raise RuntimeError(f'can\'t figure out the NIRISS filter: {filt} {pupil}')
                else:
                    filt=pupil
            else:
                # all good!
                pass
            
        fwhm_psf = self.dict_utils[instrument.upper()][filt]['psf fwhm']
        return(fwhm_psf)
        
        

    def find_stars(self, threshold=3, var_bkg=False, primaryhdr=None, scihdr=None):
        
        '''
        Parameters
        ----------
        
        threshold : float 
            The absolute image value above which to select sources.
        
        fwhm : float
            The full-width half-maximum (FWHM) of the major axis of the Gaussian kernel in units of pixels.
            
        var_bkg : bool
            Use Background2D (see description above)
            
        '''
        
        if primaryhdr is None: primaryhdr=self.primaryhdr
        if scihdr is None: scihdr=self.scihdr
        
        det = primaryhdr['DETECTOR']
        if det in ['GUIDER1','GUIDER2']:
            filt = 'NA'
            pupil = 'NA'
        else:  
            filt = primaryhdr['FILTER']
            pupil = primaryhdr['PUPIL']
        
        print('Finding stars --- Detector: {d}, Filter: {f}'.format(f=filt, d=det))
        
        #sigma_psf = self.dict_utils[filt]['psf fwhm']
        fwhm_psf = self.get_fwhm_psf(filt,pupil)
    
        print('FWHM for the filter {f}:'.format(f=filt), fwhm_psf, "px")

        zero_mask = np.where(self.data == 0.0,1,0)
        nan_mask  = np.where(np.isnan(self.data),1,0)
        full_mask = np.bitwise_or(nan_mask,zero_mask)
        if self.mask is not None: 
             full_mask = np.bitwise_or(full_mask,self.mask)
        bool_mask = np.where(full_mask>0,True,False)
       
        self.data_bkgsub, std = self.calc_bkg(self.data, mask=bool_mask, var_bkg=var_bkg)

        daofind = DAOStarFinder(threshold=threshold * std, fwhm=fwhm_psf, exclude_border=True)
        self.found_stars = daofind(self.data_bkgsub, mask=bool_mask)
        #found_stars = daofind(data_bkgsub)
                    
        print('')
        print('Number of sources found in the image:', len(self.found_stars))
        print('-------------------------------------')
        print('')
        
        return self.found_stars, self.data_bkgsub
    
    def get_radii_phot(self, filt, pupil, 
                       radii_Nfwhm = None,                       
                       radius_Nfwhm_sky_in = None, 
                       radius_Nfwhm_sky_out = None,
                       radius_Nfwhm_for_mag = None):
        if radii_Nfwhm is None:
            radii_Nfwhm = self.radii_Nfwhm
        if isinstance(radii_Nfwhm,float) or isinstance(radii_Nfwhm,int):
            radii_Nfwhm = [radii_Nfwhm]
            
        if radius_Nfwhm_for_mag is None:
            radius_Nfwhm_for_mag = radii_Nfwhm[0]

        if radius_Nfwhm_sky_in is None:
           radius_Nfwhm_sky_in  = self.radius_Nfwhm_sky_in
           
        if radius_Nfwhm_sky_out is None:
           radius_Nfwhm_sky_out  = self.radius_Nfwhm_sky_out
           
          
        
        fwhm = self.get_fwhm_psf(filt,pupil)
        
        radii = [(fwhm*x) for x in radii_Nfwhm]
        radius_for_mag = fwhm*radius_Nfwhm_for_mag
        if not(radius_for_mag in radii):
            raise RuntimeError(f'radius for mag {radius_for_mag} is not in {radii}')
        
        radius_sky_in = fwhm*radius_Nfwhm_sky_in
        radius_sky_out = fwhm*radius_Nfwhm_sky_out
        
        return(radii,radius_sky_in,radius_sky_out,radius_for_mag)

    def colname(self,basecolname,radius):
        return(f'{basecolname}_{radius:.1f}px')

    #def aperture_phot(self, radius=[3.5], sky_in=7, sky_out=10, add_radius_to_colname=False):
    def aperture_phot(self, filt=None, pupil=None, 
                      radii_Nfwhm = None,
                      radius_Nfwhm_sky_in = None, 
                      radius_Nfwhm_sky_out = None, 
                      radius_Nfwhm_for_mag =None,
                      primaryhdr=None, scihdr=None):
        
        if primaryhdr is None: primaryhdr=self.primaryhdr
        if scihdr is None: scihdr=self.scihdr

        det = primaryhdr['DETECTOR']
        if filt is None:
            if det in ['GUIDER1','GUIDER2']:
                filt = 'NA'
            else:
                filt = primaryhdr['FILTER']
        if pupil is None: 
            if det in ['GUIDER1','GUIDER2']:
                pupil = 'NA'
            else:
                pupil = primaryhdr['PUPIL']

        (self.radii_px,
         self.radius_sky_in_px,
         self.radius_sky_out_px,
         self.radius_for_mag_px) = self.get_radii_phot(filt,pupil,
                                                       radii_Nfwhm = radii_Nfwhm,
                                                       radius_Nfwhm_sky_in = radius_Nfwhm_sky_in, 
                                                       radius_Nfwhm_sky_out = radius_Nfwhm_sky_out, 
                                                       radius_Nfwhm_for_mag = radius_Nfwhm_for_mag)
        if self.verbose: print(f'radii:{self.radii_px}pixels radius_sky_in:{self.radius_sky_in_px} radius_sky_out:{self.radius_sky_out_px}  radius_for_mag:{self.radius_for_mag_px}')

        positions = np.transpose((self.found_stars['xcentroid'], self.found_stars['ycentroid']))
        
        tic = time.perf_counter()
    
        table_aper = Table()
                
        for rad in self.radii_px:
            print(f'Performing aperture photometry for radius r = {rad} px')
            aperture = CircularAperture(positions, r=rad)
            
            annulus_aperture = CircularAnnulus(positions, 
                                               r_in=self.radius_sky_in_px, 
                                               r_out=self.radius_sky_out_px)
            annulus_masks = annulus_aperture.to_mask(method='center')
    
            local_sky_median = []
            local_sky_stdev = []
            
            for annulus_mask in annulus_masks:
                
                annulus_data = annulus_mask.multiply(self.data)
                ok =np.logical_and(annulus_mask.data > 0, np.isfinite(annulus_data))
                if (np.sum(ok) >= 10):
                    annulus_data_1d = annulus_data[ok]
                    mean_sigclip, median_sigclip, stdev_sigclip = sigma_clipped_stats(annulus_data_1d, 
                                                                                     sigma=3.5, maxiters=5)
                    if mean_sigclip < 0 or median_sigclip == 0:
                        median_sigclip = -99.99
                        stdev_sigclip = -9.99
                
                else:
                    median_sigclip = -99.99
                    stdev_sigclip = -9.99
                
                local_sky_median.append(median_sigclip)
                local_sky_stdev.append(stdev_sigclip)
            
            local_sky_median = np.array(local_sky_median)
            local_sky_stdev = np.array(local_sky_stdev)
            
#            if select_dq:        
#                
#                phot = aperture_photometry(data, aperture, method='exact')
            
#            else:
                
#                zero_mask = np.where(data == 0,0,1)
#                nan_mask  = np.where(np.isnan(data),0,1)
#                zero_mask = nan_mask * zero_mask
#        
#                nan_mask = np.where(zero_mask == 0,True,False)
            boolmask = self.get_bool_mask(self.data,mask=self.mask)

            phot = aperture_photometry(self.data, aperture, method='exact', mask=boolmask)
            
            phot['annulus_median'] = local_sky_median
            phot['aper_bkg'] = local_sky_median * aperture.area
            phot['aper_sum_bkgsub'] = phot['aperture_sum'] - phot['aper_bkg']

            error_poisson = np.sqrt(phot['aperture_sum'])
            error_scatter_sky = aperture.area * local_sky_stdev**2
            error_mean_sky = local_sky_stdev**2 * aperture.area**2 / annulus_aperture.area
            fluxerr = np.sqrt(error_poisson**2 + error_scatter_sky + error_mean_sky)
            aper_table = get_jwst_ap_corr_table()
            if filt not in aper_table['filter']:
                raise RuntimeError('Do not have ap corr for filter %s'%filt)
            aper_func = scipy.interpolate.interp1d(np.array([aper_table[aper_table['filter']==filt][x] for x in aper_table.colnames if x.startswith('r')]).flatten(),
                                           np.array([aper_table[aper_table['filter']==filt][x] for x in aper_table.colnames if x.startswith('ap')]).flatten()
                                           )
    
    
            apcorr = aper_func(rad)
            phot['aper_sum_bkgsub']*=apcorr
            phot['magerr'] = 1.086 * fluxerr/phot['aper_sum_bkgsub']
            imwcs = wcs.WCS(self.im['SCI',1])
            
            pixel_scale = wcs.utils.proj_plane_pixel_scales(imwcs)[0]  * imwcs.wcs.cunit[0].to('arcsec')
            flux_units = u.MJy / u.sr * (pixel_scale * u.arcsec)**2
            flux = phot['aper_sum_bkgsub']*flux_units
            phot['mag'] = flux.to(u.ABmag).value
            table_aper.add_column(phot['aperture_sum'], name=self.colname('aper_sum',rad))
            table_aper.add_column(phot['annulus_median'], name=self.colname('annulus_median',rad))
            table_aper.add_column(phot['aper_bkg'], name=self.colname('aper_bkg',rad))
            table_aper.add_column(phot['aper_sum_bkgsub'], name=self.colname('aper_sum_bkgsub',rad))
            table_aper.add_column(fluxerr, name=self.colname('flux_err',rad))
            table_aper.add_column(phot['mag'], name='mag')
            table_aper.add_column(phot['magerr'], name='dmag')

            #if rad == self.radius_for_mag_px:
                #table_aper['mag'] = -2.5 * np.log10(table_aper[self.colname('aper_sum_bkgsub',rad)])
                #table_aper['dmag'] = 1.086 * (table_aper[self.colname('flux_err',rad)] / 
                #                              table_aper[self.colname('aper_sum_bkgsub',rad)])      

        table_aper['x'] = self.found_stars['xcentroid']
        table_aper['y'] = self.found_stars['ycentroid']
        table_aper['sharpness'] = self.found_stars['sharpness']
        table_aper['roundness1'] = self.found_stars['roundness1']
        table_aper['roundness2'] = self.found_stars['roundness2']
        table_aper = table_aper[~np.isnan(table_aper['mag'])]
        toc = time.perf_counter()
        print("Time Elapsed:", toc - tic)
    
        self.t = table_aper.to_pandas()
        return table_aper

    def clean_phottable(self,SNR_min=3.0,indices=None):
        # remove nans
        ixs = self.ix_not_null(['mag','dmag'],indices=indices)
        
        if SNR_min is not None:
            dmag_max = 1.086 * 1.0/SNR_min
            ixs = self.ix_inrange('dmag',None,dmag_max,indices=ixs)
        return(ixs)

    def xy_to_radec(self,xcol='x',ycol='y',racol='ra',deccol='dec',indices=None,
                    xshift=0.0,yshift=0.0):
        ixs = self.getindices(indices=indices)

        image_model = ImageModel(self.im)
        
        print(f'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF {xshift} {yshift}')
        ra,dec = image_model.meta.wcs(self.t.loc[ixs,xcol]+xshift, self.t.loc[ixs,ycol]+yshift)
        coord = SkyCoord(ra, dec, unit='deg')
        
        self.t.loc[ixs,racol] = coord.ra.degree
        self.t.loc[ixs,deccol] = coord.dec.degree

    def radec_to_xy(self,racol='ra',deccol='dec',xcol='x',ycol='y',indices=None):
        ixs = self.getindices(indices=indices)

        image_model = ImageModel(self.im)
        world_to_detector = image_model.meta.wcs.get_transform('world', 'detector')
        (self.t.loc[ixs,xcol], self.t.loc[ixs,ycol]) = world_to_detector(self.t.loc[ixs,racol],self.t.loc[ixs,deccol])

    def radec_to_idl(self,racol='ra', deccol='dec', xcol_idl='x_idl', ycol_idl='y_idl',
                     aperturename=None,
                     primaryhdr=None, scihdr=None, indices=None):

        if primaryhdr is None: primaryhdr=self.primaryhdr
        if scihdr is None: scihdr=self.scihdr

        if aperturename is None:
            aperturename = self.aperture

    
        indices = self.getindices()
    
        x_idl, y_idl      = radec_to_idl(self.t.loc[indices,racol], 
                                         self.t.loc[indices,deccol],
                                         aperturename, 
                                         primaryhdr, scihdr)
        self.t.loc[indices,xcol_idl]=x_idl
        self.t.loc[indices,ycol_idl]=y_idl
  
        return x_idl, y_idl
        
    def xy_to_idl(self,xcol='x', ycol='y', xcol_idl='x_idl', ycol_idl='y_idl',
                     aperturename=None,
                     primaryhdr=None, scihdr=None, indices=None):

        if primaryhdr is None: primaryhdr=self.primaryhdr
        if scihdr is None: scihdr=self.scihdr
    
        if aperturename is None:
            aperturename = self.aperture

        indices = self.getindices()
        print(self.t.loc[indices,xcol])
        x_idl, y_idl      = xy_to_idl(self.t.loc[indices,xcol], 
                                      self.t.loc[indices,ycol],
                                      aperturename, 
                                      primaryhdr, scihdr)
        self.t.loc[indices,xcol_idl]=x_idl
        self.t.loc[indices,ycol_idl]=y_idl
  
        return x_idl, y_idl
    

    def init_refcat(self,refcatname,mjd=None,
                    refcat_racol=None,refcat_deccol=None):
        self.refcat = pdastroclass()

        self.refcat.name = refcatname

        self.refcat.racol = None
        self.refcat.deccol = None
        self.refcat.short = None
        self.refcat.cols2copy = []
        self.refcat.mainfilter = None
        self.refcat.mainfilter_err = None
        self.refcat.maincolor = None
        
        if refcatname.lower()=='gaia':
            if mjd is not None:
                self.refcat.racol = 'ra1'
                self.refcat.deccol = 'dec1'
            else:
                self.refcat.racol = 'ra'
                self.refcat.deccol = 'dec'
            self.refcat.name = refcatname
            self.refcat.short = 'gaia'
            self.refcat.cols2copy = ['source_id','ra_error','dec_error','g','g_err','rp','rp_err','g_rp','g_rp_err']
            self.refcat.mainfilter = 'g'
            self.refcat.mainfilter_err = 'g_err'
            self.refcat.maincolor = 'g_rp'
        elif os.path.basename(refcatname)=='LMC_gaia_DR3.nrcposs':
            #self.refcat.load_spacesep(refcatname)
            self.refcat.racol = 'ra'
            self.refcat.deccol = 'dec'
            self.refcat.name = refcatname
            self.refcat.short = 'gaia'
            self.refcat.cols2copy = ['gaia_mag']
            self.refcat.mainfilter = 'gaia_mag'
            self.refcat.mainfilter_err = None
            self.refcat.maincolor = None
        elif os.path.basename(refcatname).lower()=='hst_lmc':
            self.refcat.racol = 'ra'
            self.refcat.deccol = 'dec'
            self.refcat.name = refcatname
            self.refcat.short = 'hst'
            self.refcat.cols2copy = ['ID','ra_error_mas','dec_error_mas','h_mag_vega','j_mag_vega','k_mag_vega','j_k']
            self.refcat.mainfilter = 'j'
            self.refcat.mainfilter_err = None
            self.refcat.maincolor = 'j_k'
        elif os.path.basename(refcatname)=='hawki':
            self.refcat.racol = 'ra'
            self.refcat.deccol = 'dec'
            self.refcat.name = refcatname
            self.refcat.short = 'hawki'
            self.refcat.cols2copy = ['ID','ra_error_mas','dec_error_mas','J2mag','K2mag','J2_K2']
            self.refcat.mainfilter = 'J2mag'
            self.refcat.mainfilter_err = None
            self.refcat.maincolor = 'J2_K2'
        else:
            if os.path.isfile(refcatname):
                # give default column names
                self.refcat.racol = 'ra'
                self.refcat.deccol = 'dec'
                self.refcat.name = os.path.basename(refcatname)
                self.refcat.short = 'reffile'
                self.refcat.cols2copy = []
                self.refcat.mainfilter = 'f814w'
                self.refcat.mainfilter_err = None
                self.refcat.maincolor = None
            else:
                raise RuntimeError(f'Don\'t know waht to do with reference catalog {refcatname}! Not a known refcat, and not a file!')
                
        if refcat_racol is not None:
            self.refcat.racol = refcat_racol
        if refcat_deccol is not None:
            self.refcat.deccol = refcat_deccol
                
        return(0)

    def load_refcat(self,refcatname,
                    ra0=None,dec0=None,radius_deg=None,
                    mjd=None,
                    pm_median=False,
                    refcat_racol=None,refcat_deccol=None,pmflag=True):

        # initialize the refcat, and set racol,deccol and other
        # parameters depending on the choice of refcat
        self.init_refcat(refcatname,mjd=mjd,
                         refcat_racol=refcat_racol,refcat_deccol=refcat_deccol)
        print('RA/Dec columns in reference catalog: ',self.refcat.racol,self.refcat.deccol)
        
        if refcatname.lower()=='gaia':
            self.refcat.t,self.refcat.racol,self.refcat.deccol = get_GAIA_sources(ra0,dec0,radius_deg,mjd=mjd,pm_median=pm_median)
            self.refcat.t['ID']=self.refcat.getindices()
        elif os.path.basename(refcatname)=='LMC_gaia_DR3.nrcposs':
            self.refcat.load_spacesep(refcatname)
            self.refcat.t['ID']=self.refcat.getindices()
        elif os.path.basename(refcatname).lower()=='hst_lmc':
            self.refcat.t,self.refcat.racol,self.refcat.deccol = get_hst_cat(ra0,dec0,radius_deg,mjd=mjd)
        elif os.path.basename(refcatname)=='hawki':
            if (mjd is not None) or pm_median:
                raise RuntimeError('Cannot correct for proper motion with catalog {refcatname}')
            self.refcat.t,self.refcat.racol,self.refcat.deccol = get_hawki_cat(ra0,dec0,radius_deg)
        else:
            if os.path.isfile(refcatname):
                print(f'LOADING refcat {refcatname}')
                self.refcat.t = get_refcat_file(refcatname,racol=self.refcat.racol,deccol=self.refcat.deccol)
                self.refcat.t['ID']=self.refcat.getindices()
            else:
                raise RuntimeError(f'Don\'t know waht to do with reference catalog {refcatname}! Not a known refcat, and not a file!')
                        
        if refcat_racol is not None:
            self.refcat.racol = refcat_racol
        if refcat_deccol is not None:
            self.refcat.deccol = refcat_deccol
        
        return(0)
    
    def getrefcatcolname(self,colname,refcatshort=None,
                         requiredflag=False,existsflag=True):
        """
        returns teh column name in the main photometry catalog
        that got copied from the reference catalog
        In general, this is refcatshort

        Parameters
        ----------
        colname : string
            column name in the reference catalog.
        refcatshort : string, optional
            short name of reference catalog that is used as prefix for the column names. The default is None.
            If None, then refcatshort is set to self.refcatshort
        requiredflag : boolean, optional
            colname cannot be None, and a RuntimeError is thrown if the table does not conain the column
        existsflag : boolean, optional
            if colname is not None, and a RuntimeError is thrown if the table does not conain the column
            
        Returns
        -------
        None.

        """
        if refcatshort is None:
            refcatshort = self.refcatshort
        if colname is not None:
           colname  = f'{refcatshort}_{colname}'
           if not(colname in self.t.columns):
               if requiredflag or existsflag:
                   raise RuntimeError(f'Column {colname} is required, but is not one of the Table columns {self.t.columns}!')
        else:
            if requiredflag:
                raise RuntimeError('This column is required and cannot be set to None!')
        return(colname)
        
    def set_important_refcatcols(self,refcatshort=None):
        """
        This routine sets the important reference cat column names 
        in the main photometric catalog (i.e., with the short name of the 
        reference cat as prefix), based on the info saved in the reference catalog
        
        Parameters
        ----------
        refcatshort : string, optional
            short name of reference catalog that is used as prefix for the column names. The default is None.
            If None, then refcatshort is set to self.refcat.short
        Raises
        ------
        RuntimeError
            DESCRIPTION.

        Returns
        -------
        None.

        """
        if refcatshort is None:
            refcatshort = self.refcat.short
            
        if refcatshort is None:
            raise RuntimeError('The short name of the reference catalog is not None, cannot define the important reference catalog columns in the matched photometric catalog!')
        
        self.refcatshort = refcatshort
        self.refcat_racol = self.getrefcatcolname('ra',requiredflag=True)
        self.refcat_deccol = self.getrefcatcolname('dec',requiredflag=True)
        self.refcat_xcol = self.getrefcatcolname('x',requiredflag=True)
        self.refcat_ycol = self.getrefcatcolname('y',requiredflag=True)

        self.refcat_mainfilter = self.getrefcatcolname(self.refcat.mainfilter,existsflag=True)
        self.refcat_mainfilter_err = self.getrefcatcolname(self.refcat.mainfilter_err,existsflag=True)
        self.refcat_maincolor = self.getrefcatcolname(self.refcat.maincolor,existsflag=True)
            
    def match_refcat(self,
                     max_sep = 1.0,
                     borderpadding=40,
                     refcatshort=None,
                     aperturename=None,
                     primaryhdr=None, 
                     scihdr=None,
                     indices=None):
        print(f'Matching reference catalog {self.refcat.name}')

        if refcatshort is None: refcatshort = self.refcat.short

        if primaryhdr is None: primaryhdr=self.primaryhdr
        if scihdr is None: scihdr=self.scihdr
        
        if aperturename is None:
            aperturename = self.aperture

        # make sure there are no NaNs        
        ixs_obj = self.ix_not_null(['ra','dec'],indices=indices)
        # get SkyCoord objects (needed for matching)
        objcoord = SkyCoord(self.t.loc[ixs_obj,'ra'],self.t.loc[ixs_obj,'dec'], unit='deg')

        
        # find the x_idl and y_idl range, so that we can cut down the objects from the outside catalog!!
        xmin = self.t.loc[ixs_obj,'x_idl'].min()
        xmax = self.t.loc[ixs_obj,'x_idl'].max()
        ymin = self.t.loc[ixs_obj,'y_idl'].min()
        ymax = self.t.loc[ixs_obj,'y_idl'].max()
        print(f'image objects are in x_idl=[{xmin:.2f},{xmax:.2f}] and y_idl=[{ymin:.2f},{ymax:.2f}] range')

        #### gaia catalog
        # get ideal coords into table
        #self.refcat.t['x_idl'], self.refcat.t['y_idl'] = radec_to_idl(self.refcat.t[self.refcat.racol], 
        #                                                              self.refcat.t[self.refcat.deccol],
        #                                                              aperturename, 
        #                                                              primaryhdr, scihdr)
        #xmin = self.refcat.t['x_idl'].min()
        #xmax = self.refcat.t['x_idl'].max()
        #ymin = self.refcat.t['y_idl'].min()
        #ymax = self.refcat.t['y_idl'].max()
        #print(f'refcat objects are in x_idl=[{xmin:.2f},{xmax:.2f}] and y_idl=[{ymin:.2f},{ymax:.2f}] range')

        # cut down to the objects that are within the image
        #ixs_cat = self.refcat.ix_inrange('x_idl',xmin-max_sep,xmax+max_sep)
        #ixs_cat = self.refcat.ix_inrange('y_idl',ymin-max_sep,ymax+max_sep,indices=ixs_cat)
        #print(f'Keeping {len(ixs_cat)} out of {len(self.refcat.getindices())} catalog objects')
        #ixs_cat = self.refcat.ix_not_null([self.refcat.racol,self.refcat.deccol],indices=ixs_cat)
        #print(f'Keeping {len(ixs_cat)}  after removing NaNs from ra/dec')
        
        image_model = ImageModel(self.im)
        world_to_detector = image_model.meta.wcs.get_transform('world', 'detector')
        (self.refcat.t['x'], self.refcat.t['y']) = world_to_detector(self.refcat.t[self.refcat.racol],self.refcat.t[self.refcat.deccol])

        #xmin = self.refcat.t['x'].min()
        #xmax = self.refcat.t['x'].max()
        #ymin = self.refcat.t['y'].min()
        #ymax = self.refcat.t['y'].max()
        #print(f'refcat objects are in x=[{xmin:.2f},{xmax:.2f}] and y=[{ymin:.2f},{ymax:.2f}] range')
        
        #sys.exit(0)

        # cut down to the objects that are within the image
        xmin = 0.0-borderpadding
        xmax = self.scihdr['NAXIS1']+borderpadding
        ymin = 0.0-borderpadding
        ymax = self.scihdr['NAXIS2']+borderpadding
        
        ixs_cat = self.refcat.ix_inrange('x',xmin,xmax)
        ixs_cat = self.refcat.ix_inrange('y',ymin,ymax,indices=ixs_cat)
        print(f'Keeping {len(ixs_cat)} out of {len(self.refcat.getindices())} catalog objects within x={xmin}-{xmax} and y={ymin}-{ymax}')
        ixs_cat = self.refcat.ix_not_null([self.refcat.racol,self.refcat.deccol],indices=ixs_cat)
        print(f'Keeping {len(ixs_cat)}  after removing NaNs from ra/dec')

        if len(ixs_cat) == 0:
            print('WARNING!!!! 0 Gaia sources from catalog within the image bounderies! skipping the rest of the steps calculating x,y of the Gaia sources etc... ')
            return(0)

        # Get the detector x,y position
        image_model = ImageModel(self.im)
        world_to_detector = image_model.meta.wcs.get_transform('world', 'detector')
        self.refcat.t.loc[ixs_cat,'x'], self.refcat.t.loc[ixs_cat,'y'] = world_to_detector(self.refcat.t.loc[ixs_cat,self.refcat.racol],self.refcat.t.loc[ixs_cat,self.refcat.deccol])

        refcatcoord = SkyCoord(self.refcat.t.loc[ixs_cat,self.refcat.racol],self.refcat.t.loc[ixs_cat,self.refcat.deccol], unit='deg')
    
        #idx, d2d, _ = match_coordinates_sky(self.t.loc[ixs_obj,'coord'], self.refcat.t.loc[ixs_cat,'coord'])
        idx, d2d, _ = match_coordinates_sky(objcoord,refcatcoord)
        # ixs_cat4obj has the same length as ixs_obj
        # for each object in ixs_obj, it contains the index to the self.refcat entry
        ixs_cat4obj = ixs_cat[idx]


        # copy over the relevant columns from refcat. The columns are preceded with '{refcatshort}_'
        cols2copy = [self.refcat.racol,self.refcat.deccol,'x','y','ID']
        if self.refcat.mainfilter is not None:
            cols2copy.append(self.refcat.mainfilter)
        if self.refcat.mainfilter_err is not None:
            cols2copy.append(self.refcat.mainfilter_err)
        if self.refcat.maincolor is not None:
            cols2copy.append(self.refcat.maincolor)
        cols2copy.extend(self.refcat.cols2copy)
        cols2copy = unique(cols2copy)

        #self.refcatshort = refcatshort
        #self.refcat_racol = None
        #self.refcat_deccol = None
        for refcat_col in cols2copy:
            if not (refcat_col in self.refcat.t.columns):
                raise RuntimeError(f'Trying to copy column {refcat_col}, but this column is not in {self.refcat.t.columns}')
            if refcat_col == self.refcat.racol:
                obj_col = f'{refcatshort}_ra'
                #self.refcat_racol = f'{refcatshort}_ra'
            elif refcat_col == self.refcat.deccol:
                obj_col = f'{refcatshort}_dec'
                #self.refcat_deccol = f'{refcatshort}_dec'
            else:
                obj_col = f'{refcatshort}_{refcat_col}'
                
            self.t.loc[ixs_obj,obj_col]=list(self.refcat.t.loc[ixs_cat4obj,refcat_col])
            if refcat_col == 'source_id':
                #print('############################################ converting source_id')
                #bla_ixs = self.ix_not_null(obj_col)
                #print(f'XXXXXXX {len(bla_ixs)} {len(self.t[obj_col])}')
                #print(self.t[obj_col])
                self.t[obj_col]=self.t[obj_col].astype(pd.Int64Dtype())
        # also add d2d
        self.t.loc[ixs_obj,f'{refcatshort}_d2d']=d2d.arcsec

        #self.refcat_xcol = f'{refcatshort}_x'
        #self.refcat_ycol = f'{refcatshort}_y'

        self.set_important_refcatcols(refcatshort=refcatshort)

        return(0)
        
    def get_photfilename(self,photfilename=None,
                         outrootdir=None,
                         outsubdir=None,
                         imagename=None):
        if photfilename is None:
            return(None)

        if photfilename.lower() == 'none':
            return(None)
        
        if photfilename.lower() == 'auto':
            if imagename is None:
                raise RuntimeError(f'could not get photfilename from {imagename}')
                
            photfilename = re.sub('\.fits$','.phot.txt',imagename)
            if photfilename == imagename:
                raise RuntimeError(f'could not get photfilename from {self.imagename}')
                
        if outrootdir is not None or outsubdir is not None:
            basename = os.path.basename(photfilename)
            if outrootdir is not None:
                outdir=outrootdir
            else:
                outdir='.'
            
            if outsubdir is not None:
                outdir+=f'/{outsubdir}'

            photfilename = f'{outdir}/{basename}'
         
            
        return(photfilename)
                    
    def get_radecinfo_image(self,im=None,nx=None,ny=None):
        if im is None: im=self.im
        image_model = ImageModel(im)
        if nx is None: nx = int(im['SCI'].header['NAXIS1'])
        if ny is None: ny = int(im['SCI'].header['NAXIS2'])
                
        ra0,dec0 = image_model.meta.wcs(nx/2.0-1,ny/2.0-1)
        coord0 = SkyCoord(ra0,dec0,unit=(u.deg, u.deg), frame='icrs')
        radius_deg = []
        for x in [0,nx-1]:        
            for y in [0,ny-1]:     
                ra,dec = image_model.meta.wcs(x,y)
                radius_deg.append(coord0.separation(SkyCoord(ra,dec,unit=(u.deg, u.deg), frame='icrs')).deg)
        radius_deg = np.amax(radius_deg)

        print(ra0,dec0,radius_deg)
        return(ra0,dec0,radius_deg)

        
    def run_phot(self,imagename, 
                 refcatname=None,
                 refcat_racol=None,
                 refcat_deccol=None,
                 pmflag = False, # apply proper motion
                 pm_median=False,# if pm_median, then the median proper motion is added instead of the individual ones
                 photfilename=None,
                 outrootdir=None,
                 outsubdir=None,
                 overwrite=False,                
                 load_photcat_if_exists=False,
                 rematch_refcat=False,
                 use_dq = False,
                 DNunits=True, 
                 SNR_min=3.0,
                 do_photometry_flag=True,
                 photcat_loaded = False,
                 Nbright4match=None,
                 xshift=0.0,# added to the x coordinate before calculating ra,dec. This can be used to correct for large shifts before matching!
                 yshift=0.0, # added to the y coordinate before calculating ra,dec. This can be used to correct for large shifts before matching!
                 ee_radius=None):
        print(f'\n### Doing photometry on {imagename}')

        # get the photfilename. photfilename='auto' removes fits from image name and replaces it with phot.txt
        self.photfilename = self.get_photfilename(photfilename,outrootdir=outrootdir,outsubdir=outsubdir,imagename=imagename)
                       
        print(f'SSSSSSSSSSSSSSSSSS {self.photfilename} xxx {self.outdir}')
        sys.exit(0)

        
        # Load photcat if wanted
        photcat_loaded = False
        if (self.photfilename is not None):
            print(f'photometry catalog filename: {self.photfilename}')
            if os.path.isfile(self.photfilename):
                if load_photcat_if_exists:
                    print(f'photcat {self.photfilename} already exists, loading it instead of recreating it')
                    self.load(self.photfilename)
                    photcat_loaded = True
                    # skip redoing the photometry
                    do_photometry_flag=False
                elif overwrite:
                    print(f'photcat {self.photfilename} already exists, but recreating it since overwrite=True')
                    rmfile(self.photfilename)
                else:
                    raise RuntimeError(f'photcat {self.photfilename} already exists, exiting! if you want to overwrite, use --overwrite option!')
        else:
            print('NO photometry catalog filename')
        
        # load the image, and prepare it. The data and mask are saved in 
        # self.data and self.mask
        self.load_image(imagename,DNunits=DNunits,use_dq=use_dq)
        # only do the photometry if not reloaded
        if do_photometry_flag:
    
            # find the stars, saved in self.found_stars
            self.find_stars()
            
            #aperture phot, saved in self.t
            self.aperture_phot()
             
        
        # get the indices of good stars
        ixs_clean = self.clean_phottable(SNR_min=SNR_min)
        print(f'{len(ixs_clean)} out of {len(self.getindices())} entries remain in photometry table')
        if Nbright4match is not None:
            ixs_sort = self.ix_sort_by_cols(['mag'],indices=ixs_clean)
            ixs_clean = ixs_sort[:Nbright4match]
            print(f'keeping the britghtest {Nbright4match} sources: {len(ixs_clean)} out of {len(self.getindices())} entries remain in photometry table')
        
        # calculate the ra,dec
        self.xy_to_radec(indices=ixs_clean,xshift=xshift,yshift=yshift)
         
        # calculate the ideal coordinates
        #self.radec_to_idl(indices=ixs_clean)
        self.xy_to_idl(indices=ixs_clean)
        # calculate the ideal coordinates
        #self.xy_to_idl(indices=ixs_clean,xcol_idl='x_idl_test',ycol_idl='y_idl_test')
        
        # get the refcat and match it, but only if either the photometry has 
        # been done (not photcat_loaded), or if rematch_refcat is requested
        if (refcatname is not None):
            if pmflag or pm_median: 
                mjd=self.scihdr['MJD-AVG']
            else: 
                mjd=None
            if (rematch_refcat or (not photcat_loaded)):
                if self.verbose: print(f'Getting {refcatname} and matching it')
                (ra0,dec0,radius_deg)=self.get_radecinfo_image()
                self.load_refcat(refcatname,
                                 ra0,dec0,radius_deg,
                                 mjd=mjd,
                                 pm_median=pm_median,
                                 refcat_racol=refcat_racol,
                                 refcat_deccol=refcat_deccol)
                self.match_refcat(indices=ixs_clean)
            else:
                # set refcat parameters like refcat.
                self.init_refcat(refcatname,mjd=mjd,
                                 refcat_racol=refcat_racol,
                                 refcat_deccol=refcat_deccol)
                self.set_important_refcatcols()
                
        
        #self.write(self.get_photfilename()+'.all')
        # save the catalog
        if self.photfilename is not None:
            print(f'Saving {self.photfilename}')
            self.write(self.photfilename,indices=ixs_clean)
        return(self.photfilename)
        
if __name__ == '__main__':
    phot = jwst_photclass()
    parser = phot.define_options()
    args = parser.parse_args()
    
    phot.verbose=args.verbose
    
    phot.run_phot(args.image,
                  refcatname='./LMC_gaia_DR3.nrcposs',
                  photfilename=args.photfilename,
                  outrootdir=args.outrootdir,
                  outsubdir=args.outsubdir,
                  overwrite=args.overwrite,
                  load_photcat_if_exists=True,
                  DNunits=True,
                  use_dq=False,
                  SNR_min=args.SNR_min,
                  xshift=args.xshift,# added to the x coordinate before calculating ra,dec. This can be used to correct for large shifts before matching!
                  yshift=args.yshift # added to the y coordinate before calculating ra,dec. This can be used to correct for large shifts before matching!
                  )
    
