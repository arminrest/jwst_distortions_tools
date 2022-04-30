#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 27 09:21:15 2022

@author: arest, mcorrenti

This is class wrapper around doing simple photometry on a single JWST image
"""

import argparse,glob,re,sys,os,time,copy

import glob as glob

import numpy as np

from astropy.io import fits
from astropy.stats import sigma_clipped_stats, SigmaClip
from astropy.table import Table

from photutils.detection import DAOStarFinder
from photutils.background import MMMBackground, MADStdBackgroundRMS, Background2D
from photutils import CircularAperture, CircularAnnulus, aperture_photometry
from astropy.visualization import simple_norm

import jwst
from jwst.datamodels import ImageModel
import pysiaf

from pdastro import pdastroclass,makepath4file,unique,AnotB,AorB,AandB
from astropy.coordinates import SkyCoord, match_coordinates_sky

def get_image_siaf_aperture(aperturename_instrument, primaryhdr, scihdr):


    instrument = primaryhdr['INSTRUME']
    apername = primaryhdr['APERNAME']
    
    ra_ref_orig = scihdr['RA_REF']
    dec_ref_orig = scihdr['DEC_REF']
    roll_ref_orig = scihdr['ROLL_REF']
    
    siaf = pysiaf.Siaf(instrument)

    aper_orig = siaf[apername]

    v2_ref_orig = aper_orig.V2Ref     
    v3_ref_orig = aper_orig.V3Ref     

    attitude = pysiaf.utils.rotations.attitude(v2_ref_orig, v3_ref_orig, ra_ref_orig, dec_ref_orig, roll_ref_orig)

    aper_orig.set_attitude_matrix(attitude)

    image_siaf_aperture = siaf[aperturename_instrument]
    image_siaf_aperture.set_attitude_matrix(attitude)
      
    return image_siaf_aperture

def radec_to_idl(ra, dec, aperturename_instrument, primaryhdr, scihdr):

    image_siaf_aperture = get_image_siaf_aperture(aperturename_instrument, primaryhdr, scihdr)

    x_idl, y_idl      = image_siaf_aperture.convert(ra, dec, 'sky', 'idl')
  
    return x_idl, y_idl

def xy_to_idl(x, y, aperturename_instrument, primaryhdr, scihdr):

    image_siaf_aperture = get_image_siaf_aperture(aperturename_instrument, primaryhdr, scihdr)

    x_idl, y_idl      = image_siaf_aperture.det(x, y, 'det', 'idl')
  
    return x_idl, y_idl
    

class jwst_photclass(pdastroclass):
    def __init__(self):
        pdastroclass.__init__(self)
        
        self.filters = ['F070W', 'F090W', 'F115W', 'F140M', 'F150W2', 'F150W', 'F162M', 'F164N', 'F182M',
                        'F187N', 'F200W', 'F210M', 'F212N', 'F250M', 'F277W', 'F300M', 'F322W2', 'F323N',
                        'F335M', 'F356W', 'F360M', 'F405N', 'F410M', 'F430M', 'F444W', 'F460M', 'F466N', 'F470N', 'F480M']
        
        self.psf_fwhm = [0.987, 1.103, 1.298, 1.553, 1.628, 1.770, 1.801, 1.494, 1.990, 2.060, 2.141, 2.304, 2.341, 1.340,
                    1.444, 1.585, 1.547, 1.711, 1.760, 1.830, 1.901, 2.165, 2.179, 2.300, 2.302, 2.459, 2.507, 2.535, 2.574]
        
        self.dict_utils = {self.filters[i]: {'psf fwhm': self.psf_fwhm[i]} for i in range(len(self.filters))}

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
        
        #self.phot=pdastroclass()

    def define_options(self,parser=None,usage=None,conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)

        parser.add_argument('image',  help='input image file')

        parser.add_argument('-v','--verbose', default=0, action='count')

        return(parser)

        
    def load_image(self, imagename, imagetype=None, DNunits=False, use_dq=False):
        self.imagename = imagename
        self.im = fits.open(imagename)
        self.primaryhdr = self.im['PRIMARY'].header
        self.scihdr = self.im['SCI'].header

        if imagetype is None:
            if re.search('cal\.fits$|tweakregstep\.fits',imagename):
                self.imagetype = 'cal'
            elif re.search('i2d\.fits$',imagename):
                self.imagetype = 'i2d'
            else:
                raise RuntimeError(f'Unknown image type for file {imagename}')
        else:
            self.imagetype = imagetype
            
        if self.imagetype == 'cal':
            #print('VVVVV',self.im.info())
            #sys.exit(0)
            dq=None
            if use_dq: dq = self.im['DQ'].data
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
    
    def get_fwhm_psf(self,filt):
        # in the future, this can be changed to get the values directly from CRDS
        fwhm_psf = self.dict_utils[filt]['psf fwhm']
        return(fwhm_psf)
        
        

    def find_stars(self, threshold=3, var_bkg=False):
        
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
        
        det = self.primaryhdr['DETECTOR']
        filt = self.primaryhdr['FILTER']
        
        print('Finding stars --- Detector: {d}, Filter: {f}'.format(f=filt, d=det))
        
        #sigma_psf = self.dict_utils[filt]['psf fwhm']
        fwhm_psf = self.get_fwhm_psf(filt)
    
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
    
    def get_radii_phot(self, filt, radii_Nfwhm = None,                       
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
        
        fwhm = self.get_fwhm_psf(filt)
        
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
    def aperture_phot(self, filt=None, radii_Nfwhm = None,
                      radius_Nfwhm_sky_in = None, 
                      radius_Nfwhm_sky_out = None, 
                      radius_Nfwhm_for_mag =None):
        
        if filt is None: filt = self.primaryhdr['FILTER']

        (self.radii_px,
         self.radius_sky_in_px,
         self.radius_sky_out_px,
         self.radius_for_mag_px) = self.get_radii_phot(filt,radii_Nfwhm = radii_Nfwhm,
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
            fluxerr = np.sqrt(error_poisson + error_scatter_sky + error_mean_sky)
      
            table_aper.add_column(phot['aperture_sum'], name=self.colname('aper_sum',rad))
            table_aper.add_column(phot['annulus_median'], name=self.colname('annulus_median',rad))
            table_aper.add_column(phot['aper_bkg'], name=self.colname('aper_bkg',rad))
            table_aper.add_column(phot['aper_sum_bkgsub'], name=self.colname('aper_sum_bkgsub',rad))
            table_aper.add_column(fluxerr, name=self.colname('flux_err',rad))
    
            if rad == self.radius_for_mag_px:
                table_aper['mag'] = -2.5 * np.log10(table_aper[self.colname('aper_sum_bkgsub',rad)])
                table_aper['dmag'] = 1.086 * (table_aper[self.colname('flux_err',rad)] / 
                                              table_aper[self.colname('aper_sum_bkgsub',rad)])      

        table_aper['x'] = self.found_stars['xcentroid']
        table_aper['y'] = self.found_stars['ycentroid']
        table_aper['sharpness'] = self.found_stars['sharpness']
        table_aper['roundness1'] = self.found_stars['roundness1']
        table_aper['roundness2'] = self.found_stars['roundness2']
    
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

    def xy_to_radec(self,xcol='x',ycol='y',racol='ra',deccol='dec',indices=None):
        ixs = self.getindices(indices=indices)

        image_model = ImageModel(self.im)
        
        ra,dec = image_model.meta.wcs(self.t.loc[ixs,xcol], self.t.loc[ixs,ycol])
        coord = SkyCoord(ra, dec, unit='deg')
        
        self.t.loc[ixs,racol] = coord.ra.degree
        self.t.loc[ixs,deccol] = coord.dec.degree

    def radec_to_xy(self,racol='ra',deccol='dec',xcol='x',ycol='y',indices=None):
        ixs = self.getindices(indices=indices)

        image_model = ImageModel(self.im)
        world_to_detector = image_model.meta.wcs.get_transform('world', 'detector')
        (self.t.loc[ixs,xcol], self.t.loc[ixs,ycol]) = world_to_detector(self.t.loc[ixs,racol],self.t.loc[ixs,deccol])

    def radec_to_idl(self,racol='ra', deccol='dec', xcol_idl='x_idl', ycol_idl='y_idl',
                     aperturename_instrument='NRCALL_FULL',
                     primaryhdr=None, scihdr=None, indices=None):

        if primaryhdr is None: primaryhdr=self.primaryhdr
        if scihdr is None: scihdr=self.scihdr
    
        indices = self.getindices()
    
        x_idl, y_idl      = radec_to_idl(self.t.loc[indices,racol], 
                                         self.t.loc[indices,deccol],
                                         aperturename_instrument, 
                                         primaryhdr, scihdr)
        self.t.loc[indices,xcol_idl]=x_idl
        self.t.loc[indices,ycol_idl]=y_idl
  
        return x_idl, y_idl
        
    def xy_to_idl(self,xcol='x', ycol='y', xcol_idl='x_idl', ycol_idl='y_idl',
                     aperturename_instrument='NRCALL_FULL',
                     primaryhdr=None, scihdr=None, indices=None):

        if primaryhdr is None: primaryhdr=self.primaryhdr
        if scihdr is None: scihdr=self.scihdr
    
        indices = self.getindices()
        print(self.t.loc[indices,xcol])
        x_idl, y_idl      = xy_to_idl(self.t.loc[indices,xcol], 
                                      self.t.loc[indices,ycol],
                                      aperturename_instrument, 
                                      primaryhdr, scihdr)
        self.t.loc[indices,xcol_idl]=x_idl
        self.t.loc[indices,ycol_idl]=y_idl
  
        return x_idl, y_idl
    
    def match_gaiacat(self,gaia_catname,
                      max_sep = 1.0,
                      aperturename_instrument='NRCALL_FULL',
                      primaryhdr=None, scihdr=None,indices=None):
        print(f'Matching Gaia catalog {gaia_catname}')

        if primaryhdr is None: primaryhdr=self.primaryhdr
        if scihdr is None: scihdr=self.scihdr

        gaiacat = pdastroclass()
        gaiacat.load_spacesep(gaia_catname)

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
        gaiacat.t['x_idl'], gaiacat.t['y_idl'] = radec_to_idl(gaiacat.t['ra'], 
                                                              gaiacat.t['dec'],
                                                              aperturename_instrument, 
                                                              primaryhdr, scihdr)

        # cut down to the objects that are within the image
        ixs_cat = gaiacat.ix_inrange('x_idl',xmin-max_sep,xmax+max_sep)
        ixs_cat = gaiacat.ix_inrange('y_idl',ymin-max_sep,ymax+max_sep,indices=ixs_cat)
        print(f'Keeping {len(ixs_cat)} out of {len(gaiacat.getindices())} catalog objects')
        ixs_cat = gaiacat.ix_not_null(['ra','dec'],indices=ixs_cat)
        print(f'Keeping {len(ixs_cat)}  after removing NaNs from ra/dec')

        # Get the detector x,y position
        image_model = ImageModel(self.im)
        world_to_detector = image_model.meta.wcs.get_transform('world', 'detector')
        gaiacat.t.loc[ixs_cat,'x'], gaiacat.t.loc[ixs_cat,'y'] = world_to_detector(gaiacat.t.loc[ixs_cat,'ra'],gaiacat.t.loc[ixs_cat,'dec'])

        gaiacoord = SkyCoord(gaiacat.t.loc[ixs_cat,'ra'],gaiacat.t.loc[ixs_cat,'dec'], unit='deg')
    
        #idx, d2d, _ = match_coordinates_sky(self.t.loc[ixs_obj,'coord'], gaiacat.t.loc[ixs_cat,'coord'])
        idx, d2d, _ = match_coordinates_sky(objcoord,gaiacoord)
        # ixs_cat4obj has the same length as ixs_obj
        # for each object in ixs_obj, it contains the index to the gaiacat entry
        ixs_cat4obj = ixs_cat[idx]

        for cat_col in ['ra','dec','x','y','x_idl','y_idl']:
            obj_col = f'cat_{cat_col}'
            vals = list(gaiacat.t.loc[ixs_cat4obj,cat_col])
            self.t.loc[ixs_obj,obj_col]=vals
        self.t.loc[ixs_obj,'gaia_mag']=list(gaiacat.t.loc[ixs_cat4obj,'gaia_mag'])
        self.t.loc[ixs_obj,'cat_index']=ixs_cat4obj
        self.t.loc[ixs_obj,'d2d']=d2d.arcsec
    
        #self.t.loc[ixs_obj,'dx_idl']=self.t.loc[ixs_obj,'cat_x_idl'] - self.t.loc[ixs_obj,'x_idl']
        #self.t.loc[ixs_obj,'dy_idl']=self.t.loc[ixs_obj,'cat_y_idl'] - self.t.loc[ixs_obj,'y_idl']
    
        #ixs_obj_small_d2d = self.ix_inrange('d2d',None,0.2)
        #print(f'{len(ixs_obj_small_d2d)} matches out of {len(ixs_obj)}')
    
        #self.t.loc[ixs_obj_small_d2d].plot.scatter('x_idl','dy_idl')
        #self.t.loc[ixs_obj_small_d2d].plot.scatter('x_idl','dx_idl')
        #self.t.loc[ixs_obj_small_d2d].plot.scatter('y_idl','dy_idl')
        #self.t.loc[ixs_obj_small_d2d].plot.scatter('y_idl','dx_idl')
        #self.t.loc[ixs_obj_small_d2d].plot.scatter('x_idl','d2d')
        
    def get_photfilename(self,photfilename=None):
        if photfilename is not None:
            return(photfilename)
        
        photfilename = re.sub('\.fits$','.phot.txt',self.imagename)
        if photfilename == self.imagename:
            raise RuntimeError(f'could not get photfilename from {self.imagename}')
        
        return(photfilename)
            
        
    def run_phot(self,imagename, gaia_catname=None, DNunits=True, SNR_min=3.0):
        if self.verbose:
            print(f'\n### Doing photometry on {imagename}')
        
        # load the image, and prepare it. The data and mask are saved in 
        # self.data and self.mask
        self.load_image(imagename,DNunits=DNunits,use_dq=True)
        
        # find the stars, saved in self.found_stars
        self.find_stars()
        
        #aperture phot, saved in self.t
        self.aperture_phot()
        
        # get the indices of good stars
        ixs_clean = self.clean_phottable(SNR_min=SNR_min)
        print(f'{len(ixs_clean)} out of {len(self.getindices())} entries remain in photometry table')
        
        # calculate the ra,dec
        self.xy_to_radec(indices=ixs_clean)
        
        # calculate the ideal coordinates
        self.radec_to_idl(indices=ixs_clean)
        # calculate the ideal coordinates
        #self.xy_to_idl(indices=ixs_clean,xcol_idl='x_idl_test',ycol_idl='y_idl_test')
        
        if gaia_catname is not None:
            self.match_gaiacat(gaia_catname,indices=ixs_clean)
        
        #self.write(self.get_photfilename()+'.all')
        # save the catalog
        photfilename = self.get_photfilename()
        print(f'Saving {photfilename}')
        self.write(photfilename,indices=ixs_clean)
        
if __name__ == '__main__':
    phot = jwst_photclass()
    parser = phot.define_options()
    args = parser.parse_args()
    
    phot.verbose=args.verbose
    
    phot.run_phot(args.image,'./LMC_gaia_DR3.nrcposs',SNR_min=10.0)
    
