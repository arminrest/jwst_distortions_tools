#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 14:32:42 2022

@author: arest, bhilbert, mcorrenti, acanipe, jpierel
"""

import os,re,sys
from pdastro import makepath,rmfile,pdastroclass,AnotB
from simple_jwst_phot import jwst_photclass
#from jwst.tweakreg import TweakRegStep
import tweakreg_hack
import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table

from jwst import datamodels

from apply_distortions_single_image import apply_distortion_singleim

# for a given catalog phot, calculate dx, dy, and make some rough cuts on dmag, d2d, and Nbright
def calc_dxdy(phot,refcatshort,
              d2d_max=None,dmag_max=None,Nbright=None, ixs=None):
    ixs = phot.getindices(ixs)
    if d2d_max is not None:
        ixs = phot.ix_inrange('d2d',None,3*d2d_max,indices=ixs)
    if dmag_max is not None:
        ixs = phot.ix_inrange('dmag',None,dmag_max,indices=ixs)
    if Nbright is not None:
        ixs_sort = phot.ix_sort_by_cols(['mag'],indices=ixs)
        ixs = ixs_sort[:Nbright]
        phot.write(columns=['mag'],indices=ixs)
        
    print(f'# of matched objects that pass initial cuts: {len(ixs)}')
    #phot.t.loc[ixs,'ddec'] = 3600.0*(phot.t.loc[ixs,'cat_dec'] - phot.t.loc[ixs,'dec'])
    #phot.t.loc[ixs,'dra'] = 3600.0*(phot.t.loc[ixs,'cat_ra'] - phot.t.loc[ixs,'ra'])*np.cos(np.deg2rad(phot.t.loc[ixs,'dec']))
    phot.t.loc[ixs,'dx'] = phot.t.loc[ixs,f'{refcatshort}_x'] - phot.t.loc[ixs,'x']
    phot.t.loc[ixs,'dy'] = phot.t.loc[ixs,f'{refcatshort}_y'] - phot.t.loc[ixs,'y']
    return(ixs)

def initplot(nrows=1, ncols=1, figsize4subplot=5, **kwargs):
    sp=[]
    xfigsize=figsize4subplot*ncols
    yfigsize=figsize4subplot*nrows
    plt.figure(figsize=(xfigsize,yfigsize))
    counter=1
    for row in range(nrows):
        for col in range(ncols):
            sp.append(plt.subplot(nrows, ncols, counter,**kwargs))
            counter+=1

    for i in range(len(sp)):
        plt.setp(sp[i].get_xticklabels(),'fontsize',12)
        plt.setp(sp[i].get_yticklabels(),'fontsize',12)
        sp[i].set_xlabel(sp[i].get_xlabel(),fontsize=14)
        sp[i].set_ylabel(sp[i].get_ylabel(),fontsize=14)
        sp[i].set_title(sp[i].get_title(),fontsize=14)

    return(sp)
     

# These are the info plots to check out dx, dy for good and bad detections
def infoplots(phot,ixs_good,ixs_bad,dy_plotlim=(-4,4),dx_plotlim=(-4,4)):

    sp = initplot(2,3)
    
    phot.t.loc[ixs_bad].plot.scatter('y','dx',ax=sp[0],ylim=dx_plotlim,color='red')
    phot.t.loc[ixs_good].plot.scatter('y','dx',ax=sp[0],ylim=dx_plotlim,ylabel='dx in pixels')
    phot.t.loc[ixs_bad].plot.scatter('x','dy',ax=sp[1],ylim=dy_plotlim,color='red')
    phot.t.loc[ixs_good].plot.scatter('x','dy',ax=sp[1],ylim=dy_plotlim,ylabel='dy in pixels')

    phot.t.loc[ixs_bad].plot.scatter('x','y',ax=sp[2],color='red')
    phot.t.loc[ixs_good].plot.scatter('x','y',ax=sp[2])
    
    phot.t.loc[ixs_bad].plot.scatter('sharpness','mag',ax=sp[3],color='red')
    phot.t.loc[ixs_good].plot.scatter('sharpness','mag',ax=sp[3])

    phot.t.loc[ixs_bad].plot.scatter('sharpness','dmag',ax=sp[4],color='red')
    phot.t.loc[ixs_good].plot.scatter('sharpness','dmag',ax=sp[4])
    
    plt.tight_layout()

    return(sp)

# plot the rotated dx or dy versus the original one
def plot_rotated(phot,ixs,d_col,col,
                 d_col_rot='__d_rot_tmp',
                 histolim=(-28,-12),
                 bins=None,
                 bin_weights_flag=False,
                 title=None):
    sp=initplot(1,2)

    phot.t.loc[ixs].plot.scatter(col,d_col,ax=sp[0],ylim=histolim,color='red',title=title)
    phot.t.loc[ixs].plot.scatter(col,d_col_rot,ax=sp[0],ylim=histolim,color='blue',ylabel=f'{d_col} in pixels')
    
    if bins is not None:
        if bin_weights_flag:
            phot.t.loc[ixs,d_col_rot].plot.hist(ax=sp[1],bins=bins,weights=phot.t.loc[ixs,'__weights'],xlim=histolim,color='blue')
        else:
            phot.t.loc[ixs,d_col_rot].plot.hist(ax=sp[1],bins=bins,xlim=histolim,color='blue')

    return(sp)

# find the maximum value in yvals, its index, and the corresponding value in xvals
# also indicate if there are multiple entries with the same maximum value (multiple_max=True)
def find_info_for_maxval(xvals,yvals,use_firstindex_if_multiple=True):
    # get the max value of the histogram, and its associated bin
    #print(histo)
    maxval = np.max(yvals)
    ixs_max = np.where(yvals==maxval)
    multiple_max = None
    if (len(ixs_max[0])==0):
        raise RuntimeError('BUUUUGGGG!!!!')
    elif (len(ixs_max[0])>1):
        if not use_firstindex_if_multiple:
            raise RuntimeError(f'More than one indices {ixs_max} for maxval={maxval}')
        #print(f'\nWARNING!! More than one indices {ixs_max} for maxval={maxval}!')
        multiple_max=True
        ix_best=ixs_max[0][0]
    else:
        multiple_max=False
        ix_best=ixs_max[0][0]
    # return ()
    
    # get the (rough) fwhm of the peak
    ix_minus = ix_best-1 
    if ix_minus<0: ix_minus=0
    ix_plus = ix_best+1
    if ix_plus>len(yvals)-1: ix_plus=len(yvals)-1
    while (ix_minus>0):
        if yvals[ix_minus]<=0.5*yvals[ix_best]:
            break
        ix_minus-=1
    while (ix_plus<len(yvals)-1):
        if yvals[ix_plus]<=0.5*yvals[ix_best]:
            break
        ix_plus+=1
    fwhm = xvals[ix_plus]-xvals[ix_minus]
    
    return(xvals[ix_best],yvals[ix_best],ix_best,fwhm,multiple_max)

# straight line. Sloppy: use global slope and intercept
def f(val,slope,intercept):
    return(val*slope+intercept)

def rotate_d_and_find_binmax(phot,ixs,d_col,col,
                             Naxis_px, # Nx or Ny, depending on col
                             d_col_rot='__d_rot_tmp',
                             binsize=0.5,
                             bin_weights_flag=True,
                             slope_min=-10.0/2048.0, # 
                             slope_max=10.0/2048.0, # 
                             slope_stepsize=1.0/2048,
                             showplots=0):
    rot_results = pdastroclass()

    if bin_weights_flag:
        phot.t.loc[ixs,'__weights']=10**(-0.4*phot.t.loc[ixs,'mag'])
    else:
        phot.t.loc[ixs,'__weights']=None
        
    slopes = np.arange(slope_min,slope_max,slope_stepsize)
    for slope in slopes:
        #slope = delta4slope_pix/Nx
        intercept = -0.5*Naxis_px * slope

        phot.t.loc[ixs,d_col_rot] = phot.t.loc[ixs,d_col] - f(phot.t.loc[ixs,col],slope,intercept)

        # get the histogram
        d_rotated = phot.t.loc[ixs,d_col_rot]
        bins = np.arange(np.min(d_rotated),np.max(d_rotated),binsize)
        if bin_weights_flag:
            histo = np.histogram(d_rotated,bins=bins,weights=phot.t.loc[ixs,'__weights'])
        else:
            histo = np.histogram(d_rotated,bins=bins)
        #sp = plt.subplot(111)
        #print(histo[1])
        #sp.plot(histo[1][1:], histo[0], 'ro')
        # get the max value of the histogram, and its associated bin center. Note that the bincenter is 
        # the value in the bins (left edge of the bin) + the half of the binsize
        (bincenter4maxval,maxval,index_maxval,fwhm,multiple_max) = find_info_for_maxval(histo[1]+0.5*binsize,histo[0])

        # Save the results
        rot_results.newrow({'slope':slope,
                            'intercept':intercept,
                            'maxval':maxval,
                            'index':index_maxval,
                            'd_bestguess':bincenter4maxval,
                            'fwhm':fwhm,
                            'multimax':multiple_max
                            })
        # plot it if wanted
        if showplots>2:
            plot_rotated(phot,ixs,
                         d_col,col,
                         d_col_rot=d_col_rot,
                         bins=bins,
                         bin_weights_flag=bin_weights_flag,
                         histolim = (bincenter4maxval-8,bincenter4maxval+8),
                         title=f'slope:{slope}')
        #sys.exit(0)

    # print the results        
    rot_results.write()

    # find the best rotation
    maxmaxval = np.max(rot_results.t['maxval'])
    ixs_maxmax = np.where(rot_results.t['maxval']==maxmaxval)
    if (len(ixs_maxmax[0])==0):
        raise RuntimeError('BUUUUGGGG!!!!')
    elif (len(ixs_maxmax[0])>1):
        #print(f'\nWARNING!! more than one bin with maxvalue={maxmaxval}!')
        multiple_max=True
        best_index=ixs_maxmax[0][0]
    else:
        multiple_max=False
        best_index=ixs_maxmax[0][0]
    print('####BEST:')
    rot_results.write(indices=[best_index])
    return(rot_results,best_index)

def sigmacut_d_rot(phot,ixs,
                   d_col,col,
                   slope,intercept,d_rot_bestguess,
                   rough_cut_px = 2.5, #This is the first rough cut:  get rid of everything d_rot_bestguess+-rough_cut_px
                   d_col_rot='__d_rot_tmp',
                   showplots=0,
                   binsize=0.5,
                   bin_weights_flag=True):

    ### recover the slope and intercept of the best binning
    phot.t.loc[ixs,d_col_rot] = phot.t.loc[ixs,d_col] - f(phot.t.loc[ixs,col],slope,intercept)
    
    # Now make the rough cut! only keep data for with dx_rotated within  d_rot_bestguess+-rough_cut_px
    ixs_roughcut = phot.ix_inrange(d_col_rot,d_rot_bestguess-rough_cut_px,d_rot_bestguess+rough_cut_px,indices=ixs)
    d_rotated = phot.t.loc[ixs,d_col_rot]
    
    if showplots>1:
        bins = np.arange(np.min(d_rotated),np.max(d_rotated),binsize)
        plot_rotated(phot,ixs_roughcut,
                     d_col,col,
                     d_col_rot=d_col_rot,
                     bins=bins,
                     bin_weights_flag=bin_weights_flag,
                     histolim = (d_rot_bestguess-3*rough_cut_px,d_rot_bestguess+3*rough_cut_px),
                     title=f'First rough cut: {d_rot_bestguess:.3f}+-{rough_cut_px:.3f} for slope={slope:.6f}')

    print('\n####################\n### d_rotated cut')
    #ixs_clean4average = phot_clear.ix_inrange(d_col,0,3,indices=ixs_clear_cut)
    phot.calcaverage_sigmacutloop(d_col_rot,verbose=3,indices=ixs_roughcut,percentile_cut_firstiteration=35)
    print(phot.statstring())
    ixs_cut = phot.statparams['ix_good']

    if showplots>1:
        #infoplots(phot,ixs_cut,ixs_roughcut,dy_plotlim=dy_plotlim,dx_plotlim=dx_plotlim)
        sp=initplot(1,2)

        phot.t.loc[ixs_roughcut].plot.scatter(col,d_col_rot,ax=sp[0],color='red')
        phot.t.loc[ixs_cut].plot.scatter(col,d_col_rot,ax=sp[0],ylabel='dx in pixels',title=f'3-sigma cut: {len(ixs_cut)} out of {len(ixs_roughcut)} left')
        
    return(ixs_cut,ixs_roughcut)

"""
# plot the rotated dx or dy versus the original one
def plot_rotated(phot,ixs,d_col,col,
                 d_col_rot='__d_rot_tmp',
                 histolim=(-28,-12),bins=None,title=None):
    sp=initplot(1,2)
    #xfigsize=10.0
    #yfigsize=5.0
    #plt.figure(figsize=(xfigsize,yfigsize))
    #sp.append(plt.subplot(121))
    #sp.append(plt.subplot(122))

    phot.t.loc[ixs].plot.scatter(col,d_col,ax=sp[0],ylim=histolim,color='red',title=title)
    phot.t.loc[ixs].plot.scatter(col,d_col_rot,ax=sp[0],ylim=histolim,color='blue',ylabel=f'{d_col} in pixels')
    
    if bins is not None:
        phot.t.loc[ixs,d_col_rot].plot.hist(ax=sp[1],bins=bins,xlim=histolim,color='blue')

    plt.tight_layout()

    return(sp)

# find the maximum value in yvals, its index, and the corresponding value in xvals
# also indicate if there are multiple entries with the same maximum value (multiple_max=True)
def find_info_for_maxval(xvals,yvals,use_firstindex_if_multiple=True):
    # get the max value of the histogram, and its associated bin
    #print(histo)
    maxval = np.max(yvals)
    ixs_max = np.where(yvals==maxval)
    multiple_max = None
    if (len(ixs_max[0])==0):
        raise RuntimeError('BUUUUGGGG!!!!')
    elif (len(ixs_max[0])>1):
        if not use_firstindex_if_multiple:
            raise RuntimeError(f'More than one indices {ixs_max} for maxval={maxval}')
        #print(f'\nWARNING!! More than one indices {ixs_max} for maxval={maxval}!')
        multiple_max=True
        ix_best=ixs_max[0][0]
    else:
        multiple_max=False
        ix_best=ixs_max[0][0]
    # return ()
    
    # get the (rough) fwhm of the peak
    ix_minus = ix_best-1 
    if ix_minus<0: ix_minus=0
    ix_plus = ix_best+1
    if ix_plus>len(yvals)-1: ix_plus=len(yvals)-1
    while (ix_minus>0):
        if yvals[ix_minus]<=0.5*yvals[ix_best]:
            break
        ix_minus-=1
    while (ix_plus<len(yvals)-1):
        if yvals[ix_plus]<=0.5*yvals[ix_best]:
            break
        ix_plus+=1
    fwhm = xvals[ix_plus]-xvals[ix_minus]
    
    return(xvals[ix_best],yvals[ix_best],ix_best,fwhm,multiple_max)

# straight line. Sloppy: use global slope and intercept
def f(val,slope,intercept):
    return(val*slope+intercept)

def rotate_d_and_find_binmax(phot,ixs,d_col,col,
                             Naxis_px, # Nx or Ny, depending on col
                             d_col_rot='__d_rot_tmp',
                             binsize=0.5,
                             slope_min=-10.0/2048.0, # 
                             slope_max=10.0/2048.0, # 
                             slope_stepsize=1.0/2048,
                             showplots=0):
    rot_results = pdastroclass()

    slopes = np.arange(slope_min,slope_max,slope_stepsize)
    for slope in slopes:
        intercept = -0.5*Naxis_px * slope

        phot.t.loc[ixs,d_col_rot] = phot.t.loc[ixs,d_col] - f(phot.t.loc[ixs,col],slope,intercept)

        # get the histogram
        d_rotated = phot.t.loc[ixs,d_col_rot]
        bins = np.arange(np.min(d_rotated),np.max(d_rotated),binsize)
        histo = np.histogram(d_rotated,bins=bins)

        # get the max value of the histogram, and its associated bin center. Note that the bincenter is 
        # the value in the bins (left edge of the bin) + the half of the binsize
        (bincenter4maxval,maxval,index_maxval,fwhm,multiple_max) = find_info_for_maxval(histo[1]+0.5*binsize,histo[0])

        # Save the results
        rot_results.newrow({'slope':slope,
                            'intercept':intercept,
                            'maxval':maxval,
                            'index':index_maxval,
                            'd_bestguess':bincenter4maxval,
                            'fwhm':fwhm,
                            'multimax':multiple_max
                            })
        # plot it if wanted
        if showplots>1:
            print('BBBB000',showplots)
            plot_rotated(phot,ixs,
                         d_col,col,
                         d_col_rot=d_col_rot,
                         bins=bins,
                         histolim = (bincenter4maxval-8,bincenter4maxval+8),
                         title=f'slope:{slope}')
            plt.show() 

    # print the results        
    rot_results.write()

    # find the best rotation
    maxmaxval = np.max(rot_results.t['maxval'])
    ixs_maxmax = np.where(rot_results.t['maxval']==maxmaxval)
    if (len(ixs_maxmax[0])==0):
        raise RuntimeError('BUUUUGGGG!!!!')
    elif (len(ixs_maxmax[0])>1):
        #print(f'\nWARNING!! more than one bin with maxvalue={maxmaxval}!')
        multiple_max=True
        best_index=ixs_maxmax[0][0]
    else:
        multiple_max=False
        best_index=ixs_maxmax[0][0]
    print('####BEST:')
    rot_results.write(indices=[best_index])
    return(rot_results,best_index)

def sigmacut_d_rot(phot,ixs,
                   d_col,col,
                   slope,intercept,d_rot_bestguess,
                   rough_cut_px = 2.5, #This is the first rough cut:  get rid of everything d_rot_bestguess+-rough_cut_px
                   d_col_rot='__d_rot_tmp',
                   showplots=0,
                   binsize=0.5):

    ### recover the slope and intercept of the best binning
    phot.t.loc[ixs,d_col_rot] = phot.t.loc[ixs,d_col] - f(phot.t.loc[ixs,col],slope,intercept)
    
    # Now make the rough cut! only keep data for with dx_rotated within  d_rot_bestguess+-rough_cut_px
    ixs_roughcut = phot.ix_inrange(d_col_rot,d_rot_bestguess-rough_cut_px,d_rot_bestguess+rough_cut_px,indices=ixs)
    d_rotated = phot.t.loc[ixs,d_col_rot]
    
    if showplots:
        bins = np.arange(np.min(d_rotated),np.max(d_rotated),binsize)
        plot_rotated(phot,ixs_roughcut,
                     d_col,col,
                     d_col_rot=d_col_rot,
                     bins=bins,
                     histolim = (d_rot_bestguess-3*rough_cut_px,d_rot_bestguess+3*rough_cut_px),
                     title=f'First rough cut: {d_rot_bestguess:.3f}+-{rough_cut_px:.3f} for slope={slope:.6f}')
        plt.show() 

    print('\n####################\n### d_rotated cut')
    #ixs_clean4average = phot_clear.ix_inrange(d_col,0,3,indices=ixs_clear_cut)
    phot.calcaverage_sigmacutloop(d_col_rot,verbose=3,indices=ixs_roughcut,percentile_cut_firstiteration=35)
    print(phot.statstring())
    ixs_cut = phot.statparams['ix_good']

    if showplots:
        #infoplots(phot,ixs_cut,ixs_roughcut,dy_plotlim=dy_plotlim,dx_plotlim=dx_plotlim)
        sp=initplot(1,2)

        phot.t.loc[ixs_roughcut].plot.scatter(col,d_col_rot,ax=sp[0],color='red')
        phot.t.loc[ixs_cut].plot.scatter(col,d_col_rot,ax=sp[0],ylabel='dx in pixels',title=f'3-sigma cut: {len(ixs_cut)} out of {len(ixs_roughcut)} left')
        
        plt.show() 

#        plot_rotated(phot,ixs_cut,
#                     d_col,col,
#                     d_col_rot=d_col_rot,
#                     histolim = (d_rot_bestguess-rough_cut_px,d_rot_bestguess+rough_cut_px),
#                     title=f'3-sigma cut: {len(ixs_cut)} out of {len(ixs_roughcut)} left')
    return(ixs_cut,ixs_roughcut)
"""

class jwst_wcs_align(apply_distortion_singleim):
    def __init__(self):
        apply_distortion_singleim.__init__(self)
        self.phot=jwst_photclass()
        
    def define_options(self,parser=None,usage=None,conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)

        parser.add_argument('cal_image',  help='cal.fits filename (_rate.fits can also be passed if distortion file is applied).')

        parser = self.default_options(parser)
        return(parser)

    def default_options(self,parser):
        parser = apply_distortion_singleim.default_options(self,parser=parser)

        parser.add_argument('--distortion_file', default=None, help='distortion file, in asdf format. If not None, distortion file is applied to the cal image')
        parser.add_argument('--skip_applydistortions_if_exists', default=False, action='store_true', help='If the output cal file already exists, skip running the level 2 pipeline to assign the distortion terms, assuming this has already been done.')

        parser.add_argument('--SNR_min', default=10.0, help='mininum SNR for object in image to be used for analysis (default=%(default)s)')

        parser.add_argument('--refcat', default='Gaia', help='reference catalog. Can be a filename or Gaia (default=%(default)s)')
        parser.add_argument('--refcat_racol', default=None, help='RA column of reference catalog. If None, then automatically determined (default=%(default)s)')
        parser.add_argument('--refcat_deccol', default=None, help='Dec column of reference catalog. If None, then automatically determined (default=%(default)s)')
        parser.add_argument('--refcat_no_pmflag', default=False, action='store_true', help='Do NOT apply the proper motion correction (only for catalogs it is applicable, e.g., gaia')
        parser.add_argument('--pm_median', default=False, action='store_true', help=' if pm_median, then the median proper motion is added instead of the individual ones')
        parser.add_argument('--photfilename', default='auto', help='photometry output filename. if "auto", the fits in the image filename is substituted with phot.txt (default=%(default)s)')
#        parser.add_argument('--photfilename', default='auto', help='photometry output filename. if "auto", the fits in the image filename is substituted with phot.txt (default=%(default)s)')

        parser.add_argument('--load_photcat_if_exists', default=False, action='store_true', help='If the photometric catalog file already exists, skip recreating it.')
        parser.add_argument('--rematch_refcat', default=False, action='store_true', help='if --load_photcat_if_exists and the photcat already exists, load the photcat, but rematch with refcat')

        parser.add_argument('--d2d_max', default=None, type=float, help='maximum distance between source in image and refcat object, in arcsec (default=%(default)s)')
        parser.add_argument('--dmag_max', default=0.1, type=float, help='maximum uncertainty of sources in image (default=%(default)s)')
        parser.add_argument('--Nbright', default=None, type=int, help='Use only Nbright brightest objects in image that are matched to refcat for alignment (default=%(default)s)')
        parser.add_argument('-p','--showplots', default=0, action='count',help='showplots=1: most important plots. showplots=2: all plots (debug/test/finetune)')
        parser.add_argument('--saveplots', default=0, action='count',help='saveplots=1: most important plots. saveplots=2: all plots (debug/test/finetune)')
        parser.add_argument('--savephottable', default=0, action='count',help='Save the final photometry table')
        

        return(parser)
        

        
    def run_align2refcat(self,imfilename,
                         phot=None,
                         ixs=None,
                         refcat_short=None,
                         racol=None,
                         deccol=None,
                         xcol='x',
                         ycol='y',
                         outdir=None,overwrite=False, 
                         skip_if_exists=False,
                         savephot=True,
                         ):
            
        if phot is None:
            phot=self.calphot
        if outdir is None:
            outdir=self.outdir
            
        if (racol is None) or (deccol is None): 
            if refcat_short is None:
                refcat_short = phot.refcat.short
                if refcat_short is None:
                    raise RuntimeError('need the refcat shortname to identify ra,dec columns!')
            racol = f'{refcat_short}_ra'
            deccol = f'{refcat_short}_dec'
            


        tweakreg = tweakreg_hack.TweakRegStep()
        cal_image = datamodels.open(imfilename)

        tweakregfilename = re.sub('_[a-zA-Z0-9]+\.fits$','_tweakregstep.fits',os.path.basename(imfilename))
        if tweakregfilename == imfilename:
            raise RuntimeError('could not determine tweakreg filename for {imfilename}')
        tweakregfilename = f'{outdir}/{tweakregfilename}'
        print(f'{tweakregfilename}')
        print(f'Setting output directory for tweakregstep.fits file to {outdir}')
        tweakreg.output_dir = outdir
        if not os.path.isdir(outdir):
            makepath(outdir)
        
        if os.path.isfile(tweakregfilename):
            if not overwrite:
                if skip_if_exists:
                    # return False means that rate2cal did not run
                    print(f'Image {tweakregfilename} already exists, skipping recreating it...')
                    return(False,tweakregfilename)
                else:
                    raise RuntimeError(f'Image {tweakregfilename} already exists! exiting. If you want to overwrite or skip rate2cal reduction, you can use "overwrite" or "skip_if_exists"')
            else:
                # make sure output filename is deleted
                rmfile(tweakregfilename)


        # It is important to set fitgeometry to rshift
        tweakreg.fitgeometry = 'rshift'
        tweakreg.align_to_gaia = False
        tweakreg.save_gaia_catalog = False
        tweakreg.save_results = True
        tweakreg.save_catalogs = False

        # the following parameters should have not impact, since these steps in tweakreg are skipped
        if 1==1:
            tweakreg.snr_threshold = 50
            tweakreg.separation = 9
            tweakreg.searchrad = 0.5
            tweakreg.minobj = 50
            tweakreg.min_gaia = 30
            tweakreg.xoffset = 0
            tweakreg.yoffset = 0
            tweakreg.brightest = 1000        
        

        tweakreg.already_matched = True
        # phot_cal.t.loc[ixs_cal_good] is the table with the good matches!
        if self.verbose: print(f'{len(ixs)} matches are passed to tweakreg {tweakreg.fitgeometry} fitting')
        t =  Table.from_pandas(phot.t.loc[ixs])
        tweakreg.refcat = t
        tweakreg.ref_racol = racol
        tweakreg.ref_deccol = deccol
        
        ### Provide your own source catalog, to be used in place of the default daofinder stuff. If you actually have a list
        ### of images, it's okay to provide a source catalog for each. 
        cal_image.source_catalog = t
        cal_data = [cal_image]
        tweakreg.source_xcol = xcol
        tweakreg.source_ycol = ycol
        
        if self.verbose: print(f'Fitting tweakreg fitgeometry={tweakreg.fitgeometry} to xy={xcol},{ycol} to ra,dec={racol},{deccol}')
        
        cal_data = [datamodels.open(cal_image)]
        tweakreg.run(cal_data)

        #make sure the image got created
        if not os.path.isfile(tweakregfilename):
            raise RuntimeError(f'Image {tweakregfilename} did not get created!!')

        # return True means that tweakrun did run
        return(True,tweakregfilename)
    
    """
    def apply_distortions(self,rate_image,
                distortion_file,
                overwrite = False,
                skip_if_exists = False,
                skip_rate2cal_if_exists = False,
                ):
        
        (runflag,calimname) = self.run_rate2cal(rate_image,
                    distortion_file,
                    overwrite = overwrite, 
                    skip_if_exists = (skip_rate2cal_if_exists |  skip_if_exists))
        print('####################',runflag,skip_if_exists)
        if not runflag: 
            if (skip_rate2cal_if_exists |  skip_if_exists):
                print(f'{calimname} already exists, skipping since skip_if_exists=True')
                return(0)
            else:
                raise RuntimeError(f'{calimname} already exists, stopping since skip_if_exists=False')
            
        return(0)
    """

    def run_all_old(self,rate_image,
                distortion_file,
                overwrite = False,
                skip_if_exists = False,
                skip_rate2cal_if_exists = False,
                skip_align2gaia_if_exists = False,
                gaia_catname_for_testing = './LMC_gaia_DR3.nrcposs',
                align_gaia_SNR_min = 10.0,
                searchrad = 3.0,
                xoffset = 0.0,
                yoffset = 0.0,
                ):
        
        print(f'###################### {skip_rate2cal_if_exists} {skip_if_exists}')
        (runflag,calimname) = self.run_rate2cal(rate_image,
                    distortion_file,
                    overwrite = overwrite, 
                    skip_if_exists = (skip_rate2cal_if_exists |  skip_if_exists))
        
        print(f'###################### XXXX1 {skip_rate2cal_if_exists} {skip_if_exists} {runflag}')
        if not runflag: 
            if (skip_rate2cal_if_exists |  skip_if_exists):
                print(f'{calimname} already exists, skipping since skip_if_exists=True')
            else:
                raise RuntimeError(f'{calimname} already exists, stopping since skip_if_exists=False')
                        
        if runflag:
            self.calphot.verbose = self.verbose
            self.calphot.run_phot(calimname,gaia_catname_for_testing,SNR_min=align_gaia_SNR_min)
        else:
            print('####### Skipping photometry on cal image!')
        
        print(f'###################### XXXX22 {skip_align2gaia_if_exists} {skip_if_exists} {runflag}')
        (runflag,tweakregfilename) = self.run_align2Gaia(calimname,
                    xoffset = xoffset,
                    yoffset = yoffset,
                    searchrad = searchrad,
                    overwrite = overwrite, 
                    skip_if_exists = (skip_align2gaia_if_exists |  skip_if_exists))
    
        if not runflag: 
            if (skip_align2gaia_if_exists |  skip_if_exists):
                print(f'{tweakregfilename} already exists, skipping since skip_if_exists=True')
            else:
                raise RuntimeError(f'{tweakregfilename} already exists, stopping since skip_if_exists=False')
    
        if runflag:
            self.gaialignphot.verbose = self.verbose
            self.gaialignphot.run_phot(tweakregfilename,gaia_catname_for_testing,SNR_min=align_gaia_SNR_min)
        else:
            print('####### Skipping photometry on tweakreg image!')
        return(0)

    def find_good_refcat_matches(self,
                                 phot=None,
                                 d2d_max = None,
                                 dmag_max = 0.1,
                                 Nbright=None,
                                 ixs=None,
                                 showplots=1,
                                 saveplots=1,
                                 savephottable=1,
                                 outbasename=None,
                                 plots_dxdy_delta_pix_ylim=7,
                                 # histo parameters
                                 binsize_px = 0.1, # this is the binsize of the dx/dy histograms
                                 bin_weights_flag=True,# If bin_weights_flag is set to True, 
                                                       #then the dx/dy bins are weighted by 
                                                       # the flux of the detection.
                                 slope_min=-10.0/2048.0, 
                                 slope_Nsteps = 200, # slope_max=-slope_min, slope_stepsize=(slope_max-slope_min)/slope_Nsteps
                                 Nfwhm = 2.0 
                                 ):
        if phot is None:
            phot=self.calphot
            
        if (saveplots or savephottable) and (outbasename is None):
            raise RuntimeError('Trying to save plots and/or phot tables, but outbasename is None!')

        Nx = phot.scihdr['NAXIS1']
        Ny = phot.scihdr['NAXIS2']
        
        # calculdate dx, dy, and do some first very rough cuts.
        ixs = calc_dxdy(phot,
                        phot.refcat.short,
                        d2d_max=d2d_max,
                        dmag_max=dmag_max,
                        Nbright=Nbright,
                        ixs=ixs)

        dx_median = phot.t.loc[ixs,'dx'].median()
        dy_median = phot.t.loc[ixs,'dy'].median()

        if self.verbose>1: print(f'Nx:{Nx} Ny:{Ny}\ndx median: {dx_median}\ndy median: {dy_median}')
        
        # these are the general limits for the y-axis for the dx/dy plots
        dy_plotlim = (dy_median-plots_dxdy_delta_pix_ylim,dy_median+plots_dxdy_delta_pix_ylim)
        dx_plotlim = (dx_median-plots_dxdy_delta_pix_ylim,dx_median+plots_dxdy_delta_pix_ylim)

        if showplots>1:
            sp = initplot(1,2)
            # plot the residuals
            phot.t.loc[ixs].plot.scatter('y','dx',ax=sp[0],ylim=dx_plotlim,color='red')
            phot.t.loc[ixs].plot.scatter('x','dy',ax=sp[1],ylim=dy_plotlim,color='red')
            plt.tight_layout()
            #if showplots: plt.show() 
            # add saveplots

        print('GGGG')        
        # Here we correct for rotation so that we can robustly remove outlier matches
        
        # A straight line with a given slope is subtracted from dx versus y, which is 
        # effectively a correction for a rotation.
        
        # Then a histogram of the resulting dx_rotated is calculated with bin size = binsize_px in pixels.
        # The maximum value of the histogram, it's associated dx_rot position, a rough upper limit estimate
        # of the fwhm of the histogram peak are then saved for each slope and intercept in dx_rot_results (pdastro table).
        
        # best_index is the index for this table for which the largest maxval of the histograms, presumably for 
        # the best rotation (or slope)
        
        # slope loops from slope_min to slope_max=-slope_min with slope_stepsize
        # slope_min=-10.0/2048.0 means that slopes with a range 
        # of +-10 pixels of 2048 pixels are probed, i.e. dx does not change
        # by more than 10 pixels over the full detector width
        #  slope_stepsize=(slope_max-slope_min)/slope_Nsteps

        slope_max=-slope_min
        slope_stepsize=(slope_max-slope_min)/slope_Nsteps

        ####BEST:
        #    slope  intercept  maxval  index  d_bestguess  fwhm  multimax
        #-0.000293        0.3     132   1919      1.91045   0.3     False
        #slope_min = -0.0005
        #slope_max = slope_min + 10*slope_stepsize
        
        (dx_rot_results,dx_best_index) = rotate_d_and_find_binmax(phot,ixs,'dx','y',
                                                                  Ny,
                                                                  binsize=binsize_px,
                                                                  bin_weights_flag=bin_weights_flag,
                                                                  slope_min=slope_min,
                                                                  slope_max=slope_max,
                                                                  slope_stepsize=slope_stepsize,
                                                                  showplots=showplots)
        print('GGGG1')        

              
        # Using the best dx_rotated, we first remove all entries with dx_rotated outside of dx_bestguess+-Nfwhm*fwhm
        # Note that FWHM ~ 2.355 stdev, so Nfwhm*fwhm should be at least 3*stdev. This is the first ROUGH cut, with 
        # which we just want to remove excessive amounts of outliers. Then a 3-sigma cut is done on the *rotated* dx
        (ixs_dx_cut,ixs_dx_roughcut) = sigmacut_d_rot(phot,ixs,'dx','y',
                                                      dx_rot_results.t.loc[dx_best_index,'slope'],
                                                      dx_rot_results.t.loc[dx_best_index,'intercept'],
                                                      dx_rot_results.t.loc[dx_best_index,'d_bestguess'],
                                                      rough_cut_px = Nfwhm*dx_rot_results.t.loc[dx_best_index,'fwhm'],
                                                      binsize=binsize_px,
                                                      bin_weights_flag=bin_weights_flag,
                                                      showplots=showplots
                                                      )
        print('GGGG2')        
        
        # Using the indices after the dx cut, we now also do a dy cut
        # slope in dy is the -slope of dx.
        slope = -dx_rot_results.t.loc[dx_best_index,'slope']
        # We only need to do it +-10*stepsize
        slope_min = slope-10*slope_stepsize
        slope_max = slope+10*slope_stepsize

        (dy_rot_results,dy_best_index) = rotate_d_and_find_binmax(phot,ixs_dx_cut,'dy','x',
                                                                  Nx,
                                                                  binsize=binsize_px,
                                                                  bin_weights_flag=bin_weights_flag,
                                                                  slope_min=slope_min,
                                                                  slope_max=slope_max,
                                                                  slope_stepsize=slope_stepsize,
                                                                  showplots=showplots)
        (ixs_dy_cut,ixs_dy_roughcut) = sigmacut_d_rot(phot,ixs_dx_cut,'dy','x',
                                                      dy_rot_results.t.loc[dy_best_index,'slope'],
                                                      dy_rot_results.t.loc[dy_best_index,'intercept'],
                                                      dy_rot_results.t.loc[dy_best_index,'d_bestguess'],
                                                      rough_cut_px = Nfwhm*dy_rot_results.t.loc[dy_best_index,'fwhm'],
                                                      binsize=binsize_px,
                                                      bin_weights_flag=bin_weights_flag,
                                                      showplots=showplots
                                                      )
        print('GGGG3')        
       

        if savephottable:
            print(f'Saving {outbasename}.good.phot.txt')
            phot.write(f'{outbasename}.good.phot.txt',indices=ixs_dy_cut)
            if savephottable>1:
                print(f'Saving {outbasename}.all.phot.txt')
                phot.write(f'{outbasename}.all.phot.txt')



        if showplots>1:
            # get the bad data points
            ixs_bad = AnotB(phot.getindices(),ixs_dy_cut)            
            infoplots(phot,ixs_dy_cut,ixs_bad,dy_plotlim=dy_plotlim,dx_plotlim=dx_plotlim)
            

        return(ixs_dy_cut)
    
    def dxdy_plot(self,phot,ixs_bestmatch, sp=None, spindex = None, title=None):
        if sp is None: sp=initplot(3,2)
        if spindex is None: spindex=range(6)
        
        dx_median = phot.t.loc[ixs_bestmatch,'dx'].median()
        dy_median = phot.t.loc[ixs_bestmatch,'dy'].median()

        dlim_big = 10
        dlim_small= 1
        dx_ylim_big = (dx_median-dlim_big,dx_median+dlim_big)
        dy_ylim_big = (dy_median-dlim_big,dy_median+dlim_big)
        dx_ylim_small = (dx_median-dlim_small,dx_median+dlim_small)
        dy_ylim_small = (dy_median-dlim_small,dy_median+dlim_small)

        phot.t.plot.scatter('y','dx',ylim=dx_ylim_big,ax=sp[spindex[0]],color='red')
        phot.t.loc[ixs_bestmatch].plot.scatter('y','dx',ylim=dx_ylim_big,ax=sp[spindex[0]],ylabel='dx in pixels',title=title)

        phot.t.plot.scatter('x','dy',ylim=dy_ylim_big,ax=sp[spindex[1]],color='red')
        phot.t.loc[ixs_bestmatch].plot.scatter('x','dy',ylim=dy_ylim_big,ax=sp[spindex[1]],ylabel='dy in pixels',title=title)

        phot.t.plot.scatter('y','dx',ylim=dx_ylim_small,ax=sp[spindex[2]],color='red')
        phot.t.loc[ixs_bestmatch].plot.scatter('y','dx',ylim=dx_ylim_small,ax=sp[spindex[2]],ylabel='dx in pixels',title=title)

        phot.t.plot.scatter('x','dx',ylim=dx_ylim_small,ax=sp[spindex[3]],color='red')
        phot.t.loc[ixs_bestmatch].plot.scatter('x','dx',ylim=dx_ylim_small,ax=sp[spindex[3]],ylabel='dx in pixels',title=title)

        phot.t.plot.scatter('x','dy',ylim=dy_ylim_small,ax=sp[spindex[4]],color='red')
        phot.t.loc[ixs_bestmatch].plot.scatter('x','dy',ylim=dy_ylim_small,ax=sp[spindex[4]],ylabel='dy in pixels',title=title)

        phot.t.plot.scatter('y','dy',ylim=dy_ylim_small,ax=sp[spindex[5]],color='red')
        phot.t.loc[ixs_bestmatch].plot.scatter('y','dy',ylim=dy_ylim_small,ax=sp[spindex[5]],ylabel='dy in pixels',title=title)

        plt.tight_layout()
        return(sp)
            
    def update_phottable_final_wcs(self,tweakregfilename,
                                   ixs_bestmatch,
                                   phot=None,
                                   showplots=1,
                                   saveplots=1,
                                   savephottable=1,
                                   overwrite=False):
        outbasename = re.sub('\.fits$','',tweakregfilename)
        if (outbasename == tweakregfilename): raise RuntimeError(f'Could not remove .fits from {tweakregfilename}')        
        if phot is None:
            phot=self.calphot

        # show or save dxdy pre WCS correction
        if showplots>=0 or saveplots:
            self.dxdy_plot(phot, ixs_bestmatch,title='pre WCS correction')
            if saveplots:
                outfilename = f'{outbasename}.phot.prewcs.png'
                if os.path.isfile(outfilename):
                    rmfile(outfilename)
                print(f'Saving {outfilename}')
                plt.savefig(outfilename)


        racol=f'{phot.refcat.short}_ra'
        deccol=f'{phot.refcat.short}_dec'
        xcol=f'{phot.refcat.short}_x'
        ycol=f'{phot.refcat.short}_y'

        phot.t.drop(columns=['dx','dy',xcol,ycol,f'{phot.refcat.short}_x_idl',f'{phot.refcat.short}_y_idl'],inplace=True)    
            
        image_model = datamodels.ImageModel(tweakregfilename)
        world_to_detector = image_model.meta.wcs.get_transform('world', 'detector')
        phot.t[xcol], phot.t[ycol] = world_to_detector(phot.t[racol],phot.t[deccol])
        
        phot.t['dx'] = phot.t[xcol] - phot.t['x']
        phot.t['dy'] = phot.t[ycol] - phot.t['y']

        if savephottable:
            outfilename = f'{outbasename}.good.phot.txt'
            if os.path.isfile(outfilename):
                if not overwrite:
                    raise RuntimeError(f'Image {outfilename} already exists! exiting. If you want to overwrite or skip rate2cal reduction, you can use "overwrite" or "skip_if_exists"')
                else:
                    # make sure output filename is deleted
                    rmfile(outfilename)
            print(f'Saving {outfilename}')
            phot.write(outfilename,indices=ixs_bestmatch)
            if savephottable>2:
                outfilename = f'{outbasename}.all.phot.txt'
                if os.path.isfile(outfilename):
                    if not overwrite:
                        raise RuntimeError(f'Image {outfilename} already exists! exiting. If you want to overwrite or skip rate2cal reduction, you can use "overwrite" or "skip_if_exists"')
                    else:
                        # make sure output filename is deleted
                        rmfile(outfilename)
                print(f'Saving {outfilename}')
                phot.write(outfilename)

        # show or save dxdy post WCS correction
        if showplots>=0 or saveplots:
            self.dxdy_plot(phot, ixs_bestmatch,title='after WCS correction')
            if saveplots:
                outfilename = f'{outbasename}.phot.finalwcs.png'
                if os.path.isfile(outfilename):
                    rmfile(outfilename)
                print(f'Saving {outfilename}')
                plt.savefig(outfilename)



    def run_all(self,input_image,
                distortion_file=None,
                overwrite = False,
                skip_if_exists = False,
                skip_applydistortions_if_exists = False,
                # refcat parameters
                refcatname = 'Gaia',
                refcat_racol='auto',
                refcat_deccol='auto',
                pmflag = False,
                pm_median=False,
                photfilename=None,
                load_photcat_if_exists=False,
                rematch_refcat=False,
                SNR_min = 10.0, # minimum S/N for photometry
                # find best matches to refcut
                d2d_max = None, # maximum distance refcat to source in image
                dmag_max = 0.1, # maximum uncertainty of source 
                Nbright=None,    # U/se only the brightest Nbright sources from image
                showplots=0,
                saveplots=0,
                savephottable=0
                ):
            
        # apply distortion coefficients if wanted.
        if distortion_file is not None:
            # apply distortion terms
            (runflag,calimname) = self.run_applydistortions(input_image,
                                                            distortion_file,
                                                            overwrite = overwrite, 
                                                            skip_if_exists = (skip_applydistortions_if_exists |  skip_if_exists))
        else:
            calimname = input_image
               
        # do the photometry
        self.calphot.verbose = self.verbose
        self.calphot.run_phot(calimname,
                              refcatname=refcatname,
                              refcat_racol=refcat_racol,
                              refcat_deccol=refcat_deccol,
                              pmflag=pmflag,
                              pm_median=pm_median,
                              photfilename=photfilename,
                              load_photcat_if_exists=load_photcat_if_exists,
                              rematch_refcat=rematch_refcat,
                              overwrite=overwrite,
                              SNR_min=SNR_min)
        
        matching_outbasename = re.sub('\.fits$','',calimname)
        if (matching_outbasename == calimname): raise RuntimeError(f'Could not remove .fits from {calimname}')        
        ixs_match= self.find_good_refcat_matches(d2d_max = d2d_max,
                                                 dmag_max = dmag_max,
                                                 Nbright = Nbright,
                                                 showplots=showplots,
                                                 saveplots=saveplots,
                                                 savephottable=savephottable,
                                                 outbasename=matching_outbasename
                                                 )        

        (runflag,tweakregfilename) = self.run_align2refcat(calimname,ixs=ixs_match,
                                                           overwrite=overwrite,skip_if_exists=skip_if_exists)
        
        self.update_phottable_final_wcs(tweakregfilename,
                                        ixs_bestmatch = ixs_match,
                                        showplots=showplots,
                                        saveplots=saveplots,
                                        savephottable=savephottable,
                                        overwrite=overwrite
                                        )

        if showplots: 
            plt.show()
        return(0)


        
if __name__ == '__main__':
    wcs_align = jwst_wcs_align()
    parser = wcs_align.define_options()
    args = parser.parse_args()
    
    wcs_align.verbose=args.verbose
    wcs_align.calphot=jwst_photclass()
    
    wcs_align.set_outdir(args.outrootdir, args.outsubdir)
    
    wcs_align.run_all(args.cal_image,
                     distortion_file = args.distortion_file,
                     overwrite = args.overwrite,
                     skip_applydistortions_if_exists=args.skip_applydistortions_if_exists,
                     refcatname = args.refcat,
                     refcat_racol = args.refcat_racol,
                     refcat_deccol = args.refcat_deccol,
                     pmflag = not args.refcat_no_pmflag,
                     pm_median = args.pm_median,
                     photfilename = args.photfilename,
                     load_photcat_if_exists=args.load_photcat_if_exists,
                     rematch_refcat=args.rematch_refcat,
                     SNR_min = args.SNR_min,
                     d2d_max = args.d2d_max, # maximum distance refcat to source in image
                     dmag_max = args.dmag_max, # maximum uncertainty of source 
                     Nbright=args.Nbright,    # U/se only the brightest Nbright sources from image
                     showplots=args.showplots,
                     saveplots=args.saveplots,# 
                     savephottable=args.savephottable
                     )
    
