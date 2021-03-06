#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 14:32:42 2022

@author: arest, bhilbert, mcorrenti, acanipe
"""

import os,re
from pdastro import makepath,rmfile
from simple_jwst_phot import jwst_photclass
#from jwst.tweakreg import TweakRegStep
import tweakreg_hack

from jwst import datamodels
from apply_distortions_single_image import apply_distortion_singleim


class test_distortion_singleim(apply_distortion_singleim):
    def __init__(self):
        apply_distortion_singleim.__init__(self)
        self.calphot=jwst_photclass()
        self.gaialignphot=jwst_photclass()
        
    def default_options(self,parser):
        parser = apply_distortion_singleim.default_options(self,parser=parser)

        parser.add_argument('--gaia_catname_for_testing', default='./LMC_gaia_DR3.nrcposs', help='Gaia catalog used for TESTING, not for tweakreg! (default=%(default)s)')

        parser.add_argument('--align_gaia_SNR_min', default=10.0, help='when aligning with Gaia: mininum SNR above noise to trigger object in image (default=%(default)s)')

        parser.add_argument('--skip_rate2cal_if_exists', default=False, action='store_true', help='If the output cal file already exists, skip running the level 2 pipeline to assign the new distortion terms, assuming this has already been done, but still do the photometry.')
        parser.add_argument('--skip_align2gaia_if_exists', default=False, action='store_true', help='If the output cal file already exists, skip running the level 2 pipeline to assign the new distortion terms, assuming this has already been done, but still do the photometry.')

        parser.add_argument('--xoffset', type=float, default=0.0, help='Initial guess for X offset in arcsec for align_gaia step. (default=%(default)s)')
        parser.add_argument('--yoffset', type=float, default=0.0, help='Initial guess for Y offset in arcsec for align_gaia step. (default=%(default)s)')
        parser.add_argument('--searchrad', type=float, default=3.0, help='The search radius in arcsec for a match for align_gaia step. Note: currently, this valiue is divided by 10 when the actual search is done! (default=%(default)s)')

        return(parser)
        
        
    def run_align2Gaia(self,cal_image,tweakreg=None,
                       kernel_fwhm=None, #float(default=2.5) # Gaussian kernel FWHM in pixels
                       snr_threshold = 50.0, #float(default=10.0) # SNR threshold above the bkg
                       separation=9, # minimum separation between reference stars in arcsec
                       searchrad=0.3, # (default=1.0): The search radius in arcsec for a match
                       minobj = 50, #integer(default=15) # Minimum number of objects acceptable for matching
                       min_gaia = 30, #integer(min=0, default=5) # Min number of GAIA sources needed
                       xoffset = 0, # Initial guess for X offset in arcsec. (Default=0.0)
                       yoffset = 0, # Initial guess for Y offset in arcsec. (Default=0.0)
                       brightest = 1000,
                       outdir=None,overwrite=False, skip_if_exists=False):
        if tweakreg is None:
            #tweakreg = TweakRegStep()
            tweakreg = tweakreg_hack.TweakRegStep()

        tweakreg.align_to_gaia = True
        tweakreg.save_results = True
        tweakreg.save_catalogs = True
        tweakreg.snr_threshold = snr_threshold
        tweakreg.searchrad = searchrad
        tweakreg.separation = separation
        tweakreg.minobj = minobj
        tweakreg.min_gaia = min_gaia
        tweakreg.brightest= brightest
        tweakreg.xoffset = xoffset
        tweakreg.yoffset = yoffset

        # only correct for rotation and shift
        tweakreg.fitgeometry = 'rshift'
        
        if outdir is None:
            outdir=self.outdir
        #if outdir is None:
        #    outdir=os.path.dirname(cal_image)
        tweakregfilename = f'{outdir}/'+re.sub('cal\.fits$','tweakregstep.fits',os.path.basename(cal_image))
        if self.verbose: print(f'Setting output directory for tweakregstep.fits file to {outdir}')
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
                # make sure cal frame is deleted
                rmfile(tweakregfilename)

        cal_data = [datamodels.open(cal_image)]
        tweakreg.run(cal_data)

        #make sure the image got created
        if not os.path.isfile(tweakregfilename):
            raise RuntimeError(f'Image {tweakregfilename} did not get created!!')

        # return True means that rate2cal did run
        return(True,tweakregfilename)

    def run_all(self,rate_image,
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
        (runflag,calimname) = self.run_applydistortions_rate2cal(rate_image,
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
        
if __name__ == '__main__':
    testdist = test_distortion_singleim()
    parser = testdist.define_options()
    args = parser.parse_args()
    
    testdist.verbose=args.verbose
    testdist.calphot=jwst_photclass()
    testdist.gaialignphot=jwst_photclass()
    
    testdist.set_outdir(args.outrootdir, args.outsubdir)
    
    testdist.run_all(args.rate_image,
                     args.distortion_file,
                     overwrite = args.overwrite,
                     skip_if_exists = args.skip_if_exists,
                     skip_rate2cal_if_exists = args.skip_rate2cal_if_exists,
                     skip_align2gaia_if_exists = args.skip_align2gaia_if_exists,
                     gaia_catname_for_testing = args.gaia_catname_for_testing,
                     align_gaia_SNR_min = args.align_gaia_SNR_min,
                     searchrad  = args.searchrad,
                     xoffset = args.xoffset,
                     yoffset = args.yoffset)
    
    """
    sys.exit(0)
    
    (runflag,calimname) = testdist.run_rate2cal(args.rate_image,
                args.distortion_file,
                overwrite = args.overwrite, 
                skip_if_exists = (args.skip_rate2cal_if_exists |  args.skip_if_exists))
    
    if not runflag and args.skip_if_exists:
        print(f'{calimname} already exists, stopping since skip_if_exists=True')
        sys.exit(0)
        
    if runflag:
        testdist.calphot.verbose = testdist.verbose
        testdist.calphot.run_phot(calimname,args.gaia_catname_for_testing,SNR_min=args.align_gaia_SNR_min)
    else:
        print('####### Skipping photometry on cal image!')
    
    (runflag,tweakregfilename) = testdist.run_align2Gaia(calimname,
                xoffset = args.xoffset,
                yoffset = args.yoffset,
                searchrad = args.searchrad,
                overwrite = args.overwrite, 
                skip_if_exists = (args.skip_align2gaia_if_exists |  args.skip_if_exists))

    if not runflag and args.skip_if_exists:
        print(f'{tweakregfilename} already exists, stopping since skip_if_exists=True')
        sys.exit(0)

    if runflag:
        testdist.gaialignphot.verbose = testdist.verbose
        testdist.gaialignphot.run_phot(tweakregfilename,args.gaia_catname_for_testing,SNR_min=args.align_gaia_SNR_min)
    else:
        print('####### Skipping photometry on tweakreg image!')
    """