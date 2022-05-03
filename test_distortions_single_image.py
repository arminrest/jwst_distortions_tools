#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 14:32:42 2022

@author: arest, bhilbert, mcorrenti, acanipe
"""

from jwst.pipeline.calwebb_image2 import Image2Pipeline
import glob,argparse,os,re,sys
from pdastro import pdastroclass,makepath4file,makepath,unique,AnotB,AorB,AandB,rmfile
from simple_jwst_phot import jwst_photclass
from jwst.tweakreg import TweakRegStep
from jwst import datamodels

class test_distortion_singleim:
    def __init__(self):
        self.calphot=jwst_photclass()
        self.gaialignphot=jwst_photclass()

    def define_options(self,parser=None,usage=None,conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)

        # default for token is $MAST_API_TOKEN
        if 'JWST_DISTORTION_OUTROOTDIR' in os.environ:
            outrootdir = os.environ['JWST_DISTORTION_OUTROOTDIR']
        else:
            outrootdir = None


        parser.add_argument('rate_image',  help='rate fits file')
        parser.add_argument('distortion_file',  help='distortion file, in asdf format')

        parser.add_argument('--gaia_catname_for_testing', default='./LMC_gaia_DR3.nrcposs', help='Gaia catalog used for TESTING, not for tweakreg! (default=%(default)s)')

        parser.add_argument('--align_gaia_SNR_min', default=10.0, help='when aligning with Gaia: mininum SNR above noise to trigger object in image (default=%(default)s)')

        parser.add_argument('--outrootdir', default=outrootdir, help='Directory in which the cal images are located, which will be used to test the distortions. (default=%(default)s)')
        parser.add_argument('--outsubdir', default=None, help='outsubdir added to output root directory (default=%(default)s)')
        parser.add_argument('--overwrite', default=False, action='store_true', help='overwrite files if they exist.')

        parser.add_argument('--skip_rate2cal_if_exists', default=False, action='store_true', help='If the output cal file already exists, skip running the level 2 pipeline to assign the new distortion terms, assuming this has already been done, but still do the photometry.')
        parser.add_argument('--skip_align2gaia_if_exists', default=False, action='store_true', help='If the output cal file already exists, skip running the level 2 pipeline to assign the new distortion terms, assuming this has already been done, but still do the photometry.')
        parser.add_argument('--skip_if_exists', default=False, action='store_true', help='Skip doing the analysis of a given input image if the cal file already exists, assuming the full analysis has been already done')

        parser.add_argument('-v','--verbose', default=0, action='count')

        return(parser)

    def set_outdir(self,outrootdir,outsubdir):
        self.outdir = outrootdir
        if self.outdir is None: self.outdir = '.'
        
        if outsubdir is not None:
            self.outdir+=f'/{outsubdir}'
        
        return(self.outdir)


    def run_rate2cal(self,rate_image,distortion_file,outdir=None,
                     overwrite=False, skip_if_exists=False):
        """
        Calwebb_image2 - does flat fielding, attaches WCS from distortion reffile, flux calibration
                         Outputs *cal.fits files
        checking format and if input files exist

        Parameters
        ----------
        rate_image : TYPE
            DESCRIPTION.
        distortion_file : TYPE
            DESCRIPTION.
        outdir : TYPE, optional
            DESCRIPTION. The default is None.
        overwrite : TYPE, optional
            DESCRIPTION. The default is False.
        skip_if_exists : TYPE, optional
            DESCRIPTION. The default is False.

        Raises
        ------
        RuntimeError
            DESCRIPTION.

        Returns
        -------
        (True/False,calimagename)
        If True, Calwebb_image2 has been run and the cal image has been created
        If False, then calimagename already existed, and re-creation is skipped since skip_if_exists=True

        """

        print(f'reducing file {rate_image} with distortion file {distortion_file}')
        stage2_img = Image2Pipeline()

        if re.search('_rate\.fits$',rate_image) is None:
            raise RuntimeError(f'rate image {rate_image} does not have _rate.fits suffix')
        if not os.path.isfile(rate_image):
            raise RuntimeError(f'rate image {rate_image} does not exist')

        if distortion_file.lower() == 'none':
            print('WARNING!! not applying any distortion file!!',stage2_img.assign_wcs.override_distortion)
        else:
            if re.search('\.asdf$',distortion_file) is None:
                raise RuntimeError(f'distortion file {distortion_file} does not have .asdf suffix. asdf format required.')
            if not os.path.isfile(distortion_file):
                raise RuntimeError(f'distortion file {distortion_file} does not exist')
            stage2_img.assign_wcs.override_distortion = distortion_file
            
        if outdir is None:
            outdir=self.outdir
            
        stage2_img.save_results = True
 
        #if outdir is None:
        #    outdir=os.path.dirname(rate_image)
        calfilename = f'{outdir}/'+re.sub('rate\.fits$','cal.fits',os.path.basename(rate_image))
        if self.verbose: print(f'Setting output directory for cal.fits file to {outdir}')
        stage2_img.output_dir = outdir
        if not os.path.isdir(outdir):
            makepath(outdir)

        if os.path.isfile(calfilename):
            if not overwrite:
                if skip_if_exists:
                    # return False means that rate2cal did not run
                    print(f'Image {calfilename} already exists, skipping recreating it...')
                    return(False,calfilename)
                else:
                    raise RuntimeError(f'Image {calfilename} already exists! exiting. If you want to overwrite or skip rate2cal reduction, you can use "overwrite" or "skip_if_exists"')
            else:
                # make sure cal frame is deleted
                rmfile(calfilename)
                
        print(f'Creating {calfilename}')
        stage2_img.run(rate_image)
        
        #make sure the image got created
        if not os.path.isfile(calfilename):
            raise RuntimeError(f'Image {calfilename} did not get created!!')
            
        # return True means that rate2cal did run
        return(True,calfilename)
        
    def run_align2Gaia(self,cal_image,tweakreg=None,
                       kernel_fwhm=None, #float(default=2.5) # Gaussian kernel FWHM in pixels
                       snr_threshold = 20.0, #float(default=10.0) # SNR threshold above the bkg
                       searchrad=3.0, # (default=1.0): The search radius in arcsec for a match
                       minobj = 50, #integer(default=15) # Minimum number of objects acceptable for matching
                       min_gaia = 30, #integer(min=0, default=5) # Min number of GAIA sources needed
                       outdir=None,overwrite=False, skip_if_exists=False):
        
        if tweakreg is None:
            tweakreg = TweakRegStep()
        tweakreg.align_to_gaia = True
        tweakreg.save_results = True

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
        
if __name__ == '__main__':
    testdist = test_distortion_singleim()
    parser = testdist.define_options()
    args = parser.parse_args()
    
    testdist.verbose=args.verbose
    
    testdist.set_outdir(args.outrootdir, args.outsubdir)
    
    (runflag,calimname) = testdist.run_rate2cal(args.rate_image,
                args.distortion_file,
                overwrite = args.overwrite, 
                skip_if_exists = (args.skip_rate2cal_if_exists |  args.skip_if_exists))
    
    if not runflag and args.skip_if_exists:
        print(f'{calimname} already exists, stopping since skip_if_exists=True')
        sys.exit(0)
        
    testdist.calphot.verbose = testdist.verbose
    testdist.calphot.run_phot(calimname,args.gaia_catname_for_testing,SNR_min=args.align_gaia_SNR_min)
    
    (runflag,tweakregfilename) = testdist.run_align2Gaia(calimname,
                overwrite = args.overwrite, 
                skip_if_exists = (args.skip_align2gaia_if_exists |  args.skip_if_exists))

    if not runflag and args.skip_if_exists:
        print(f'{tweakregfilename} already exists, stopping since skip_if_exists=True')
        sys.exit(0)

    testdist.gaialignphot.verbose = testdist.verbose
    testdist.gaialignphot.run_phot(tweakregfilename,args.gaia_catname_for_testing,SNR_min=args.align_gaia_SNR_min)
