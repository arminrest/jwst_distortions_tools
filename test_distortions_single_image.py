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

class test_distortion_singleim(jwst_photclass):
    def __init__(self):
        jwst_photclass.__init__(self)

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

        parser.add_argument('--outrootdir', default=outrootdir, help='Directory in which the cal images are located, which will be used to test the distortions. (default=%(default)s)')
        parser.add_argument('--outsubdir', default=None, help='outsubdir added to output root directory (default=%(default)s)')
        parser.add_argument('--overwrite', default=False, action='store_true', help='overwrite files if they exist.')

        parser.add_argument('--skip_rate2cal_if_exists', default=False, action='store_true', help='If the output cal file already exists, skip running the level 2 pipeline to assign the new distortion terms, assuming this has already been done, but still do the photometry.')
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
        if re.search('_rate\.fits$',rate_image) is None:
            raise RuntimeError(f'rate image {rate_image} does not have _rate.fits suffix')
        if re.search('\.asdf$',distortion_file) is None:
            raise RuntimeError(f'distortion file {distortion_file} does not have .asdf suffix. asdf format required.')
        if not os.path.isfile(rate_image):
            raise RuntimeError(f'rate image {rate_image} does not exist')
        if not os.path.isfile(distortion_file):
            raise RuntimeError(f'distortion file {distortion_file} does not exist')
            
        if outdir is None:
            outdir=self.outdir
            
        print(f'reducing file {rate_image} with distortion file {distortion_file}')
        stage2_img = Image2Pipeline()
        stage2_img.save_results = True
        stage2_img.override_distortion = distortion_file
        
        calfilename = re.sub('rate\.fits$','cal.fits',rate_image)
        if outdir is not None:
            calfilename = f'{outdir}/{os.path.basename(calfilename)}'
            if self.verbose: print(f'Setting output directory for cal.fits file to {outdir}')
            stage2_img.output_dir = outdir
            if not os.path.isdir(outdir):
                makepath(outdir)

        print(calfilename)
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
        
        

if __name__ == '__main__':
    testdist = test_distortion_singleim()
    parser = testdist.define_options()
    args = parser.parse_args()
    
    testdist.verbose=args.verbose
    
    testdist.set_outdir(args.outrootdir, args.outsubdir)
    
    (runflag,calimname) = testdist.run_rate2cal(args.rate_image,
                args.distortion_file,
                overwrite = args.overwrite, 
                skip_exists = (args.skip_rate2cal_if_exists |  args.skip_if_exists))
    
    if not runflag and args.skip_if_exists:
        print(f'{calimname} already exists, stopping since skip_if_exists=True')
        sys.exit(0)
    