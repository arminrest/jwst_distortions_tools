#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 16 10:50:49 2022

@author: arest
"""

from jwst.pipeline.calwebb_image2 import Image2Pipeline
from jwst.assign_wcs import AssignWcsStep

import argparse,re,os
from pdastro import makepath,rmfile

class apply_distortion_singleim:
    def __init__(self):
        self.verbose=0
        self.outdir = None
                
    def define_options(self,parser=None,usage=None,conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)

        parser.add_argument('cal_image',  help='rate fits file.')
        parser.add_argument('distortion_file',  help='distortion file, in asdf format')

        parser = self.default_options(parser)
        return(parser)

    def default_options(self,parser):

        # default directory for input images
        #if 'JWST_DISTORTION_IMAGEDIR' in os.environ:
        #    ratedir = os.environ['JWST_DISTORTION_IMAGEDIR']
        #else:
        #    ratedir = None

        # default directory for output
        if 'JWST_OUTROOTDIR' in os.environ:
            outrootdir = os.environ['JWST_OUTROOTDIR']
        else:
            outrootdir = None

        #parser.add_argument('--rate_dir', default=ratedir, help='Directory in which the rate images are located, which will be used to test the distortions. (default=%(default)s)')

        parser.add_argument('--outrootdir', default=outrootdir, help='output root directory. The output directoy is the output root directory + the outsubdir if not None (default=%(default)s)')
        parser.add_argument('--outsubdir', default=None, help='outsubdir added to output root directory (default=%(default)s)')
        parser.add_argument('--overwrite', default=False, action='store_true', help='overwrite files if they exist.')

        parser.add_argument('--skip_if_exists', default=False, action='store_true', help='Skip doing the analysis of a given input image if the cal file already exists, assuming the full analysis has been already done')

        parser.add_argument('-v','--verbose', default=0, action='count')

        return(parser)
    
    def set_outdir(self,outrootdir=None,outsubdir=None):
        self.outdir = outrootdir
        if self.outdir is None: self.outdir = '.'
        
        if outsubdir is not None:
            self.outdir+=f'/{outsubdir}'
        
        return(self.outdir)

    def run_applydistortions_assignwcs(self,cal_image,distortion_file,outdir=None,
                     overwrite=False, skip_if_exists=False):
        """
        Calwebb_image2 - does flat fielding, attaches WCS from distortion reffile, flux calibration
                         Outputs *cal.fits files
        checking format and if input files exist

        Parameters
        ----------
        cal_image : TYPE
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

        print(f'assigning WCS to file {cal_image} using distortion file {distortion_file}')
        step = AssignWcsStep()

        if not os.path.isfile(cal_image):
            raise RuntimeError(f'image {cal_image} does not exist')

        if distortion_file.lower() == 'none':
            print('WARNING!! not applying any distortion file!!')
        else:
            if re.search('\.asdf$',distortion_file) is None:
                raise RuntimeError(f'distortion file {distortion_file} does not have .asdf suffix. asdf format required.')
            if not os.path.isfile(distortion_file):
                raise RuntimeError(f'distortion file {distortion_file} does not exist')
            step.override_distortion = distortion_file
            
        if outdir is None:
            outdir=self.outdir
            
        step.save_results = True
 
        #if outdir is None:
        #    outdir=os.path.dirname(cal_image)
        assignwcsfilename = re.sub('\_[a-zA-Z0-9]+\.fits$','_assignwcsstep.fits',os.path.basename(cal_image))
        if assignwcsfilename == os.path.basename(cal_image):
            raise RuntimeError('Could not get assignwcsstep filename from {os.path.basename(cal_image)}!')
        assignwcsfilename = f'{outdir}/{assignwcsfilename}'
        
        if self.verbose: print(f'Setting output directory for assignwcsstep.fits file to {outdir}')
        step.output_dir = outdir
        if not os.path.isdir(outdir):
            makepath(outdir)

        if os.path.isfile(assignwcsfilename):
            if skip_if_exists:
                # return False means that rate2cal did not run
                print(f'Image {assignwcsfilename} already exists, skipping recreating it...')
                return(False,assignwcsfilename)
            else:
                if overwrite:
                    print(f'WARNING! {assignwcsfilename} exists, deleting it since "overwrite" is set!')
                    # make sure cal frame is deleted
                    rmfile(assignwcsfilename)
                else:
                    raise RuntimeError(f'Image {assignwcsfilename} already exists! exiting. If you want to overwrite or skip, you can use "overwrite" or "skip_if_exists"')
                
        print(f'Creating {assignwcsfilename}')
        step.run(cal_image)
        
        #make sure the image got created
        if not os.path.isfile(assignwcsfilename):
            raise RuntimeError(f'Image {assignwcsfilename} did not get created!!')
        else:
            print(f'distortions applied to {assignwcsfilename}!!')
            
        # return True means that rate2cal did run
        return(True,assignwcsfilename)

    def run_applydistortions_rate2cal(self,rate_image,distortion_file,outdir=None,
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
                print(f'WARNING! {calfilename} exists, deleting it since "overwrite" is set!')
                # make sure cal frame is deleted
                rmfile(calfilename)
                
        print(f'Creating {calfilename}')
        stage2_img.run(rate_image)
        
        #make sure the image got created
        if not os.path.isfile(calfilename):
            raise RuntimeError(f'Image {calfilename} did not get created!!')
            
        # return True means that rate2cal did run
        return(True,calfilename)

    def run_applydistortions(self,input_image,distortion_file,outdir=None,
                             overwrite=False, skip_if_exists=False):
        if re.search('rate\.fits$',input_image) is not None:
            # if it is a rate file, do the full level2b steps
            if self.verbose: print('This is a rate.fits file! doing the fill level2b part, and applying new distortions!')
            (runflag,calimname) = self.run_applydistortions_rate2cal(input_image,
                                                                     distortion_file,
                                                                     outdir=outdir,
                                                                     overwrite = overwrite, 
                                                                     skip_if_exists = skip_if_exists)
            
        else:
            if self.verbose: print('This is NOT a rate.fits file! Therefore just applying the distortions with AssignWcsStep')
            (runflag,calimname) = self.run_applydistortions_assignwcs(input_image,
                                                                      distortion_file,
                                                                      outdir=outdir,
                                                                      overwrite = overwrite, 
                                                                      skip_if_exists = skip_if_exists)
        return(runflag,calimname)

if __name__ == '__main__':

    applydist = apply_distortion_singleim()
    parser = applydist.define_options()
    args = parser.parse_args()
    
    applydist.verbose=args.verbose
    
    applydist.set_outdir(args.outrootdir, args.outsubdir)

    applydist.run_applydistortions(args.cal_image,
                                   args.distortion_file,
                                   overwrite = args.overwrite,
                                   skip_if_exists = args.skip_if_exists)
