#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 09:39:07 2022

@author: arest
"""

import argparse,sys,os
from test_distortions_single_image import test_distortion_singleim
from apply_distortions import apply_distortions

class test_distortions(apply_distortions):
    def __init__(self):
        apply_distortions.__init__(self)
        
        self.test_distortion_singleim = test_distortion_singleim()
        
    def define_options(self,parser=None,usage=None,conflict_handler='resolve'):
        parser = apply_distortions.define_options(self,parser=parser,usage=usage,conflict_handler=conflict_handler)
        parser = self.test_distortion_singleim.default_options(parser)
 
        return(parser)

        ###### DELETE
        
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)

        # default for token is $MAST_API_TOKEN
        if 'JWST_DISTORTION_IMAGEDIR' in os.environ:
            ratedir = os.environ['JWST_DISTORTION_IMAGEDIR']
        else:
            ratedir = None


        parser.add_argument('--rate_dir', default=ratedir, help='Directory in which the rate images are located, which will be used to test the distortions. (default=%(default)s)')
        parser.add_argument('--rate_files', nargs='+', default=['*_rate.fits'], help='list of rate file(pattern)s to which the distortion files are applied to. "rate_dir" is used if not None (default=%(default)s)')
        parser.add_argument('--distortion_files', nargs='+', default=['*.asdf'], help='list of the distortion file(pattern)s to be applied. (default=%(default)s)')

        parser.add_argument('--outdir', default='same_as_distortionfiles', help='output directory. If "same_as_distortionfiles", then the cal/photometry images are saved in the same directory as the distortion coefficient files (default=%(default)s)')
        #parser.add_argument('--add2basename', default=None, help='This is added to the basename. (default=%(default)s)')
        parser.add_argument('--overwrite', default=False, action='store_true', help='overwrite files if they exist.')
        parser.add_argument('--skip_rate2cal_if_exists', default=False, action='store_true', help='If the output cal file already exists, skip running the level 2 pipeline to assign the new distortion terms, assuming this has already been done, but still do the photometry.')

        parser.add_argument('--ignore_filters', default=False, action='store_true', help='distortions are grouped by aperture/filter/pupil. Use this option if you want to create distortion files independent of filter.')
        parser.add_argument('--ignore_pupils', default=False, action='store_true', help='distortions are grouped by aperture/filter/pupil. Use this option if you want to create distortion files independent of pupil.')

        parser.add_argument('--apertures', nargs='+', default=None, help='constrain the rate file list to these apertures (default=%(default)s)')
        parser.add_argument('--filters', nargs='+', default=None, help='constrain the rate file list to these filters (default=%(default)s)')
        parser.add_argument('--pupils', nargs='+', default=None, help='constrain the rate file list to these pupils (default=%(default)s)')

        parser.add_argument('-v','--verbose', default=0, action='count')

        return(parser)
    
    def run_images(self, ixs, 
                   outdir=None,
                   outsubdir=None,
                   overwrite = False,
                   skip_if_exists = False,
                   skip_rate2cal_if_exists = False,
                   skip_align2gaia_if_exists = False,
                   gaia_catname_for_testing = './LMC_gaia_DR3.nrcposs',
                   align_gaia_SNR_min = 10.0,
                   searchrad = 3.0,
                   xoffset = 0.0,
                   yoffset = 0.0):
        
        self.t.loc[ixs,'errorflag'] = None
        for ix in ixs:
            distfile = self.t.loc[ix,'distortion_match']
            ratefile = self.t.loc[ix,'filename']
            testdist_singleim = test_distortion_singleim()
            testdist_singleim.verbose = self.verbose
            
            if outdir == 'same_as_distortionfiles':
                outdir = os.path.dirname(distfile)
                
            self.test_distortion_singleim.set_outdir(outdir,outsubdir)

            try:
                self.test_distortion_singleim.run_all(ratefile,distfile,                     
                                                      overwrite = overwrite,
                                                      skip_if_exists = skip_if_exists,
                                                      skip_rate2cal_if_exists = skip_rate2cal_if_exists,
                                                      skip_align2gaia_if_exists = skip_align2gaia_if_exists,
                                                      gaia_catname_for_testing = gaia_catname_for_testing,
                                                      align_gaia_SNR_min = align_gaia_SNR_min,
                                                      searchrad  = searchrad,
                                                      xoffset = xoffset,
                                                      yoffset = yoffset
                                          )
                self.t.loc[ix,'errorflag']=False
            except:
                self.t.loc[ix,'errorflag']=True

if __name__ == '__main__':

    testdist = test_distortions()
    parser = testdist.define_options()
    args = parser.parse_args()
    
    testdist.verbose=args.verbose
    
    testdist.get_rate_files(args.rate_files,directory=args.rate_dir)
    testdist.get_distortion_files(args.distortion_files,directory=None)
    
    ixs_matches,ixs_not_matches = testdist.match_distortion4ratefile(require_filter=not args.ignore_filters, 
                                                                     require_pupil=not args.ignore_pupils,
                                                                     apertures=args.apertures, 
                                                                     filts=args.filters, 
                                                                     pupils=args.pupils)
    
    if len(ixs_matches)==0:
        print('NO IMES FOUND!! exiting...')
        sys.exit(0)

    testdist.run_images(ixs_matches,
                        outdir=args.outrootdir,
                        outsubdir=args.outsubdir,
                        overwrite = args.overwrite,
                        skip_if_exists = args.skip_if_exists,
                        skip_rate2cal_if_exists = args.skip_rate2cal_if_exists,
                        skip_align2gaia_if_exists = args.skip_align2gaia_if_exists,
                        gaia_catname_for_testing = args.gaia_catname_for_testing,
                        align_gaia_SNR_min = args.align_gaia_SNR_min,
                        searchrad  = args.searchrad,
                        xoffset = args.xoffset,
                        yoffset = args.yoffset)
    testdist.write()