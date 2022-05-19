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
        
        self.debug = False
        
        self.test_distortion_singleim = test_distortion_singleim()
        
    def define_options(self,parser=None,usage=None,conflict_handler='resolve'):
        parser = apply_distortions.define_options(self,parser=parser,usage=usage,conflict_handler=conflict_handler)
        parser = self.test_distortion_singleim.default_options(parser)
        #parser.add_argument('-d','--debug', default=False, action='store_true', help='debug mode: throw exceptions instead of try except block')

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
                                                      yoffset = yoffset)
                self.t.loc[ix,'errorflag']=False
            except Exception as e:
                print(f'ERROR while running {ratefile} {distfile}: {e}')
                self.t.loc[ix,'errorflag']=True

if __name__ == '__main__':

    testdist = test_distortions()
    parser = testdist.define_options()
    args = parser.parse_args()
    
    testdist.verbose=args.verbose
    #testdist.debug=args.debug
    
    testdist.get_rate_files(args.rate_files,directory=args.rate_dir)
    testdist.get_distortion_files(args.distortion_files,directory=None)
    
    ixs_matches,ixs_not_matches = testdist.match_distortion4ratefile(apertures=args.apertures, 
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