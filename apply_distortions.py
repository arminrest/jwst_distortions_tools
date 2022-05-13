#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 25 09:39:07 2022

@author: arest
"""

import argparse,glob,re,sys,os
from pdastro import pdastroclass,makepath4file,unique,AnotB,AorB,AandB
import pandas as pd
from astropy.io import fits
import numpy as np
from test_distortions_single_image import test_distortion_singleim

class apply_distortions(pdastroclass):
    def __init__(self):
        pdastroclass.__init__(self)
        
        self.verbose=0
        
        self.distortionfiles = pdastroclass()
        
        self.aperture_col = 'AperName'
        self.filter_col = 'filter'
        self.pupil_col = 'pupil'
        

    def define_options(self,parser=None,usage=None,conflict_handler='resolve'):
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
    
    def get_files(self,filepatterns,directory=None):
        filenames=[]
        for filepattern in filepatterns:
            #(tmpdir,basename) = os.path.split(filepattern)
            #print(f'{tmpdir},{basename}')
            #if tmpdir =='' and (directory is not None):
            if directory is not None:
                filepattern=os.path.join(directory,filepattern)
            if self.verbose>2: print(f'Looking for filepattern {filepattern}')
            filenames.extend(glob.glob(filepattern))
        
        for i in range(len(filenames)):
            filenames[i] = os.path.abspath(filenames[i])
        if self.verbose: print(f'Found {len(filenames)} files matching filepatterns {filepatterns}')
        return(filenames)
    
    def get_rate_files(self,filepatterns,directory=None):
        self.t['filename'] = self.get_files(filepatterns,directory=directory)
        for ix in self.getindices():
            hdr = fits.getheader(self.t.loc[ix,'filename'])
            detector = re.sub('long$','5',hdr['DETECTOR'].lower())
            self.t.loc[ix,self.aperture_col]=f'{detector}_{hdr["SUBARRAY"].lower()}'
            self.t.loc[ix,'detector']=f'{hdr["DETECTOR"].lower()}'
            self.t.loc[ix,'subarray']=f'{hdr["SUBARRAY"].lower()}'
            self.t.loc[ix,self.filter_col]=f'{hdr["FILTER"].lower()}'
            self.t.loc[ix,self.pupil_col]=f'{hdr["PUPIL"].lower()}'
        if self.verbose:
            print('##################\n### Rate files:')
            self.write()
        
    def get_distortion_files(self,filepatterns,directory=None):
        self.distortionfiles.t['filename'] = self.get_files(filepatterns,directory=directory)
        for ix in self.distortionfiles.getindices():
            
            # get the filter and save it in the 'filter' column
            #m = re.search('^([a-zA-Z0-9]+_[a-zA-Z0-9]+)_(f[a-zA-Z0-9]+)_([a-zA-Z0-9]+)',os.path.basename(self.distortionfiles.t.loc[ix,'filename']))
            #if m is None:
            #    raise RuntimeError(f'could not parse filename {os.path.basename(self.distortionfiles.t.loc[ix,"filename"])} for aperture, filter and/or pupil!')
            #aperture,filt,pupil = m.groups()

            m1 = re.search('distortion_coeffs_([a-zA-Z0-9]+_[a-zA-Z0-9]+)_([a-zA-Z0-9]+)_([a-zA-Z0-9]+).*\.asdf',os.path.basename(self.distortionfiles.t.loc[ix,'filename']))
            m2 = re.search('^([a-zA-Z0-9]+_[a-zA-Z0-9]+).*_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)\..*distcoeff\.asdf',os.path.basename(self.distortionfiles.t.loc[ix,'filename']))
            if m1 is not None:
                aperture,filt,pupil = m1.groups()
            elif m2 is not None:
                aperture,filt,pupil = m2.groups()        
            else:
                raise RuntimeError(f'could not parse filename {os.path.basename(self.distortionfiles.t.loc[ix,"filename"])} for aperture, filter and/or pupil!')


            self.distortionfiles.t.loc[ix,self.aperture_col]=f'{aperture}'
            self.distortionfiles.t.loc[ix,self.filter_col]=f'{filt}'
            self.distortionfiles.t.loc[ix,self.pupil_col]=f'{pupil}'
        if self.verbose:
            print('##################\n### Distortion files:')
            self.distortionfiles.write()

    def match_distortion4ratefile(self, ixs_rate=None, require_filter=True, require_pupil=True,
                                  apertures=None, filts=None, pupils=None):
        
        ixs_matches = []
        ixs_not_matches = []
        
        ixs_rate = self.getindices(indices=ixs_rate)

        # if wanted, only use rate files with the specified apertures
        if apertures is not None:
            tmp_ixs = []
            for aperture in apertures:
                tmp_ixs.extend(self.ix_equal(self.aperture_col,aperture,indices=ixs_rate))
            ixs_rate = unique(tmp_ixs)
            ixs_rate.sort()
            
        # if wanted, only use rate files with the specified filters
        if filts is not None:
            tmp_ixs = []
            for filt in filts:
                tmp_ixs.extend(self.ix_equal(self.filter_col,filt,indices=ixs_rate))
            ixs_rate = unique(tmp_ixs)
            ixs_rate.sort()
        
        # if wanted, only use rate files with the specified pupils
        if pupils is not None:
            tmp_ixs = []
            for pupil in pupils:
                tmp_ixs.extend(self.ix_equal(self.pupil_col,pupil,indices=ixs_rate))
            ixs_rate = unique(tmp_ixs)
            ixs_rate.sort()
            
        
        for ix_rate in ixs_rate:
            ixs_coeff = self.distortionfiles.ix_equal(self.aperture_col,self.t.loc[ix_rate,self.aperture_col])
            #self.distortionfiles.write(indices=ixs_coeff)

            if require_filter:
                #print('filter cut')
                ixs_coeff = self.distortionfiles.ix_equal(self.filter_col,self.t.loc[ix_rate,self.filter_col],indices=ixs_coeff)
                #self.distortionfiles.write(indices=ixs_coeff)

            if require_pupil:
                #print('pupil cut')
                ixs_coeff = self.distortionfiles.ix_equal(self.pupil_col,self.t.loc[ix_rate,self.pupil_col],indices=ixs_coeff)
                #self.distortionfiles.write(indices=ixs_coeff)
            
            # check the matches!
            if len(ixs_coeff)==1:
                self.t.loc[ix_rate,'distortion_match']=self.distortionfiles.t.loc[ixs_coeff[0],'filename']
                ixs_matches.append(ix_rate)
            elif len(ixs_coeff)==0:
                print(f'WARNING! could not find match for {self.t.loc[ix_rate,"filename"]}')
                self.t.loc[ix_rate,'distortion_match']=np.nan
                ixs_not_matches.append(ix_rate)
            elif len(ixs_coeff)>1:
                print(f'ERROR: more than one match for {self.t.loc[ix_rate,"filename"]}!')
                self.distortionfiles.write(indices=ixs_coeff)
                raise RuntimeError(f'more than one match for {self.t.loc[ix_rate,"filename"]}!')
                
        print(f'{len(ixs_matches)} out of {len(ixs_rate)} matched!')
        if self.verbose:
            print('Matches:')
            self.write(indices=ixs_matches)
            
        if len(ixs_not_matches)>0:
            print(f'\n**********************************\n*** WARNING *** {len(ixs_not_matches)} out of {len(ixs_rate)} did not matched:')
            self.write(indices=ixs_not_matches)
            
        return(ixs_matches,ixs_not_matches)

    def apply_distortions(self, ixs, outdir='same_as_distortionfiles', 
                          skip_rate2cal_if_exists=None, outsubdir=None):
        self.t.loc[ixs,'errorflag'] = None
        for ix in ixs:
            distfile = self.t.loc[ix,'distortion_match']
            ratefile = self.t.loc[ix,'filename']
            applydist_singleim = test_distortion_singleim()
            applydist_singleim.verbose = self.verbose
            
            if outdir == 'same_as_distortionfiles':
                outdir = os.path.dirname(distfile)
                
            applydist_singleim.set_outdir(outdir)

            try:
                applydist_singleim.apply_distortions(ratefile,distfile,skip_rate2cal_if_exists=skip_rate2cal_if_exists)
                self.t.loc[ix,'errorflag']=False
            except:
                self.t.loc[ix,'errorflag']=True

if __name__ == '__main__':

    applydist = apply_distortions()
    parser = applydist.define_options()
    args = parser.parse_args()
    
    applydist.verbose=args.verbose
    
    applydist.get_rate_files(args.rate_files,directory=args.rate_dir)
    applydist.get_distortion_files(args.distortion_files,directory=None)
    
    ixs_matches,ixs_not_matches = applydist.match_distortion4ratefile(require_filter=not args.ignore_filters, 
                                                                     require_pupil=not args.ignore_pupils,
                                                                     apertures=args.apertures, 
                                                                     filts=args.filters, 
                                                                     pupils=args.pupils)
    
    if len(ixs_matches)==0:
        print('NO IMES FOUND!! exiting...')
        sys.exit(0)

    applydist.apply_distortions(ixs_matches,outdir=args.outdir,
                               skip_rate2cal_if_exists=args.skip_rate2cal_if_exists)
    applydist.write()