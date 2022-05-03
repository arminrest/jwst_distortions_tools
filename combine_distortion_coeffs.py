#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 10:35:28 2022

@author: arest
"""
import re,os,sys
import argparse,glob

# pdastroclass is wrapper around pandas.
from pdastro import pdastroclass,pdastrostatsclass,makepath4file,unique
import pandas as pd
from pandas.core.dtypes.common import is_string_dtype
from distortion2asdf import coeffs2asdf
from plot_distortion_diffs import plot_distortionfiles_diffs

class combine_coeffs(pdastrostatsclass):
    def __init__(self):
        pdastrostatsclass.__init__(self)
        
        self.verbose=0
        
        self.coeffs = ['Sci2IdlX','Sci2IdlY','Idl2SciX','Idl2SciY']
        self.aperture_col = 'AperName'
        self.siaf_index_col = 'siaf_index'
        self.cols2copy = [self.aperture_col,self.siaf_index_col,'exponent_x','exponent_y']
        
        self.format_coeff = '{:.10e}'
        
        self.results = pdastroclass()

    def define_options(self,parser=None,usage=None,conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)

        # default for token is $MAST_API_TOKEN
        if 'JWST_DISTORTION_OUTROOTDIR' in os.environ:
            outrootdir = os.environ['JWST_DISTORTION_OUTROOTDIR']
        else:
            outrootdir = None

        parser.add_argument('coeff_filepatterns', nargs='+', type=str, default=None, help='list of coefficient file(pattern)s')

        parser.add_argument('-v','--verbose', default=0, action='count')

        parser.add_argument('--apertures', default=None, nargs='+', choices=['NRCA1_FULL','NRCA2_FULL','NRCA3_FULL','NRCA4_FULL','NRCA5_FULL','NRCB1_FULL','NRCB2_FULL','NRCB3_FULL','NRCB4_FULL','NRCB5_FULL'], help='output root directory (default=%(default)s)')

        parser.add_argument('-s','--save_coefficients', action='store_true', default=False, help='Save the coefficients, using outrootdir, outsubdir, and outbasename')
        parser.add_argument('--outrootdir', default=outrootdir, help='output root directory (default=%(default)s)')
        parser.add_argument('--outsubdir', default=None, help='outsubdir added to output root directory (default=%(default)s)')
        parser.add_argument('--add2basename', default=None, help='This is added to the basename. (default=%(default)s)')
        parser.add_argument('--overwrite', default=False, action='store_true', help='overwrite files if they exist.')

        parser.add_argument('--ignore_filters', default=False, action='store_true', help='distortions are grouped by aperture/filter/pupil. Use this option if you want to create distortion files independent of filter.')
        parser.add_argument('--ignore_pupils', default=False, action='store_true', help='distortions are grouped by aperture/filter/pupil. Use this option if you want to create distortion files independent of pupil.')

        parser.add_argument('--showplot', default=False, action='store_true', help='show the difference of the input coefficients with respect to the combined distortion file')
        parser.add_argument('--saveplot', default=False, action='store_true', help='save the difference of the input coefficients with respect to the combined distortion file as pdf file')
        return(parser)
    
    def load_coeff_files(self,coeff_filepatterns,require_filter=True,require_pupil=True):
        counter = 0
        frames={}
        for filepattern in coeff_filepatterns:
            filenames = glob.glob(filepattern)
            for filename in filenames:
                if self.verbose: print(f'Loading {filename}')
                # read the file
                frames[counter] = pd.read_csv(filename,skipinitialspace=True,comment='#')

                # some of the column names have spaces, removed them!!
                mapper={}
                for col in frames[counter].columns:
                    if re.search('\s+',col):
                        mapper[col]=re.sub('\s+','',col)
                    if is_string_dtype(frames[counter][col]):
                        frames[counter][col] = frames[counter][col].str.strip()
                frames[counter] = frames[counter].rename(columns=mapper)
                
                # save the filename 
                frames[counter]['filename']=filename
                
                # get the filter and save it in the 'filter' column
                m = re.search('distortion_coeffs_[a-zA-Z0-9]+_[a-zA-Z0-9]+_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_',filename)
                if m is None:
                    if require_filter or require_pupil:
                        raise RuntimeError(f'could not parse filename {filename} for filter and/or pupil!')
                    else: 
                        print(f'WARNING! could not parse filename {filename} for filter and/or pupil!')
                    filt=pupil=None
                else:
                    filt,pupil = m.groups()

                frames[counter]['filter']=filt
                frames[counter]['pupil']=pupil

                
                if len(frames[counter])<1:
                    raise RuntimeError(f'file {filename} has no data!')
                                    
                counter+=1
        if self.verbose: print(f'Loaded {counter} coeff files')
        self.t = pd.concat(frames,ignore_index=True)

        for col in ['siaf_index','exponent_x','exponent_y']:       
            if col in self.t.columns:
                self.t[col]=self.t[col].astype('int')

        
    def calc_average_coefficients(self, aperture, indices=None, filt=None, pupil=None,
                                  saveflag=False,outbasename='test',overwrite=False,
                                  showplot=False,saveplot=False):
        print(f'\n######################################################\n### Determining coefficients for {aperture}, filter={filt}, pupil={pupil}\n######################################################')


        # if None, use all indices
        indices=self.getindices(indices=indices)

        siaf_indexs = unique(self.t.loc[indices,self.siaf_index_col])
        siaf_indexs.sort()
        
        inputfilenames = unique(self.t.loc[indices,'filename'])
        
        if self.verbose: print(f'siaf_indexs: {siaf_indexs}')
        
        self.set_statstring_format(format_floats=self.format_coeff)
        
        ixs_result = []
        ixs_input = []
        resultsstats = pdastrostatsclass()
        resultsstats.param2columnmapping = resultsstats.intializecols4statparams(format4outvals=self.format_coeff)

        
        for siaf_index in siaf_indexs:
            ixs_siaf_index = self.ix_equal(self.siaf_index_col, siaf_index, indices=indices)
            ixs_input.extend(ixs_siaf_index)
            if self.verbose>2:
                print(f'# input coefficients for siaf_index {siaf_index}:')
                self.write(indices=ixs_siaf_index)
            
            colvals={}
            for col in self.cols2copy:
                vals = unique(self.t.loc[ixs_siaf_index,col])
                if len(vals)!=1:
                    raise RuntimeError(f'BUG! more than one value in col {col}: {vals}')
                colvals[col]=vals[0]
            colvals['filter']=filt
            colvals['pupil']=pupil            
            #print(colvals)
            ix_result = self.results.newrow(colvals)
            ixs_result.append(ix_result)
            
            for coeff_col in self.coeffs:
                self.calcaverage_sigmacutloop(coeff_col,indices=ixs_siaf_index)
                self.results.t.loc[ix_result,coeff_col]=self.statparams['mean']
                
                if self.verbose>3: 
                    print(coeff_col,self.statstring())

                ix_resultstats = resultsstats.newrow(colvals)
                resultsstats.t.loc[ix_resultstats,'coeff']=coeff_col
                resultsstats.t.loc[ix_resultstats,'filter']=filt
                resultsstats.t.loc[ix_resultstats,'pupil']=pupil
                resultsstats.statresults2table(self.statparams,
                                                    resultsstats.param2columnmapping,
                                                    destindex=ix_resultstats)
                
        #make sure these columns are integer
        for col in ['siaf_index','exponent_x','exponent_y']:       
            if col in self.results.t.columns: self.results.t[col]=self.results.t[col].astype('int')
            if col in resultsstats.t.columns: resultsstats.t[col]=resultsstats.t[col].astype('int')
        
        if self.verbose:
            resultsstats.write()
            
        if saveflag:
            if not overwrite:
                if os.path.exists(f'{outbasename}.distcoeff.txt') or os.path.exists(f'{outbasename}.distcoeff.asdf'):
                    raise RuntimeError(f'File(s) {outbasename}.distcoeff.txt and/or {outbasename}.distcoeff.asdf already exist! Exiting. Use --overwrite if you want to overwrite them')

            formatters = {}
            for col in self.coeffs: formatters[col]=self.format_coeff.format


            # save the results
            makepath4file(outbasename)
            print(f'##### Saving coefficients for {outbasename}')
            print(f'Saving {outbasename}.distcoeff.txt')
            self.results.write(outbasename+'.distcoeff.txt',indices=ixs_result,formatters=formatters)

            

            coeffs = coeffs2asdf()
            print(f'Saving {outbasename}.distcoeff.asdf')
            coeffs.coefffile2adfs(outbasename+'.distcoeff.txt',outname=f'{outbasename}.distcoeff.asdf')

            print(f'Saving {outbasename}.singlefile.txt')
            self.write(outbasename+'.singlefile.txt',indices=ixs_input,formatters=formatters)

            print(f'Saving {outbasename}.statsinfo.txt')
            resultsstats.write(outbasename+'.statsinfo.txt',formatters=formatters)

            if showplot or saveplot:
                if saveplot:
                    plotfilename = f'{outbasename}.distcoeff.pdf'
                else:
                    plotfilename = None
                plot_distortionfiles_diffs(f'{outbasename}.distcoeff.txt',inputfilenames,
                                           showplot=showplot,output_plot_name=plotfilename)        


    def get_outdir(self,outrootdir,outsubdir):
        outdir = outrootdir
        if outdir is None: outdir = '.'
        
        if outsubdir is not None:
            outdir+=f'/{outsubdir}'
        
        return(outdir)
        
    def get_outbasename(self,outdir,aperture,filt=None,pupil=None,add2basename=None,Nfiles=None):
        outname = f'{outdir}'
        outname += f'/{aperture}'.lower()
        
        if filt is not None:
            outname += f'_{filt}'.lower()
        else:
            outname += '_na'
        
        if pupil is not None:
            outname += f'_{pupil}'.lower()
        else:
            outname += '_na'
        
        if add2basename is not None:
            outname += f'_{add2basename}'
            
        if Nfiles is not None:
            outname += f'.N{Nfiles:02}'
       
        return(outname)
        
 

    def calc_average_coefficients_all(self, indices=None, apertures=None, 
                                      require_filter=True, require_pupil=True,
                                      saveflag=True, overwrite=False,
                                      outrootdir=None,outsubdir=None,add2basename=None,
                                      showplot=False,saveplot=False
                                      ):
        # if None, use all indices
        indices=self.getindices(indices=indices)
        
        if apertures is None:
            apertures = unique(self.t.loc[indices,self.aperture_col])
        apertures.sort()
        
        outdir = self.get_outdir(outrootdir=outrootdir,outsubdir=outsubdir)
        
        print(f'Determining average coeff values for apertures {apertures}')
        
        for aperture in apertures:
            # get the indices for the given aperture
            ixs2use = self.ix_equal(self.aperture_col, aperture, indices=indices)
            if require_filter: 
                # get the indices for the given filter
                filts = unique(self.t.loc[ixs2use,'filter'])
                filts.sort()

                for filt in filts:
                    #print(f'Determining average coeff values for aperture {aperture}, filter {filt}')
                    ixs2use_filt = self.ix_equal('filter', filt, indices=ixs2use)
                    
                    if require_pupil: 
                        # get the indices for the given pupil
                        pupils = unique(self.t.loc[ixs2use_filt,'pupil'])
                        pupils.sort()

                        for pupil in pupils:
                            #print(f'Determining average coeff values for aperture {aperture}, filter {filt}, pupil {pupil}')
                            ixs2use_pupil = self.ix_equal('pupil', pupil, indices=ixs2use_filt)
                            outbasename = self.get_outbasename(outdir,aperture,filt=filt,pupil=pupil,add2basename=add2basename)
                            self.calc_average_coefficients(aperture,indices=ixs2use_pupil, 
                                                           filt=filt, pupil=pupil,
                                                           saveflag=saveflag,outbasename=outbasename,
                                                           overwrite=overwrite,
                                                           showplot=showplot,saveplot=saveplot)
                    else:
                        pupil=None
                        outbasename = self.get_outbasename(outdir,aperture,filt=filt,pupil=pupil,add2basename=add2basename)
                        self.calc_average_coefficients(aperture,indices=ixs2use_filt, 
                                                       filt=filt, pupil=pupil,
                                                       saveflag=saveflag,outbasename=outbasename,
                                                       overwrite=overwrite,
                                                       showplot=showplot,saveplot=saveplot)
                        
                    
            else:
                filt=pupil=None
                outbasename = self.get_outbasename(outdir,aperture,filt=filt,pupil=pupil,add2basename=add2basename)
                self.calc_average_coefficients(aperture,indices=ixs2use, 
                                               filt=filt, pupil=pupil,
                                               saveflag=saveflag,outbasename=outbasename,
                                               overwrite=overwrite,
                                               showplot=showplot,saveplot=saveplot)

        for col in ['siaf_index','exponent_x','exponent_y']:       
            if col in self.results.t.columns:
                self.results.t[col]=self.results.t[col].astype('int')

        self.results.write()

if __name__ == '__main__':

    coeffs = combine_coeffs()
    parser = coeffs.define_options()
    args = parser.parse_args()
    
    coeffs.verbose=args.verbose

    coeffs.load_coeff_files(args.coeff_filepatterns, 
                            require_filter = not args.ignore_filters,
                            require_pupil = not args.ignore_pupils)
    
    coeffs.calc_average_coefficients_all(apertures=args.apertures, 
                                         require_filter = not args.ignore_filters,
                                         require_pupil = not args.ignore_pupils,
                                         saveflag = args.save_coefficients,
                                         overwrite = args.overwrite,
                                         outrootdir = args.outrootdir, 
                                         outsubdir = args.outsubdir,
                                         add2basename = args.add2basename,
                                         showplot=args.showplot,
                                         saveplot=args.saveplot
                                         )
    
#    if args.save_coefficients:
#        coeffs.save_coefficients_all(args.outrootdir, 
#                                     args.outsubdir,
#                                     args.add2basename, 
#                                     apertures=None, 
#                                     require_filter = not args.ignore_filters,
#                                     require_pupil = not args.ignore_pupils,
#                                     overwrite=args.overwrite)