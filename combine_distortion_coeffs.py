#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 10:35:28 2022

@author: arest
"""
import re,os,sys
import argparse,glob,copy

# pdastroclass is wrapper around pandas.
from pdastro import pdastroclass,pdastrostatsclass,makepath4file,unique
import pandas as pd
from pandas.core.dtypes.common import is_string_dtype
from distortion2asdf import coeffs2asdf
from plot_distortion_diffs import plot_distortionfiles_diffs
import pysiaf
import numpy as np

pid_info = {}
pid_info[1069] = {'description':'Distortion reference file created from NRC-21 (PID=1069) data during commissioning',
                  'author':'V. Platais, E. Egami, A. Rest',
                  'pedigree':'INFLIGHT 2022-04-27 2022-04-28'}
pid_info[1070] = {'description':'Distortion reference file created from NRC-21b (PID=1070) data during commissioning',
                  'author':'E. Egami, J. Girard, A. Rest',
                  'pedigree':'INFLIGHT 2022-04-27 2022-04-28'}


class combine_coeffs(pdastrostatsclass):
    def __init__(self):
        pdastrostatsclass.__init__(self)
        
        self.verbose=0
        
        self.coeffs = ['Sci2IdlX','Sci2IdlY','Idl2SciX','Idl2SciY']
        self.aperture_col = 'AperName'
        self.siaf_index_col = 'siaf_index'
        self.cols2copy = [self.aperture_col,self.siaf_index_col,'exponent_x','exponent_y']
        
        self.format_coeff = '{:.10e}'
        
        self.results = coeffs2asdf()
        
        self.overview = pdastroclass()

    def define_options(self,parser=None,usage=None,conflict_handler='resolve'):
        if parser is None:
            parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)

        # default for token is $MAST_API_TOKEN
        if 'JWST_DISTORTION_OUTROOTDIR' in os.environ:
            outrootdir = os.environ['JWST_DISTORTION_OUTROOTDIR']
        else:
            outrootdir = None

        parser.add_argument('coeff_filepatterns', nargs='+', type=str, default=None, help='list of coefficient file(pattern)s. This can also be a *.singlefile.txt (output from this script) or a *goodfiles.txt with a column "filename", and the files in that list are used')

        parser.add_argument('--skip_if_file_not_exists', default=False, action='store_true', help='Do not throw an error if an input file does not exist, just skip it.')


        parser.add_argument('-v','--verbose', default=0, action='count')

        parser.add_argument('--apertures', default=None, nargs='+', choices=['NRCA1_FULL','NRCA2_FULL','NRCA3_FULL','NRCA4_FULL','NRCA5_FULL','NRCB1_FULL','NRCB2_FULL','NRCB3_FULL','NRCB4_FULL','NRCB5_FULL'], help='output root directory (default=%(default)s)')

        parser.add_argument('-s','--save_coefficients', action='store_true', default=False, help='Save the coefficients, using outrootdir, outsubdir, and outbasename')
        parser.add_argument('--outrootdir', default=outrootdir, help='output root directory (default=%(default)s)')
        parser.add_argument('--outsubdir', default=None, help='outsubdir added to output root directory (default=%(default)s)')
        parser.add_argument('--add2basename', default=None, help='This is added to the basename. (default=%(default)s)')
        parser.add_argument('--overwrite', default=False, action='store_true', help='overwrite files if they exist.')

        parser.add_argument('--ignore_filters', default=False, action='store_true', help='distortions are grouped by aperture/filter/pupil. Use this option if you want to create distortion files independent of filter.')
        parser.add_argument('--ignore_pupils', default=False, action='store_true', help='distortions are grouped by aperture/filter/pupil. Use this option if you want to create distortion files independent of pupil.')

        parser.add_argument('--vecmax_limits_mas', nargs='+', type=float, help='The cuts on vec_max are applied in order')

        parser.add_argument('--showplot', default=False, action='store_true', help='show the difference of the input coefficients with respect to the combined distortion file')
        parser.add_argument('--saveplot', default=False, action='store_true', help='save the difference of the input coefficients with respect to the combined distortion file as pdf file')
        parser.add_argument('--coron_region', type=str, default='all', choices=['top','topcore','bottom','all','full'], help='for coronography: specify the region of interest to be plotted (default=%(default)s)')

        parser.add_argument('--save_overview', default=None, help='Save the overview table (default=%(default)s)')


        return(parser)
    
    def load_coeff_files(self,coeff_filepatterns,
                         skip_if_file_not_exists = False,
                         require_filter=True,require_pupil=True):
        counter = 0
        frames={}
        filenames = []
        for filepattern in coeff_filepatterns:
            if re.search('singlefile\.txt$',filepattern):
                # get the input filenames from the singlefile files!
                infofiles=glob.glob(filepattern)
                if len(infofiles)==0:
                    raise RuntimeError(f'Could not find any files that match {filepattern}')
                for infofile in infofiles:
                    info = pdastroclass()
                    info.load(infofile)
                    filenames.extend(unique(info.t['filename']))
            elif re.search('goodfiles\.txt$',filepattern):
                # get the input filenames from the singlefile files!
                infofiles=glob.glob(filepattern)
                if len(infofiles)==0:
                    raise RuntimeError(f'Could not find any files that match {filepattern}')
                for infofile in infofiles:
                    info = pdastroclass()
                    info.load(infofile,comment='#')
                    filenames.extend(unique(info.t['filename']))
            else:
                newfilenames = glob.glob(filepattern)
                if len(newfilenames)==0:
                    raise RuntimeError(f'Could not find any files that match {filepattern}')
                filenames.extend(newfilenames)            
        filenames = unique(filenames)
        filenames.sort()

        for filename in filenames:
            if not os.path.isfile(filename):
                if skip_if_file_not_exists:
                    print(f'\n*** WARNING ****\n file{filename} does not exist! skipping')
                    continue
                else:
                    raise RuntimeError(f'file{filename} does not exist!')
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
            
            #print(frames[counter]['AperName'],unique(frames[counter]['AperName']))
            aperture = unique(frames[counter]['AperName'])[0]
            
            # get the filter and save it in the 'filter' column
            m1 = re.search(f'distortion_coeffs_{aperture.lower()}_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_jw',os.path.basename(filename))
            m2 = re.search(f'^{aperture.lower()}_([a-zA-Z0-9]+)_([a-zA-Z0-9]+).*\.distcoeff\.txt',os.path.basename(filename))
            if m1 is not None:
                filt,pupil = m1.groups()
            elif m2 is not None:
                filt,pupil = m2.groups()        
            else:
                if require_filter or require_pupil:
                    raise RuntimeError(f'could not parse filename {os.path.basename(filename)} for filter and/or pupil!')
                else: 
                    print(f'WARNING! could not parse filename {os.path.basename(filename)} for filter and/or pupil!')
                filt=pupil=None
            
            #m = re.search('distortion_coeffs_[a-zA-Z0-9]+_[a-zA-Z0-9]+_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_',os.path.basename(filename))
            #if m is None:
            #    if require_filter or require_pupil:
            #        raise RuntimeError(f'could not parse filename {filename} for filter and/or pupil!')
            #    else: 
            #        print(f'WARNING! could not parse filename {filename} for filter and/or pupil!')
            #    filt=pupil=None
            #else:
            #    filt,pupil = m.groups()

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

    def get_mesh(self,detector,aperref,subarr,nx=25,ny=25,coron_region='all'):
        x0 = aperref.XSciRef
        y0 = aperref.YSciRef
        if subarr == 'FULL' or coron_region=='full' or (detector in ['NRCA1','NRCA3']):
            x = np.linspace(1, aperref.XSciSize, nx)
            y = np.linspace(1, aperref.YSciSize, ny)
        elif detector=='NRCA5':
            print(f'coron_region={coron_region}')
            x = np.linspace(171, 1800, nx)
            if coron_region=='all':
                y = np.linspace(1, 1820, ny)
            elif coron_region=='top':
                y = np.linspace(1470, 1820, ny)
            elif coron_region=='topcore':
                x = np.linspace(300, 1650, nx)
                y = np.linspace(1520, 1750, ny)
            elif coron_region=='bottom':
                y = np.linspace(1, 1470, ny)
        elif detector=='NRCA2':
            print(f'coron_region={coron_region}')
            x = np.linspace(401, 2048, nx)
            if coron_region=='all':
                y = np.linspace(1, 1800, ny)
            elif coron_region=='top':
                y = np.linspace(1151, 1800, ny)
            elif coron_region=='topcore':
                x = np.linspace(500, 1900, nx)
                y = np.linspace(1250, 1700, ny)
            elif coron_region=='bottom':
                y = np.linspace(1, 1150, ny)
            else:
                raise RuntimeError(f'coron_region={coron_region} not known!')
        elif detector=='NRCA4':
            print(f'coron_region={coron_region}')
            x = np.linspace(1, 1500, nx)
            if coron_region=='all':
                y = np.linspace(1, 1700, ny)
            elif coron_region=='top':
                y = np.linspace(1151, 1800, ny)
            elif coron_region=='topcore':
                x = np.linspace(100, 1400, nx)
                y = np.linspace(1250, 1700, ny)
            elif coron_region=='bottom':
                y = np.linspace(1, 1150, ny)
            else:
                raise RuntimeError(f'coron_region={coron_region} not known!')
        else:
            raise RuntimeError(f'subarr={subarr} and/or detector {detector} not known!')
        xg, yg = np.meshgrid(x-x0, y-y0)
        return(xg,yg)


    def distortion_diffs_vecmax(self, coeffref, t2, coron_region='all'):
        
        print(f'TESTTTTTT: {coeffref.instrument} {coeffref.aperture} {coeffref.detector} {coeffref.subarr}') 
        #print('t2:',t2)
        #print('t2Sci:',t2['Sci2IdlX'])
        #sys.exit(0)
        siafref = pysiaf.Siaf(coeffref.instrument)
        aperref = siafref[coeffref.aperture]
        (xg,yg) = self.get_mesh(coeffref.detector, aperref,coeffref.subarr,coron_region=coron_region)
        
        
        number_of_coefficients = len(coeffref.t['Sci2IdlX'])
        poly_degree = pysiaf.utils.polynomial.polynomial_degree(number_of_coefficients)
    
    
        xg_idl1 = pysiaf.utils.polynomial.poly(list(coeffref.t['Sci2IdlX']), xg, yg, order=poly_degree)
        yg_idl1 = pysiaf.utils.polynomial.poly(list(coeffref.t['Sci2IdlY']), xg, yg, order=poly_degree)
        xg_idl2 = pysiaf.utils.polynomial.poly(list(t2['Sci2IdlX']), xg, yg, order=poly_degree)
        yg_idl2 = pysiaf.utils.polynomial.poly(list(t2['Sci2IdlY']), xg, yg, order=poly_degree)
        
        dx = xg_idl2 - xg_idl1
        dy = yg_idl2 - yg_idl1
    
        vec = np.sqrt(dx**2+dy**2)
        vec_max = np.max(vec)
        
        return(vec_max,dx,dy)
        
    def calc_average_coefficients(self, aperture, indices=None, filt=None, pupil=None,
                                  saveflag=False,outbasename='test',overwrite=False,
                                  showplot=False,saveplot=False,coron_region='all'):


        # if None, use all indices
        indices=self.getindices(indices=indices)

        print(f'\n######################################################\n### Determining coefficients for {aperture}, filter={filt}, pupil={pupil}: {len(indices)} input files\n######################################################')

        siaf_indexs = unique(self.t.loc[indices,self.siaf_index_col])
        siaf_indexs.sort()
        
        inputfilenames = unique(self.t.loc[indices,'filename'])
        inputfilenames.sort()
        
        if self.verbose: print(f'siaf_indexs: {siaf_indexs}')
        
        self.set_statstring_format(format_floats=self.format_coeff)
        
        ixs_result = []
        ixs_input = []
        
        self.results = coeffs2asdf()
        
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
            
        #### Calculate vec_max!!!
        self.vecmax = pdastroclass()
        self.results.get_instrument_info()
        
        # This is the slower, but cleaner way
        #for inputfilename in inputfilenames:
            #continue
        #    ixs_file = self.ix_equal('filename',inputfilename,indices=indices)
        #    (vec_max,dx,dy) = self.distortion_diffs_vecmax_(self.results, self.t.loc[ixs_file])

        siafref = pysiaf.Siaf(self.results.instrument)
        aperref = siafref[self.results.aperture]
        (xg,yg) = self.get_mesh(self.results.detector, aperref,self.results.subarr,coron_region=coron_region)
        number_of_coefficients = len(self.results.t['Sci2IdlX'])
        poly_degree = pysiaf.utils.polynomial.polynomial_degree(number_of_coefficients)
        xg_idl1 = pysiaf.utils.polynomial.poly(list(self.results.t['Sci2IdlX']), xg, yg, order=poly_degree)
        yg_idl1 = pysiaf.utils.polynomial.poly(list(self.results.t['Sci2IdlY']), xg, yg, order=poly_degree)

        for inputfilename in inputfilenames:
            ixs_file = self.ix_equal('filename',inputfilename,indices=indices)
            xg_idl2 = pysiaf.utils.polynomial.poly(list(self.t.loc[ixs_file,'Sci2IdlX']), xg, yg, order=poly_degree)
            yg_idl2 = pysiaf.utils.polynomial.poly(list(self.t.loc[ixs_file,'Sci2IdlY']), xg, yg, order=poly_degree)

            dx = xg_idl2 - xg_idl1
            dy = yg_idl2 - yg_idl1
        
            vec = np.sqrt(dx**2+dy**2)
            vec_max = np.max(vec)
            
            self.vecmax.newrow({self.aperture_col:self.results.aperture,
                                'filter':filt,
                                'pupil':pupil,
                                'vec_max_mas':vec_max*1000.0,
                                'filename':inputfilename,
                                'outbasename':outbasename
                                })

        self.vecmax.write()
        
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

            # convert to asdf file: get history etc
            history = ['Files used to create this reference file:']
            filenames = unique(self.t.loc[ixs_input,'filename'])
            filenames.sort()
            
            # get the PID in order to get author etc defined at the top of the script
            pids=[]
            for filename in filenames:
                m = re.search('_jw(\d\d\d\d\d)\d\d\d\d\d\d_',filename)
                if m is not None:    
                    pids.append(m.groups()[0])
                history.append(os.path.basename(filename))
            pids = unique(pids)
            if len(pids)==1 and int(pids[0]) in pid_info:
                pid=int(pids[0])
                author = pid_info[pid]['author']
                descrip = pid_info[pid]['description']
                pedigree = pid_info[pid]['pedigree']
                print(f'pid={pid}, setting author={author}, pedigree={pedigree}, and description={descrip}!')
            else: 
                print(f'\n!!!!!!!!!!!!!!!!\n!!!! WARNING!!\n!!!! pids={pids} not in pid_info, therefore not setting author, pedigree, and description!')
                author = descrip = pedigree = None
            print(f'Saving {outbasename}.distcoeff.asdf')
            #print(history,pids)
            #sys.exit(0)
            self.results.coefffile2adfs(outbasename+'.distcoeff.txt', filt=filt, 
                                        outname=f'{outbasename}.distcoeff.asdf',
                                        history=history,
                                        author=author, 
                                        descrip=descrip, 
                                        pedigree=pedigree)

            print(f'Saving {outbasename}.singlefile.txt')
            self.write(outbasename+'.singlefile.txt',indices=ixs_input,formatters=formatters)

            print(f'Saving {outbasename}.statsinfo.txt')
            resultsstats.write(outbasename+'.statsinfo.txt',formatters=formatters)

            print(f'Saving {outbasename}.vec_max.txt')
            self.vecmax.write(outbasename+'.vec_max.txt')

            if showplot or saveplot:
                if saveplot:
                    plotfilename = f'{outbasename}.distcoeff.pdf'
                else:
                    plotfilename = None
                plot_distortionfiles_diffs(f'{outbasename}.distcoeff.txt',inputfilenames,
                                           coron_region=args.coron_region,
                                           #save_vec_max= args.save_vec_max,
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
            outname += f'.{add2basename}'
            
        if Nfiles is not None:
            outname += f'.N{Nfiles:02}'
       
        return(outname)
        
    def calc_average_coefficients_vecmax_cut(self, aperture, vecmax_limits_mas=None, 
                                             indices=None, filt=None, pupil=None,
                                             **kwargs):

        print(f'\n######################################################\n### Determining coefficients for {aperture}, filter={filt}, pupil={pupil}\n######################################################')


        # if None, use all indices
        indices=self.getindices(indices=indices)
        counter=1
        
        kwargs0 = copy.deepcopy(kwargs)
        if vecmax_limits_mas is not None:
            kwargs0['saveplot']=False
            kwargs0['showplot']=False
            kwargs0['saveflag']=False
        self.calc_average_coefficients(aperture,indices=indices, filt=filt, pupil=pupil, **kwargs0)
        if vecmax_limits_mas is None or (vecmax_limits_mas==[]):
            return(0)
        
        for counter,vecmax_limit in enumerate(vecmax_limits_mas):
            
            ixs_vecmax_cut = self.vecmax.ix_inrange('vec_max_mas',None,vecmax_limit)
            print(f'!!!!!!!!!!!!!!!!!!!!!!!! ### vec_max cut: keeping {len(ixs_vecmax_cut)} out of {len(self.vecmax.getindices())} with vec_max<={vecmax_limit}mas')
            ixs_keep = []
            for ix_vecmax in ixs_vecmax_cut:
                filename = self.vecmax.t.loc[ix_vecmax,'filename']
                ixs_keep.extend(self.ix_equal('filename',filename,indices=indices))
            
            if counter==len(vecmax_limits_mas)-1:
                # final calculation! save the results, plots etc
                self.calc_average_coefficients(aperture,indices=ixs_keep, filt=filt, pupil=pupil, **kwargs)
                self.vecmax.write()
                Nin = len(unique(self.t.loc[indices,'filename']))
                Ngood = len(unique(self.t.loc[ixs_keep,'filename']))
                self.overview.newrow({self.aperture_col:self.results.aperture,
                                      'filter':filt,
                                      'pupil':pupil,
                                      'Nin':Nin,
                                      'Ngood':Ngood,
                                      'Ncut':Nin-Ngood,
                                      'median_vec_max_mas':np.median(self.vecmax.t['vec_max_mas'])
                                      })
            else:
                self.calc_average_coefficients(aperture,indices=ixs_keep, filt=filt, pupil=pupil, **kwargs0)
                

 

    def calc_average_coefficients_all(self,  vecmax_limits_mas=None, 
                                      indices=None, apertures=None, 
                                      require_filter=True, require_pupil=True,
                                      saveflag=True, overwrite=False,
                                      outrootdir=None,outsubdir=None,add2basename=None,
                                      showplot=False,saveplot=False,
                                      save_overview=None
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
                            self.calc_average_coefficients_vecmax_cut(aperture,
                                                           vecmax_limits_mas=vecmax_limits_mas,
                                                           indices=ixs2use_pupil, 
                                                           filt=filt, pupil=pupil,
                                                           saveflag=saveflag,outbasename=outbasename,
                                                           overwrite=overwrite,
                                                           showplot=showplot,saveplot=saveplot)
                    else:
                        pupil=None
                        outbasename = self.get_outbasename(outdir,aperture,filt=filt,pupil=pupil,add2basename=add2basename)
                        self.calc_average_coefficients_vecmax_cut(aperture,
                                                       vecmax_limits_mas=vecmax_limits_mas,
                                                       indices=ixs2use_filt, 
                                                       filt=filt, pupil=pupil,
                                                       saveflag=saveflag,outbasename=outbasename,
                                                       overwrite=overwrite,
                                                       showplot=showplot,saveplot=saveplot)
                        
                    
            else:
                filt=pupil=None
                outbasename = self.get_outbasename(outdir,aperture,filt=filt,pupil=pupil,add2basename=add2basename)
                self.calc_average_coefficients_vecmax_cut(aperture,
                                               vecmax_limits_mas=vecmax_limits_mas,
                                               indices=ixs2use, 
                                               filt=filt, pupil=pupil,
                                               saveflag=saveflag,outbasename=outbasename,
                                               overwrite=overwrite,
                                               showplot=showplot,saveplot=saveplot)

        for col in ['siaf_index','exponent_x','exponent_y']:       
            if col in self.results.t.columns:
                self.results.t[col]=self.results.t[col].astype('int')

        print('HELLO',showplot)
        self.overview.write()
        if save_overview is not None:
            overview_filename=f'{outdir}/{os.path.basename(save_overview)}'
            print(f'Saving overview to {overview_filename}')
            self.overview.write(overview_filename)
           
            

if __name__ == '__main__':

    coeffs = combine_coeffs()
    parser = coeffs.define_options()
    args = parser.parse_args()
    
    coeffs.verbose=args.verbose

    coeffs.load_coeff_files(args.coeff_filepatterns, 
                            skip_if_file_not_exists = args.skip_if_file_not_exists,
                            require_filter = not args.ignore_filters,
                            require_pupil = not args.ignore_pupils)
    
    coeffs.calc_average_coefficients_all(apertures=args.apertures, 
                                         require_filter = not args.ignore_filters,
                                         require_pupil = not args.ignore_pupils,
                                         vecmax_limits_mas = args.vecmax_limits_mas,
                                         saveflag = args.save_coefficients,
                                         overwrite = args.overwrite,
                                         outrootdir = args.outrootdir, 
                                         outsubdir = args.outsubdir,
                                         add2basename = args.add2basename,
                                         showplot=args.showplot,
                                         saveplot=args.saveplot,
                                         save_overview=args.save_overview
                                         )
    
#    if args.save_coefficients:
#        coeffs.save_coefficients_all(args.outrootdir, 
#                                     args.outsubdir,
#                                     args.add2basename, 
#                                     apertures=None, 
#                                     require_filter = not args.ignore_filters,
#                                     require_pupil = not args.ignore_pupils,
#                                     overwrite=args.overwrite)