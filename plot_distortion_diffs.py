#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, argparse,glob,re,os
import pysiaf
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from distortion2asdf import coeffs2asdf
from jwst import datamodels
from pdastro import pdastroclass,unique

font1 = {'family': 'helvetica', 'color': 'black', 'weight': 'normal', 'size': '12'}
font2 = {'family': 'helvetica', 'color': 'black', 'weight': 'normal', 'size': '20'}




def overplot_distortion_diffs(t1,t2, xg, yg, plotref=True):
    number_of_coefficients = len(t1['Sci2IdlX'])
    poly_degree = pysiaf.utils.polynomial.polynomial_degree(number_of_coefficients)


    xg_idl1 = pysiaf.utils.polynomial.poly(t1['Sci2IdlX'], xg, yg, order=poly_degree)
    yg_idl1 = pysiaf.utils.polynomial.poly(t1['Sci2IdlY'], xg, yg, order=poly_degree)
    xg_idl2 = pysiaf.utils.polynomial.poly(t2['Sci2IdlX'], xg, yg, order=poly_degree)
    yg_idl2 = pysiaf.utils.polynomial.poly(t2['Sci2IdlY'], xg, yg, order=poly_degree)
    
    dx = xg_idl2 - xg_idl1
    dy = yg_idl2 - yg_idl1

    vec = np.sqrt(dx**2+dy**2)
    vec_max = np.max(vec)
    
    plt.quiver(xg_idl1, yg_idl1, dx,dy, color='blue')
    return(vec_max)

def get_mesh(detector,aperref,subarr,coron_region='all'):
    nx, ny = (25, 25)
    if subarr == 'FULL' or coron_region=='full':
        x = np.linspace(1, aperref.XSciSize, nx)
        y = np.linspace(1, aperref.YSciSize, ny)
        x0 = aperref.XSciRef
        y0 = aperref.YSciRef
        #x0 = 0.0
        #y0 = 0.0
        xg, yg = np.meshgrid(x-x0, y-y0)
        
        print(f'XSciRef:{aperref.XSciRef} YSciRef:{aperref.YSciRef} :{aperref.XSciSize} YSciSize:{aperref.YSciSize}')
    elif detector=='NRCA5':
        print(f'coron_region={coron_region}')
        x = np.linspace(171, 1800, nx)
        if coron_region=='all':
            #x = np.linspace(150, 1700, nx)
            #y = np.linspace(1, 1800, ny)
            y = np.linspace(1, 1820, ny)
            x0 = aperref.XSciRef
            y0 = aperref.YSciRef
            #x0 = 0.0
            #y0 = 0.0
            print(x,y)
            xg, yg = np.meshgrid(x-x0, y-y0)
        elif coron_region=='top':
            #x = np.linspace(150, 1700, nx)
            y = np.linspace(1470, 1820, ny)
            x0 = aperref.XSciRef
            y0 = aperref.YSciRef
            #x0 = 0.0
            #y0 = 0.0
            xg, yg = np.meshgrid(x-x0, y-y0)
        elif coron_region=='topcore':
            x = np.linspace(300, 1650, nx)
            y = np.linspace(1520, 1750, ny)
            #y = np.linspace(1, 1750, ny)
            x0 = aperref.XSciRef
            y0 = aperref.YSciRef
            #x0 = 0.0
            #y0 = 0.0
            xg, yg = np.meshgrid(x-x0, y-y0)
        elif coron_region=='bottom':
            #x = np.linspace(150, 1700, nx)
            y = np.linspace(1, 1470, ny)
            x0 = aperref.XSciRef
            y0 = aperref.YSciRef
            #x0 = 0.0
            #y0 = 0.0
            xg, yg = np.meshgrid(x-x0, y-y0)
    elif detector=='NRCA2':
        print(f'coron_region={coron_region}')
        x = np.linspace(401, 2048, nx)
        if coron_region=='all':
            #x = np.linspace(150, 1700, nx)
            #y = np.linspace(1, 1800, ny)
            y = np.linspace(1, 1800, ny)
            x0 = aperref.XSciRef
            y0 = aperref.YSciRef
            #x0 = 0.0
            #y0 = 0.0
            print(x,y)
            xg, yg = np.meshgrid(x-x0, y-y0)
        elif coron_region=='top':
            #x = np.linspace(150, 1700, nx)
            y = np.linspace(1151, 1800, ny)
            x0 = aperref.XSciRef
            y0 = aperref.YSciRef
            #x0 = 0.0
            #y0 = 0.0
            xg, yg = np.meshgrid(x-x0, y-y0)
        elif coron_region=='topcore':
            x = np.linspace(500, 1900, nx)
            y = np.linspace(1250, 1700, ny)
            #y = np.linspace(1, 1750, ny)
            x0 = aperref.XSciRef
            y0 = aperref.YSciRef
            #x0 = 0.0
            #y0 = 0.0
            xg, yg = np.meshgrid(x-x0, y-y0)
        elif coron_region=='bottom':
            #x = np.linspace(150, 1700, nx)
            y = np.linspace(1, 1150, ny)
            x0 = aperref.XSciRef
            y0 = aperref.YSciRef
            #x0 = 0.0
            #y0 = 0.0
            xg, yg = np.meshgrid(x-x0, y-y0)
        else:
            raise RuntimeError(f'coron_region={coron_region} not known!')
    elif detector=='NRCA4':
        print(f'coron_region={coron_region}')
        x = np.linspace(1, 1500, nx)
        if coron_region=='all':
            #x = np.linspace(150, 1700, nx)
            #y = np.linspace(1, 1800, ny)
            y = np.linspace(1, 1700, ny)
            x0 = aperref.XSciRef
            y0 = aperref.YSciRef
            #x0 = 0.0
            #y0 = 0.0
            print(x,y)
            xg, yg = np.meshgrid(x-x0, y-y0)
        elif coron_region=='top':
            #x = np.linspace(150, 1700, nx)
            y = np.linspace(1151, 1800, ny)
            x0 = aperref.XSciRef
            y0 = aperref.YSciRef
            #x0 = 0.0
            #y0 = 0.0
            xg, yg = np.meshgrid(x-x0, y-y0)
        elif coron_region=='topcore':
            x = np.linspace(100, 1400, nx)
            y = np.linspace(1250, 1700, ny)
            #y = np.linspace(1, 1750, ny)
            x0 = aperref.XSciRef
            y0 = aperref.YSciRef
            #x0 = 0.0
            #y0 = 0.0
            xg, yg = np.meshgrid(x-x0, y-y0)
        elif coron_region=='bottom':
            #x = np.linspace(150, 1700, nx)
            y = np.linspace(1, 1150, ny)
            x0 = aperref.XSciRef
            y0 = aperref.YSciRef
            #x0 = 0.0
            #y0 = 0.0
            xg, yg = np.meshgrid(x-x0, y-y0)
        else:
            raise RuntimeError(f'coron_region={coron_region} not known!')
    else:
        raise RuntimeError(f'subarr={subarr} and/or detector {detector} not known!')
    return(xg,yg)

def plot_distortion_diffs(coeffref,coefflist, coron_region='all', output_plot_name=None,showplot=False):
    plt.clf()
    plt.rc('font', family='serif')
    plt.figure(figsize=(12,12))
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)
    plt.xlabel("x_idl [arcsec]", fontdict=font2)
    plt.ylabel("y_idl [arcsec]", fontdict=font2)
    #plt.rcParams['axes.titlesize'] = plt.rcParams['axes.labelsize'] = 30
    #plt.rcParams['xtick.labelsize'] = plt.rcParams['ytick.labelsize'] = 20

    coeffref.get_instrument_info()
    print(f'Instrument: {coeffref.instrument} Aperture:{coeffref.aperture} subarray:{coeffref.subarr}')
    siafref = pysiaf.Siaf(coeffref.instrument)
    aperref = siafref[coeffref.aperture]
    
    (xg,yg) = get_mesh(coeffref.detector, aperref,coeffref.subarr,coron_region=coron_region)

    vec_maxs = []
    for i,coeff in enumerate(coefflist):
        coeff.get_instrument_info()
        # error checking!
        if coeff.instrument!=coeffref.instrument:
            raise RuntimeError(f'inconsistent instruments: {coeff.instrument}!={coeffref.instrument}')
        if coeff.aperture!=coeffref.aperture:
            raise RuntimeError(f'inconsistent aperture: {coeff.aperture}!={coeffref.aperture}')

        vec_max = overplot_distortion_diffs(coeffref.t,coeff.t,xg,yg,plotref=(i==0))
        print(f'{coeff.filename} {vec_max} maxvec_arcsec')
        vec_maxs.append(float(f'{vec_max*1000:.1f}'))
    max_vec_maxs = np.amax(np.array(vec_maxs))
    #m1 = re.search(f'distortion_coeffs_[a-zA-Z0-9]+_[a-zA-Z0-9]+_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_',os.path.basename(coeffref.filename))
    #m2 = re.search('^[a-zA-Z0-9]+_[a-zA-Z0-9]+_([a-zA-Z0-9]+)_([a-zA-Z0-9]+).*\.distcoeff\.txt',os.path.basename(coeffref.filename))
    m1 = re.search(f'distortion_coeffs_{coeffref.aperture.lower()}_([a-zA-Z0-9]+)_([a-zA-Z0-9]+)_jw',os.path.basename(coeffref.filename))
    m2 = re.search(f'^{coeffref.aperture.lower()}_([a-zA-Z0-9]+)_([a-zA-Z0-9]+).*\.distcoeff\.txt',os.path.basename(coeffref.filename))
    if m1 is not None:
        filt,pupil = m1.groups()
    elif m2 is not None:
        filt,pupil = m2.groups()        
    else:
        print('Could not determine filter/pupil!!! {}')
        filt,pupil = (None,None)
    title = f'{coeffref.instrument} {coeffref.aperture} filt/pupil={filt}/{pupil} max(vec_maxs)={max_vec_maxs}mas\n vec_maxs={vec_maxs}mas'
    plt.title(title, fontdict=font2)

    plt.tight_layout()

    if output_plot_name is not None:
        print(f'Saving plot in {output_plot_name}')
        plt.savefig(output_plot_name)
        
    if showplot:
        plt.show()

def plot_distortionfiles_diffs(coefffileref,coefffilelist,
                               coron_region = 'all',
                               output_plot_name=None,
                               showplot=False):
    coeffref = coeffs2asdf()
    print(f'Loading reference coefficient file {coefffileref}')
    coeffref.load_coeff_file(coefffileref)
    
    coefflist = []
    for filename in coefffilelist:
        coeffs = coeffs2asdf()
        print(f'Loading {filename}')
        coeffs.load_coeff_file(filename)
        coefflist.append(coeffs)
    plot_distortion_diffs(coeffref,coefflist,coron_region=coron_region,
                          output_plot_name=output_plot_name,showplot=showplot)
    
def define_options(parser=None,usage=None,conflict_handler='resolve'):
    if parser is None:
        parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)

    parser.add_argument('coeff_filepatterns', nargs='+', type=str, default=None, help='list of coefficient file(pattern)s. The first one will be the reference coefficients.')

    parser.add_argument('--coron_region', type=str, default='all', choices=['top','topcore','bottom','all','full'], help='for coronography: specify the region of interest (default=%(default)s)')
    

    parser.add_argument('-v','--verbose', default=0, action='count')
    parser.add_argument('-p','--showplot', action='store_true', default=False, help='show the residual plot (default=%(default)s)')
    parser.add_argument('-s','--saveplot', default=None, help='Save the plot in the given filename')

    return(parser)

if __name__ == '__main__':

    parser = define_options()
    args = parser.parse_args()
    
    filenames=[]
    for filepattern in args.coeff_filepatterns:
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

    if len(filenames)<2:
        raise RuntimeError(f'At least 2 files required, only {filenames}')
    
    refile = filenames[0]
    filenames=filenames[1:]
    filenames.sort()

    plot_distortionfiles_diffs(refile,filenames,
                               coron_region = args.coron_region,
                               showplot=args.showplot,
                               output_plot_name=args.saveplot)
