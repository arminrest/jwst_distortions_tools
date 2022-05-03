#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, argparse,glob,re
import pysiaf
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import ascii
from distortion2asdf import coeffs2asdf
from jwst import datamodels

font1 = {'family': 'helvetica', 'color': 'black', 'weight': 'normal', 'size': '12'}
font2 = {'family': 'helvetica', 'color': 'black', 'weight': 'normal', 'size': '20'}


def overplot_distortion_diffs(t1,t2,plotref=True,
                              XSciRef=1024.5,YSciRef=1024.5,
                              XSciSize=2048, YSciSize=2048):
    number_of_coefficients = len(t1['Sci2IdlX'])
    poly_degree = pysiaf.utils.polynomial.polynomial_degree(number_of_coefficients)

    nx, ny = (25, 25)
    x = np.linspace(1, XSciSize, nx)
    y = np.linspace(1, YSciSize, ny)
    x0 = XSciRef
    y0 = YSciRef
    #x0 = 0.0
    #y0 = 0.0
    xg, yg = np.meshgrid(x-x0, y-y0)

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

def plot_distortion_diffs(coeffref,coefflist,output_plot_name=None,showplot=False):
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
    print(f'Instrument: {coeffref.instrument} Aperture:{coeffref.aperture}')
    siafref = pysiaf.Siaf(coeffref.instrument)
    aperref = siafref[coeffref.aperture]
    
    print(f'XSciRef:{aperref.XSciRef} YSciRef:{aperref.YSciRef} :{aperref.XSciSize} YSciSize:{aperref.YSciSize}')
    
    vec_maxs = []
    for i,coeff in enumerate(coefflist):
        coeff.get_instrument_info()
        # error checking!
        if coeff.instrument!=coeffref.instrument:
            raise RuntimeError(f'inconsistent instruments: {coeff.instrument}!={coeffref.instrument}')
        if coeff.aperture!=coeffref.aperture:
            raise RuntimeError(f'inconsistent aperture: {coeff.aperture}!={coeffref.aperture}')

        vec_max = overplot_distortion_diffs(coeffref.t,coeff.t,plotref=(i==0),
                                            XSciRef=aperref.XSciRef,YSciRef=aperref.YSciRef,
                                            XSciSize=aperref.XSciSize, YSciSize=aperref.YSciSize)
        print(f'{coeff.filename} max vec: {vec_max}')
        vec_maxs.append(f'{vec_max:.4f}')
    plt.title(f'{coeffref.instrument} {coeffref.aperture}\n vec_maxs={vec_maxs}arcsec', fontdict=font2)
    plt.tight_layout()

    if output_plot_name is not None:
        print(f'Saving plot in {output_plot_name}')
        plt.savefig(output_plot_name)
        
    if showplot:
        plt.show()

def plot_distortionfiles_diffs(coefffileref,coefffilelist,output_plot_name=None,showplot=False):
    coeffref = coeffs2asdf()
    print(f'Loading reference coefficient file {coefffileref}')
    coeffref.load_coeff_file(coefffileref)
    
    coefflist = []
    for filename in coefffilelist:
        coeffs = coeffs2asdf()
        print(f'Loading {filename}')
        coeffs.load_coeff_file(filename)
        coefflist.append(coeffs)
   
    plot_distortion_diffs(coeffref,coefflist,output_plot_name=output_plot_name,showplot=showplot)
    
def define_options(parser=None,usage=None,conflict_handler='resolve'):
    if parser is None:
        parser = argparse.ArgumentParser(usage=usage,conflict_handler=conflict_handler)

    parser.add_argument('coeff_filepatterns', nargs='+', type=str, default=None, help='list of coefficient file(pattern)s. The first one will be the reference coefficients.')

    parser.add_argument('-v','--verbose', default=0, action='count')
    parser.add_argument('-p','--showplot', action='store_true', default=False, help='show the residual plot (default=%(default)s)')
    parser.add_argument('-s','--saveplot', default=None, help='Save the plot in the given filename')

    return(parser)

if __name__ == '__main__':

    parser = define_options()
    args = parser.parse_args()
    
    filenames=[]
    for filepattern in args.coeff_filepatterns:
        newfiles = glob.glob(filepattern)
        if len(newfiles)==0:
            raise RuntimeError(f'No files for file(pattern) {filepattern}')
        filenames.extend(newfiles)
    
    if len(filenames)<2:
        raise RuntimeError(f'At least 2 files required, only {filenames}')
    
    refile = filenames[0]
    filenames=filenames[1:]
    filenames.sort()

    plot_distortionfiles_diffs(refile,filenames,showplot=args.showplot,output_plot_name=args.saveplot)
