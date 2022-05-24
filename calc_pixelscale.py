#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys, argparse,glob,re,os
import pysiaf
import matplotlib.pyplot as plt
import numpy as np
from distortion2asdf import coeffs2asdf
from astropy.io import fits


class calc_pixelscale(coeffs2asdf):
    def __init__(self):
        coeffs2asdf.__init__(self)
        
        self.verbose=0
        
        self.x0=self.y0=1024.5
        self.Nx=self.Ny=2048
        self.nx_mesh=self.ny_mesh=16
        
    def get_mesh(self,Nx=None,Ny=None,nx_mesh=None,ny_mesh=None,x0=None,y0=None):
        if Nx is None: Nx=self.Nx
        if Ny is None: Ny=self.Ny
        if nx_mesh is None: nx_mesh=self.nx_mesh
        if ny_mesh is None: ny_mesh=self.ny_mesh
        if x0 is None: x0=self.x0
        if y0 is None: y0=self.y0
        
        
        x = np.linspace(1, Nx, nx_mesh)
        y = np.linspace(1, Ny, ny_mesh)
        xg, yg = np.meshgrid(x-x0, y-x0)
        return(xg,yg)

    def get_mesh_fullframe(self,Nx=None,Ny=None,x0=None,y0=None):
        if Nx is None: Nx=self.Nx
        if Ny is None: Ny=self.Ny
        if x0 is None: x0=self.x0
        if y0 is None: y0=self.y0

        x = np.arange(1, Nx+1)
        y = np.arange(1, Ny+1)
        xg, yg = np.meshgrid(x-x0, y-y0)
        return(xg,yg)

    def calc_pxscale_mesh(self):
        self.get_instrument_info()
        print(f'Instrument: {self.instrument} Detector:{self.detector} Aperture:{self.aperture} subarray:{self.subarr}')
        
        (xg,yg) = self.get_mesh_fullframe()
        print(xg,yg)

        x_pixscl = pysiaf.utils.polynomial.dpdx(self.t['Sci2IdlX'], xg, yg)
        y_pixscl = pysiaf.utils.polynomial.dpdy(self.t['Sci2IdlY'], xg, yg)
        return(x_pixscl,y_pixscl)

    def save_pxscale_mesh(self,x_pixscl,y_pixscl,outfilebasename=None,overwrite=False):
        if outfilebasename is None:
            outfilebasename=re.sub('\.txt','',self.filename)
        
        outfilename=f'{outfilebasename}.x_pixscl.fits'
        if self.verbose: print(f'Saving {outfilename}')
        if os.path.isfile(outfilename) and (not overwrite):
            raise RuntimeError(f'{outfilename} already exists, exiting. Use overwrite option to overwrite!')
        fits.writeto(outfilename,x_pixscl,overwrite=True)

        outfilename=f'{outfilebasename}.y_pixscl.fits'
        if self.verbose: print(f'Saving {outfilename}')
        if os.path.isfile(outfilename) and (not overwrite):
            raise RuntimeError(f'{outfilename} already exists, exiting. Use overwrite option to overwrite!')
        fits.writeto(outfilename,y_pixscl,overwrite=True)
        
        outfilename=f'{outfilebasename}.xy_pixscl.fits'
        if self.verbose: print(f'Saving {outfilename}')
        if os.path.isfile(outfilename) and (not overwrite):
            raise RuntimeError(f'{outfilename} already exists, exiting. Use overwrite option to overwrite!')
        fits.writeto(outfilename,y_pixscl*y_pixscl,overwrite=True)
        
if __name__ == '__main__':

    pixelscale = calc_pixelscale()
    parser = pixelscale.define_options()
    parser.add_argument('--overwrite', default=False, action='store_true', help='overwrite files if they exist.')

    args = parser.parse_args()
    
    pixelscale.verbose=args.verbose
    
    filenames=[]
    for filepattern in args.coeff_filepatterns:
        newfilenames = glob.glob(filepattern)
        if len(newfilenames)==0:
            raise RuntimeError(f'Could not find any files that match {filepattern}')
        filenames.extend(newfilenames)            

    filenames.sort()

    for filename in filenames:
        if pixelscale.verbose: print(f'#### {filename}')
        pixelscale.load_coeff_file(filename)
        (x_pixscl,y_pixscl) = pixelscale.calc_pxscale_mesh()
        pixelscale.save_pxscale_mesh(x_pixscl,y_pixscl,overwrite=args.overwrite)
