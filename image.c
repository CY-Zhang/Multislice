/*          *** image.c ***

------------------------------------------------------------------------
Copyright 1998,2008 Earl J. Kirkland

This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.

---------------------- NO WARRANTY ------------------
THIS PROGRAM IS PROVIDED AS-IS WITH ABSOLUTELY NO WARRANTY
OR GUARANTEE OF ANY KIND, EITHER EXPRESSED OR IMPLIED,
INCLUDING BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
MERCHANABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
IN NO EVENT SHALL THE AUTHOR BE LIABLE
FOR DAMAGES RESULTING FROM THE USE OR INABILITY TO USE THIS
PROGRAM (INCLUDING BUT NOT LIMITED TO LOSS OF DATA OR DATA
BEING RENDERED INACCURATE OR LOSSES SUSTAINED BY YOU OR
THIRD PARTIES OR A FAILURE OF THE PROGRAM TO OPERATE WITH
ANY OTHER PROGRAM). 

------------------------------------------------------------------------

  ANSI C and TIFF version

  rewritten in ANSI-C 29-sep-1995 E. Kirkland
  convert to TIFF file format 20-may-1996 ejk
  remove commas in input format 11-july-1997 ejk
  fixed sign error 23-jan-1998 ejk
  add 2-fold and 3-fold astigmatism to coherent mode 24-jan-1998 ejk
  add "\n" to one printf() format 28-jan-1998 ejk
  changed scaling in diffraction pattern mode 2-feb-1998 ejk
  fixed sqrt() in diff patt scaling 19-feb-1998 ejk
  changed variable PI in tcross() to pi because Linux gcc
      predefines PI  9-feb-1999 ejk
  update memory allocation routines 13-nov-1999 ejk
  change void main() to int main() for better portability
             22-jan-2000 ejk
  change data type of nxl,nyl to long32 for compatibility with new
      tiffsubs.c  17-jul-2007 ejk
  convert to GPL 4-jul-2008 ejk

  This file is formatted for a tab size of 8 characters

*/
#include <stdio.h>      /* ANSI C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "fft2dc.h"     /* FFT routines */
#include "slicelib.h"   /* misc. routines for multislice */
#include "tiffsubs.h"   /* file I/O routines in TIFF format */

#define NCMAX   132     /* max number of characters per line */
#define NPARAM  64      /* number of parameters */

#define MCPIX   0       /* coherent image mode */
#define MPCPIX  1       /* partial coherent mode */
#define MDIFFP  2       /* diffraction mode */

int main()
{
        char filein[NCMAX], fileout[NCMAX];
        char datetime[20];
        const char version[] = "3-jul-2008 (ejk)";

        int ix, iy, nx, ny, ixmid, iymid, npix,
                itens, mode, nsum, nbits[3], samples;
        int lcenter=0, lapert=0;
        long32 nxl, nyl;   /*  tiffsubs.h data type = 32 bit integer */

        float kx,ky, ky2, k2, k2max, v0, wavlen, scale, pixc,
                ax, by, cz, rx, ry, pi, rmin, rmax, aimin, aimax,
                Cs,df,  alpha0, ddf, objlx, objly,
                tr, ti, wr, wi, dfa2, dfa2phi, dfa3, dfa3phi;
        float *param;
        float  **pixr, **pixi, **pix2r, **pix2i;
        double sum, time, chi, chi1, chi2, phi, clog;

        void tcross( float **rpixin, float **ipixin, float **rpixo,
               float **ipixo, int nx, int ny, float *p );

/*  echo version date */

        printf( "image version dated %s\n", version );
        printf("Copyright (C) 1998, 2008 Earl J. Kirkland\n" );
        printf( "This program is provided AS-IS with ABSOLUTELY NO WARRANTY\n "
            " under the GNU general public license\n\n" );

        printf( "calculate TEM images with defocus\n\n");

        pi = (float) (4.0 * atan( 1.0 ));

/*  get input file name */

        printf("Name of file with input multislice result:\n");
        scanf("%s", filein );

/*  open input file and check sizes etc. */

        if( topenFloat( filein ) != 1 ) {
                printf( "Cannot open input file %s\n", filein );
                exit( 0 );
        }
        tsize( &nxl, &nyl, nbits, &samples );
        nx = (int) nxl;
        nx = nx/2;              /* should be complex */
        ny = (int) nyl;
        ixmid = nx/2;
        iymid = ny/2;

        if( (nx != powerof2(nx) ) || ( ny != powerof2(ny) )) {
                printf("Nx, Ny = %d, %d\n", nx,ny );
                printf("   not a power of 2, try again.\n");
                tclose();
                exit( 0 );
        }

/*  ask mode */

        printf("Type %d for coherent real space image,\n"
                "  or %d for partially coherent real space image,\n"
                "  or %d for diffraction pattern output:\n",
                MCPIX, MPCPIX, MDIFFP );
        scanf("%d", &mode );

/*  get imaging parameters if required */

        Cs = 0.0F;
        df = 0.0F;
        k2max = 0.0F;
        objlx = 0.0F;
        objly = 0.0F;

        if( (mode == MCPIX) || (mode == MPCPIX) ) {
                printf("Name of file to get defocused output:\n");
                scanf("%s", fileout );

                printf("Spherical aberration in mm.:\n");
                scanf("%f", &Cs);
                Cs = Cs * 1.0e7F;

                printf("Defocus in Angstroms:\n");
                scanf("%f", &df );

                printf("Objective aperture size in mrad:\n");
                scanf("%f", &k2max );

                if( mode == MPCPIX ) {
                        printf("Illumination semi-angle in mrad:\n");
                        scanf("%f", &alpha0 );
                        alpha0 = alpha0 * 0.001F;
                        printf("Defocus spread in Angstroms:\n");
                        scanf("%f", &ddf );
                } else if( mode == MCPIX ) {
                        printf( "Magnitude and angle of two-fold astigmatism"
                                " (in Angst. and degrees):\n");
                        scanf( "%f %f", &dfa2, &dfa2phi);
                        dfa2phi = dfa2phi * pi /180.0F;
                        printf( "Magnitude and angle of three-fold astigmatism"
                                " (in Angst. and degrees):\n");
                        scanf( "%f %f", &dfa3, &dfa3phi);
                        dfa3phi = dfa3phi * pi /180.0F;
                        printf("Objective lens and aperture"
                                " center x,y in mrad\n"
                                " (i.e. non-zero for dark field):\n");
                        scanf("%f %f", &objlx, &objly);
                }

/*  get diffraction pattern parameters if required  */

        } else {
                printf("Name of file to get diffraction pattern:\n");
                scanf("%s", fileout );

                lcenter = askYN("Do you want to include central beam");
                lapert = askYN("Do you want to impose the aperture");
                if ( lapert == 1 ) {
                        printf("Aperture size in mrad:\n");
                        scanf("%f", &k2max );
                        printf("Objective lens and aperture"
                                " center x,y in mrad\n"
                                "  (i.e. non-zero for dark field):\n");
                        scanf("%f %f", &objlx, &objly );
                }

                printf( "Type 0 for linear scale,\n"
                        "  or 1 to do logarithmic intensity scale:\n"
                        "  or 2 to do log(1+c*pixel) scale:\n");
                scanf("%d", &itens );
                if( itens == 2 ) {
                        printf( "Type scaling constant c:\n");
                        scanf( "%lf", &clog );
                }
        }

/*  read in specimen parameters and multislice result  */

        time = cputim();        /* get elapsed time for comparison */

        param = (float*) malloc1D( NPARAM, sizeof(float), "param" );
        for( ix=0; ix<NPARAM; ix++) param[ix] = 0.0F;
        pixr = (float**) malloc2D( 2*nx, ny, sizeof(float), "pixr" );
        pixi = pixr + nx;

        if( treadFloatPix( pixr, nxl, nyl, &npix, datetime, param )
                 != 1 ) {
                printf("Cannot read input file %s.\n", filein);
                exit(0);
        }
        tclose();
        if( npix != 2 ) {
                printf("Input file %s must be complex, can't continue.\n",
                        filein );
                exit( 0 );
        }

        ax = param[pDX] * ((float) nx);
        by = param[pDY] * ((float) ny);
        if( param[pDY] <= 0.0F ) by = param[pDX] * ((float) ny);
        cz = param[pC];
        rx  = 1.0F / ax;
        ry  = 1.0F / by;

        v0 = param[pENERGY];
        wavlen = (float) wavelength( v0 );
        printf("Starting pix energy = %8.2f keV\n", v0);

        rmax  = param[pRMAX];
        aimax = param[pIMAX];
        rmin  = param[pRMIN];
        aimin = param[pIMIN];
        printf("Starting pix range %g %g real\n"
                "                   %g %g imag\n", rmin,rmax,aimin,aimax);

        param[pDEFOCUS] = df;
        param[pASTIG] = 0.0F;
        param[pTHETA] = 0.0F;
        param[pOAPERT] = k2max * 0.001F;
        param[pCS] = Cs;
        param[pWAVEL] = wavlen;

        k2max = k2max*0.001F/wavlen;
        k2max = k2max * k2max;

        objlx = objlx * 0.001F / wavlen;   /* convert to spatial freq. */
        objly = objly * 0.001F / wavlen;

/*   coherent real space image mode 

     Convolute multislice result with objective lens
        aberration function (coherent case)
        with shifted objective lens for dark field

   NOTE zero freq. is in the bottom left corner and
     expands into all other corners - not in the center
     this is required for FFT
*/

        if( mode == MCPIX ) {
           chi1  = pi * wavlen;
           chi2  = 0.5 * Cs * wavlen * wavlen;

           fft2d( pixr, pixi, nx, ny, +1);

                for( iy=0; iy<ny; iy++) {
                   ky = (float) iy;
                   if( iy > iymid ) ky = (float)(iy-ny);
                   ky = ky*ry - objly;
                   ky2 = ky*ky;
                   for( ix=0; ix<nx; ix++) {
                        kx = (float) ix;
                        if( ix > ixmid ) kx = (float) (ix-nx);
                        kx = kx*rx - objlx;
                        k2 = kx*kx + ky2;
                        if ( k2 < k2max ) {
                                phi = atan2( ky, kx );
                                chi = chi1*k2* ( chi2*k2 - df 
                                   + dfa2*sin( 2.0*(phi-dfa2phi) ) 
                                   + 2.0F*dfa3*wavlen*sqrt(k2)*
                                        sin( 3.0*(phi-dfa3phi) )/3.0 );
                                tr = (float)  cos( chi );
                                ti = (float) -sin( chi );
                                wr = pixr[ix][iy];
                                wi = pixi[ix][iy];
                                pixr[ix][iy] = tr*wr - ti*wi;
                                pixi[ix][iy] = tr*wi + ti*wr;
                        } else {
                                pixr[ix][iy] = pixi[ix][iy] = 0.0F;
                        }
                   }
                }

                fft2d( pixr, pixi, nx, ny, -1);

/*  calculate partially coherent image here  
    NOTE this assumes that the parameter offsets in txcoef()
        match those in slicelib.h !!!! (this is bad practice but....)
*/

        } else if( mode == MPCPIX) {

                param[pCAPERT] = alpha0;
                param[pDDF] = ddf;
                param[ 24 ] = 0.0F;

                pix2r = (float**) malloc2D( nx, ny, sizeof(float), "pix2r" );
                pix2i = (float**) malloc2D( nx, ny, sizeof(float), "pix2i" );
                fft2d( pixr, pixi, nx, ny, +1);
                invert2D( pixr, nx, ny );
                invert2D( pixi, nx, ny );
                tcross( pixr, pixi, pix2r, pix2i, nx, ny, param );

                for( ix=0; ix<nx; ix++)
                for( iy=0; iy<ny; iy++) {
                        pixr[ix][iy] = pix2r[ix][iy];
                        pixi[ix][iy] = pix2i[ix][iy];
                }

                invert2D( pixr, nx, ny );
                invert2D( pixi, nx, ny );
                fft2d( pixr, pixi, nx, ny, -1);

/*  do diffraction pattern here  */

        } else {

                fft2d ( pixr, pixi, nx, ny, +1);
                invert2D( pixr, nx, ny);
                invert2D( pixi, nx, ny);
                if( lcenter == 0 ) 
                        pixr[ixmid][iymid] = pixi[ixmid][iymid] = 0.0F;

                if ( lapert ) {
                        for( iy=0; iy<ny; iy++) {
                                ky = (iy-iymid)*ry - objly;
                                ky2 = ky*ky;
                                for( ix=0; ix<nx; ix++) {
                                        kx = (ix-ixmid)*rx - objlx;
                                        k2 = kx*kx + ky2;
                                        if ( k2 > k2max )
                                                pixr[ix][iy] = pixi[ix][iy] = 0.0F;
                                } /* end for ix */
                        } /* end for iy... */
                }

        }

/*  Output results and find min and max to echo */

        scale = 1.0F / ( ((float)nx) * ((float)ny) );

        sum = 0.0;
        nsum = 0;
        for( iy=0; iy<ny; iy++) {
           for( ix=0; ix<nx; ix++) {
              tr = pixr[ix][iy];
              ti = pixi[ix][iy];
              if( mode == MPCPIX ) pixc = tr;
              else  pixc = tr*tr + ti*ti;
              if ( (mode == MDIFFP) &&  (itens == 1)) {
                if( pixc > 1.e-30F)  pixc = (float) log( (double) pixc );
                        else pixc = -30.0F;
              } else if ( (mode == MDIFFP) && (itens == 2)) {
                 pixc = (float) log( 1.0 + clog*sqrt((double)pixc) );
              } 

              if( (ix == 0) && (iy == 0) ) {
                rmin = pixc;
                rmax = rmin;
              } else if( (ix != ixmid) && (iy != iymid) ) {
                if( pixc < rmin ) rmin = pixc;
                if( pixc > rmax ) rmax = pixc;
              }
              if( (ix>(3*nx)/8) && (ix<(5*nx)/8) &&
                  (iy>(3*ny)/8) && (iy<(5*ny)/8) ) {
                 sum = sum + pixc;
                 nsum += 1;
              }
              pixr[ix][iy] = pixc;
           }  /* end for ix... */
        } /* end for iy... */

        param[pRMAX] = rmax;
        param[pIMAX] = 0.0F;
        param[pRMIN] = rmin;
        param[pIMIN] = 0.0F;
        if ( mode == MDIFFP ) {
                param[pRMIN] = (float) (0.05*rmin + 0.95*sum/nsum);
                param[pDX] = 1.0F / (nx*param[pDX]);
                param[pDY] = 1.0F / (ny*param[pDY]);
        }
        tcreateFloatPixFile( fileout, pixr, (long) nx,
                 (long) ny, 1, param );

        printf("Pix range %f to %f\n",  param[pRMIN], rmax );

        time = cputim() - time;
        printf("Elapsed time = %f sec\n", time );

        return 0;

} /* end main() */

/*------------------- tcross ------------------------------*/
/*
 Subroutines to do exact nonlinear partial coherence transfer as in
      [1]  M.A. O'Keefe, 37th EMSA (1979) p.556
      [2]  K. Ishizuka, Ultramicroscopy 5 (1980) p.55-65

  both input and output are in Fourier space with zero frequency
  in the center (NX/2 , NY/2)

 perform 2D weighted convolution with transmission cross coefficient
  TXCOEF to calculate partially coherent CTEM images
    NOTE: this version uses Friedels law

  rpixin, ipixin  : real,imag. input pix array dimensioned as [ix][iy]
                     range used = P(17) = aperture
  rpixo, ipixo    : real,imag. output pix array dimensioned as [ix][iy]
                     range output= 2.0*P(17) = 2 x aperture
  nx, ny : actual image size to work with
  p[43]  : parameter array ( dx,dy, defocus etc.) explained in TXCOEF

  started 14-NOV-85 on CIPRES2 Earl J. Kirkland
      started by modifying XFERFN.FTN (form RSX CIPRES)
  changed TXCOEF to be more general (i.e. include alignment and leading
      factors so that it is complete - little or no speed lose)
         28-NOV-85 EJK
  fixed typo (DT3 to D3T) and sign (-D4) in TXCOEF 17-jun-1986 ejk
  changed sign convention to be consistent with forward propagation
        of electrons (see Self et.al. UM 11 (1983) p.35-52 
          - only TXCOEF effected   EJK 2-JUL-86
  converted to C 8-nov-1995 ejk
  fixed sign error 23-jan-1998 ejk

*/

void tcross( float **rpixin, float **ipixin, float **rpixo, float **ipixo,
        int nx, int ny, float *p )
{
        int j, ix, iy, ixmid, iymid, ixmin, ixmax, iymin, iymax,
                ix2min, ix2max, iy0, iy2min, iy2max, ix1, ix2, iy1, iy2,
                *ixi, *ixf, ixi0, ixf0;
        float scale, k2p, k2pp, rx, ry, txcoefr, txcoefi,
                xr, xi, yr, yi, zr, zi, pi, *kxp, *kyp, *kxp2, *kyp2;

        void txcoef( float kxp, float kyp, float kxpp, float kypp,
                 float p[], float *txcoefr, float *txcoefi );


/*  get scratch arrays and init params */

        ixi = (int*) malloc1D( nx, sizeof(int), "ixi" );
        ixf = (int*) malloc1D( nx, sizeof(int), "ixf" );
        kxp = (float*) malloc1D( nx, sizeof(float), "kxp" );
        kyp = (float*) malloc1D( ny, sizeof(float), "kyp" );
        kxp2 = (float*) malloc1D( nx, sizeof(float), "kxp" );
        kyp2 = (float*) malloc1D( ny, sizeof(float), "kyp2" );

        pi = (float) (4.0 * atan( 1.0 ));
        xr = p[17]/p[19];
        p[31] = xr*xr;
        xi = p[19];
        p[32] = 0.5F*pi* p[18]* xi*xi*xi;
        p[33] = pi* p[19] * p[11];

        xr = pi*p[21]* p[22];
        p[34] = xr*xr;
        p[35] = p[18]* p[19]* p[19];
        xr = pi*p[21];
        p[36] = xr*xr;
        xr =  pi* p[19]* p[22];
        p[37] = xr*xr/4.0F;
        xr = pi*pi*p[21]*p[21]*p[22] ;
        p[38] = xr*xr;
        xr =  pi*p[21]*p[22] ;
        p[39] = pi*( xr*xr )*p[19];

/* initialize  */

        rx = p[14]* nx;
        ry = p[15]* ny;
        if( ry <= 0.0F) ry = rx;
        rx = 1.0F/rx;
        ry = 1.0F/ry;

        scale=1.0F/( ((float)nx) * ((float)ny) );

        ixmid = nx/2;
        iymid = ny/2;

        /* find range of convolution */

        j = (int) ( p[17] / (rx*p[19]) +1.0F );
        ixmin  = ixmid - j;
        ixmax  = ixmid + j;
        ix2min = ixmid - 2*j;
        ix2max = ixmid + 2*j;
        j = (int) ( p[17] / (ry*p[19]) +1.0F );
        iymin  = iymid - j;
        iymax  = iymid + j;
        iy2min = iymid - 2*j;
        iy2max = iymid + 2*j;
        
        if( ixmin < 0 ) ixmin = 0;
        if( ixmax > (nx-1) ) ixmax = (nx-1);
        if( ix2min < 0 ) ix2min = 0;
        if( ix2max > (nx-1) ) ix2max = (nx-1);
        if( iymin < 0 ) iymin = 0;
        if( iymax > (ny-1) ) iymax = (ny-1);
        if( iy2min < 0 ) iy2min = 0;
        if( iy2max > (ny-1) ) iy2max = (ny-1);

        for ( ix=ix2min; ix<=ix2max; ix++) {
           ix1= ixmin - ix + ixmid;
           if( ixmin > ix1 ) ixi[ix] = ixmin;  else  ixi[ix] = ix1;
           ix1= ixmax - ix + ixmid;
           if( ixmax < ix1 ) ixf[ix] = ixmax;   else  ixf[ix] = ix1;
        }

        for( ix=ix2min; ix<=ix2max; ix++) {
                kxp[ix] = (ix-ixmid) * rx;
                kxp2[ix]= kxp[ix] * kxp[ix];
        }
        for( iy=iy2min; iy<=iy2max; iy++) {
                kyp[iy] = (iy-iymid) * ry;
                kyp2[iy]= kyp[iy] * kyp[iy];
        }

        for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
                rpixo[ix][iy] = ipixo[ix][iy] = 0.0F;

/*  do convolution in bottom half plane 
        the origin is at (ixmid,iymid)
        (ix,iy) = output point
        (ix1,iy1) = input point from psi
        (ix2,iy2) = input point from psi*
*/

        for( iy=iy2min; iy<=iymid; iy++ )
           for( iy1=iymin; iy1<=iymax; iy1++ ) {
                iy2 = iy1 + (iy - iymid);
                if( (iy2>=iymin) && (iy2<=iymax)) {
                   for( ix=ix2min; ix<=ix2max; ix++ ) {
                        ixi0 = ixi[ix];
                        ixf0 = ixf[ix];
                        for( ix1=ixi0; ix1<=ixf0; ix1++ ) {
                           ix2= ix1 + (ix - ixmid);
                           k2p = kxp2[ix1] + kyp2[iy1];
                           k2pp= kxp2[ix2] + kyp2[iy2];
                           if ( (k2p<=p[31]) && (k2pp<=p[31]) ) {
                                txcoef( kxp[ix1],kyp[iy1],kxp[ix2],kyp[iy2],
                                   p, &txcoefr, &txcoefi);
                                xr = rpixin[ix1][iy1];
                                xi = ipixin[ix1][iy1];
                                yr =  rpixin[ix2][iy2];
                                yi = -ipixin[ix2][iy2];
                                zr = xr*yr - xi*yi;
                                zi = xr*yi + xi*yr;
                                rpixo[ix][iy] += zr*txcoefr - zi*txcoefi;
                                ipixo[ix][iy] += zr*txcoefi + zi*txcoefr;
                           }  /* end if( ( k2p... */
                        }  /* end for( ix1=... */
                   }  /* end for( ix=... */
                }  /* end if( (iy2... */
        }  /* end for iy1... */

/*  scale result */

        for( ix= ix2min; ix<=ix2max; ix++) 
        for( iy= iy2min; iy<=iymid; iy++) {
                rpixo[ix][iy] *= scale;
                ipixo[ix][iy] *= scale;
        }

/*  Invoke Friedel's law to get top half plane  */

        iy0 = iymid+1;
        for( iy=iy0; iy<=iy2max; iy++) {
                iy1= iymid - (iy-iymid);
                for( ix=ix2min; ix<=ix2max; ix++) {
                        ix1= ixmid - (ix-ixmid);
                        rpixo[ix][iy]=  rpixo[ix1][iy1];
                        ipixo[ix][iy]= -ipixo[ix1][iy1];
                }
           rpixo[0][iy] = ipixo[0][iy] = 0.0F;  /* ??? */
        }

        free( ixi );
        free( ixf );
        free( kxp );
        free( kyp );
        free( kxp2 );
        free( kyp2 );

        return;
}     /*   end tcross() */

/*------------------------ txcoef -------------------------*/
/*
   The cross spectral transfer function as in 
     M.A. O'Keefe, 37th EMSA (1979) p.556 EJK 28-SEP-83

   switched to my derivation 5-DEC-83 EJK
   changed sign convention to be consistent with forward propagation
     of electrons (see Self et.al. UM 11 (1983) p.35-52 
       EJK 2-JUL-86
   converted to C  6-nov-1995 ejk
   fixed sign error 23-jan-1998 ejk

  the following are the array index definitions
  (the order is historical - don't ask why! )

   P(11) = defocus
   P(12) = defocus astigmatism
   P(13) = angle of astigmatism
   P(14) = dx in A pixel dimension
   P(15) = dy in A pixel dimension
   P(17) = aperture in radians
   P(18) = Cs in A
   P(19) = wavelength in A
   P(21) = illumination semiangle in rad
   P(22) = defocus spread in A
   P(24) = Debye Waller temp factor/4 in A
   P(31) = P(17)/P(19) **2 maximum k^2 value
 
   P(32) = PI*P(18)* (P(19)**3) /2        ; coherent part
   P(33) = PI* P(19) * P(11)
 
   P(34) = (PI* P(21)* P(22))**2
   P(35) = P(18)* P(19)* P(19)
   P(36) = (PI* P(21)) **2
   P(37) = (PI* P(19)* P(22))**2 /4.
   P(38) = PI^4 * P(21)^4 *P(22)^2
   P(39) = PI^3 * P(22)^2 * P(21)^2 *P(19)
      
     KP = kprime
     KPP= kprime + k
 
   txcoefr, txcoefi = returned real and imag. parts of function
                (i.e. C can't return a complex value )
*/

void txcoef( float kxp, float kyp, float kxpp, float kypp,
                 float p[], float *txcoefr, float *txcoefi )
{
        double kx, ky, k2 ,kp2, kpp2, chip, chipp,
                d1, d2, d3, d3t, d4, d4t, vx, vy, w, u, xr;

        kp2  = kxp*kxp + kyp*kyp;
        kpp2 = kxpp*kxpp + kypp*kypp;

        if( ( kp2 > p[31] ) || ( kpp2 > p[31] ) ) {
          *txcoefr = *txcoefi = 0.0F;
          return;
        }

        kx = kxpp - kxp;
        ky = kypp - kyp;
        k2 = kx*kx + ky*ky;

        /*      v = Wc1/(2*pi*lambda)
                w = Wc2/(-pi*lambda)
                u = 1 + pi^2 * beta^2 * delta0^2 k^2
        */
        vx = p[35]*( kp2*kxp - kpp2*kxpp ) + p[11]*kx;
        vy = p[35]*( kp2*kyp - kpp2*kypp ) + p[11]*ky;
        w  = kp2 - kpp2;
        u  = 1.0 + p[34]*k2;

        d1  = p[36] * (vx*vx+vy*vy);
        d2  = p[37] *w*w /u;
        d3t = vx*kx + vy*ky;
        d3  = p[38]*d3t*d3t/u;
        d4t = p[39]*w/u;
        d4  = d4t*d3t;

        chip  = kp2*  (p[32]*kp2  - p[33] );
        chipp = kpp2* (p[32]*kpp2 - p[33] );

        chip = -chip + chipp + d4;
        xr = exp( -d1 -d2 +d3 -p[24]*(kp2+kpp2) ) / sqrt(u);
        *txcoefr = (float) ( cos(chip) * xr );
        *txcoefi = (float) ( sin(chip) * xr );

        return;

}  /* end txcoef() */
