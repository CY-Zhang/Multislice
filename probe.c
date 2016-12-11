/*          *** probe.c ***

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

        ANSI-C version

        Calculate a focused probe wavefunction in real space

        this file is formatted for  a tab size of 8 characters

        rewritten in C 6-dec-1995 ejk
        fixed sign error in aberration function 1-mar-1997 ejk
        removed commas from keyboard input 3-oct-1997 ejk
        updated memory allocation routines 20-nov-1999 ejk
        change void main() to int main() for better portability
             22-jan-2000 ejk
        add Cs5=5th order aberration and astigmatism  19-jul-2005 ejk
        small cosmetic changes 18-jul-2007 ejk
        convert to GPL 3-jul-2008 ejk
	modify chi() to support CEOS aberration measurements
*/

#include <stdio.h>      /*  ANSI-C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "fft2dc.h"     /*  FFT's */
#include "slicelib.h"   /* define parameter offsets */
#include "tiffsubs.h"   /* file I/O libraries */

#define NCMAX   132     /* characters per line to read featom.tab */
#define NPARAM  64      /* number of parameters */


/*------------------------ chi() ---------------------*/
/*
        calculate the aberration function chi
        with 3rd and 5th order Cs plus astigmatism plus offset

        put in a subroutine so I only have to get it right once

  input values:

  Aberration parameters are:
	C1, defocus, in nm
	A1, 2-fold astigmatism, in nm
	A2, 3-fold astigmatism, in nm
	B2, axial coma, in nm
	C3, primary spherical aberration, in um
	A3, 4-fold astigmasitm, in um
	S3, star aberration, in um
	A4, 5-fold astigmatism, in um
	D4, 3-lobe aberration, in um
	B4, axial coma, in um
	C5, 5th order spherical aberration, in mm
	A5, 5th order spherical aberration, in mm
  
In the aberrations array, the first column contains the aberration coefficient, the
second column contains the rotation angle for that aberration, the third column
contains the sine of the angle, and the fourth column contains the cosine.

        wl = wavelength (in Ang.)
        kx, ky = spatial freq (in 1/Ang)
        dx,dy = offset (in Ang)

*/
double chi( double **aber, double wl, double kx, double ky, double dx, double dy)
{
    double c, phi;
    double kxr, kyr;

    kxr = kx;
    kyr = ky;

    c = 0.0;
      
    /* Evaluate the phase shifts from the various aberrations */
    c += -1.0*(dx*kx + dy*ky);                                                              /* probe center shift */
    c += 0.5*wl*aber[0][0]*(kx*kx + ky*ky);                                                     /* C1, defocus */

    kxr = kx*aber[1][3] + ky*aber[1][2];  kyr = ky*aber[1][3] - kx*aber[1][2];
    c += 0.5*wl*aber[1][0]*(kxr*kxr - kyr*kyr);                                                 /* A1, 2-fold astigmatism */
    kxr = kx*aber[2][3] + ky*aber[2][2];  kyr = ky*aber[2][3] - kx*aber[2][2];
    c += (1.0/3.0)*wl*wl*aber[2][0]*(pow(kxr,3) - 3*kxr*kyr*kyr);                               /* A2, 3-fold astigmatism */
    kxr = kx*aber[3][3] + ky*aber[3][2];  kyr = ky*aber[3][3] - kx*aber[3][2];
    c += wl*wl*aber[3][0]*(pow(kxr, 3) + kxr*kyr*kyr);                                          /* B2, axial coma */

    c += 0.25*pow(wl,3)*aber[4][0]*(pow(kx, 4) + 2*pow(kx,2)*pow(ky,2) + pow(ky, 4));           /* C3, primary spherical aberration */

    kxr = kx*aber[5][3] + ky*aber[5][2];  kyr = ky*aber[5][3] - kx*aber[5][2];
    c += 0.25*pow(wl,3)*aber[5][0]*(pow(kxr, 4) - 6*pow(kxr,2)*pow(kyr, 2) + pow(kyr,4));	/* A3, 4-fold astigmasitm */
    kxr = kx*aber[6][3] + ky*aber[6][2];  kyr = ky*aber[6][3] - kx*aber[6][2];
    c += pow(wl,3)*aber[6][0]*(pow(kxr,4) - pow(kyr,4));                                        /* S3, star aberration */ 
    kxr = kx*aber[7][3] + ky*aber[7][2];  kyr = ky*aber[7][3] - kx*aber[7][2];
    c += 0.2*pow(wl,4)*aber[7][0]*(pow(kxr, 5) - 10*pow(kxr,3)*pow(kyr,2) + 5*kxr*pow(kyr,4));	/* A4, 5-fold astigmatism */
    kxr = kx*aber[8][3] + ky*aber[8][2];  kyr = ky*aber[8][3] - kx*aber[8][2];
    c += pow(wl,4)*aber[8][0]*(pow(kxr, 5) - 2*pow(kxr,3)*pow(kyr,2) - 3*kxr*pow(kyr, 4));	/* D4, 3-lobe aberration */
    kxr = kx*aber[9][3] + ky*aber[9][2];  kyr = ky*aber[9][3] - kx*aber[9][2];
    c += pow(wl,4)*aber[9][0]*(pow(kxr,5) + 2*pow(kxr,3)*pow(kyr,2) + kxr*pow(kyr,4));	        /* B4, axial coma */

    c += (1.0/6.0)*pow(wl,5)*aber[10][0]*(pow(kx,6) + 3*pow(kx,4)*pow(ky,2) 
					  + 3*pow(kx,2)*pow(ky,4) + pow(ky,6));                 /* C5, 5th order spherical aberration */

    kxr = kx*aber[11][3] + ky*aber[11][2];  kyr = ky*aber[11][3] - kx*aber[11][2];
    c += (1.0/6.0)*pow(wl,5)*aber[11][0]*(pow(kxr,6) - 15*pow(kxr,4)*pow(kyr,2) 
					  + 15*pow(kxr,2)*pow(kyr,4) - pow(kyr,6));             /* A5, 5th order spherical aberration */

    c *= 2.0*3.1415927;

    /*phi = atan2( ky, kx );
    
    //        c = chi1*k2* ( (chi2 + chi3*k2)*k2 - df 
    //               + dfa2*sin( 2.0*(phi-dfa2phi) ) 
    //                + 2.0*dfa3*wavlen*sqrt(k2)* sin( 3.0*(phi-dfa3phi) )/3.0 )
    //                        - ( (dx2p*kx) + (dy2p*ky) ); */
    return( c );
}

int main()
{
        char fileout[NCMAX];
        int ix, iy, nx, ny, ixmid, iymid, i, ismoth, npixels;
        float rmin, rmax, aimin, aimax;
        float *param, **pixr, **pixi;
        double  kx, ky, ky2, k2, k2max, keV, wavlen, ax, by, rx, ry,
                rx2, ry2, pi, dx, dy, scale, pixel,
                sum, chi0, time, apert;
	double **aber;

	aber = (double**)malloc2D(12, 4, sizeof(double), "Aberrations array");

/*  Echo version date etc.  */
        printf( "probe version dated 3-jul-2008 ejk\n");
        printf("Copyright (C) 1998, 2008 Earl J. Kirkland\n" );
        printf( "This program is provided AS-IS with ABSOLUTELY NO WARRANTY\n "
            " under the GNU general public license\n\n" );

        printf( "calculate focused probe wave function\n\n");

        pi = 4.0 * atan( 1.0 );

/*  Get desired image size, parameters etc. */

        printf("Name of file to get focused probe wave function:\n");
        scanf("%s", fileout );

        printf("Desired size of output image in pixels Nx,Ny:\n");
        scanf("%d %d", &nx, &ny );

        if( (nx != powerof2(nx)) || (ny != powerof2(ny)) ) {
                printf("Nx=%d, Ny=%d must be a power of 2,\n"
                        "try again.\n", nx, ny);
                exit( 0 );
        }

        printf("Size of output image in Angstroms ax,by:\n");
        scanf("%lf %lf", &ax, &by );

        printf("Probe parameters: V0(kv),  aperture (mrad):\n");
        scanf("%lg %lg", &keV, &apert );

	printf("First order aberrations: C1 (nm), A1 (nm, deg): \n");
	scanf("%lg %lg %lg", &aber[0][0], &aber[1][0], &aber[1][1]);
	aber[0][1] = 0.0; /* no angle for C1 */
	aber[0][0] *= 10.0; /* nm to Angstroms */
	aber[1][0] *= 10.0;

	printf("Second order aberrations: A2 (nm, deg), B2 (nm, deg): \n");
	scanf("%lg %lg %lg %lg", &aber[2][0], &aber[2][1], &aber[3][0], &aber[3][1]);
	aber[2][0] *= 10.0;
	aber[3][0] *= 10.0;

	printf("Third order aberrations: C3 (um), A3 (um, deg), S3 (nm, deg): \n");
	scanf("%lg %lg %lg %lg %lg", &aber[4][0], &aber[5][0], &aber[5][1], &aber[6][0], &aber[6][1]);
	aber[4][1] = 0.0; /* no angle for C3 */
	aber[4][0] *= 1.0e4; /* um to Angstroms */
	aber[5][0] *= 1.0e4;
	aber[6][0] *= 1.0e4;

	printf("Fourth order aberrations: A4 (um, deg), D4 (um, deg), B4 (um, deg): \n");
	scanf("%lg %lg %lg %lg %lg %lg", &aber[7][0], &aber[7][1], &aber[8][0], &aber[8][1], 
	      &aber[9][0], &aber[9][1]);
	aber[7][0] *= 1.0e4; /* um to Angstroms */
	aber[8][0] *= 1.0e4;
	aber[9][0] *= 1.0e4;

	printf("Fifth order aberrations: C5 (mm), A5 (mm, deg): \n");
	scanf("%lg %lg %lg", &aber[10][0], &aber[11][0], &aber[11][1]);
	aber[10][1] = 0.0; /* no angle for C5 */
	aber[10][0] *= 1.0e7; /* mm to Angstroms */
	aber[11][0] *= 1.0e7;

	/* Convert angles from degress to radians, pre-compute the angle sin and cosine */
	for(i=0; i<12; i++) {
	  aber[i][1] *= 180.0/pi;
	  aber[i][2] = sin(aber[i][1]);
	  aber[i][3] = cos(aber[i][1]);
	}

        printf("Type 1 for smooth aperture:\n");
        scanf("%d", &ismoth );

        printf("Probe position in Angstroms:\n");
        scanf("%lf %lf", &dx, &dy );

/*  Calculate misc constants  */

	printf("Aberrations:\n");
	printf("C1 = %g Ang, A1 = %g Ang, %g rad\n", aber[0][0], aber[1][0], aber[1][1]);


        time = cputim( );
        
        rx  = 1.0/ax;
        rx2 = rx * rx;
        ry  = 1.0/by;
        ry2 = ry * ry;
        
        ixmid = nx/2;
        iymid = ny/2;
        
        wavlen = wavelength( keV );
        printf("electron wavelength = %g Angstroms\n", wavlen);

        k2max = apert*0.001/wavlen;
        k2max = k2max * k2max;

        param = (float*) malloc1D( NPARAM, sizeof(float), "probe-param" );
        for( i=0; i<NPARAM; i++) param[i] = 0.0F;
        pixr = (float**) malloc2D( 2*nx, ny, sizeof(float), "pixr" );
        pixi = pixr + nx;

/*   Calculate MTF 
        NOTE zero freq. is in the bottom left corner and
        expands into all other corners - not in the center
        this is required for FFT

        PIXEL = diagonal width of pixel squared
        if a pixel is on the aperture boundary give it a weight
        of 1/2 otherwise 1 or 0
*/
        pixel = ( rx2 + ry2 );
        npixels = 0;


        for( iy=0; iy<ny; iy++) {
           ky = (double) iy;
           if( iy > iymid ) ky = (double) (iy-ny);
           ky = ky*ry;
           ky2 = ky*ky;
           for( ix=0; ix<nx; ix++) {
                kx = (double) ix;
                if( ix > ixmid ) kx = (double) (ix-nx);
                kx = kx*rx;
                k2 = kx*kx + ky2;
                if ( ( ismoth != 0) && 
                        ( fabs(k2-k2max) <= pixel) ) {
                   /*  old chi= chi1*k2* (chi2*k2-df) -
                        2.0F*pi*( (dx*kx/ax) + (dy*ky/by) ); */
		  chi0 = chi( aber, wavlen, kx, ky, dx, dy );
                   pixr[ix][iy]= (float) ( 0.5 * cos(chi0));
                   pixi[ix][iy]= (float) (-0.5 * sin(chi0));
                   printf("smooth by 0.5 at ix=%d, iy=%d\n", ix, iy );p
                } else if ( k2 <= k2max ) {
		  chi0 = chi( aber, wavlen, kx, ky, dx, dy );
                   pixr[ix][iy]= (float)  cos(chi0);
                   pixi[ix][iy]= (float) -sin(chi0);
                   npixels++;
                } else {
                   pixr[ix][iy] = pixi[ix][iy] = 0.0F;
                }
           }
        }

        printf("there were %d pixels inside the aperture\n", npixels );

        fft2d( pixr, pixi, nx, ny, -1);

/*  Normalize probe intensity to unity  */

        sum = 0.0;
        for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) 
                sum +=  pixr[ix][iy]*pixr[ix][iy]
                                + pixi[ix][iy]*pixi[ix][iy];

        scale = 1.0 / sum;
        scale = scale * ((double)nx) * ((double)ny);
        scale = (double) sqrt( scale );

        for( ix=0; ix<nx; ix++) 
           for( iy=0; iy<ny; iy++) {
                pixr[ix][iy] *= (float) scale;
                pixi[ix][iy] *= (float) scale;
        }

/*  Output results and find min and max to echo
    remember that complex pix are stored in the file in FORTRAN
                order for compatibility
*/

        rmin = pixr[0][0];
        rmax = rmin;
        aimin = pixi[0][0];
        aimax = aimin;
        for( iy=0; iy<ny; iy++) {
           for( ix=0; ix<nx; ix++) {
                if( pixr[ix][iy] < rmin ) rmin = pixr[ix][iy];
                if( pixr[ix][iy] > rmax ) rmax = pixr[ix][iy];
                if( pixi[ix][iy] < aimin ) aimin = pixi[ix][iy];
                if( pixi[ix][iy] > aimax ) aimax = pixi[ix][iy];
           }
        }

        param[pRMAX] = rmax;
        param[pIMAX] = aimax;
        param[pRMIN] = rmin;
        param[pIMIN] = aimin;
        param[pDEFOCUS]= (float) aber[0][0];
        param[pDX]= (float) (ax / nx);
        param[pDY]= (float) (by / ny);
        param[pENERGY]= (float) keV;
        param[pWAVEL]= (float) ( sqrt(k2max) * wavlen);
        param[pCS]= (float) aber[4][0];
        param[27]= (float) dx;
        param[28]= (float) dy;

        if( tcreateFloatPixFile( fileout, pixr, (long) (2*nx), (long) ny,
                         2, param ) != 1 )
                printf( "probe cannot write an output file.\n");

        printf( "Pix range %15.7g to %15.7g real,\n"
                "      and %15.7g to %15.7g imaginary\n",
                rmin, rmax, aimin, aimax );
        time = cputim() - time;
        printf("\nCPU time = %f sec\n", time );

	/* free2D(aber); */

        return 0;

}  /* end main() */
