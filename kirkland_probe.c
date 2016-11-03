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

        chi1 = pi * wavlen;
        chi2 = 0.5 * Cs3 * 1.0e7*wavlen*wavlen;
        chi3 = Cs5 * 1.0e7 * wavlen*wavlen*wavlen*wavlen /3.0;

        df = defocus in Ang.
        wavlen = wavelength (in Ang.)
        kx, ky = spatial freq (in 1/Ang)
        k2 = kx*kx + ky*ky
        dx2p,dy2p = 2*pi*offset (in Ang)
        dfa2, dfa2phi = magnitude and angle of 2nd order astigmatism
        dfa3, dfa3phi = magnitude and angle of 3rd order astigmatism

        - should premultiply dfa3 a little more for efficiency
*/
double chi( double chi1, double chi2, double chi3, double df, double wavlen,
           double kx, double ky, double k2, double dx2p, double dy2p,
           double dfa2, double dfa2phi, double dfa3, double dfa3phi)
{
    double c, phi;

        phi = atan2( ky, kx );

        c = chi1*k2* ( (chi2 + chi3*k2)*k2 - df 
                + dfa2*sin( 2.0*(phi-dfa2phi) ) 
                + 2.0*dfa3*wavlen*sqrt(k2)* sin( 3.0*(phi-dfa3phi) )/3.0 )
                        - ( (dx2p*kx) + (dy2p*ky) );
        return( c );
}

int main()
{
        char fileout[NCMAX];
        int ix, iy, nx, ny, ixmid, iymid, i, ismoth, npixels;
        float rmin, rmax, aimin, aimax;
        float *param, **pixr, **pixi;
        double  kx, ky, ky2, k2, k2max, keV, wavlen, ax, by, rx, ry,
                rx2, ry2, pi, dx, dy, dx2p, dy2p, scale, pixel,
                Cs3, Cs5, df, chi1, chi2, chi3, sum, chi0, time,
                apert, dfa2, dfa2phi, dfa3, dfa3phi;

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

        printf("Probe parameters, V0(kv), Cs3(mm), Cs5(mm),"
               " df(Angstroms), apert(mrad):\n");
        scanf("%lg %lg %lg %lg %lg",
              &keV, &Cs3, &Cs5, &df, &apert );

        printf( "Magnitude and angle of 2-fold astigmatism"
                " (in Ang. and degrees):\n");
        scanf( "%lf %lf", &dfa2, &dfa2phi);
        dfa2phi = dfa2phi * pi /180.0;
        printf( "Magnitude and angle of 3-fold astigmatism"
                " (in Ang. and degrees):\n");
        scanf( "%lf %lf", &dfa3, &dfa3phi);
        dfa3phi = dfa3phi * pi /180.0;

        printf("Type 1 for smooth aperture:\n");
        scanf("%d", &ismoth );

        printf("Probe position in Angstroms:\n");
        scanf("%lf %lf", &dx, &dy );
        dx2p = 2.0*pi*dx;
        dy2p = 2.0*pi*dy;

/*  Calculate misc constants  */

        time = cputim( );
        
        rx  = 1.0/ax;
        rx2 = rx * rx;
        ry  = 1.0/by;
        ry2 = ry * ry;
        
        ixmid = nx/2;
        iymid = ny/2;
        
        wavlen = wavelength( keV );
        printf("electron wavelength = %g Angstroms\n", wavlen);

        chi1 = pi * wavlen;
        chi2 = 0.5 * Cs3*1.0e7 *wavlen*wavlen;
        chi3 = Cs5*1.0e7 * wavlen*wavlen*wavlen*wavlen /3.0;

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

        printf( "dfa2= %f\n", dfa2 );

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
                    chi0 = chi( chi1, chi2, chi3, df, wavlen, kx, ky,
                            k2, dx2p, dy2p, dfa2, dfa2phi, dfa3, dfa3phi);
                   pixr[ix][iy]= (float) ( 0.5 * cos(chi0));
                   pixi[ix][iy]= (float) (-0.5 * sin(chi0));
                   printf("smooth by 0.5 at ix=%d, iy=%d\n", ix, iy );
                } else if ( k2 <= k2max ) {
                    chi0 = chi( chi1, chi2, chi3, df, wavlen, kx, ky,
                            k2, dx2p, dy2p, dfa2, dfa2phi, dfa3, dfa3phi);
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
        param[pDEFOCUS]= (float) df;
        param[pDX]= (float) (ax / nx);
        param[pDY]= (float) (by / ny);
        param[pENERGY]= (float) keV;
        param[pWAVEL]= (float) ( sqrt(k2max) * wavlen);
        param[pCS]= (float) Cs3;
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

        return 0;

}  /* end main() */
