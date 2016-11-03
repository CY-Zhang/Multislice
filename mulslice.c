/*
          *** mulslice.c ***

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

  rewritten in C 26-july-1995 E. Kirkland
  in working form 5-Oct-1995 ejk
  switch to TIFF file format 5-may-1996 ejk
  remove slice expansion (i.e. put in atompot) 4-aug-1996 ejk
  add dx,dy to output parameters 19-feb-1997 ejk
  remove commas from input formats 11-july-1997 ejk
  fixed small problem with anti-aliasing 5-jan-1998 ejk
  added astigmatism in pc mode and inc. beam tilt 28-jan-1998 ejk
  fixed format of error message 16-feb-1998 ejk
  update memory allocation routines 13-nov-1999 ejk
  change void main() to int main() for better portability
             22-jan-2000 ejk
  change data type of nxl,nyl to long32 for compatibility with new
      tiffsubs.c  17-jul-2007 ejk
  convert to GPL 3-jul-2008 ejk

    PIX   = final pix for partial coherence mode
    WAVE  = current specimen transmitted wavefunction
    TRANS = single slice transmission function
    PROPX,PROPY  = propagator function factored as two 1D arrays
    TEMP  = scratch array

  This program calls  subroutines from slicelib.c,
        tiffsubs.c, fft2dc.c

  ax,by  = unit cell size in x,y
  BW     = Antialiasing bandwidth limit factor
  acmin  = minimum illumination angle
  acmax  = maximum illumination angle
  Cs     = spherical aberration coefficient
  DF0    = defocus (mean value)
  SIGMAF = defocus spread (standard deviation)
  DFDELT = sampling interval for defocus integration
  
  this file is formatted for a TAB size of 8 characters 
  
*/

#include <stdio.h>      /* ANSI C libraries */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "fft2dc.h"     /* FFT routines */
#include "slicelib.h"   /* misc. routines for multislice */
#include "tiffsubs.h"   /* file I/O routines in TIFF format */

#define BW (2.0F/3.0F)  /* bandwidth limit */
#define ABERR 1.0e-4    /* max error for a,b */

#define NSMAX 1000      /* max number of slices */
#define NPARAM  64      /* number of parameters */
#define NLMAX   52      /* maximum number of layers */

#define NCMAX 256       /* max characters in file names */
#define NCINMAX  500    /* max number of characters in stacking spec */

int main()
{
        int lstart=0, lpartl=0, lbeams=0;
        int ix,iy, nx,ny, nx2, ny2, ixmid,iymid, i, islice, nslice,
                nacx,nacy, iqx, iqy, ilayer, nlayer, npix,
                nslic0, nslic1, ndf, idf, nbout, ib, nbits[3], samples;
        int *layer, *hbeam, *kbeam;
        long nbeams, nillum;
        long32 nxl, nyl, nxl2, nyl2;   /*  data type for tiffsubs */

        float *kx, *ky, *kx2, *ky2, *xpos, *ypos, *cz, *param, *sparam;
        float k2, k2max, scale, v0, vz, mm0, wavlen, rx, ry,
                ax, by, pi, rmin, rmax, aimin, aimax, x,y,
                ax2,by2, rx2,ry2, cztot, ctiltx, ctilty, tctx, tcty,
                acmin, acmax, Cs, df, df0, sigmaf, dfdelt, aobj,
                qx, qy, qy2, q2, q2min, q2max, sumdf, pdf, k2maxo,
                dfa2, dfa2phi, dfa3, dfa3phi, btiltx, btilty;

        float tr, ti, wr, wi;

        char **filein, fileout[NCMAX], filestart[NCMAX], filebeam[NCMAX];
        char datetime[20];

        char *cin, *cin2;
        double sum, timer, xdf, chi, chi1, chi2, phi, t;

        float **waver, **wavei, ***transr, ***transi,
                **propxr, **propxi, **propyr, **propyi,
                **tempr, **tempi, **pix;
                
        FILE *fp1;

        /*  set up symbolic mapping
                this must be the same as in parlay */

        char cname[] = 
                "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

        /*  echo version date */

        printf("mulslice version dated 3-jul-2008 ejk\n");
        printf("Copyright (C) 1998, 2008 Earl J. Kirkland\n" );
        printf( "This program is provided AS-IS with ABSOLUTELY NO WARRANTY\n "
            " under the GNU general public license\n\n" );

        printf("perform traditional multislice calculation\n\n");

        pi = (float) (4.0 * atan( 1.0 ));
        param = (float*) malloc1D( NPARAM, sizeof(float), "param" );
        sparam = (float*) malloc1D( NPARAM, sizeof(float), "sparam" );

        /*  read in the layer stacking sequence and parse it
            multiple line continuation is signified with a  trailing '-'
            a trailing '/echo' displays the results of parsing
        */

        cin2 = cin = (char*) malloc( NCINMAX * sizeof( char ) );
        if( cin == NULL ) {
                printf("cannot allocate stacking sequence storage.\n");
                exit(0);
        }
        for(i=0; i<NCINMAX; i++) cin[i] = 0;

        printf("Type in the stacking sequence :\n");
        do {
           scanf("%s", cin2 );
        }while( ( (cin2=strchr(cin,'-')) != NULL  )
          && ( strlen(cin) < (NCINMAX-80) ) );

        layer = (int*) malloc1D( NSMAX, sizeof(int), "layer" );

        if( parlay( cin, layer, NSMAX, NLMAX, &nslice, 1)
           < 0 ) exit( 0 );

        /*  Find total number of layers  */

        nlayer = 0;
        for( i=0; i<nslice; i++) 
                if( layer[i] > nlayer ) nlayer = layer[i];

        nlayer += 1;

        /*  Get input file name etc. */

        printf("\nType in the name of %d atomic potential layers :\n\n",
                nlayer);

        filein = (char**) malloc2D( nlayer, NCMAX, sizeof(char), "filein" );
        for( i=0; i<nlayer; i++) {
           printf("Name of file with input atomic potential %c :\n",
                         cname[i]);
           scanf("%s", filein[i] );
        }

/*  get more file names etc. */

        printf("Name of file to get binary output of multislice result:\n");
        scanf("%s", fileout );

        lpartl = askYN("Do you want to include partial coherence");

        if( lpartl == 1 ) {
                printf("Illumination angle min, max in mrad:\n");
                scanf("%f %f", &acmin, &acmax);
                acmin  = acmin  * 0.001F;
                acmax  = acmax  * 0.001F;
                printf("Spherical aberration (in mm.):\n");
                scanf("%g", &Cs);
                Cs = Cs * 1.0e7F;
                printf("Defocus, mean, standard deviation, and"
                       " sampling size (in Ang.)=\n");
                scanf("%f %f %f", &df0, &sigmaf, &dfdelt);
                printf("Objective aperture (in mrad) =\n");
                scanf("%f", &aobj);
                aobj = aobj * 0.001F;
                printf( "Magnitude and angle of 2-fold astigmatism"
                        " (in Ang. and degrees):\n");
                scanf( "%f %f", &dfa2, &dfa2phi);
                        dfa2phi = dfa2phi * pi /180.0F;
                printf( "Magnitude and angle of 3-fold astigmatism"
                        " (in Ang. and degrees):\n");
                scanf( "%f %f", &dfa3, &dfa3phi);
                dfa3phi = dfa3phi * pi /180.0F;
                lstart = 0;
        } else {
                printf("NOTE, the program image must also be run.\n");
                lstart = askYN("Do you want to start from previous result");
        }

        if ( lstart == 1 ) {
                printf("Name of file to start from:\n");
                scanf("%s", filestart);
        } else {
                printf("Incident beam energy in kev:\n");
                scanf("%g", &v0);
        }

        printf("Crystal tilt x,y in mrad.:\n");
        scanf("%f %f", &ctiltx, &ctilty);
        ctiltx = ctiltx * 0.001F;
        ctilty = ctilty * 0.001F;

        if( lpartl == 0 ) {
                printf("Incident beam tilt x,y in mrad.:\n");
                scanf("%f %f", &btiltx, &btilty);
                btiltx = btiltx * 0.001F;
                btilty = btilty * 0.001F;

                lbeams = askYN("Do you want to record the (real,imag) value\n"
                        " of selected beams vs. thickness");
                if( lbeams == 1 ) {
                        printf("Name of file for beams info:\n");
                        scanf("%s", filebeam );
                        printf("Number of beams:\n");
                        scanf("%d", &nbout);
                        if( nbout<1 ) nbout = 1;
                        hbeam = (int*) malloc1D( nbout, sizeof(int), "hbeam" );
                        kbeam = (int*) malloc1D( nbout, sizeof(int), "kbeam" );
                        for( ib=0; ib<nbout; ib++) {
                                printf("Beam %d, h,k=\n", ib+1);
                                scanf("%d %d", &hbeam[ib], &kbeam[ib] );
                        }
                }

        }
        timer = cputim();

/*  get starting value of transmitted wavefunction if required
   (this can only be used in coherent mode)
    remember to save params for final output pix  */

        if ( lstart == 1 ) {
                if( topenFloat( filestart ) != 1 ) {
                        printf("Cannot open input file: %s .\n", filestart ); 
                        exit( 0 );
                }
                tsize( &nxl, &nyl, nbits, &samples );
                nx = (int) nxl;
                ny = (int) nyl;

                waver = (float**) malloc2D( nx, ny, sizeof(float), "waver" );
                nx = nx/2;
                wavei = waver + nx;
                if( treadFloatPix( waver, nxl, nyl, &npix, datetime, sparam )
                         != 1 ) {
                        printf("Cannot read input file %s.\n", filestart);
                        exit(0);
                }
                tclose();
                if( npix != 2 ) {
                   printf("Input file %s must be complex, can't continue.\n",
                        filestart );
                   exit( 0 );
                }
                if( (nx != powerof2(nx)) || (ny != powerof2(ny)) ) {
                        printf("Nx=%d, Ny=%d must be a power of 2\n"
                           "problem in starting image, try again.\n", nx, ny);
                        exit( 0 );
                }

                ax = sparam[pDX] * nx;
                by = sparam[pDY] * ny;
                v0     = sparam[pENERGY];
                nslic0 = (int) sparam[pNSLICES];
                printf("Starting pix range %g to %g real\n"
                       "                   %g to %g imag\n",
                       sparam[pRMIN], sparam[pRMAX], sparam[pIMIN],
                       sparam[pIMAX] );
                printf("Beam voltage = %g kV\n", v0);
                printf("Old crystal tilt x,y = %g, %g mrad\n",
                       1000.*sparam[pXCTILT],
                        1000.*sparam[pYCTILT]);

        } else {

          nslic0 = 0;
        }

/*  calculate relativistic factor and electron wavelength */

        mm0 = 1.0F + v0/511.0F;
        wavlen = (float) wavelength( v0 );
        printf("Wavelength = %f Angstroms\n", wavlen );

/*  read in atomic potential and specimen parameters
    and calculate specimen transmission function
    for a single slice in transr,i */

        cz = (float*) malloc1D( nlayer, sizeof(float), "cz" );
        transr = (float***) malloc( nlayer * sizeof(float**) );
        transi = (float***) malloc( nlayer * sizeof(float**) );
        if( (transr==NULL) || (transi==NULL) ) {
                printf("Cannot allocate memory for transr,i.\n");
                exit(0);
        }

        for( ilayer=0;  ilayer<nlayer; ilayer++ ) {

                if( topenFloat( filein[ilayer] ) != 1 ) {
                        printf("Cannot open file %s.\n", filein[ilayer] );
                        exit( 0 );
                }
                tsize( &nxl2, &nyl2, nbits, &samples );
                nx2 = (int) nxl2;
                ny2 = (int) nyl2;
                transr[ilayer] = (float**) malloc2D( 2*nx2, ny2,
                                sizeof(float), "transr" );
                transi[ilayer] = transr[ilayer] + nx2;
                if( treadFloatPix( transr[ilayer], nxl2, nyl2, &npix,
                        datetime, param ) != 1 ) {
                        printf("Cannot read input file %s.\n", filein[ilayer]);
                        exit(0);
                }
                tclose();
                if( npix != 1 ) {
                        printf("Input potential file %s is not real.\n",
                                filein[ilayer] );
                        exit( 0 );
                }
 
                cz[ilayer] = param[pC];
                printf("layer %c, cz = %f\n", cname[ilayer],
                        cz[ilayer]);

                if ( ( lstart==1 ) || ( ilayer != 0) ) {
                        if ( (nx!=nx2) || (ny!=ny2) ) {
                                printf("pix size incompatible.\n");
                                printf("old size = %d, %d\n", nx,ny);
                                printf("new size = %d, %d\n", nx2, ny2 );
                                printf("layer = %1c\n", cname[ilayer]);
                                exit( 0 );
                        }
                        ax2 = param[pDX] * nx;
                        by2 = param[pDY] * ny;
                        if( ( fabs( ax-ax2 ) > fabs(ABERR*ax) ) ||
                            ( fabs( by-by2 ) > fabs(ABERR*by) ) ) {
                                printf("incompatible lattice constants\n");
                                printf("potential    a,b,c = %g, %g, %g\n",
                                       ax2,by2,cz[ilayer] );
                                printf("starting pix a,b,c = %g, %g\n",
                                       ax,by);
                                printf("   layer = %1c\n", cname[ilayer] );
                                exit( 0 );
                        }
                } else {
                        nx = nx2;
                        ny = ny2;
                        ax = param[pDX] * nx;
                        by = param[pDY] * ny;
                        if( (nx != powerof2( nxl2)) ||
                            (ny != powerof2( nyl2))) {
                                printf("Nx, Ny must be a power"
                                       "of 2, try again.\n");
                                printf("layer = %c, nx=%d, ny=%d\n",
                                         cname[ilayer], nxl, nyl  );
                                exit( 0 );
                        }
                }

                scale = wavlen * mm0;
                for( iy=0; iy<ny; iy++) {
                        for( ix=0; ix<nx; ix++) {
                                vz= scale * transr[ilayer][ix][iy];
                                transr[ilayer][ix][iy] = (float) cos(vz);
                                transi[ilayer][ix][iy] = (float) sin(vz);
                        }
                }

        }  /* end for(ilayer=... */

        printf("Size in pixels Nx x Ny= %d x %d = %d beams\n",
               nx,ny, nx*ny);
        printf("Lattice constant a = %12.4f, b = %12.4f\n", ax,by);

/*  calculate the total specimen thickness and echo */

        cztot = 0.0F;
        for( islice=0; islice<nslice; islice++) 
                cztot += cz[ layer[islice] ];
        printf("Total specimen thickness = %g Angstroms\n", cztot);
        
/*  calculate spatial frequencies and positions for future use */

        rx = 1.0F/ax;
        rx2= rx*rx;
        ry = 1.0F/by;
        ry2= ry*ry;
        ixmid = nx/2;
        iymid = ny/2;

        kx   = (float*) malloc1D( nx, sizeof(float), "kx" );
        kx2  = (float*) malloc1D( nx, sizeof(float), "kx2" );
        xpos = (float*) malloc1D( nx, sizeof(float), "xpos" );
        freqn( kx, kx2, xpos, nx, ax );

        ky   = (float*) malloc1D( ny, sizeof(float), "ky" );
        ky2  = (float*) malloc1D( ny, sizeof(float), "ky2" );
        ypos = (float*) malloc1D( ny, sizeof(float), "ypos" );
        freqn( ky, ky2, ypos, ny, by );

        if( lstart == 0 ) {
                waver = (float**) malloc2D( 2*nx, ny, sizeof(float), "waver" );
                wavei = waver + nx;
                qx = btiltx / wavlen;   /* add incident beam tilt */
                qy = btilty / wavlen;
                for( ix=0; ix<nx; ix++)
                for( iy=0; iy<ny; iy++) {
                        t = 2.0*pi*( qx*xpos[ix] + qy*ypos[iy] );
                        waver[ix][iy] = (float) cos( t );
                        wavei[ix][iy] = (float) sin( t );
                }
        }

/*  calculate propagator function, and bandwidth limit the transmission
   function for anti-aliasing.
   set to zero outside sampling circle
*/
        k2max = nx/(2.0F*ax);
        tctx = ny/(2.0F*by);
        if( tctx < k2max ) k2max = tctx;
        k2max = BW * k2max;
        printf("Bandwidth limited to a real space resolution of %f Angstroms\n",
                                         1.0F/k2max);
        printf("   (= %.2f mrad)  for symmetrical anti-aliasing.\n",
                 wavlen*k2max*1000.0F);
        k2max = k2max*k2max;

        tctx = (float) (2.0 * tan(ctiltx));
        tcty = (float) (2.0 * tan(ctilty));
        
        propxr = (float**) malloc2D( nlayer, nx, sizeof(float), "propxr" );
        propxi = (float**) malloc2D( nlayer, nx, sizeof(float), "propxi" );
        propyr = (float**) malloc2D( nlayer, ny, sizeof(float), "propyr" );
        propyi = (float**) malloc2D( nlayer, ny, sizeof(float), "propyi" );

        for( ilayer=0; ilayer<nlayer; ilayer++) {

                scale = pi * cz[ilayer];
                nbeams = 0;

                for( ix=0; ix<nx; ix++) {
                        t = scale * ( kx2[ix]*wavlen - kx[ix]*tctx );
                        propxr[ilayer][ix] = (float)  cos(t);
                        propxi[ilayer][ix] = (float) -sin(t);
                }
                for( iy=0; iy<ny; iy++) {
                        t = scale * ( ky2[iy]*wavlen - ky[iy]*tcty );
                        propyr[ilayer][iy] = (float)  cos(t);
                        propyi[ilayer][iy] = (float) -sin(t);
                }

                fft2d( transr[ilayer], transi[ilayer], nx, ny, +1);
                for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
                        k2= ky2[iy] + kx2[ix];
                        if (k2 < k2max) nbeams++;
                        else transr[ilayer][ix][iy] =
                                         transi[ilayer][ix][iy] = 0.0F;
                }
                fft2d( transr[ilayer], transi[ilayer], nx, ny, -1);

        } /* end for(ilayer... */

        printf("Number of symmetrical non-aliasing beams = %ld\n", nbeams);

/*  iterate the multislice algorithm proper

   NOTE: zero freq. is in the bottom left corner and
     expande into all other corners - not in the center
     this is required for the FFT - don't waste time rearranging

-------- partial coherence method ----------------------------
   force the integrals to include the origin and to be symmetric
   about the origin and to have the same periodic boundary
   conditions as the sampling grid
*/
        if( lpartl == 1 ) {

                printf("Illumination angle sampling (in mrad) = %f, %f\n\n",
                        1000.*rx*wavlen, 1000.*ry*wavlen);

                pix = (float**) malloc2D( nx, ny, sizeof(float), "pix" );
                for( ix=0; ix<nx; ix++)
                        for( iy=0; iy<ny; iy++) pix[ix][iy] = 0.0F;
                
                tempr = (float**) malloc2D( nx, ny, sizeof(float), "tempr" );
                tempi = (float**) malloc2D( nx, ny, sizeof(float), "tempi" );

                scale = 1.0F / ( ((float)nx) * ((float)ny) );
                ndf = (int) ( ( 2.5F * sigmaf ) / dfdelt );

                nacx = (int) ( ( acmax / ( wavlen * rx ) ) + 1.5F );
                nacy = (int) ( ( acmax / ( wavlen * ry ) ) + 1.5F );

                q2max = acmax / wavlen;
                q2max = q2max*q2max;

                q2min = acmin / wavlen;
                q2min = q2min*q2min;

                k2maxo = aobj / wavlen;
                k2maxo = k2maxo*k2maxo;

                chi1 = pi * wavlen;
                chi2 = 0.5 * Cs * wavlen *wavlen;
                nillum = 0;

                /*  integrate over the illumination angles */

                for( iqy= -nacy; iqy<=nacy; iqy++) {
                   qy = iqy * ry;
                   qy2 = qy * qy;

                   for( iqx= -nacx; iqx<=nacx; iqx++) {
                        qx = iqx * rx;
                        q2 = qx*qx + qy2;

                        if( (q2 <= q2max) && (q2 >= q2min) ) {
                           nillum += 1;
                           for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
                                t = 2.0*pi*( qx*xpos[ix] + qy*ypos[iy] );
                                waver[ix][iy] = (float) cos(t);
                                wavei[ix][iy] = (float) sin(t);
                           }
                           for( islice=0; islice<nslice; islice++) {
                                ilayer = layer[islice];
                                transmit( waver,  wavei, 
                                        transr[ilayer], transi[ilayer],
                                         nx, ny );
                                fft2d( waver, wavei, nx, ny, +1);
                                propagate( waver,  wavei, 
                                        propxr[ilayer], propxi[ilayer],
                                        propyr[ilayer], propyi[ilayer],
                                        kx2,  ky2,  k2max, nx, ny );
                                fft2d( waver, wavei, nx, ny, -1);
                           }

                           sum = 0.0;
                           for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
                                sum += waver[ix][iy]*waver[ix][iy]
                                        + wavei[ix][iy]*wavei[ix][iy];
                           sum = sum * scale;

                           printf("Illum. angle = %7.3f, %7.3f mrad",
                                1000.*qx*wavlen, 1000.*qy*wavlen);
                           printf(", total intensity= %f\n", sum );

                        /*  integrate over +/- 2.5 sigma of defocus */

                           fft2d( waver, wavei, nx, ny, +1);
                           sumdf = 0.0F;

                           for( idf= -ndf; idf<=ndf; idf++) {
                                df = df0 + idf*dfdelt;

                                for( ix=0; ix<nx; ix++)
                                for( iy=0; iy<ny; iy++) {
                                   k2 = kx2[ix] + ky2[iy];
                                   if( k2 <= k2maxo ) {
                                        phi = atan2( ky[iy], kx[ix] );
                                        chi = chi1*k2* ( chi2*k2 - df 
                                           + dfa2*sin( 2.0*(phi-dfa2phi) ) 
                                           + 2.0F*dfa3*wavlen*sqrt(k2)*
                                           sin( 3.0*(phi-dfa3phi) )/3.0 );
                                        tr = (float)  cos(chi);
                                        ti = (float) -sin(chi);
                                        wr = waver[ix][iy];
                                        wi = wavei[ix][iy];
                                        tempr[ix][iy] = wr*tr - wi*ti;
                                        tempi[ix][iy] = wr*ti + wi*tr;
                                   } else {
                                        tempr[ix][iy] = 0.0F;
                                        tempi[ix][iy] = 0.0F;
                                   }
                                }

                                fft2d( tempr, tempi, nx, ny, -1);

                                xdf = (double) ( (df - df0) /sigmaf );
                                pdf = (float) exp( -0.5F * xdf*xdf );
                                sumdf += pdf;

                                for( ix=0; ix<nx; ix++)
                                for( iy=0; iy<ny; iy++) {
                                        x = tempr[ix][iy];
                                        y = tempi[ix][iy];
                                        pix[ix][iy] += pdf* ( x*x + y*y );
                                }

                           }/* end for(idf..) */
                        }/* end if( q2...) */

                   } /* end for( iqx..) */
                } /* end for( iqy..) */

                printf("Total number of illumination angle = %ld\n",
                                nillum);
                printf("Total number of defocus values = %d\n", 2*ndf+1);
                scale = 1.0F / ( ((float)nillum) * sumdf );
                rmin  = pix[0][0] * scale;
                rmax  = rmin;
                aimin = 0.0F;
                aimax = 0.0F;

                for( ix=0; ix<nx; ix++)
                for( iy=0; iy<ny; iy++) {
                        pix[ix][iy] = pix[ix][iy] * scale;
                        if( pix[ix][iy] < rmin ) rmin = pix[ix][iy];
                        if( pix[ix][iy] > rmax ) rmax = pix[ix][iy];
                }

/* ------------------- coherent method -------------------- 
        remember that waver[][] was initialized above */

        } else {
                if( lbeams == 1 ) {
                        fp1 = fopen( filebeam, "w" );
                        if( fp1==NULL) {
                                printf("can't open file %s\n", filebeam);
                                exit(0);
                        }
                        fprintf( fp1, " (h,k) = " );
                        for(ib=0; ib<nbout; ib++)
                                fprintf(fp1," (%d,%d)", hbeam[ib], kbeam[ib]);
                        fprintf( fp1, "\n" );
                        fprintf( fp1, "nslice, (real,imag) (real,imag) ...\n\n");
                        for( ib=0; ib<nbout; ib++) {
                                if( hbeam[ib] < 0 ) hbeam[ib] = nx + hbeam[ib];
                                if( kbeam[ib] < 0 ) kbeam[ib] = ny + kbeam[ib];
                                if( hbeam[ib] < 0 ) hbeam[ib] = 0;
                                if( kbeam[ib] < 0 ) kbeam[ib] = 0;
                                if( hbeam[ib] > nx-1 ) hbeam[ib] = nx-1;
                                if( kbeam[ib] > ny-1 ) kbeam[ib] = ny-1;
                        }
                }
        
                nslic1 = nslic0;
                scale = 1.0F / ( ((float)nx) * ((float)ny) );

                for( islice=0; islice<nslice; islice++ ) {

                   ilayer =layer[islice];
                   transmit( waver,  wavei, 
                        transr[ilayer], transi[ilayer], nx, ny );

                   fft2d( waver, wavei, nx, ny, +1);
                   /* remember: prop must be here to anti-alias */
                   propagate( waver,  wavei, 
                                propxr[ilayer], propxi[ilayer],
                                propyr[ilayer], propyi[ilayer],
                                kx2,  ky2,  k2max, nx, ny );
                   if( lbeams == 1 )  {
                        fprintf( fp1, "%5d", nslic0+1);
                        for( ib=0; ib<nbout; ib++) 
                                fprintf(fp1, "%10.6f %10.6f",
                                   scale*waver[hbeam[ib]][kbeam[ib]],
                                   scale*wavei[hbeam[ib]][kbeam[ib]]);
                        fprintf( fp1, "\n");
                   }
                   fft2d( waver, wavei, nx, ny, -1);

                   sum = 0.0;
                   for( ix=0; ix<nx; ix++)
                   for( iy=0; iy<ny; iy++)
                        sum += waver[ix][iy]*waver[ix][iy] +
                                wavei[ix][iy]*wavei[ix][iy];
                   sum = sum * scale;

                   nslic0 +=  1;
                   printf("slice %4d, layer = %c, integrated intensity = %f\n",
                        nslic0, cname[ilayer], sum );

                } /* end for(islice...) */
        
                rmin  = waver[0][0];
                rmax  = rmin;
                aimin = wavei[0][0];
                aimax = aimin;

                for( ix=0; ix<nx; ix++)
                for( iy=0; iy<ny; iy++) {
                        x = waver[ix][iy];
                        y = wavei[ix][iy];
                        if( x < rmin ) rmin = x;
                        if( x > rmax ) rmax = x;
                        if( y < aimin ) aimin = y;
                        if( y > aimax ) aimax = y;
                }

        } /* end else .. coherent section */

/*  output results and find min and max to echo
    remember that complex pix are stored in the file in FORTRAN
                order for compatibility */
        if( lstart == 1 )
                for( ix=0; ix<NPARAM; ix++ ) param[ix] = sparam[ix];
        else
                for( ix=0; ix<NPARAM; ix++ ) param[ix] = 0.0F;
        param[pRMAX]  = rmax;
        param[pIMAX]  = aimax;
        param[pRMIN]  = rmin;
        param[pIMIN]  = aimin;
        param[pXCTILT]  = ctiltx;
        param[pYCTILT] = ctilty;
        param[pENERGY] = v0;
        param[pDX] = (float) ( ax/((float)nx) );
        param[pDY] = (float) ( by/((float)ny) );
        param[pWAVEL] = wavlen;
        param[pNSLICES] = (float) nslic0;
        if ( lpartl == 1 ) {
          param[pDEFOCUS] = df0;
          param[pOAPERT] = aobj;
          param[pCS] = Cs;
          param[pCAPERT] = acmax;
          param[pDDF] = sigmaf;
        }
        if ( lpartl == 1 )
                tcreateFloatPixFile( fileout, pix, (long) nx,
                         (long) ny, 1, param );
        else 
                tcreateFloatPixFile( fileout, waver, (long) (2*nx),
                         (long) ny, 2, param );
        printf( "pix range %g to %g real,\n"
                "          %g to %g imag\n",  rmin,rmax,aimin,aimax );

        printf("Total CPU time = %f sec.\n", cputim()-timer );

        return 0;

} /* end main() */
