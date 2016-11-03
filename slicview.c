/*          *** slic(e)view.c ***

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

  ANSI-C an TIFF version

  Generate a 3D perspective view of a multislice
  specimen description with a hard sphere atom model.

  This file is formatted for a tab size of 8 characters.

  started from model3d.f, atompot.f and multislice.f
       12-July-1990 Earl J. Kirkland
  started conversion to C 5-Dec-1995 ejk
  converted to TIFF 30-jun-1996 ejk
  misc debugging 30-jun-1997 ejk
  add EPS postscript output 1-jul-1997 ejk
  add autoslic xyz format 1-feb-1998 ejk
  update memory allocation routines 20-nov-1999 ejk
  change void main() to int main() for better portability
             22-jan-2000 ejk
  change lpix to long32 to work with new tiffsubs 18-jul-2007 ejk
  convert to GPL 3-jul-2008 ejk

*/

#include <stdio.h>      /* ANSI-C routines */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "slicelib.h"   /* misc. routines for multislice */
#include "tiffsubs.h"   /* file I/O routines in TIFF format */

#define NCMAX   256     /* max number of characters per line */
#define NCINMAX 500     /* max number of characters in stacking spec */

#define NAMAX   20000   /* max number of atoms */
#define MAXLAY  26      /* max number of different types of layers  */
#define MAXSLI  1000    /* maximum number of slices  */

#define TIFFmode 1      /* output a TIFF image file */
#define EPSmode 2       /* output an encapsulated postscript file */
#define ATOMPOTinput 1  /* atompot input data format mode */
#define XYZinput 2      /* autoslic xzy input format mode */

/* define subroutines in this file */

void center( float *x, float ax );
void model3( float x[], float y[], float z[], float s[],
                int icolor[], int npts, float d, float size0,
                double rotat, double tilt,
                long32 **pix, int nxout, int nyout, float *scalepix,
                int mode, float EPSxsize, float EPSysize );
void iswap( int *i, int *j );
void fswap( float *a, float *b );
void sphere( double x0, double y0, double size, int icolor, 
                long32 **l2pix, int nxout, int nyout );

/*  encapsulated postscript routine and data */
int EPSheader( const char filname[], double xsize, double ysize );
int EPSsphere( double xp, double yp, double spsize );
int EPSend(  );
#define xpos0   1.0     /*  offset for postscript coord. */
#define ypos0   1.0     /*       = lower left corner in inches */
FILE *fpEPS;  /* for postscript file routines */

int main() 
{
        int nx, ny, ixc, iyc, i, nsym, jj,  k, islice,
                nlayer, natoms, mode, inputformat, ncx, ncy, ncz;

        int *layer, nslice, iz, *ncellx, *ncelly, *icolor,*Znum;
        char **filein, fileout[NCMAX], cline[NCMAX], *cin, *cin2,
                xyzfile[NCMAX];
        
        long32 **lpix;   /*  tiffsubs data type (32 bit integer) */
        float *x, *y, *z, *s, *cz, *symx1, *symx2, *symy1, *symy2,
                scale, ax, by, cz2, total, occ, x0, y0, cztot,
                *occ2, *wobble, d, size, EPSxsize, EPSysize;
        double rotat, tilt, pi, time;

        FILE *fp;

        /*  set up symbolic mapping
                this must be the same as in parlay */

        char cname[] = 
                "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

        /* ------- echo version date ------------------- */
        printf("sliceview C version dated 3-jul-2008 ejk\n");
        printf("Copyright (C) 1998, 2008 Earl J. Kirkland\n" );
        printf( "This program is provided AS-IS with ABSOLUTELY NO WARRANTY\n "
            " under the GNU general public license\n\n" );

        printf("Create a 3d perspective view of a "
                "multislice specimen\n\n");

        printf("Type %d for atompot input format,\n"
                "  or %d for XYZ autoslic input format:\n",
                ATOMPOTinput, XYZinput );
        scanf( "%d", &inputformat );

        /*  get output file name  */

        printf("\nName of output file to get "
                "3D perspective view of multislice specimen:\n");
        scanf( "%s", fileout );

        printf( "Type %d for TIFF output or %d for EPS output.\n",
                TIFFmode, EPSmode);
        scanf( "%d", &mode );

        if( mode == TIFFmode ) {
                printf("Image dimensions in (integer) pixels Nx Ny :\n");
                scanf( "%d %d", &nx, &ny );
        } else if( mode == EPSmode ) {
                printf("EPS size in inch. (real) xsize ysize :\n");
                scanf( "%f %f", &EPSxsize, &EPSysize );
        } else {
                printf( "Bad output mode, try again.\n" );
                exit( 0 );
        }

        /* remember that sphere size is relative to the full scale dimension */
        printf("Type viewing distance and sphere size:\n");
        scanf( "%f %f", &d, &size);

        printf("Type rotation, tilt angle in degrees:\n");
        scanf( "%lf %lf", &rotat, &tilt );
        pi = 4.0 * atan( 1.0 );
        rotat = pi * rotat/180.0;
        tilt  = pi * tilt /180.0;

/* ----------  atompot input format ---------------------- */

    if( inputformat == ATOMPOTinput ) {
        printf("maximum number of slices = %d\n", MAXSLI);   /*???*/
        /*  read in the layer stacking sequence and parse it
           multiple line continuation is signified with a 
                        trailing '-'
            a trailing '/echo' displays the results of parsing
        */
        cin2 = cin = (char*) malloc1D( NCINMAX, sizeof(char), "cin" );
        for(i=0; i<NCINMAX; i++) cin[i] = 0;

        printf("Type in the stacking sequence :\n");
        do {
           scanf("%s", cin2 );
        }while( ( (cin2=strchr(cin,'-')) != NULL  )
                && ( strlen(cin) < (NCINMAX-80) ) );

        layer = (int*) malloc1D( MAXSLI, sizeof(int), "layer" );
        if( (cin2=strstr( cin, "/echo") ) != NULL ) (*cin2) = '\0';

        if( parlay( cin, layer, MAXSLI, MAXLAY, &nslice, 1)
                < 0 ) exit( 0 );

        if( cin2 != NULL ) {
                for( i=0; i<nslice; i++) printf(" %d ", layer[i] );
                printf("\n");
        }

        /*  Find total number of layers  */
        nlayer = 0;
        for( i=0; i<nslice; i++) 
        if( layer[i] > nlayer ) nlayer = layer[i];

        nlayer += 1;
        /*  Get input file names for atompot data  */
        filein = (char**) malloc2D( nlayer, NCMAX, sizeof(char), "filein" );
        ncellx = (int*) malloc1D( nlayer, sizeof(int), "ncellx" );
        ncelly = (int*) malloc1D( nlayer, sizeof(int), "ncelly" );
        cz = (float*) malloc1D( nlayer, sizeof(float), "cz" );
        printf("Type in the names of %d atompot description files.\n",
                nlayer );
        for( i=0; i<nlayer; i++) {
                printf("Name of file with input crystal data for layer"
                        " %c :\n", cname[i]);
                scanf( "%s", filein[i] );

                printf( "replicate layer %c unit cell by ncellx,ncelly:\n",
                        cname[i] );
                scanf( "%d %d",  &ncellx[i], &ncelly[i] );
                if( ncellx[i] < 1 ) ncellx[i] = 1;
                if( ncelly[i] < 1 ) ncelly[i] = 1;

                /* test that file name is correct before continuing */
                fp = fopen( filein[i], "r" );
                if( fp == NULL ) {
                        perror(" can't open crystal data input file");
                        exit(0 );
                }
                fscanf( fp, "%f %f %f", &ax, &by, &cz[i] );
                fclose( fp );
        }

        time = cputim();  /* just for fun */

/*  read in specimen parameters and expand  */

        x = (float*) malloc1D( NAMAX, sizeof(float), "x" );
        y = (float*) malloc1D( NAMAX, sizeof(float), "y" );
        z = (float*) malloc1D( NAMAX, sizeof(float), "y" );
        s = (float*) malloc1D( NAMAX, sizeof(float), "s" );
        icolor = (int*) malloc1D( NAMAX, sizeof(int), "icolor" );

        symx1 = symx2 = symy1 = symy2 = NULL;

        total = 0.0F;
        natoms = 0;
        for( k=0; k<nlayer; k++) {

           fp = fopen( filein[k], "r");
           if( fp == NULL ) {
                perror(" can't open crystal data input file");
                exit(0 );
           }
           ReadLine( fp, cline, NCMAX, filein[k] );
           sscanf( cline, "%f %f %f", &ax, &by, &cz[k] );
           printf("Lattice constants = %f  %f  %f Angstroms\n",
                   ax, by, cz[k]);

           ReadLine( fp, cline, NCMAX, filein[k] );
           sscanf( cline, "%d", &nsym);
           if (nsym > 0) {
                symx1 = (float*) realloc( symx1, nsym * sizeof(float) );
                symx2 = (float*) realloc( symx2, nsym * sizeof(float) );
                symy1 = (float*) realloc( symy1, nsym * sizeof(float) );
                symy2 = (float*) realloc( symy2, nsym * sizeof(float) );
                if( (symx1==NULL) || (symx2==NULL) ||
                    (symy1==NULL) || (symy2==NULL) ) {
                        printf("cannot allocate symmetry storage\n");
                        exit( 0 );
                }
                for( i=0; i<nsym; i++)  {
                   ReadLine( fp, cline, NCMAX, filein[k] );
                   sscanf( cline, "%f %f %f %f",
                     &symx1[i], &symx2[i], &symy1[i], &symy2[i]);
                }
           }

           /*  remember that this must be expanded on-the-fly so it
                can't be the same as in atompot  */

L60:       ReadLine( fp, cline, NCMAX, filein[k] );
           sscanf( cline, "%d", &iz);
           if( (strlen(cline) >1) && (iz >= 1) ) {

L70:            ReadLine( fp, cline, NCMAX, filein[k] );
                sscanf( cline, "%f %f %f", &occ, &x0, &y0);
                if( (strlen(cline) >1) && (occ > 0.0) ) {

                /*  calculate the total z as a negative number for
                        compatibility with model3D()  */

                cztot = 0.0F;
                for( islice=0; islice<nslice; islice++) {
                   cztot = cztot - cz[ layer[islice] ];      
          
                   /*  put in the 'raw' coordinates  */
                   if( layer[islice] == k ) {
                        for( iyc=0; iyc<ncelly[k]; iyc++)
                        for( ixc=0; ixc<ncellx[k]; ixc++) {
                           if( natoms >= NAMAX ) {
                                printf("too many atoms\n");
                                printf("  maximum allowed = %d\n",NAMAX);
                                fclose( fp );
                                exit( 0 );
                           }
                           x[natoms] = (x0 + ixc) * ax;
                           y[natoms] = (y0 + iyc) * by;
                           z[natoms] = cztot;
                           s[natoms] = 1.0F;
                           icolor[natoms] = 1;
                           total = total + occ;
                           natoms = natoms + 1;

                           /*  Now apply all symmetry operations.  */

                           if( nsym > 0 ) {
                                for( jj=0; jj<nsym; jj++) {
                                   if( natoms >= NAMAX ) {
                                        printf("too many atoms\n");
                                        printf("  maximum allowed = %d\n", NAMAX);
                                        fclose( fp );
                                        exit( 0 );
                                   }
                                   x[natoms] = symx1[jj]*x0 + symx2[jj];
                                   center( &x[natoms], 1.0F );
                                   x[natoms] = ( x[natoms] + ixc ) * ax;

                                   y[natoms] = symy1[jj]*y0 + symy2[jj];
                                   center( &y[natoms], 1.0F );
                                   y[natoms] = ( y[natoms] + iyc ) * by;

                                   z[natoms] = cztot;
                                   s[natoms] = 1.0F;
                                   icolor[natoms] = 1;
                                   total = total + occ;
                                   natoms = natoms + 1;
                                }  /* end for(jj...) */
                           }  /* end if( nsym>0) */

                        } /* end for( ixc..) */
                   }  /* end if( layer[]] == k ) */

                }  /* end for( islice...) */

                goto L70;
             }  /*  end if( strlen()... */

             goto L60;
          }  /*  end if( strlen()... */

        } /* end for(k=1,nlayer) */

        printf("Total specimen thickness = %f\n", fabs(cztot));

/* -------- XYZ autoslic input format --------------- */

    } else if( inputformat == XYZinput ) {

        printf( "Name of file with specimen structure"
                " in autoslic xyz format:\n");
        scanf( "%s", xyzfile );
        printf("Replicate by Nx, Ny, Nz unit cells:\n");
        scanf( "%d %d %d", &ncx, &ncy, &ncz );

        time = cputim();  /* just for fun */

        natoms = ReadXYZcoord( xyzfile, ncx, ncy, ncz, &ax, &by, &cz2,
                &Znum, &x, &y, &z, &occ2, &wobble, cline, NCMAX );

        /* fill in required data */
        s = (float*) malloc1D( natoms, sizeof(float), "s" );
        icolor = (int*) malloc1D( natoms, sizeof(int), "icolor" );
        total = 0.0F;
        for( i=0; i<natoms; i++) {
                s[i] = 1.0F;
                icolor[i] = 1;
                total += occ2[i];
        }
    }
        printf("Total number of atoms = %d\n", natoms );
        printf("with a total occupancy of %f\n", total);

        /*  draw the perspective image  */

        if( mode == TIFFmode ) {
                /* use a long32** for compatibility with tiffsubs
                        although an unsigned char would work */
                lpix = (long32**) malloc2D( ny, nx, sizeof(long32), "lpix" );
                model3( x, y, z, s, icolor, natoms, d, size, rotat, tilt,
                        lpix, nx, ny, &scale, mode, EPSxsize, EPSysize );

                /*  Output TIFF results */
                if( tcreatePixFile( fileout, lpix, (long) nx, (long) ny,
                        0, 0, 8, 0, 0, 1.0, 1.0 ) != 1 )
                                printf("Cannot write output file.\n");
        } else if( mode == EPSmode ) {
                EPSheader( fileout, EPSxsize, EPSysize );
                nx = ny = 1;
                /* note really used */
                lpix = (long32**) malloc2D( ny, nx, sizeof(long), "lpix" );
                model3( x, y, z, s, icolor, natoms, d, size, rotat, tilt,
                        lpix, nx, ny, &scale, mode, EPSxsize, EPSysize );
                EPSend();
        }

        time = cputim() - time;
        printf("CPU time = %f sec\n", time );

}/*  end main() */


/*--------------------- center() -----------------*/
/*
  Subroutine to make sure that a coordinate is within the proper
  unit cell boundaries. If it isn't then translate by the
  appropriate number of lattice constants.
*/
void center( float *x, float ax )
{
        int i;

        if( *x < 0.0 ) {
                i = (int) fabs( (double) (*x / ax) );
                *x = *x + ( i + 1 ) * ax;
        } else if ( *x > ax ) 
                *x = (float) fmod( (double)*x , (double)ax );

        return;

}  /* end center() */


/*-------------------- model3() -----------------*/
/*
  Subroutine to draw 3D perspective of molecules
  into a 2D image with hidden surfaces.
  Each atom will be represented by a shaded solid sphere.

  The CRT is in xy plane with rotation and tilt
  instead of azimuthal and polar angles.

  NOTE: Color is not currently supported. However,
    leave the color code in for compatibility with
    old data files (from when there was color) and
    in case I ever figure out how to do it here.

   rewritten in C 22-feb-1996 E. Kirkland
   added EPS mode (kind of a kluge) 1-july-1997 ejk

  This program calls iswap(), fswap() and sphere()

* x[npts], y[npts], y[npts] = real valued x,y,z coordinates
                                of each atom
* y[npts]      = real size of  each atom
* icolor[npts] = integer color code for each atom
                        (not currently supported)
  npts  = integer number of atoms (or points)
  d     = real viewing distance
  size0 = real initial size of each sphere
          (relative to full size of the 2D image)
  rotat = real rotation angle (in radians)
  tilt  = real tilt angle (in radians)
* pix[][] = long int array to get 2d perspective image
                (remember that x,y indexes are reversed as
                pix[iy][ix] for tiffsubs convention)
  nxout, nyout   = portion of array to use
* scalepix       = real value to get the lateral scale of pix[][]
  mode = 1 for TIFF and 2 for EPS format

  NOTE 1: nxout,nyout and pix[][] are not used in EPS mode
        and EPSxsize and EPSysize are ignored in TIFF mode

  * = these arguments are changed

  NOTE 2: on exit X,Y,Z,S,ICOLOR will have the scaled and sorted
      coordinates (the original values are overwritten)
*/
void model3( float x[], float y[], float z[], float s[],
                int icolor[], int npts, float d, float size0,
                double rotat, double tilt,
                long32 **pix, int nxout, int nyout, float *scalepix,
                int mode, float EPSxsize, float EPSysize )
{
        int  ix, iy, i, j, k, m;
        float xmin, xmax, ymin, ymax,zmin,zmax, smax,
                cr, sr, st, ct, x0,y0,z0, scale, fs, fsx, fsy;
        double  size, pi, size2;

/*  initialize graphics image buffer  */

        if( mode == TIFFmode ) {
           for( ix=0; ix<nxout; ix++)
                for( iy=0; iy<nyout; iy++) pix[iy][ix] = 0;
        }

/*  Calculate misc constants  */

        cr = (float) cos( rotat );
        sr = (float) sin( rotat );
        ct = (float) cos( tilt );
        st = (float) sin( tilt );
        pi = 4.0 * atan( 1.0 );

/*  find range of atomic coordinates  */

        xmin = ymin = zmin = 1.0e30F;
        xmax = ymax = zmax = smax = -xmin;

        for( i=0; i<npts; i++) {
           if( s[i] > smax ) smax = s[i];
           if( x[i] < xmin ) xmin = x[i];
           if( x[i] > xmax ) xmax = x[i];
           if( y[i] < ymin ) ymin = y[i];
           if( y[i] > ymax ) ymax = y[i];
           if( z[i] < zmin ) zmin = z[i];
           if( z[i] > zmax ) zmax = z[i];
        }

/*  Move to center of molecule and rotate  */

        xmin = 0.5F*( xmax + xmin );
        ymin = 0.5F*( ymax + ymin );
        zmin = 0.5F*( zmax + zmin );

/*  Translate  
        - don't translate in xy version
           rotate about the axis
*/
        for( i=0; i<npts; i++) {
                x[i] = x[i] - xmin;
                y[i] = y[i] - ymin;
                z[i] = z[i] - zmin;

                x0 = x[i];      /*  Rotation  */
                y0 = y[i];
                x[i] =  cr*x0 - sr*y0;
                y[i] =  sr*x0 + cr*y0;

                y0 = y[i];      /*  Tilt  */
                z0 = z[i];
                y[i] =  ct*y0 + st*z0;
                z[i] = -st*y0 + ct*z0;
        }

/*  Sort by depth  ( Shell sort )  */

        printf( "Sorting atoms by depth...\n");

        m = 1;
        j = 2;
        do{ m+=1; j*=2;} while( j <= npts );
        m = j / 2;

        do{             /* ??? not most efficient stride but it works */
                k = npts - m;
                for( j=0; j<k; j++)
                for( i=j; i>=0; i-=m) {
                        if( z[i+m] < z[i] ) {
                                fswap( &x[i], &x[i+m] );
                                fswap( &y[i], &y[i+m] );
                                fswap( &z[i], &z[i+m] );
                                fswap( &s[i], &s[i+m] );
                                iswap( &icolor[i], &icolor[i+m] );
                        }
                }
                m = m/2;
        } while( m > 0 );

/*  Test sort routine
   DELETE this after awhile
*/
        for( i=1; i<npts; i++) 
                if( z[i-1] > z[i] ) printf("Bad sort in model3d !\n");

/*  Form 2D perspective projection and find scale of x,y */

        xmin = ymin = 1.0e30F;
        xmax = ymax = -xmin;

        for( i=0; i<npts; i++) {
                x[i] = x[i] /( d-z[i] );
                y[i] = y[i] /( d-z[i] );
                if( x[i] < xmin ) xmin = x[i];
                if( x[i] > xmax ) xmax = x[i];
                if( y[i] < ymin ) ymin = y[i];
                if( y[i] > ymax ) ymax = y[i];
        }
        scale = xmax - xmin;
        if( (ymax-ymin) > scale ) scale = ymax-ymin;
        scale = 1.10F * scale;
        scale = ( 1.0F - size0 ) / scale;

/*
   Scale coord. to ( 0.0 - 1.0 ) and display molecule
   SIZE0 must be in the 0.0-1.0 range
*/
        printf("Drawing atoms...\n");
        size = size0*( d - z[npts-1] );
        if( mode == TIFFmode ) {
                if( nxout > nyout ) fs = (float) nxout; 
                        else fs = (float) nyout;
                fsx = ((float)nxout) / fs;
                fsy = ((float)nyout) / fs;
                x0 = 0.5F* ( fsx - (xmax-xmin)* scale );
                y0 = 0.5F* ( fsy - (ymax-ymin)* scale );
                for( i=0; i<npts; i++) {
                        x[i] = ( x[i] - xmin ) * scale  + x0;
                        y[i] = ( y[i] - ymin ) * scale  + y0;
                        size2 = (s[i]/smax) * size / ( d - z[i] );
                        sphere( x[i], y[i], size2, icolor[i], 
                                pix, nxout, nyout );
                }
        } else {
                if( EPSxsize > EPSysize ) fs = EPSxsize; 
                        else fs = EPSysize;
                fsx = EPSxsize / fs;
                fsy = EPSysize / fs;
                x0 = 0.5F* ( fsx - (xmax-xmin)* scale );
                y0 = 0.5F* ( fsy - (ymax-ymin)* scale );
                for( i=0; i<npts; i++){
                        x[i] = ( x[i] - xmin ) * scale  + x0;
                        y[i] = ( y[i] - ymin ) * scale  + y0;
                        size2 = (s[i]/smax) * size / ( d - z[i] );
                        EPSsphere( fs*x[i], fs*y[i], fs*size2 );
                }
        }

/*  Save the scale  */

        if( (xmax-xmin) > (ymax-ymin) ) scale = (xmax-xmin);
                else scale = (ymax-ymin);
        *scalepix = scale;
        return;

}  /* end model3d() */

/*------------------------ iswap() ------------------------*/
/*
  Swap 2 int's (for sorting)
*/
void iswap( int *i, int *j )
{
        int it;
        it = *i;
        *i = *j;
        *j = it;
        return;
}

/*------------------------ fswap() ------------------------*/
/*
  Swap 2 float's (for sorting)
*/
void fswap( float *a, float *b )
{
        float t;
         t = *a;
        *a = *b;
        *b = t;
        return;
}

/*------------------------ sphere() ------------------------*/
/*
   Generate shaded sphere at IX0,IY0 with
   radius IRAD  (in units of pixels)

  x0, y0        = coordinates of sphere range= 0.0 to 1.0
  size          = diameter of sphere range 0.0 to 1.0
                  (this is a percentage of the full scale)
  icolor        =  0 red,  1 green,   2 blue
                        (not currently implemented!)
  l2pix[][]     = integer*2 array to get image
                  (remember that it ix,iy indexes are reversed )
  nxout, nyout  = size of output image in pixels
*/
void sphere( double x0, double y0, double size, int icolor, 
                long32 **l2pix, int nxout, int nyout )
{
        int ix, iy, ixd, iyd, ir2, ir, iy1, iy2,
                iyd2, ix1, ix2, ival1, ix0, iy0, irad, nm1;
        double val2, anout;

/*  scale coord. to max screen size of 0-NOUT  */

        nm1 = nxout;
        if( nyout > nm1 ) nm1 = nyout;
        nm1 -= 1;
        anout = (double) nm1 ;
        ix0 = (int) (x0 * anout + 0.5);
        iy0 = (int) (y0 * anout + 0.5);
        irad = (int) (0.5 * anout * size + 0.5);

/*  set up scale appropriate for unsigned char pixels  */

        ival1 = 255;
        val2  = 150.0;

/*  Calculate sphere size  */

        iy1 = iy0 - irad;
        iy2 = iy0 + irad;

        if( (nyout-1) < iy1 ) iy1 = (nyout-1);
        if( 0 > iy1 ) iy1 = 0;

        if( (nyout-1) < iy2 ) iy2 = (nyout-1);
        if( 0 > iy2 ) iy2 = 0;

        ir2 = irad*irad;

/*  Generate a shaded circle that looks like a 3D sphere */

        for( iy=iy1; iy<=iy2; iy++) {
           iyd  = iy - iy0;
           iyd2 = iyd*iyd;

           if( iyd2 < ir2 ) {
                ixd  = (int) (sqrt( ( (double)(ir2-iyd2) ) ) + 0.5);
                ix1  = ix0 - ixd;
                ix2  = ix0 + ixd;

                if( (nxout-1) < ix1 ) ix1 = (nxout-1);
                if( 0 > ix1 ) ix1 = 0;
                if( (nxout-1) < ix2 ) ix2 = (nxout-1);
                if( 0 > ix2 ) ix2 = 0;

                for( ix=ix1; ix<=ix2; ix++) {
                 ir = ix - ix0;
                 ir = iyd2 + ir*ir;
                  if( ir < ir2 )
                        l2pix[iy][ix] = ival1 - 
                                (long) floor( (val2 * ir) / ir2 + 0.5);
                }
           }
        }

        return;

}  /* end sphere */


/*--------------- EPSheader() ----------------------------*/
/*
  started from EPSheader() in publish.c 1-jul-1997 Earl J.Kirkland

  Output an image in postscript format

  filnam      : output file name
  xsize,ysize : (real valued) size in inches of whole image

  returned value is the error code (+1 for success )

*/ 
int EPSheader( const char filname[], double xsize, double ysize )
{
        int ix, iy, i, i1;
        
        time_t caltime;
        struct tm *mytime;
        char ctemp[32];

/*  Open output file  */

        fpEPS = fopen( filname, "w+" );
        if( fpEPS == NULL ) {
                printf( "Cannot open file %s in pspix()\n", filname );
                return( -1 );
        }

/*  Write EPSF postscript header info  */

        fprintf( fpEPS, "%s\n", "%!PS-Adobe-3.0 EPSF-3.0");
        fprintf( fpEPS, "%s\n", "%%Creator: publish.c");
        caltime = time( NULL );
        mytime = localtime( &caltime );
        strftime( ctemp, 34, "%H:%M, %d-%b-%Y",  mytime );
        fprintf( fpEPS, "%s %s\n", "%%CreationDate:", ctemp);
        fprintf( fpEPS, "%s %s\n", "%%Title:", filname);
        i  = (int) (xpos0 * 72);
        i1 = (int) (ypos0 * 72);
        ix = (int) ( (xpos0 + xsize) * 72);
        iy = (int) ( (ypos0 + ysize) * 72);
        fprintf( fpEPS, "%s %d %d %d %d\n", "%%BoundingBox: ",
                i, i1, ix, iy);
        fprintf( fpEPS, "%s\n", "%%EndComments");

/*  define functions  */

        fprintf( fpEPS, "%s\n", "%%BeginProlog");
        fprintf( fpEPS, "/inch { 72 mul } def\n");

        fprintf( fpEPS, "%%\n");
        fprintf( fpEPS, "%% macro to make a unit circle at (0,0)\n");
        fprintf( fpEPS, "%%\n /circle {newpath 0 0 1 0 360 arc\n" 
                        " closepath fill} def\n\n");

        fprintf( fpEPS, "%%\n");
        fprintf( fpEPS, "%%  macro to make a shaded sphere\n");
        fprintf( fpEPS, "%%  call as-->  xscale yscale xpos ypos sp\n");
        fprintf( fpEPS, "%%\n");
        fprintf( fpEPS, "/sp { gsave translate scale \n");
        fprintf( fpEPS, "   0.0 0.04 1 { sqrt 1 exch sub setgray circle \n");
        fprintf( fpEPS, "   0.98 0.98 scale  } for grestore } def\n\n");
         
        fprintf( fpEPS, "%s\n", "%%EndProlog");

        return( +1 );

}  /* end EPSheader() */

/*--------------- EPSsphere() ----------------------------*/
/*
  started  1-jul-1997 Earl J.Kirkland

  Output a shaded sphere in postscript format
  EPSheader(0 must be called before this routine and
  EPSend() after this routine

  xp,yp = (real valued) position of sphere in inches
  spsize = size of sphere in inches

  (remember 1 inch = 72 points )

  returned value is the error code (+1 for success )

*/ 
int EPSsphere( double xp, double yp, double spsize )
{

        fprintf( fpEPS, "%9.2f %9.2f %9.2f %9.2f sp\n", 
                0.5*spsize*72, 0.5*spsize*72,
                (xpos0 + xp)*72, (ypos0 + yp)*72 );
        return( +1 );

}  /* end EPSsphere() */

/*--------------- EPSend() ----------------------------*/
/*
  started from EPSend() in publish.c  1-jul-1997 Earl J.Kirkland

  write the closing statements for EPS file

  returned value is the error code (+1 for success )

*/ 
int EPSend(  )
{
        /*  postcript trailer  */

        fprintf( fpEPS, "showpage\n");
        fprintf( fpEPS, "%s\n", "%%EOF" );

        fclose( fpEPS );
        return( +1 );

}  /* end EPSend() */


