/*         *** display.c ***

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

   ANSI C version

  Display multislice binary data files

  The currently supported image output devices are:
        1  Print out single (1D) line of data as text
        2  Print 2D image as ASCII text file
        3  (EPS) Postscript greyscale in file
        4  (EPS) Postscript contour plot in file

   This file is formatted for a tab size of 8 characters

   started  FORTRAN version 18-JUL-84 EJK
    :
    :
   started conversion to C and removed a lot of old output device
      options  24-feb-1996 ejk
   convert to TIFF image file format 19-apr-1997 ejk
   added EPS headers 23-apr-1997 ejk
   add CreationDate to EPS header 4-may-1997 ejk
   fixed 1D mode  24-jun-1997 ejk
   removed commas from input format 16-july-1997 ejk
   update memory allocation routines 13-nov-1999 ejk
   add error checking on file open in 1D mode 14-nov-1999 ejk
   fixed several errors with interpolation in 1D line scan mode
             20-nov-1999 ejk
    change void main() to int main() for better portability
             21-jan-2000 ejk
    change data type of nxl,nyl to long32 to be compatible with
        new tiffsubs  17-jul-2007 ejk
    convert to GPL 3-jul-2008 ejk

*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "slicelib.h"   /* misc. routines for multislice */
#include "tiffsubs.h"   /* file I/O routines in TIFF format */

/*  define mode codes for each output device */
#define maxMode         4       /* maximum number of output modes */
#define lineMode        1       /* Print out single line of data */
#define textMode        2       /* ASCII text file */
#define psGreyMode      3       /* Postscript greyscale in file pspix.eps */
#define psContourMode   4       /* Postscript contour plot in file pspix.eps */

#define NCMAX 256       /* max characters in file names */
#define NPARAM  64      /* number of parameters */

/* define functions at the end of this file */
double pinter( float **pix, int nx, int ny, double x, double y, int linter );
int pspix( const char filnam[], float **rpix, int nx, int ny,
        float rmin, float rmax, double xsize, double ysize,
        const char ctitle[], int lines, int lchar );
int pscontour( const char filnam[], float **rpix, int nx, int ny,
        float zmin, float zmax, int nclevl, double xsize, double ysize,
        const char ctitle[], int lines, int lchar );

int main()
{
        char datetime[20], infile[NCMAX], ctitle[300], outfile[NCMAX];

        int lrescale, lcomp, ix, iy, nx, ny, mode, i, nclev, npout,
            nbits[3], samples, npix, imag;
                
        long32 nxl, nyl, **lpix;   /* tiffsubs.h data types */

        float **pix, **impix, *param, **thepix;
        float r, rmin, rmax, rmin0, rmax0, imin, imax, imin0, imax0,
             amin, amax;

        double xout, yout, a, b, scalex, scaley, xsize, ysize, 
                xinit, xfinal, yinit, yfinal, x, y, zipixel, zpixel;

        FILE *fp;
        time_t caltime;
        struct tm *mytime;

/*  get output device type  */

        printf("display version dated 3-jul-2008 (ejk)\n");
        printf("Copyright (C) 1998, 2008 Earl J. Kirkland\n" );
        printf( "This program is provided AS-IS with ABSOLUTELY NO WARRANTY\n "
            " under the GNU general public license\n\n" );

        printf("generate display of TEMSIM images\n");

        do{
           printf( "The available image output modes are:\n"
                   "    code  mode\n"
                   "     %d    print single (1D) line of image data as text\n"
                   "     %d    ASCII text file of 2D image\n"
                   "     %d    postscript (EPS) greyscale in file\n"
                   "     %d    postscript (EPS) contour plot in file\n"
                   "Enter code number:\n", lineMode, textMode, psGreyMode, 
                                        psContourMode );
           scanf("%d", &mode );
        }while( (mode<1) || ( mode > maxMode ) );

/*  Get input/output file name  */

        printf("Name of file that has binary data to display:\n");
        scanf("%s", infile );

        printf("Name of output file:\n");
        scanf("%s", outfile );

        param = (float*) malloc1D( NPARAM, sizeof(float), "param" );

/* First try to read it as an extended TIFF floating point image.
        1. Remember that there is a plain 8 bit image in the front
            of the file so don't test in the other order.
        2. Store all forms in a float2D pix so the output modes have
            a common format. */
                        
        if( tFloatTest( infile ) == 1 ) {
           if( topenFloat( infile ) != 1 ) {
                printf("Cannot open file %s\n", infile );
                exit( 0 );
           }
           tsize( &nxl, &nyl, nbits, &samples );
           nx = (int) nxl;
           ny = (int) nyl;
           pix = (float**) malloc2D( nx, ny, sizeof(float), "pix" );
           if(treadFloatPix( pix, nxl,nyl, &npix, datetime, param) != 1)
           {    printf("Cannot read input file %s.\n", infile );
                exit(0);
           }
           tclose();
           if( (npix<1) || (npix>2) ) {
                printf( "bad npix = %d in TIFF file %s\n", npix, infile );
                exit( 0 );
           }

           /* real and imag. pix stacked together */
           if( npix > 1 ) {
                nx = nx / npix;
                impix = pix + nx;       /* imag pix */
           }
           a = param[pDX] * ((float)nx);
           b = param[pDY] * ((float)ny);
           rmin = param[pRMIN];
           rmax = param[pRMAX];
           imin = param[pIMIN];
           imax = param[pIMAX];
     
        /*  Next try to read an 8/16 bit integer pix 
            In principle the conversion of long lpix to float
                pix could be done in place to save memory but
                you can't always guarantee that a float and a long
                are 32 bits on every CPU, so use two arrays   */

        } else if( tifftest( infile ) == 1 ) {

           if( topen( infile ) != 1 ) {
                printf("Cannot open file %s\n", infile );
                exit( 0 );
           }
           tsize( &nxl, &nyl, nbits, &samples );
           nx = (int) nxl;
           ny = (int) nyl;
           pix = (float**) malloc2D( nx, ny, sizeof(float), "pix" );
           lpix = (long32**) malloc2D( ny, nx, sizeof(long32), "lpix" );
           if( treadPix( lpix ) != 1 ) {
                printf("Cannot read input file %s.\n", infile );
                exit(0);
           }
           tclose();
           npix = 1;
           a = (float) nx;
           b = (float) ny;
           strcpy( datetime, " " );
           rmin = rmax = (float) lpix[0][0];
           for( ix=0; ix<nx; ix++) 
              for( iy=0; iy<ny; iy++) {
                  r = pix[ix][iy] = (float) lpix[iy][ix];
                  if( r < rmin ) rmin = r;
                  if( r > rmax ) rmax = r;
            }
        
        } else{
                printf( "%s is not a valid TIFF file.\n", infile );
                exit( 0 );
        }

        printf("Size of old image in pixels, Nx, Ny = %d, %d\n"
                "  created %s\n",  nx, ny, datetime);
 
/*  echo range of image  */

        if( npix == 1 ) {                       /* integer or real image */
                printf("Pix goes from %f to %f\n", rmin, rmax );
        } else if( npix == 2 ) {                /* complex image */
                printf("Pix goes from %f to %f (real)\n", rmin, rmax );
                printf("     and from %f to %f (imag)\n", imin, imax );
        }
        printf("Lattice constants are: ax, by = %f, %f\n", a, b );

/*  ask for options if appropriate */

        if( (mode != lineMode) && (npix==2) ){
            printf("Type 1 for real part and 2 for imaginary part:\n");
            scanf( "%d", &imag );
            if( (imag <1) || (imag>2) ) imag = 1;
        } else imag = 1;

        if( (mode == psContourMode) || (mode == psGreyMode) ) {
                lrescale = askYN("Do you want to rescale greyscale range");
                lcomp = askYN("Do you want to complement image");
        } else {
                lcomp = 0;
                lrescale = 0;
        }

        if( lrescale == 1 ) {
            if( npix == 1 ) {           /*  real valued image */
                printf("New pix display range, min and max =\n");
                scanf("%f %f", &rmin0, &rmax0);
            } else if( npix == 2 ) {    /* complex valued image */
                printf("New real pix display range, min and max =\n");
                scanf("%f %f", &rmin0, &rmax0);
                printf("New imaginary pix display range, min and max =\n");
                scanf("%f %f", &imin0, &imax0);
            }
        }

        if( lrescale == 1) {
          rmin = rmin0;
          rmax = rmax0;
        }

        if(  mode == psContourMode ) {
                printf("Number of contour levels:\n");
                scanf("%d", &nclev);
                if( nclev <= 0) nclev = 7;
        }

        if( mode == lineMode)  {
                printf("Type initial x y position:\n");
                scanf("%lf %lf", &xinit, &yinit );
                if( (xinit<0) || (xinit>a) || (yinit<0) || (yinit>b) ) {
                        printf( "x and/or y out of range\n" );
                        exit( 0 );
                }
                printf("Type final x,y position:\n");
                scanf("%lf %lf", &xfinal, &yfinal);
                if( (xfinal<0) || (xfinal>a) || (yfinal<0) || (yfinal>b) ) {
                        printf( "x and/or y out of range\n" );
                        exit( 0 );
                }
                printf("Type number of points for output:\n");
                scanf("%d", &npout);
                if( npout <= 0 ) {
                        printf( "bad range\n" );
                        exit( 0 );
                }
        }

/*  Complement if requested  */

        if ( lcomp == 1 ) {
                for( ix=0; ix<nx; ix++) 
                for( iy=0; iy<ny; iy++) 
                        pix[ix][iy] = ( rmax - pix[ix][iy] ) + rmin;
            if( npix == 2 ) {
                for( ix=0; ix<nx; ix++) 
                for( iy=0; iy<ny; iy++) 
                        impix[ix][iy] = ( imax - impix[ix][iy] ) + imin;
            }
        }

/*  Print out single (1D) line as text here  */

        if ( mode == lineMode ) {
           fp = fopen( outfile, "w+" );
           if( fp == NULL ) {
                   printf( "cannot open %s for writing\n", outfile );
                   exit( 0 );
           }
           fprintf(fp, "C\n");
           fprintf(fp, "C  Image Data File : \n");
           fprintf(fp, "C   %s\n", infile);
           fprintf(fp, "C\n");
           fprintf(fp, "C  npix : %d\n", npix) ;
           fprintf(fp, "C  Created : %s\n", datetime );
           fprintf(fp, "C  initial X,Y : %f, %f\n", xinit, yinit);
           fprintf(fp, "C    final X,Y : %f, %f\n", xfinal, yfinal);
           fprintf(fp, "C\n");
           fprintf(fp, "C        x              y              value\n");
           fprintf(fp, "C\n");

           scalex = ( xfinal - xinit ) / ( npout - 1);
           scaley = ( yfinal - yinit ) / ( npout - 1);
           for( i=0; i<npout; i++) {
                xout = scalex * i + xinit;
                yout = scaley * i + yinit;
                x = (nx-1) * xout / a;
                y = (ny-1) * yout / b;
                zpixel = pinter( pix, nx, ny, x, y, 1 );
                if( npix == 2 ) {
                   zipixel = pinter( impix, nx, ny, x, y, 1 );
                   fprintf( fp, "%15.7f %15.7f %15.7f %15.7f\n",
                         xout, yout, zpixel, zipixel );
                } else
                   fprintf( fp, "%15.7f %15.7f %15.7f\n", xout, yout, zpixel );
           }

          fclose( fp );

/*  Print out 2D image as text here  */

        } else if ( mode == textMode ) {
           fp = fopen( outfile, "w+" );
           if( fp == NULL ) {
                   printf( "cannot open file %s for writing\n", outfile );
                   exit(0);
           }

           for( iy=(ny-1); iy>=0; iy--) {
                for( ix=0; ix<nx; ix++)
                    if( imag == 1 ) fprintf( fp, "%12.4g ", pix[ix][iy]);
                    else fprintf( fp, "%12.4g ", impix[ix][iy]);
                fprintf( fp, "\n" );
           }

          fclose( fp );

/*  Postscript laserprinter greyscale or contour plot
     let postscript do all the scaling etc
*/
        } else if( (mode == psGreyMode) || (mode == psContourMode ) ) {

           strcpy( ctitle, "Pix: " );
           strncpy( (ctitle+5), infile, 55 );   /* 1st line */
           
           if( imag == 1 ) {
               amin = rmin;
               amax = rmax;
               thepix = pix;
           } else {
               amin = imin;
               amax = imax;
               thepix = impix;
           }
           i = sprintf( (ctitle+70), "range: %12.6g to %12.6g;",
               amin, amax );
           caltime = time( NULL );
           mytime = localtime( &caltime );
           strftime( (ctitle+70+i), 34, " printed at: %H:%M, %d-%b-%Y",
               mytime );
           
           sprintf( (ctitle+140), "Nx, Ny = %d, %d ; a,b = %f, %f",
                         nx, ny, a, b );

           if( a > b ) {
                xsize = 6.0;
                ysize = xsize * ( b/a );
           } else {
                ysize = 6.0;
                xsize = ysize * ( a/b );
           }         
           if( mode == psGreyMode ) {
                pspix( outfile,  thepix, nx, ny, 
                        amin, amax, xsize, ysize, ctitle, 3, 70 );
           } else if( mode == psContourMode ) {
                pscontour( outfile,  thepix, nx, ny, amin,
                        amax, nclev, xsize, ysize, ctitle, 3, 70 );
           }
           printf("Please print file %s.\n", outfile);

        }  /* end if( mode == ...) loop */

        return 0;

}  /* end main() */


/*--------------------- pinter() ----------------------------*/
/*
  Bilinear interpolation from pix array

  pix[][]  = real input array with the image
  nx,ny    = integer dimension of the image
  x,y      = real coord of image point to interpolate at
  linter   = logical - if true then interpolate
                  otherwise just take nearest value
*/
double pinter( float **pix, int nx, int ny, double x, double y,
              int linter )
{
        int ix, iy, ix2, iy2;
        double  a, b, c, d, x1, y1, ans;

        if ( linter == 1 ) {
                ix = (int) x;
                iy = (int) y;
                if( ix > (nx-2) ) ix = (nx-2);
                if( ix <    0 )   ix = 0;
                if( iy > (ny-2) ) iy = (ny-2);
                if( iy <    0 )   iy = 0;
                ix2 = ix + 1;
                iy2 = iy + 1;
                x1 = (float) ix;
                y1 = (float) iy;

                d = pix[ix][iy]  - pix[ix][iy2] - pix[ix2][iy]
                                 + pix[ix2][iy2];
                c = pix[ix][iy2] - pix[ix][iy]  - d*x1;
                b = pix[ix2][iy] - pix[ix][iy]  - d*y1;
                a = pix[ix][iy]  - b*x1 - c*y1 - d*x1*y1;

                ans = a + b*x + c*y + d*x*y;

        } else {
                ix = (int) (x + 0.5);
                iy = (int) (y + 0.5);
                if( ix > (nx-1) ) ix = (nx-1);
                if( ix <    0 )   ix = 0;
                if( iy > (ny-1) ) iy = (ny-1);
                if( iy <    0 )   iy = 0;

                ans = pix[ix][iy];
        }
        
        return( ans );

}  /* end pinter() */


/*--------------- pspix() ----------------------------*/
/*
  started 18-may-1988 Earl J.Kirkland
  fixed round-off err problem with intensity scaling
         23-may-1988 ejk
  converted to C 3-mar-1996 ejk
  converted to EPS 23-apr-1997 ejk
  added CreationDate 4-may-1997 ejk

  Output a sampled image in postscript format

  filnam      : output file name
  rpix[][]    : real array holding image to output
  nx,ny       : (integer valued) dimension of the image
  rmin, rmax  : (real valued) range of image values to output
  xsize,ysize : (real valued) size in inches of output pix
                  if <= 0 then default to 7x7 inch
  ctitle      : character array with title
                   substring split on to different lines
  lines       : (integer valued) number of to split 
                  ctitle into
  lchar       : number of character per line

  returned value is the error code (+1 for success )

*/ 
int pspix( const char filnam[], float **rpix, int nx, int ny,
        float rmin, float rmax, double xsize, double ysize,
        const char ctitle[], int lines, int lchar )
{
        int ix, iy, i, i1, ip;
        long nbytes;
        double scale, xsize2, ysize2, xpos0=1.0, ypos0=1.0, ypos;
        
        FILE *fp;
        time_t caltime;
        struct tm *mytime;
        char ctemp[32];

/*  Open output file  */

        fp = fopen( filnam, "w+" );
        if( fp == NULL ) {
                printf( "Cannot open file %s in pspix()\n", filnam );
                return( -1 );
        }

/*  Write EPSF postscript header info  */

        fprintf( fp, "%s\n", "%!PS-Adobe-3.0 EPSF-3.0");
        fprintf( fp, "%s\n", "%%Creator: pspix()");
        caltime = time( NULL );
        mytime = localtime( &caltime );
        strftime( ctemp, 34, "%H:%M, %d-%b-%Y",  mytime );
        fprintf( fp, "%s %s\n", "%%CreationDate:", ctemp);
        fprintf( fp, "%s %s\n", "%%Title:", filnam);
        i  = (int) (xpos0 * 72) - 2;
        i1 = (int) (ypos0 * 72) - 2;
        ix = (int) ( (xpos0 + xsize) * 72) + 2;
        iy = (int) ( (ypos0 + ysize) * 72) + 2;
        if( (lines > 0) && (lchar > 0) ) iy += 18*(lines+2);
        fprintf( fp, "%s %d %d %d %d\n", "%%BoundingBox: ", i, i1, ix, iy);
        fprintf( fp, "%s\n", "%%EndComments");

/*  define functions  */

        fprintf( fp, "%s\n", "%%BeginProlog");
        fprintf( fp, "/picstr %d string def\n", nx );
        fprintf( fp, "/inch { 72 mul } def\n");
        fprintf( fp, "/makeimage { %d %d 8 [ %d 0 0 %d 0 0]\n",
                 nx, ny, nx, ny);
        fprintf( fp, "{currentfile picstr readhexstring pop} image} def\n");
        fprintf( fp, "%s\n", "%%EndProlog");

/*  Display title if requested  */

        if( (lines > 0) && (lchar > 0) ) {
                fprintf( fp, "/Times-Roman findfont 12 scalefont setfont\n");

                for( i=0; i<lines; i++) {
                        ypos = ypos0 + ysize + 
                                (18.0/72.0)*((double)(lines-i)) + 0.25;
                        fprintf( fp, "%f inch %f inch moveto\n", xpos0, ypos);
                        fprintf( fp, "(%.*s) show\n", lchar, &ctitle[i*lchar]);
                }
                fprintf( fp, "\n" );
        }

/*  define size and position of the ps image  */

        fprintf( fp, "%f inch %f inch translate\n", xpos0, ypos0 );

        if ( xsize <= 0.0)  xsize2 = 7.0;
                else xsize2 = xsize;

        if ( ysize <= 0.0)  ysize2 = 7.0;
                else ysize2 = ysize;

        fprintf( fp, "%f inch %f inch scale\n", xsize2, ysize2);
        nbytes = (long) ( ny*ceil( (2.0*nx)/70.0) + 1);
        fprintf( fp, "%s %ld Hex Lines\n", "%%BeginData:", nbytes);
        fprintf( fp, "makeimage\n");

/*  Output image
        normalize to 8 bits/pixels with saturation
*/
        scale = 255.0 / ( rmax - rmin );

        for( iy=0; iy<ny; iy++ ) {
                i = 0;
                for( ix=0; ix<nx; ix++) {
                        ip = (int) ( scale * ( rpix[ix][iy] - rmin ) + 0.5 );
                        if( ip < 0 ) ip = 0;
                        if( ip > 255 ) ip = 255;
                        fprintf( fp, "%2.2x", ip );
                        i = i + 2;
                        if( i > 70  ) {
                                fprintf( fp, "\n");  i = 0;
                        }
                }
                if( i != 0 ) fprintf( fp, "\n");
        }  /* end for(iy...) */
        fprintf( fp, "%s\n", "%%EndData");

/*  postcript trailer  */

        fprintf( fp, "showpage\n");
        fprintf( fp, "%s\n", "%%EOF" );

        fclose( fp );
        return( +1 );

}  /* end pspix() */


/*--------------- pscontour() -----------------------------*/
/*
  started from pspix 21-june-1992 Earl J. Kirkland
  converted to C 3-mar-1996 ejk
  add EPS format 22-apr-1997 ejk

  output a sampled image in postscript format as a contour plot

  filnam[]    : output file name
  rpix[][]    : real array holding image to output
  nx,ny       : (integer valued) dimension of the image
  zmin, zmax  : (real valued) range of image values to output
  nclev       : (integer) number of levels to contour
  xsize,ysize : (real valued) size in inches of output pix
                  if .le. 0 then default to 8x8 inch
  ctitle      : character*n array with title
                   substring split on to different lines
  lines       : (integer valued) number of to split 
                  ctitle into
  lchar       : number of character per line

        returned value is the error code ( +1 for success )
*/ 
int pscontour( const char filnam[], float **rpix, int nx, int ny,
        float zmin, float zmax, int nclevl, double xsize, double ysize,
        const char ctitle[], int lines, int lchar )
{
        char cline[83];
        int ix, iy, i1, i, npt, nclev, ictr, icorn;
                
        double xsize2, ysize2, xpos, ypos, x1, y1, ctr, xm, ym,
                d0,d1, dmax, dmin, delta, scaled, scalex, scaley, dx,dy;

        double xpos0=1.0, ypos0=1.0;  /* origin of page */

        int idx[]  = {-1, -1,  0, 0, -1};
        int iddx[] = { 0, 1, 0, -1};
        int idy[]  = { 0, -1, -1, 0,  0 };
        int iddy[] = {-1, 0, 1,  0};

        FILE *fp;
        time_t caltime;
        struct tm *mytime;
        char ctemp[32];

/*  Open output file  */

        fp = fopen( filnam, "w+" );
        if( fp == NULL ) {
                printf( "Cannot open file in pscontour()\n" );
                return( -1 );
        }

/*  Write EPSF postscript header info  */

        fprintf( fp, "%s\n", "%!PS-Adobe-3.0 EPSF-3.0");
        fprintf( fp, "%s\n", "%%Creator: pscontour()");
        caltime = time( NULL );
        mytime = localtime( &caltime );
        strftime( ctemp, 34, "%H:%M, %d-%b-%Y",  mytime );
        fprintf( fp, "%s %s\n", "%%CreationDate:", ctemp);
        fprintf( fp, "%s %s\n", "%%Title: ", filnam);
        i  = (int) (xpos0 * 72) - 2;
        i1 = (int) (ypos0 * 72) - 2;
        ix = (int) ( (xpos0 + xsize) * 72) + 5;
        iy = (int) ( (ypos0 + ysize) * 72) + 5;
        if( (lines > 0) && (lchar > 0) ) iy += 18*(lines+2);
        fprintf( fp, "%s %d %d %d %d\n", "%%BoundingBox: ", i, i1, ix, iy);
        fprintf( fp, "%s\n", "%%EndComments");

/* define functions */

        fprintf( fp, "%s\n", "%%BeginProlog");
        fprintf( fp, "/inch { 72 mul } def\n");
        fprintf( fp, "/ml { moveto lineto stroke } def\n");
        fprintf( fp, "%s\n", "%%EndProlog");

/*  Display title if requested  */

        if( (lines > 0) && (lchar > 0) ) {
                fprintf( fp, "/Times-Roman findfont 12 scalefont setfont\n");

                for( i=0; i<lines; i++) {
                        ypos = ypos0 + ysize + 
                                (18.0/72.0)*((double)(lines-i)) + 0.25;
                        fprintf( fp, "%f inch %f inch moveto\n", xpos0, ypos);
                        fprintf( fp, "(%.*s) show\n", lchar, &ctitle[i*lchar]);
                }

                fprintf( fp, "\n" );
        }

/*  calculate and display contour levels  */

        if( nclevl > 1 ) {
                scaled= (zmax-zmin) / ((float)(nclevl+1));
                nclev = nclevl;
        } else {
                scaled = 1.0F;
                nclev = 1;
        }

        ypos = ypos0 + ysize + (18.0/72.0)*(-0.05) + 0.25;
        fprintf( fp, "%f inch %f inch moveto\n", xpos0, ypos);
        fprintf( fp, "(%d Contour Levels = ) show\n", nclev);

        for( ictr=0; ictr<nclev; ictr++) {
                ctr= zmin + ((float)(ictr+1)) * scaled;
                i1 = ictr*9 ;
                if( i1 < (80-9) ) {
                        sprintf( &cline[i1],"%8.2f,",  ctr);
                }
        }
        if( i1 > 80 ) {
            strcpy( &cline[i1+9], "...");
            i1 = i1+3;
        }
        cline[i1+9] = '\0';
        ypos = ypos0 + ysize + (18.0/72.0)*(-0.75) + 0.25;
        fprintf( fp, "%f inch %f inch moveto\n", xpos0, ypos);
        fprintf( fp, "(%s) show\n",  cline);

/*  define size and position of the ps image  */

        fprintf( fp, "%f inch %f inch translate\n",  xpos0, ypos0 );

        if ( xsize <= 0.0)  xsize2 = 8.0;
                else  xsize2 = xsize;

        if ( ysize <= 0.0) ysize2 = 8.0;
                else ysize2 = ysize;

/*  Convert inches to points  */

        xsize2 = xsize2 * 72.0;
        ysize2 = ysize2 * 72.0;

/*  Calculate scale points/pixel  */

        scalex = xsize2 / ((float)(nx-1));
        scaley = ysize2 / ((float)(ny-1));

/*    draw a box  */

        fprintf( fp, "stroke newpath 2 setlinewidth\n");
        xpos = 0.0;
        ypos = 0.0;

        fprintf( fp, "%.2f %.2f moveto\n",  xpos, ypos);

        fprintf( fp, "%.2f %.2f lineto\n", xpos+xsize2, ypos);
        fprintf( fp, "%.2f %.2f lineto\n", xpos+xsize2, ypos+ysize2);
        fprintf( fp, "%.2f %.2f lineto\n", xpos, ypos+ysize2);
        fprintf( fp, "%.2f %.2f lineto\n", xpos, ypos);

        fprintf( fp, "stroke newpath 0.5 setlinewidth\n");

/*  Plot contours  (do NOT contour rmin and rmax)  */

        for( iy=1; iy<ny; iy++)     /*  loop over data points  */
        for( ix=1; ix<nx; ix++) {

        /*  loop over contour levels  */

           for( ictr=1; ictr<=nclev; ictr++) {
                ctr = zmin + ((float)ictr) * scaled;
                npt = 0;

/*  loop over four corners of current grid point  */

                for( icorn=0; icorn<4; icorn++) {
                   d0= rpix[ix +idx[icorn] ][ iy + idy[icorn] ];
                   d1= rpix[ix +idx[icorn+1] ][ iy + idy[icorn+1] ];
                   if( d0 > d1 ) dmax = d0; else dmax = d1;
                   if( d0 < d1 ) dmin = d0; else dmin = d1;
                   if( (ctr <= dmax) && (ctr > dmin) ) {
                        npt   = npt + 1;
                        delta = (ctr - d0) / (d1 - d0);

                        dx  = delta * ((float)iddx[icorn]) +
                                         ((float)(ix-1+idx[icorn]));
                        x1 = xpos + dx * scalex;
                        if( xpos > x1 ) x1 = xpos;
                        if( (xsize2+xpos) < x1) x1 = xsize2+xpos;

                        dy  = delta * ((float)iddy[icorn]) +
                                ((float)(iy-1+idy[icorn]));
                        y1 = ypos + dy * scaley;
                        if( ypos > y1 ) y1 = ypos;
                        if( (ysize2+ypos) < y1 ) y1 = ysize2+ypos;

                        if ( npt <= 1 ) {
                                xm = x1; ym = y1;
                        } else {
                                fprintf( fp, "%.2f %.2f %.2f %.2f ml\n",
                                        xm, ym, x1, y1);
                                npt = 0;
                        }
                   }
                }  /* end for(icorn...) */
           }  /* end for(ictr...) */
        }  /* end for(ix... ) */

/*  postscript trailer  */

        fprintf( fp, "showpage\n");
        fprintf( fp, "%s\n", "%%EOF" );

        fclose( fp );

        return( +1 );

}  /* end pscontour() */
