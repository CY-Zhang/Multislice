/*      *** autostem.c ***

modified 02-03-10 to use CEOS probe aberrations up to A5.  pmv
02-04-10 switched iseed to millisecond instead of second timer.  pmv
02-23-10 switched iseed to user input instead of millisecond timer.  pmv
03-08-10 switched all file open from mode w+ to mode w for Condor compability, pmv
03-18-10 added Monte Carlo integration of dE / Cc
06-10-16 added sample rotation option, cz
11-03-16 minor bug fixed for using more than 2 detector range, cz


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

-----------------------------------------------------------------------------
   ANSI-C and TIFF version

  Calculate STEM images and line scans from a non-periodic
  distribution of (x,y,z) atomic coordinates using repetitive multislice
  
  The transmission function for each slice can take a lot of
  computer time.  This program propagates the STEM probe for many
  beam position through the specimen at the same time to avoid recalculating
  the specimen transmission function for each probe position. This
  requires a huge amount of memory.  In 1D line scan mode a 2D probe
  wave function is stored for all positions.  In 2D image mode the 2D
  probe wave functions for a whole line are stored simulataneously.

  this file is formatted for a tab size of 4 characters

  multithreaded code using openMP may be turned on/off with
  symbol USE_OPENMP (ignore undefined pragma warning if off)

  Ref:
  [1] Zhiheng Yu, Philip Batson, John Silcox, "Artifacts in Aberration
      Corrected ADF-STEM Imaging", Ultramicroscopy 96 (2003) p.275-284

  started from stemslic.c and autoslic.c 19-jul-1998 E. Kirkland
  convert to simultaneous transmission of many probes at the same
    time to reuse transmission functions  28-jul-1998 ejk
  finished 2D image mode (1D mode works OK) 29-jul-1998 ejk
  last modified 29-jul-1998 ejk
  add Mote-Carlo integration of source size 24-aug-1998 ejk
  fixed typo in random dither in y direction of probe position
    1-feb-1999 ejk
  fixed error in nprobes in image mode (should be nyout
      but it was nxout- at top question) 3-july-2001 ejk
  update memory allocation routines,
     change void main() to int main() for better portability,
     and add 5th order spherical aberration  3-mar-2004 ejk
  start modification to keep multiple thickness's 21-jul-2005 ejk
      - in working form 26-jul-2005 ejk
  put in faster sorting routine sortByZ() 6-feb-2006 ejk
  add some openMP multithreading 23-may-2006 ejk
  move openMP stuff into conditional compilation
     so I only need one version of the code,
     and move sortbyz() to sliclib 23-aug-2006 ejk
  add periodic() to put probe position back into the supercell
     and recal. x,y position without source size wobble
     for output in 1D mode 31-aug-2006 ejk
  change propagation range to be whole unit cell not just
     range of atoms to treat sparsely populated spec.
     better 14-sep-2006 ejk
  change range to start at -0.25*deltaz 26-sep-2006 ejk
  minor cosmetic changes 11-jul-2007 ejk
  minor fixes for openMP stuff 27-feb-2008 ejk
  move vzatomLUT() to slicelib.c and reformat for TAB size of 4 char
               11-jun-2008 ejk
  convert to GPL 3-jul-2008 ejk
  fix bug in multithreading (add more private var.) 
      5-nov-2008 ejk
*/

#include <stdio.h>  /* ANSI C libraries used */
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "fft2dc.h" /* FFT routines */
#include "slicelib.h"   /* misc. routines for multislice */
#include "tiffsubs.h"   /* file I/O routines in TIFF format */


/*#define USE_OPENMP*/      /* define to use openMP */

#ifdef USE_OPENMP
#include <omp.h>
/*  get wall time for benchmarking openMP */
#define walltim() ( omp_get_wtime() )
double walltimer;
#endif

#define BW (2.0F/3.0F)  /* antialiasing bandwidth limit factor */
#define ABERR 1.0e-5    /* max error for a,b */

#define NSMAX   1000    /* maximum number of total slices */
#define NCMAX   512 /* max number of characters per line */
#define NCINMAX 500 /* max number of characters in stacking spec */
#define NLMAX   52  /* maximum number of layers */
#define NPARAM  64  /* number of parameters */

#define NZMAX   103 /* max atomic number Z */
#define NRMAX   100 /* number of in look-up-table in vzatomLUT */
#define RMIN    0.01    /* r (in Ang) range of LUT for vzatomLUT() */
#define RMAX    5.0

/* global specimen and probe parameters to share */

float ***prober, ***probei; /* complex probe wave functions */
float **transr, **transi;   /* complex transmission functions */
float *propxr, *propxi;     /* complex propagator vs x */
float *propyr, *propyi;     /* complex propagator vs y */
float zmin, zmax;
float **sensmap, ***sensmap_probe, ***sensmap_save, **sensmap_rotated, **sensmap_rotated_save;

int nx, ny, nxprobe, nyprobe, nslice;
int natom;
float *kx, *ky, *kx2, *ky2, *kxp, *kyp, *kxp2, *kyp2, *mradx, *mrady, *mradx_sens, *mrady_sens;
float  **rmin, **rmax;
float *xa, *ya, *za, *occ, *wobble;
float *xa2, *ya2, *za2, *occ2;
int natom, *Znum, *Znum2;
double wavlen, k2maxp, apert1, apert2, pi, keV;
float ax, by, cz;                   /*  specimen dimensions */
double **aber; /* aberration parameter list */
double  Df0, cDf_sigma;  /* parameters for chromatic aberration defocus spread */
double *almin, *almax, *k2max, *k2min, deltaz;
long nbeamt;
unsigned long iseed;  /* seed for random number gen.  shared with STEMsignals, so global */

/* define functions at end of this file (i.e. so main can be 1st) */
double periodic( double pos, double size );
void STEMsignals( double x[], double y[], int npos,
         double ***detect, int ndetect,
         double ThickSave[], int nThick, double sum[], float ***sensmap_probe );

/* spline interpolation coeff. */
int splineInit=0, *nspline;
double  *splinx, **spliny, **splinb, **splinc, **splind;

void trlayer( const float x[], const float y[], const float occ[],
    const int Znum[], const int natom, 
    const float ax, const float by, const float kev,
    float **transr, float **transi, const long nx, const long ny,
    double *phirms, long *nbeams, const float k2max );


/* evaluate the aberration phase plate using CEOS aberration definitions */
/*
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
double aphase( double **aber, double wl, double kx, double ky, double dx, double dy)
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

    return( c );
}



int main()
{
    char filein[NCMAX], fileout[NCMAX], fileoutpre[NCMAX],
        description[NCMAX], datetime[24], ch[NCMAX];
  
    const char version[] = "16-jun-2016 (pmv)";

    int ix, iy, nx2, ny2, i, idetect, nout, nxout, nyout,
        ncellx, ncelly, ncellz, iwobble, nwobble,
        ndetect, nprobes, ip, nThick, it, senx, seny, npixx = 1,
        minIdx, minIdy, nbits[3], j, samples;

    int  l1d, lwobble, lxzimage, OutCut;
    int px_in_1, px_in_2, centerx, centery;

    long nbeamp;

    float *x2, *y2;
    float *param, ***pixr, **pixout, temp, pmin, pmax;
    float wmin, wmax, xmin,xmax, ymin, ymax, temperature;
    float *param_detector, reference_x, reference_y, mindifference, difference, reference;
    float px_in_1f, px_in_2f, px_out_1f, px_out_2f, temp_sens, *ratio;
    float length, ori_angle, new_angle, mradx_rot, mrady_rot;

    double scale, *x, *y, sum, *sums, w, ***detect, ***detect2,
       tctx, tcty, xi,xf, yi,yf, dx, dy, totmin, totmax,
       ctiltx, ctilty, timer, sourcesize, sourceFWHM, *ThickSave;


    double Cc, dE;
    double rot_angle, val, i_new, j_new;
    char **Sensitivity; 
    char **probe_file;

    int i_new_int, j_new_int;


    FILE *fp;

/* start by announcing version etc */

    printf("autostem version dated %s\n", version );
    printf("Copyright (C) 1998, 2008 Earl J. Kirkland\n" );
    printf( "This program is provided AS-IS with ABSOLUTELY NO WARRANTY\n "
            " under the GNU general public license\n\n" );
#ifdef USE_OPENMP
    printf( "multithreaded using openMP\n");
#endif
    printf( "calculate STEM images\n\n");

  /* read in detector sensitivity map, CZ 3-28-16*/

      Sensitivity = (char**)malloc2D(1,NCMAX, sizeof(char),"Sensitivity");
      printf("Name of detector file is : %s\n",Sensitivity[0]);
      scanf("%s", Sensitivity[0]);
      topenFloat(Sensitivity[0]);
      tsize( &senx, &seny, nbits, &samples);
      senx = (int) senx;
      seny = (int) seny;
      sensmap = (float**)malloc2D(senx, seny, sizeof(float), "sensmap");
      sensmap_rotated = (float**)malloc2D(senx*2, seny*2, sizeof(float),"sensmap_rotated");
      sensmap_rotated_save = (float**)malloc2D(senx*2, seny*2, sizeof(float),"sensmap_rotated_save");
      param_detector  = (float*) malloc1D( 5, sizeof(float), "param_detector" );
      param = (float*) malloc1D( NPARAM, sizeof(float), "param" );   /* ??? ejk 19-jul-98 */
      treadFloatPix(sensmap, senx, seny, &npixx, datetime, param_detector);


  /*read in center of detector in the input map, image center used as default if
     input number is out of range CZ 4-15-16*/
      printf("Type in the center of detector in the input file in px\n");
      scanf("%d %d",&centerx, &centery);
      if ((centerx<0)||(centerx>senx)||(centery<0)||(centery>seny))
      {
        printf("Input center position out of range, image center used as default center\n");
        centerx = senx/2;
        centery = seny/2;
      }

      /*read in the rotation angle of sample compared with detector
      CZ 5-3-16 */
      pi = 4.0 * atan( 1.0 );
      val = pi/180;
      printf("Type in the sample rotation (counterclockwise as positive)\n");
      scanf("%lf",&rot_angle); 
      /*do rotation part when rotation angle is not zero CZ 6-15-16*/
      if (rot_angle!=0)
      {
        /*printf("val = %lf, rot_angle = %lf\n",val, rot_angle );
        printf("cos = %lf, sin = %lf\n", cos(rot_angle*val), sin(rot_angle*val));*/
        for (i = 0; i < senx*2; i++)
        {
          for (j = 0; j < seny*2; j++)
          {
            sensmap_rotated[i][j]=0;
          }
        }
        for (i = 0; i < senx; i++)
        {
          for (j = 0; j < seny; j++)
          {
            i_new_int = i - centerx;
            j_new_int = j - centery;
            i_new = i_new_int*cos(rot_angle*val)+j_new_int*sin(rot_angle*val);
            j_new = j_new_int*cos(rot_angle*val)-i_new_int*sin(rot_angle*val);
            /*printf("for i = %d, j = %d, i_new = %lf, j_new = %lf\n",i, j, i_new, j_new);*/
            i_new = i_new + centerx + senx/2;
            j_new = j_new + centery + seny/2;
            i_new_int = (int)floor(i_new);
            j_new_int = (int)floor(j_new);
            if (sensmap_rotated[i_new_int][j_new_int] == 0)
            {
              sensmap_rotated[i_new_int][j_new_int] = sensmap[i][j];
              sensmap_rotated_save[i_new_int][j_new_int] = sensmap[i][j]/1.5;
            }
            else
            {
              sensmap_rotated[i_new_int][j_new_int] = (sensmap_rotated[i_new_int][j_new_int]+sensmap[i][j])/2;
              sensmap_rotated_save[i_new_int][j_new_int] = sensmap_rotated[i_new_int][j_new_int]/1.5;
            }

            /*printf("rotated[%d][%d] = sensmap[%d][%d] = %f\n", i_new_int, j_new_int, i, j, sensmap[i][j]);*/
          }
        } 
        for (i = 1; i < (senx*2-1); i++)
        {
          for (j = 1; j < (seny*2-1); j++)
          {
            if (sensmap_rotated[i][j] == 0)
            {
              sensmap_rotated[i][j] = (sensmap_rotated[i-1][j-1]+sensmap_rotated[i-1][j+1]+sensmap_rotated[i+1][j-1]+sensmap_rotated[i+1][j+1])/4;
              sensmap_rotated_save[i][j] = (sensmap_rotated[i-1][j-1]+sensmap_rotated[i-1][j+1]+sensmap_rotated[i+1][j-1]+sensmap_rotated[i+1][j+1])/6;
            }
          }
        }
        for (i = 0; i < senx; i++)
        {
          for (j = 0; j < seny; j++)
          {
            sensmap[i][j] = sensmap_rotated[i+senx/2][j+seny/2];
          }
        }
      }
      /*printf("sensmap[256][256] = %f\n",sensmap[256][256] );*/

/*  get simulation options */

    aber = (double**)malloc2D(12, 4, sizeof(double), "Aberration array");

    printf("Name of file with input atomic "
           "potential in x,y,z format:\n");
    scanf("%500s", filein );

    printf("Replicate unit cell by NCELLX,NCELLY,NCELLZ :\n");
    scanf("%d %d %d", &ncellx, &ncelly, &ncellz);
    if( ncellx < 1 ) ncellx = 1;
    if( ncelly < 1 ) ncelly = 1;
    if( ncellz < 1 ) ncellz = 1;

/*  get more parameter etc */

    printf("STEM probe parameters, V0(kv), apert1,2(mrad) :\n");
    scanf("%lg %lg %lg", &keV, &apert1, &apert2);

    wavlen = wavelength( keV );
    printf("wavelength = %f Angstroms\n", wavlen);
    if( apert1 > apert2 ) {
       printf("Bad probe aperture specification.\n");
       printf("apert1 must be less than apert2.\n");
       printf("apert1=%f, apert2 = %f\n", apert1, apert2);
       exit( 0 );
    }

    printf("First order aberrations: C1 (nm), A1 (nm, deg): \n");
    scanf("%lg %lg %lg", &aber[0][0], &aber[1][0], &aber[1][1]);
    aber[0][1] = 0.0; /* no angle for C1 */
    aber[0][0] *= 10.0; /* nm to Angstroms */
    aber[1][0] *= 10.0;
    Df0 = aber[0][0]; /* store initial defocus for chromatic aberration spread later */

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


    do { 
       printf("Size of specimen transmission function"
          " Nx,Ny in pixels : \n");
       scanf("%d %d", &nx, &ny);
       nx2 = (int) powerof2( (long) nx );
       ny2 = (int) powerof2( (long) ny );
       if( (nx != nx2) || (ny != ny2) ) {
        printf(" Nx=%d, Ny=%d must be a power of 2, try again.\n",
               nx, ny);
        exit(1);
        }
    } while(  (nx != nx2) || (ny != ny2) );

    do { 
       printf("Size of probe wave function"
          " Nx,Ny in pixels : \n");
       scanf("%d %d", &nxprobe, &nyprobe);
       nx2 = (int) powerof2( (long) nxprobe );
       ny2 = (int) powerof2( (long) nyprobe );
       if( (nxprobe != nx2) || (nyprobe != ny2) ) {
        printf(" Nx=%d, Ny=%d must be a power of 2, try again.\n",
               nxprobe, nyprobe);
        exit(1);
        }
    } while(  (nxprobe != nx2) || (nyprobe != ny2) );


    printf("Crystal tilt x,y in mrad. :\n");
    scanf( "%lf %lf", &ctiltx, &ctilty );
    ctiltx = ctiltx * 0.001;
    ctilty = ctilty * 0.001;

    l1d = askYN("Do you want to calculate a 1D line scan");

    if( l1d == 1 ) {
        lxzimage = askYN("Do you want to save all depth information as xz image");
        nThick = 1;
    } else {
        do { printf("Number of thickness levels to save, including"
                " the end(>=1):\n");
         scanf( "%d", &nThick );
        } while (nThick <= 0);
        ThickSave = (double*) malloc1D( nThick, sizeof(double), "ThickSave");
        if( nThick > 1 ) {
           printf( "type thickness (in Ang.) of %d intermediate layers"
            " :\n", (nThick-1) ); 
           for( it=0; it<(nThick-1); it++) scanf( "%lf", &ThickSave[it] );
        }
    }

    printf("File name prefix to get output of STEM multislice result "
        "(no extension):\n");
    scanf("%500s", fileoutpre);

    do {    printf("Number of detector geometries (>=1):\n");
            scanf( "%d", &ndetect );
    } while (ndetect <= 0);

    almin  = (double*)  malloc1D( ndetect, sizeof(double), "almin" );
    almax  = (double*)  malloc1D( ndetect, sizeof(double), "almax" );
    probe_file = (char**) malloc2D(ndetect+1, NCMAX, sizeof( char ),"probe_file" );
    ratio = (float*) malloc1D( ndetect,sizeof(float), "ratio");

/*Calculate mrad/px ratio based on the input inner and outer angle, sens = 0.5 used as 
    detector border threshold*/

    for( idetect=0; idetect<ndetect; idetect++) 
    {
        printf("Detector %3d: Type, min,max angles(mrad)"
            " of collector :\n", idetect+1);
        scanf("%lg %lg",
            &almin[idetect], &almax[idetect] );
        OutCut = askYN("Is the outer collection angle limited by aperture?");
        if (OutCut != 1)
        {
          almax[idetect]==10000; /* if the outer collection angle is not determined by aperture, use a large enough outer angle */
        }
        printf("Name of file to save adapted sensitivity map for detector %d: \n",idetect+1);
        scanf("%s", probe_file[idetect+1]);
        probe_file[0] = "detector_rotated.tif";
        px_in_1 = 0;
        px_in_2 = 0;
        for (i = centery; i >0; i--)
        {
          temp_sens = sensmap[centerx][i];
          /*printf("sensmap[%d][%d] = %f \n", centerx, i, temp_sens);*/
          if ((temp_sens>0.05) && (px_in_1 == 0))
          {
            px_in_1 = i;
          }
        }

        for (i = centery; i < seny; i++)
        {
          temp_sens = sensmap[centerx][i];
          if ((temp_sens>0.05) && (px_in_2 == 0))
          {
            px_in_2 = i;
          }
        }
        px_in_1f = (float)px_in_1;
        px_in_2f = (float)px_in_2;
        /* temp scheme for ratio, use only inner angle input to determine the mrad/px ratio */
        ratio[idetect] = 2*almin[idetect]/(fabs(px_in_1f-(float)centery)+fabs(px_in_2f-(float)centery));
        
        printf("ratio of detector %d is %f mrad/px\n", idetect+1, ratio[idetect]);
        printf("px_in_if = %f , px_in_2f = %f\n", px_in_1f, px_in_2f);
    }  /* end for(idetect=.. */

    /*  calculate spatial frequencies for future use
        (one set for transmission function and one for probe
        wavefunction)
    NOTE: zero freg is in the bottom left corner and
        expands into all other corners - not in the center
        this is required for FFT - don't waste time rearranging

    remember : the x values are the same for both sets
    
    x2, y2 are not used anywhere else - just for freqn()
    
    */

    kxp  = (float*) malloc1D( nxprobe, sizeof(float), "kxp" );
    kyp  = (float*) malloc1D( nyprobe, sizeof(float), "kyp" );
    kxp2 = (float*) malloc1D( nxprobe, sizeof(float), "kxp2" );
    kyp2 = (float*) malloc1D( nyprobe, sizeof(float), "kyp2" );
    x2  = (float*) malloc1D( nx, sizeof(float), "x2" );
    y2  = (float*) malloc1D( ny, sizeof(float), "y2" );
    mradx = (float*)malloc1D(nxprobe, sizeof(float), "mradx");
    mrady = (float*)malloc1D(nyprobe, sizeof(float), "mrady"); 
    mradx_sens = (float*)malloc1D(senx, sizeof(float), "mradx_sens"); 
    mrady_sens = (float*)malloc1D(seny, sizeof(float), "mrady_sens"); 
    


    if( l1d == 1 ) {
        printf("xi, xf, yi, yf, nout :\n");
        scanf("%lg %lg %lg %lg %d", &xi, &xf, &yi, &yf, &nout);
        nprobes = nout;
    }else {
        printf("xi,xf,yi,yf, nxout,nyout :\n");
        scanf("%lg %lg %lg %lg %d %d",
             &xi, &xf, &yi, &yf, &nxout, &nyout);
        nprobes = nyout;
    }

    /*  remember that the slice thickness must be > atom size
        to use projected atomic potential */
    printf("Slice thickness (in Angstroms):\n");
    scanf("%lf", &deltaz );
    if( deltaz < 1.0 ) {
        printf("WARNING: this slice thickness is probably too thin"
            " for autostem to work properly.\n");
    }

    lwobble = askYN("Do you want to include thermal vibrations");
    if( lwobble == 1 ) {
        printf( "Type the temperature in degrees K:\n");
        scanf( "%g", &temperature );
        /* get random number seed from time if available 
            otherwise ask for a seed */
        printf( "Type number of configurations to average over:\n");
        scanf( "%d", &nwobble );
        if( nwobble < 1 ) nwobble = 1;
	printf("Type initial seed for random number generator:\n");
	scanf("%ld", &iseed);
	printf( "Random number seed initialized to %ld\n", iseed );
        printf( "Type source size (FWHM in Ang.):\n" );
        scanf( "%lf", &sourceFWHM );
	printf( "Type chromatic aberration coefficient (mm) and energy FWHM (eV): \n");
	scanf( "%lf  %lf", &Cc, &dE );
    } else {
        temperature = 0.0F;
        nwobble = 1;
        sourceFWHM = 0.0;
	Cc = 0.0;
	dE = 0.0;
    }
    /* convert FWHM to standard deviation 
            by dividing by 2*sqrt(2*ln(2)) */
    sourcesize = sourceFWHM / 2.354820045;

    /* convert Cc and dE into defocus spread, following Reimer */
    cDf_sigma = (Cc*1.0e7)*(dE/(keV*1.0e3))*( (1.0 + keV / 511.0)/(1 + keV/1022.0) ) / 2.354820045;
    if(lwobble == 1) {
      printf( "Chromatic defocus standard deviation is %g Ang. \n", cDf_sigma);
    }

    

    /*  read in atomic potential and specimen parameters
        and calculate specimen transmission function
        for a single slice in transr,i
    */

    timer = cputim();   /* get initial CPU time */
#ifdef USE_OPENMP
    walltimer = walltim();  /* wall time for opneMP */
#endif

    for( i=0; i<NPARAM; i++) param[i] = 0.0F;

/*  calculate relativistic factor and electron wavelength */

    wavlen = (float) wavelength( keV );
    printf("electron wavelength = %g Angstroms\n", wavlen);

/*  read in specimen coordinates and scattering factors */

    natom = ReadXYZcoord( filein, ncellx, ncelly, ncellz,
        &ax, &by, &cz, &Znum, &xa, &ya, &za, &occ, &wobble,
        description, NCMAX );

    printf("%d atomic coordinates read in\n", natom );
    printf("%s", description );

    printf("Lattice constant a,b,c = %12.4f, %12.4f, %12.4f\n", ax,by,cz);

    freqn( kxp, kxp2, x2, nxprobe, ax*((double)nxprobe)/nx );
    freqn( kyp, kyp2, y2, nyprobe, by*((double)nyprobe)/ny );

/*Initialize value in mrad list, both for detector sens and beam*/


        for (i = 0; i < nxprobe; i++)
        {
            mradx[i] = kxp[i]/0.001*wavlen;
            mrady[i] = kyp[i]/0.001*wavlen;
            /*printf("OK, mradx[%d] = %f, kxp[%d] = %f\n", i, mradx[i], i, kxp[i]);*/
        }/*mradx and mrady value same for all detectors*/
        
        sensmap_probe = (float***)malloc(ndetect* sizeof(float**));
        sensmap_save = (float***)malloc(ndetect* sizeof(float**));

        for (idetect = 0; idetect < ndetect ; idetect++)
        {
            sensmap_probe[idetect] = (float**)malloc2D(nxprobe, nyprobe, sizeof(float),"sensmap_probe");
            sensmap_save[idetect] = (float**)malloc2D(nxprobe, nyprobe, sizeof(float),"sensmap_save");
        }

        /* find value for each pixel in sensmap_probe, which is used to account for detector sensitivity */

        printf("senx = %d seny = %d\n",senx, seny );
        printf("nxprobe = %d nyprobe = %d\n",nxprobe, nyprobe );
      for (idetect = 0; idetect<ndetect; idetect++)
      {
        for (i = 0; i < senx; i++)
        {
            mradx_sens[i] = (i-centerx)*ratio[idetect];
            mrady_sens[i] = (i-centery)*ratio[idetect];
        }
        /*Find the closest value for each pixel on sensmap_probe*/
        for (ix = 0; ix < nxprobe; ix++)
        {
            for (iy = 0; iy < nyprobe; iy++)
            {
                mindifference = 100000;
                minIdx = 0;
                reference = mradx[ix];

                for (j = 0; j < senx; j++)
                {
                    difference = fabs(reference - mradx_sens[j]);
                    if (difference < mindifference)
                    {
                        mindifference = difference;
                        minIdx = j;
                    }
                }

                mindifference = 100000;
                minIdy = 0;
                reference = mrady[iy];

                for (j = 0; j < seny; j++)
                {
                    difference = fabs(reference - mrady_sens[j]);
                    if (difference < mindifference)
                    {
                        mindifference = difference;
                        minIdy = j;
                    }
                }  
            sensmap_probe[idetect][ix][iy] = sensmap[minIdx][minIdy]; 
            sensmap_save[idetect][ix][iy] = sensmap_probe[idetect][ix][iy]/1.5;
            }
        }
        param[1] = 1;
        param[3] = 0;
        tcreateFloatPixFile( probe_file[idetect+1], sensmap_save[idetect], (long)nxprobe, (long)nyprobe, 1, param);
      }
      if(rot_angle!=0)
      {
          tcreateFloatPixFile( probe_file[0], sensmap_rotated_save, senx*2, seny*2, 1, param);       
      }
        /*param[1] = 1;
        param[3] = 0;*/
        /*printf("param[1] = %d param[2] = %d param[3] = %d param[4] = %d\n",param[1], param[2], param[3], param[4] );*/
        /*tcreateFloatPixFile( probe_file, sensmap_probe, (long)nxprobe,
                    (long) nyprobe, 1, param);*/  /*won't create adapted sensmap here*/

    /* calculate thickness levels to save (1D mode) or check range (2D mode) */
    if( l1d == 1 ) {
        if( lxzimage == 1 ) {
            /* save all thickness levels  */
            nThick = (int) ( cz/deltaz + 0.5 );
            ThickSave = (double*) malloc1D( nThick, sizeof(double), "ThickSave");
            for( it=0; it<nThick; it++) {
                ThickSave[it] = deltaz*(it+1);
            }
        } else {
            nThick = 1;
            ThickSave = (double*) malloc1D( nThick, sizeof(double), "ThickSave");
            ThickSave[0] = cz;
        }
        printf( "save up to %d thickness levels\n", nThick );  /* diagnostic */
    } else {
        ThickSave[nThick-1] = cz;  /*  always save the last level */
        for( it=0; it<(nThick-1); it++) 
        if( (ThickSave[it] < 0.0) || (ThickSave[it] > cz) ) {
            printf("Bad thickness level = %g A, allowed range= "
                "0.0 to %f A\n", ThickSave[it], cz );
            exit( 0 );
        }
    }  /* end if( l1d == ... */

    if( lwobble == 0 ) {
        printf( "Sorting atoms by depth...\n");
        sortByZ( xa, ya, za, occ, Znum, natom );
    }
    /* to add random offsets */
    xa2 = (float*) malloc1D( natom, sizeof(float),  "xa2" );
    ya2 = (float*) malloc1D( natom, sizeof(float), "ya2" );
    za2 = (float*) malloc1D( natom, sizeof(float), "za2" );
    Znum2 = (int*) malloc1D( natom, sizeof(int), "Znum2" );
    occ2 = (float*) malloc1D( natom, sizeof(float), "occ2" );

    /*  calculate the total specimen volume and echo */
    xmin = xmax = xa[0];
    ymin = ymax = ya[0];
    zmin = zmax = za[0];
    wmin = wmax = wobble[0];

    for( i=0; i<natom; i++) {
        if( xa[i] < xmin ) xmin = xa[i];
        if( xa[i] > xmax ) xmax = xa[i];
        if( ya[i] < ymin ) ymin = ya[i];
        if( ya[i] > ymax ) ymax = ya[i];
        if( za[i] < zmin ) zmin = za[i];
        if( za[i] > zmax ) zmax = za[i];
        if( wobble[i] < wmin ) wmin = wobble[i];
        if( wobble[i] > wmax ) wmax = wobble[i];
    }
    printf("Total specimen range is\n %g to %g in x\n"
           " %g to %g in y\n %g to %g in z\n", xmin, xmax,
           ymin, ymax, zmin, zmax );
    if( lwobble == 1 )
        printf("Range of thermal rms displacements (300K) = %g to %g\n",
            wmin, wmax );

    /*  check for valid scan coordinates  */

    if( (xi < 0.0) || (xi > ax) ||
        (xf < 0.0) || (xf > ax) ||
        (yi < 0.0) || (yi > by) ||
        (yf < 0.0) || (yf > by) ) {
            printf("WARNING: Coordinates out of range; will be made periodic.\n");
            printf("xi,xf,yi,yf= %f, %f, %f, %f\n", xi, xf, yi, yf );
    }

    /*  check that requested probe size is not bigger 
        than transmission function size (or too small)
    */
    if( (nxprobe > nx) || (nxprobe < 2) ) {
        nxprobe = nx;
        printf("Probe size reset to nx = %d\n", nxprobe);
    }

    if( (nyprobe > ny) || (nyprobe < 2) ) {
        nyprobe = ny;
        printf("probe size reset to ny = %d\n", nyprobe);
    }

    kx  = (float*) malloc1D( nx, sizeof(float), "kx" );
    ky  = (float*) malloc1D( ny, sizeof(float), "ky" );
    kx2 = (float*) malloc1D( nx, sizeof(float), "kx2" );
    ky2 = (float*) malloc1D( ny, sizeof(float), "ky2" );

    freqn( kx, kx2, x2, nx, ax );
    freqn( ky, ky2, y2, ny, by );



    free( x2 );
    free( y2 );

/* impose anti-aliasing bandwidth limit on transmission functions */

    sum = ((double)nx)/(2.0*ax);
    k2maxp = ((double)ny)/(2.0*by);
    if( sum < k2maxp ) k2maxp = sum;
    k2maxp= BW * k2maxp;
    printf("Bandwidth limited to a real space resolution of %f Angstroms\n",
                     1.0F/k2maxp);
    printf("   (= %.2f mrad)  for symmetrical anti-aliasing.\n",
         wavlen*k2maxp*1000.0F);
    k2maxp = k2maxp * k2maxp;

/*  allocate some more arrays and initialize wavefunction */

    transr = (float**) malloc2D( nx, ny, sizeof(float), "transr" );
    transi = (float**) malloc2D( nx, ny, sizeof(float), "transi" );

    /* calculate propagator functions with probe sample size
        impose anti-aliasing bandwidth limit
    */
    propxr = (float*) malloc1D( nxprobe, sizeof(float), "propxr" );
    propxi = (float*) malloc1D( nxprobe, sizeof(float), "propxi" );
    propyr = (float*) malloc1D( nyprobe, sizeof(float), "propyr" );
    propyi = (float*) malloc1D( nyprobe, sizeof(float), "propyi" );

    tctx = 2.0 * tan(ctiltx);
    tcty = 2.0 * tan(ctilty);

    scale = pi * deltaz;
    for( ix=0; ix<nxprobe; ix++) {
        w = scale * ( kxp2[ix] * wavlen - kxp[ix]*tctx );
        propxr[ix]= (float)  cos(w);
        propxi[ix]= (float) -sin(w);
    }

    for( iy=0; iy<nyprobe; iy++) {
        w = scale * ( kyp2[iy] * wavlen - kyp[iy]*tcty );
        propyr[iy]= (float)  cos(w);
        propyi[iy]= (float) -sin(w);
    }

    nbeamp = 0;
    for( iy=0; iy<nyprobe; iy++)
    for( ix=0; ix<nxprobe; ix++) {
        if( (kyp2[iy] + kxp2[ix]) < k2maxp ) nbeamp++;
    }

    printf("Number of symmetrical anti-aliasing "
           "beams in probe = %ld\n", nbeamp);

    /*  convert aperture dimensions */

    k2min = (double*) malloc1D( ndetect, sizeof(double), "k2min" );
    k2max = (double*) malloc1D( ndetect, sizeof(double), "k2max" );

    for( idetect=0; idetect<ndetect; idetect++) {
        k2max[idetect] = 0.001 * almax[idetect]/wavlen;
        k2max[idetect] = k2max[idetect] * k2max[idetect];
        k2min[idetect] = 0.001 * almin[idetect]/wavlen;
        k2min[idetect] = k2min[idetect] * k2min[idetect];
    }

    /*  init the min/max record of total integrated intensity */

    totmin =  10.0;
    totmax = -10.0;
    detect  = (double***) malloc3D( nThick, ndetect, nprobes,
        sizeof(double), "detect" );
    detect2 = (double***) malloc3D( nThick, ndetect, nprobes,
        sizeof(double), "detect2" );
    sums = (double*) malloc1D( nprobes, sizeof(double), "sums" ); 
    rmin = (float**) malloc2D( nThick, ndetect, sizeof(float), "rmin" );
    rmax = (float**) malloc2D( nThick, ndetect, sizeof(float), "rmax" );

    /* allocate probe wave function arrays*/
    prober = (float***) malloc3D( nprobes, nxprobe, nyprobe,
                sizeof(float), "prober" );
    probei = (float***) malloc3D( nprobes, nxprobe, nyprobe,
                sizeof(float), "probei" );

    fflush(NULL);

/* ------------- start here for a full image output -------------- */
/*
  do one whole line at once NOT the whole image (which may be huge)
*/
    if( l1d == 0 ) {
       printf("output file size in pixels is %d x %d\n",
          nxout, nyout );
       if( nprobes != nyout ) {
           printf( "Error, nprobes=%d must be the same as"
               "nyout=%d, in image mode.\n", nprobes, nyout );
           exit( 0 );
       }
       /* double up first index to mimic a 4D array */
       pixr = (float***) malloc3D( ndetect*nThick, nxout, nyout,
           sizeof(float), "pixr"  );
       for( i=0; i<(nThick*ndetect); i++) {
           for( ix=0; ix<nxout; ix++)
           for( iy=0; iy<nyout; iy++)
            pixr[i][ix][iy] = 0.0F;
       }

       /*  iterate the multislice algorithm proper for each
        position of the focused probe */

       dx = (xf-xi)/((double)(nxout-1));
       dy = (yf-yi)/((double)(nyout-1));
       x = (double*) malloc1D( nprobes, sizeof(double), "x" );
       y = (double*) malloc1D( nprobes, sizeof(double), "y" );

        /*  add random thermal displacements 
            scaled by temperature if requested 
        remember that initial wobble is at 300K for
            each direction */

        for( iwobble=0; iwobble<nwobble; iwobble++) {
            if( lwobble == 1 ){
                scale = (float) sqrt(temperature/300.0) ;
                for( i=0; i<natom; i++) {
                xa2[i] = xa[i] + 
                    (float)(wobble[i]*rangauss(&iseed)*scale);
                ya2[i] = ya[i] + 
                    (float)(wobble[i]*rangauss(&iseed)*scale);
                za2[i] = za[i] + 
                        (float)(wobble[i]*rangauss(&iseed)*scale);
                occ2[i] = occ[i];
                Znum2[i] = Znum[i];
                }
                sortByZ( xa2, ya2, za2, occ2, Znum2, natom );
                printf("\n\n\nconfiguration # %d\n", iwobble+1 );
                printf( "The new range of z is %g to %g\n",
                    za2[0], za2[natom-1] );
		fflush(NULL);
            } else for( i=0; i<natom; i++) {
                xa2[i] = xa[i];
                ya2[i] = ya[i];
                za2[i] = za[i];
                occ2[i] = occ[i];
                Znum2[i] = Znum[i];
            }
            zmin = za2[0];  /* reset zmin/max after wobble */
            zmax = za2[natom-1];
    
            for( ix=0; ix<nxout; ix++) {
    
                for( iy=0; iy<nyout; iy++) {
                    x[iy] = xi + dx * ((double) ix)
                        + sourcesize * rangauss(&iseed);
                    y[iy] = yi + dy * ((double) iy)
                        + sourcesize * rangauss(&iseed);
                    x[iy] = periodic( x[iy], ax );   /* put back in supercell */
                    y[iy] = periodic( y[iy], by );  /* if necessary */
                }
    
                STEMsignals( x, y, nyout, detect, ndetect, 
                ThickSave, nThick, sums, sensmap_probe );
                for( iy=0; iy<nyout; iy++) {
                    if( sums[iy] < totmin ) totmin = sums[iy];
                    if( sums[iy] > totmax ) totmax = sums[iy];
                    for( it=0; it<nThick; it++){
                        for( idetect=0; idetect<ndetect; idetect++)
                        pixr[idetect + it*ndetect][ix][iy] += (float)
                            (detect[it][idetect][iy]/((double)nwobble));
                    }
                    if( sums[iy] < 0.9) 
                        printf("Warning integrated intensity too small, = "
                        "%g at x,y= %g, %g\n", sums[iy], x[iy], y[iy] );
                    if( sums[iy] > 1.1) 
                        printf("Warning integrated intensity too large, = "
                        "%g at x,y= %g, %g\n", sums[iy], x[iy], y[iy] );
                }
            
            } /* end for(ix...) */
    
        } /* end for(iwobble... ) */

        /*  output data files  */
        for( it=0; it<nThick; it++)
        for( i=0; i<ndetect; i++) {
            rmin[it][i] = rmax[it][i] = pixr[i+it*ndetect][0][0];
            for( ix=0; ix<nxout; ix++)
            for( iy=0; iy<nyout; iy++) {
                temp = pixr[i+it*ndetect][ix][iy];
                if( temp < rmin[it][i] )rmin[it][i] = (float) temp;
                if( temp > rmax[it][i] )rmax[it][i] = (float) temp;
            }
        }

        /*  directory file listing parameters for each image file */
        sprintf( fileout, "%s.txt", fileoutpre );
        fp = fopen( fileout, "w" );  /* changed from "w+" for Condor compatbility */
        if( fp == NULL ) {
            printf("Cannot open output file %s.\n", fileout );
            exit( 0 );
        }
     
       fprintf(fp, "C\n");
       fprintf(fp, "C   output of autostem version %s\n", version);
       fprintf(fp, "C\n");
       fprintf(fp, "C   nslice= %d\n", nslice);
       fprintf(fp,"deltaz= %g, file in= %s\n", deltaz, filein );
       fprintf(fp, "V0= %g, Apert= %g mrad to %g mrad\n", keV, apert1, apert2 );
       fprintf(fp, "Aberrations:\n");
       fprintf(fp, "C1 = %g nm\n", aber[0][0]/10.0);
       fprintf(fp, "A1 = %g nm, %g deg\n", aber[1][0]/10.0, aber[1][1]*180.0/pi);
       fprintf(fp, "A2 = %g nm, %g deg\n", aber[2][0]/10.0, aber[2][1]*180.0/pi);
       fprintf(fp, "B2 = %g nm, %g deg\n", aber[3][0]/10.0, aber[3][1]*180.0/pi);
       fprintf(fp, "C3 = %g um\n", aber[4][0]/1.0e4);
       fprintf(fp, "A3 = %g um, %g deg\n", aber[5][0]/1.0e4, aber[5][1]*180.0/pi);
       fprintf(fp, "S3 = %g um, %g deg\n", aber[6][0]/1.0e4, aber[6][1]*180.0/pi);
       fprintf(fp, "A4 = %g um, %g deg\n", aber[7][0]/1.0e4, aber[7][1]*180.0/pi);
       fprintf(fp, "D4 = %g um, %g deg\n", aber[8][0]/1.0e4, aber[8][1]*180.0/pi);
       fprintf(fp, "B4 = %g um, %g deg\n", aber[9][0]/1.0e4, aber[9][1]*180.0/pi);
       fprintf(fp, "C5 = %g mm\n", aber[10][0]/1.0e7);
       fprintf(fp, "A5 = %g mm, %g deg\n", aber[11][0]/1.0e7, aber[11][1]*180.0/pi);
       fprintf(fp, "Crystal tilt x,y= %lg, %lg\n", ctiltx,ctilty);

       for(  idetect=0; idetect<ndetect; idetect++) {
            fprintf(fp, "Detector %d, Almin= %g mrad, Almax= %g mrad\n",
                idetect, almin[idetect], almax[idetect] );
       }

       fprintf(fp, "ax= %g A, by= %g A, cz= %g A\n", ax,by,cz);
       fprintf(fp, "Number of symmetrical anti-aliasing "
           "beams in probe wave function= %ld\n", nbeamp );
       fprintf(fp, "with a resolution (in Angstroms) = %g\n",
           1.0/sqrt(k2maxp) );
       if( lwobble == 1 ) {
          fprintf( fp,
           "Number of thermal configurations = %d\n", nwobble );
              fprintf( fp, "Source size = %g Ang. (FWHM) \n", sourceFWHM );
	      fprintf( fp, "Cc = %g mm.  dE = %g eV (FWHM).  Defocus standard dev. = %g Ang. \n", 
		       Cc, dE, cDf_sigma);
       }
       fprintf( fp, "\n" );

        /*  store params plus min and max */
        param[pIMAX]    = 0.0F;
        param[pIMIN]    = 0.0F;
        param[pXCTILT]  = (float) ctiltx;
        param[pYCTILT]  = (float) ctilty;
        param[pDEFOCUS] = (float) aber[0][0];
        param[pDX]      = (float) dx;
        param[pDY]      = (float) dy;
        param[pENERGY]  = (float) keV;
        param[pOAPERT]  = (float) apert2;
        param[pCS]      = (float) aber[4][0];
        param[pWAVEL]   = (float) wavlen;
        param[pNSLICES] = (float) -1.0; /* ??? ejk 19-jul-1998 */
        param[35]       = (float) aber[10][0];  /* ??? ejk 3-mar-2004 */

        for( it=0; it<nThick; it++)
        for( i=0; i<ndetect; i++) {
            sprintf( fileout, "%s%03d%03d.tif", fileoutpre, i, it );
            printf("%s: output pix range : %g to %g\n", fileout, 
                rmin[it][i], rmax[it][i]);
            param[pRMAX] = rmax[it][i];
            param[pRMIN] = rmin[it][i];
            param[pMINDET] = (float) ( almin[i] * 0.001 );
            param[pMAXDET] = (float) ( almax[i] * 0.001 );
            if( tcreateFloatPixFile( fileout, pixr[i+it*ndetect],
                (long) nxout, (long) nyout, 1, param ) != 1 ) {
                printf("Cannot write output file %s.\n", fileout );
            }
            fprintf(fp,"file: %s, detector= %g to %g mrad, "
                "thicknes= %g A, range= %g to %g\n", fileout,
                almin[i], almax[i], ThickSave[it], rmin[it][i], rmax[it][i]);
        }
        fclose( fp );

/* ------------- start here for 1d line scan ---------------- */

    } else if ( l1d == 1 ) {

       dx = (xf-xi)/((double)(nout-1));
       dy = (yf-yi)/((double)(nout-1));
       x = (double*) malloc1D( nprobes, sizeof(double), "x" );
       y = (double*) malloc1D( nprobes, sizeof(double), "y" );
            for( ip=0; ip<nout; ip++) {
                for( it=0; it<nThick; it++)
                for( idetect=0; idetect<ndetect; idetect++)
                    detect[it][idetect][ip] = 0.0;
       }

       /*  add random thermal displacements scaled by temperature
            if requested 
        remember that initial wobble is at 300K for each direction */

       for( iwobble=0; iwobble<nwobble; iwobble++) {
    
            if( lwobble == 1 ){
                scale = (float) sqrt(temperature/300.0) ;
                for( i=0; i<natom; i++) {
                    xa2[i] = xa[i] + 
                        (float)(wobble[i]*rangauss(&iseed)*scale);
                    ya2[i] = ya[i] + 
                        (float)(wobble[i]*rangauss(&iseed)*scale);
                    za2[i] = za[i] + 
                        (float)(wobble[i]*rangauss(&iseed)*scale);
                    occ2[i] = occ[i];
                    Znum2[i] = Znum[i];
                }
                printf("configuration # %d\n", iwobble+1 );
                sortByZ( xa2, ya2, za2, occ2, Znum2, natom );
                printf( "The new range of z is %g to %g\n",
                    za2[0], za2[natom-1] );
		fflush(NULL);
            } else for( i=0; i<natom; i++) {
                xa2[i] = xa[i];
                ya2[i] = ya[i];
                za2[i] = za[i];
                occ2[i] = occ[i];
                Znum2[i] = Znum[i];
            }
            zmin = za2[0];      /* reset zmin/max after wobble */
            zmax = za2[natom-1];
            for( ip=0; ip<nout; ip++) {
                x[ip] = xi + dx * ((double)ip)
                            + sourcesize * rangauss(&iseed);
                y[ip] = yi + dy * ((double)ip)
                            + sourcesize * rangauss(&iseed);
                x[ip] = periodic( x[ip], ax );   /* put back in supercell */
                y[ip] = periodic( y[ip], by );  /* if necessary */
            }
           
            STEMsignals( x, y, nprobes, detect2, ndetect, 
                ThickSave, nThick, sums, sensmap_probe );
            for( ip=0; ip<nprobes; ip++) {
                if( sums[ip] < totmin ) totmin = sums[ip];
                if( sums[ip] > totmax ) totmax = sums[ip];
                for( it=0; it<nThick; it++){
                for( idetect=0; idetect<ndetect; idetect++)
                   detect[it][idetect][ip] += 
                        detect2[it][idetect][ip]/((double)nwobble);
                }
                if( sums[ip] < 0.9) 
                printf("Warning integrated intensity too small, = %g"
                    " at x,y= %g, %g\n", sums[ip], x[ip], y[ip] );
                if( sums[ip] > 1.1) 
                printf("Warning integrated intensity too large, = %g"
                    " at x,y= %g, %g\n", sums[ip], x[ip], y[ip] );
            }

       }  /* end for(iwobble... */

       /* ------ first output text data ---------------- */
       sprintf( fileout, "%s.dat", fileoutpre );
       printf("output file= %s\n", fileout);

       fp = fopen( fileout, "w" ); /* changed from "w+" for Condor compatibility 03/08/10 pmv */
       if( fp == NULL ) {
           printf("Cannot open output file %s.\n", fileout );
             exit( 0 );
       }

       fprintf(fp, "C\n");
       fprintf(fp, "C   output of autostem version %s\n", version);
       fprintf(fp, "C\n");
       fprintf(fp, "C   nslice= %d\n", nslice);
       fprintf(fp,"deltaz= %g, file in= %s\n", deltaz, filein );
       fprintf(fp, "V0= %g, Apert= %g mrad to %g mrad\n", keV, apert1, apert2 );
       fprintf(fp, "Aberrations:\n");
       fprintf(fp, "C1 = %g nm\n", aber[0][0]/10.0);
       fprintf(fp, "A1 = %g nm, %g deg\n", aber[1][0]/10.0, aber[1][1]*180.0/pi);
       fprintf(fp, "A2 = %g nm, %g deg\n", aber[2][0]/10.0, aber[2][1]*180.0/pi);
       fprintf(fp, "B2 = %g nm, %g deg\n", aber[3][0]/10.0, aber[3][1]*180.0/pi);
       fprintf(fp, "C3 = %g um\n", aber[4][0]/1.0e4);
       fprintf(fp, "A3 = %g um, %g deg\n", aber[5][0]/1.0e4, aber[5][1]*180.0/pi);
       fprintf(fp, "S3 = %g um, %g deg\n", aber[6][0]/1.0e4, aber[6][1]*180.0/pi);
       fprintf(fp, "A4 = %g um, %g deg\n", aber[7][0]/1.0e4, aber[7][1]*180.0/pi);
       fprintf(fp, "D4 = %g um, %g deg\n", aber[8][0]/1.0e4, aber[8][1]*180.0/pi);
       fprintf(fp, "B4 = %g um, %g deg\n", aber[9][0]/1.0e4, aber[9][1]*180.0/pi);
       fprintf(fp, "C5 = %g mm\n", aber[10][0]/1.0e7);
       fprintf(fp, "A5 = %g mm, %g deg\n", aber[11][0]/1.0e7, aber[11][1]*180.0/pi);
       fprintf(fp, "Crystal tilt x,y= %lg, %lg\n", ctiltx,ctilty);

       for(  idetect=0; idetect<ndetect; idetect++) {
            fprintf(fp, "Detector %d, Almin= %g mrad, Almax= %g mrad\n",
                idetect, almin[idetect], almax[idetect] );
       }

       fprintf(fp, "ax= %g A, by= %g A, cz= %g A\n", ax,by,cz);
       fprintf(fp, "Number of symmetrical anti-aliasing "
           "beams in probe wave function= %ld\n", nbeamp );
       fprintf(fp, "with a resolution (in Angstroms) = %g\n",
           1.0/sqrt(k2maxp) );
       if( lwobble == 1 ) {
              fprintf( fp,
               "Number of thermal configurations = %d\n", nwobble );
                  fprintf( fp, "Source size = %g Ang. (FWHM) \n", sourceFWHM );
		  fprintf( fp, "Cc = %g mm.  dE = %g eV (FWHM).  Defocus standard dev. = %g Ang. \n", 
			   Cc, dE, cDf_sigma);
       }
       fprintf(fp, "C     x      y     signal\n");

       for( ip=0; ip<nprobes; ip++) {
           /* recalculate mean x,y without source size wobble */
           x[ip] = xi + dx * ((double)ip);
           y[ip] = yi + dy * ((double)ip);
           fprintf(fp, "%14.7g %14.7g", x[ip], y[ip]);
           for(i=0; i<ndetect; i++) 
               fprintf(fp, "%14.7g", detect[nThick-1][i][ip] );
           fprintf(fp, "\n");
       }

       fclose( fp );

       /* ------ next output xz image data ---------------- */

       if( lxzimage == 1 ) {
    
            /*  directory file listing parameters for each image file */
            sprintf( fileout, "%s.txt", fileoutpre );
            fp = fopen( fileout, "w" ); /* changed from "w+" for Condor compatibility 03/08/10 pmv*/
            if( fp == NULL ) {
                printf("Cannot open output file %s.\n", fileout );
                exit( 0 );
            }
    
            /*  store params plus min and max */
            param[pIMAX]    = 0.0F;
            param[pIMIN]    = 0.0F;
            param[pXCTILT]  = (float) ctiltx;
            param[pYCTILT]  = (float) ctilty;
            param[pDEFOCUS] = (float) aber[0][0];
            param[pDX]  = (float) dx;
            param[pDY]  = (float) dy;
            param[pENERGY]  = (float) keV;
            param[pOAPERT]  = (float) apert2;
            param[pCS]  = (float) aber[4][0];
            param[pWAVEL]   = (float) wavlen;
            param[pNSLICES] = (float) -1.0; /* ??? ejk 19-jul-1998 */
            param[35]   = (float) aber[10][0];  /* ??? ejk 3-mar-2004 */
    
            pixout = (float**) malloc2D( nprobes, nThick, 
                    sizeof(float), "pixout"  );
    
            for( idetect=0; idetect<ndetect; idetect++){
                sprintf( fileout, "%s%03d.tif", fileoutpre, idetect );
                printf("output file= %s\n", fileout);
    
                /* convert to float and fix pixel order */
                pmin = pmax = (float) detect[0][idetect][0];
                for( ix=0; ix<nprobes; ix++)
                for( iy=0; iy<nThick; iy++) {
                    temp = pixout[ix][iy] = (float) detect[iy][idetect][ix];
                    if( temp < pmin )pmin = temp;
                    if( temp > pmax )pmax = temp;
                }
    
                printf("%s: output pix range : %g to %g\n", fileout, pmin, pmax);
                param[pRMAX] = pmax;
                param[pRMIN] = pmin;
                param[pMINDET] = (float) ( almin[i] * 0.001 );
                param[pMAXDET] = (float) ( almax[i] * 0.001 );
                if( tcreateFloatPixFile( fileout, pixout,
                    (long) nprobes, (long) nThick, 1, param ) != 1 ) {
                        printf("Cannot write output file %s.\n", fileout );
                }
                fprintf(fp,"file: %s, detector= %g to %g mrad, range= %g to %g\n",
                fileout, almin[idetect], almax[idetect], pmin, pmax);
            }
            fclose( fp );

       }  /* end if( lxzimage==1... */

    } /* end if( l1d.. ) */

    printf("Number of symmetrical anti-aliasing "
           "beams in trans. function = %ld\n", nbeamt);

    /*  echo min/max of total integrated intensity */
    printf("The total integrated intensity range was:\n");
    printf("   %g to %g\n\n",  totmin, totmax );

    printf("CPU time = %g sec.\n", cputim()-timer);
#ifdef USE_OPENMP
    printf("wall time = %g sec.\n", walltim() - walltimer);
#endif

    return( 0 );

}  /* end main() */

/*------------------------ periodic() ---------------------*/
/*
     make probe positions periodic in the supercell
     in case some wobble off the edge with source size of user excess

    pos = input position (x or y);
    size = supercell size ( 0 to size)

    return positive value  0 <= x < size
*/
double periodic( double pos, double size )
{
    double x=pos;
    while( x < 0 ) x += size;
    x = fmod( x, size );
    return( x );
}

/*------------------------ STEMsignals() ---------------------*/
/*

  NOTE: this is NOT the same as STEMsignal() in stemslice.c

  subroutine to calculate the stem signal arising from a given
  probe position

  iterate the multislice algorithm proper for each position of
  the focused probe

  This version uses massive amounts of memory to avoid
  recalculating the transmission functions more than necessary

   note zero freg is in the bottom left corner and
     expands into all other corners - not in the center
     this is required for fft - don't waste time rearranging

     change to propagate thru whole unit cell not just
     range of atoms on 14-sep-2006 ejk

  x[],y[]     = real positions of the incident probe
  npos        = int number of positions
  detect[][][]= real array to get signal into each detector
            for each probe position and thickness
  ndetect     = number of detector geometries
  ThickSave[] = thicknesses at which to save data (other than the last)
  nThick      = number of thickness levels (including the last)
  sum         = real total integrated intensity
  
  the assumed global variables are:
  
  nxprobe,nyprobe   = int size of probe wavefunction in pixels
  nx,ny         = int size of transmission function in pixels
  layer[]       = int array with slice layer indecies
  prober[][], probei[][] = float real,image probe wavefunction
  transr[][], transi[][] = float real,imag transmission function
  propxr[][], propxi[][] = float real,imag propagator vs x
  propyr[][], propyi[][] = float real,imag propagator vs y
  ax,by,cz      = float unit cell size in Angs
  kxp[], kyp[]      = float spatial frequencies vs x, y
  kxp2[], kyp2[]    = float square of kxp[], kyp[]
  apert1, apert2    = double min,max objective aperture (mrad)
  k2maxp            = double max spatial freq of probe squared
  pi                = double constant PI
  wavlen            = double electron wavelength in Angs
  aber              = array of all aberrations and angles
  Df0               = mean defocus, entered by user
  cDf_sigma         = chromatic aberration defocus distribtion standard deviation

  xa[],ya[],za[]    = atom coordinates
  occ[]         = atomic occupancy
  Znum[]        = atomic numbers
  natom         = number of atoms
  deltaz        = slice thickness
  v0            = beam energy
  nbeamt        = number of beams in transmission function
  zmin, zmax        = range of z coord. of the atoms
  nslice        = number of slices
  
    NOTE:  too many thing come in as globals, but...

*/
void STEMsignals( double x[], double y[], int npos,
         double ***detect, int ndetect,
         double ThickSave[], int nThick, double sum[], float ***sensmap_probe )
{
    int ix, iy, ixt, iyt, idetect, *ixoff, *iyoff, ixmid, iymid;
    int istart, na, ip, i, it;

    long nxprobel, nyprobel, nxl, nyl;

    float scale, prr, pri, tr, ti, *trr, *tri;

    double xoff,yoff, k2maxa, k2maxb, chi, w, k2, phirms;
    double sum0, delta, zslice, totalz;

/*   make sure x,y are ok */

    for( ip=0; ip<npos; ip++) {
        if( (x[ip] < 0.0) || (x[ip] > ax) ||
            (y[ip] < 0.0) || (y[ip] > by) ) {
            sum[ip] = -1.2345;
            printf("bad x=%f,y=%f in STEMsignals()\n", x[ip], y[ip]);
            return;
        }
    }

    /*  generate a probe at position x,y which must be inside the
        0->ax and 0->by region
        normalize to have unit integrated intensity
        probe position (x,y) = (0,0) = lower left corner

    NOTE: The probe wavefunction is nxprobe x nyprobe which may
        be smaller than the transmission function.
        Position the probe near the center of the nxprobe x nyprobe
        region because it does not wrap around properly
        (i.e. manually wrap around in the nx x ny region if
         the probe is near a boundary of the nx x ny region))
    */

    ixmid = nxprobe/2;
    iymid = nyprobe/2;
    k2maxa = apert1 * 0.001/wavlen;
    k2maxa = k2maxa *k2maxa;
    k2maxb = apert2 * 0.001/wavlen;
    k2maxb = k2maxb * k2maxb;

    ixoff = (int*) malloc1D( npos, sizeof(int), "ixoff" );
    iyoff = (int*) malloc1D( npos, sizeof(int), "iyoff" );

    /*  calculate all of the probe wave functions at once
        to reuse the transmission functions which take a long
        time to calculate*/

/*  paralleling this loop has little effect */
/*#pragma omp parallel for private(ix,iy,xoff,yoff,sum0,k2,w,phi,chi,scale,tr,ti) */
    for( ip=0; ip<npos; ip++) {
        ixoff[ip] = (int) floor( x[ip]*((double)nx) / ax ) - ixmid;
        xoff  = x[ip] - ax*((double)ixoff[ip])/((double)nx);

        iyoff[ip] = (int) floor( y[ip]*((double)ny) / by ) - iymid;
        yoff  = y[ip] - by*((double)iyoff[ip])/((double)ny);

	aber[0][0] = Df0 + cDf_sigma*rangauss(&iseed);  /* pmv chromatic aberration defocus */
	/*printf( "Chromatic defocus = %g Angstrom. \n", aber[0][0]);*/ 
      /*stop printing Cc due to large amount of output cz 7-15-16*/

        sum0 = 0.0;
        for( ix=0; ix<nxprobe; ix++)
        for( iy=0; iy<nyprobe; iy++) {
            k2 = kxp2[ix] + kyp2[iy];
           if( (k2 >= k2maxa) && (k2 <= k2maxb) ) {
	      /*chi = aphase(aber, wavlen, kxp[ix], kyp[iy], x[ip], y[ip]);*/
	      chi = aphase(aber, wavlen, kxp[ix], kyp[iy], xoff, yoff);
                prober[ip][ix][iy] = tr = (float) cos( chi );
                probei[ip][ix][iy] = ti = (float) sin( chi );
                sum0 += (double) (tr*tr + ti*ti);
            } else {
                 prober[ip][ix][iy] = 0.0F;
                 probei[ip][ix][iy] = 0.0F;
            }
        }

        scale = (float) ( 1.0/sqrt(sum0) );
        for( ix=0; ix<nxprobe; ix++)
        for( iy=0; iy<nyprobe; iy++) {
            prober[ip][ix][iy] *= scale;
            probei[ip][ix][iy] *= scale;
        }

    }  /* end for( ip...) */

    /*  transmit thru nslice layers
        beware ixoff,iyoff must be with one nx,ny
            of the array bounds
    */

    nxprobel = (long) nxprobe;
    nyprobel = (long) nyprobe;

    nxl = (long) nx;
    nyl = (long) ny;
    
    scale = 1.0F / ( ((float)nx) * ((float)ny) );

    zslice = 0.75*deltaz;  /*  start a little before top of unit cell */
    istart = 0;
    nslice = 0;
    it = 0;     /* thickness level index */

    if( zmax > cz ) totalz = zmax;
        else totalz = cz;
    printf( "specimen range is 0 to %g Ang.\n", totalz );
    
    /* range of unit cell */
    while(  (zslice < (totalz+0.25*deltaz)) || (istart<natom) ) {

       /* find range of atoms for current slice */
       na = 0;
       for(i=istart; i<natom; i++) {
          if( za2[i] < zslice ) na++; else break;
       }

       printf( "slice ending at z= %g Ang. with %d atoms\n", zslice, na );

       /* calculate transmission function and bandwidth limit */
       if( na > 0 ) trlayer( &xa2[istart], &ya2[istart], &occ2[istart],
            &Znum2[istart], na, (float)ax, (float)by, (float)keV,
            transr, transi,
            nxl, nyl, &phirms, &nbeamt, (float) k2maxp );

/*#pragma omp parallel for private(ix,iy,ixt,iyt,trr,tri,prr,pri) cannot be recognized by ODIE*/ 
       for( ip=0; ip<npos; ip++) {
           /* apply transmission function if there are atoms in this slice */
           if( na > 0 ) {
                fft2d( prober[ip], probei[ip], nxprobel, nyprobel, -1);
    
                for( ix=0; ix<nxprobe; ix++) {
                    ixt = ix + ixoff[ip];
                    if( ixt >= nx ) ixt = ixt - nx;
                    else if( ixt < 0 ) ixt = ixt + nx;
                    trr = transr[ixt];
                    tri = transi[ixt];
                    for( iy=0; iy<nyprobe; iy++) {
                        iyt = iy + iyoff[ip];
                        if( iyt >= ny ) iyt = iyt - ny;
                        else if( iyt < 0 ) iyt = iyt + ny;
                        prr = prober[ip][ix][iy];
                        pri = probei[ip][ix][iy];
                        prober[ip][ix][iy] = prr*trr[iyt] - pri*tri[iyt];
                        probei[ip][ix][iy] = prr*tri[iyt] + pri*trr[iyt];
                    } /* end for(iy...) */
                }  /* end for(ix...) */
                fft2d( prober[ip], probei[ip], nxprobel, nyprobel, +1); //+1 or -1 used to mark fft or ifft
           }

    
            /*  multiplied by the propagator function */
            propagate( prober[ip], probei[ip], 
                propxr, propxi, propyr, propyi,
                kxp2, kyp2, (float)k2maxp, nxprobe, nyprobe );

        }  /* end  for( ip=... */

       /*  if this is a good thickness level then save the ADF signals
           - remember that the lest level may be off by one layer with
              thermal displacements so special case it
       */
       
       /*  look at all values because they may not be in order */
       for( it = 0; it<nThick; it++ ) 
        if( fabs(ThickSave[it]-zslice)<fabs(0.5*deltaz)) {
            
            printf( "save ADF signals, thickness level %d\n", it);  /* diagnostic */
    
            for( ip=0; ip<npos; ip++) {
                /*  zero sum count */
                sum[ip] = 0.0;
                for(ix=0; ix<ndetect; ix++) detect[it][ix][ip] = 0.0;
    
                /*  sum intensity incident on the detector
                and calculate total integrated intensity
                */
                for( ix=0; ix<nxprobe; ix++)
                for( iy=0; iy<nyprobe; iy++) {
                    prr = prober[ip][ix][iy];
                    pri = probei[ip][ix][iy];
                    delta = prr*prr + pri*pri;
                    sum[ip] += delta;
                    k2 = kxp2[ix] + kyp2[iy];
                    /*for( idetect=0; idetect<ndetect; idetect++) {
                        if( (k2 >= k2min[idetect] ) &&
                        (k2 <= k2max[idetect] ) )
                        detect[it][idetect][ip] += delta*sensmap_probe[idetect][ix][iy];
                    }*/
                    for( idetect=0; idetect<ndetect; idetect++) {
                        if(k2 <= k2max[idetect] ) {
                          detect[it][idetect][ip] += delta*sensmap_probe[idetect][ix][iy];  
                        }
                    }/* New scheme for calculate intensity, apply k2max only as aperture and ignore k2min*/
                } /* end for(iy..) */
            }  /* end for( ip.. */
    
       }  /* end if( ((it...*/

       nslice++;
       zslice += deltaz;
       istart += na;

    }  /* end while( istart...) */

    free( ixoff );
    free( iyoff );
    return;

}/* end STEMsignals() */


/*--------------------- trlayer() -----------------------*/
/*
  Calculate complex specimen transmission function
  for one layer using real space projected atomic potentials

  x[],y[] = real array of atomic coordinates
  occ[]   = real array of occupancies
  Znum[]  = array of atomic numbers
  natom   = number of atoms
  ax, by  = size of transmission function in Angstroms
  kev     = beam energy in keV
  transr  = 2D array to get real part of specimen
        transmission function
  transi  = 2D array to get imag part of specimen
        transmission function
  nx, ny  = dimensions of transmission functions
  *phirms = average phase shift of projected atomic potential
  *nbeams = will get number of Fourier coefficients
  k2max   = square of max k = bandwidth limit

*/
void trlayer( const float x[], const float y[], const float occ[],
    const int Znum[], const int natom, 
    const float ax, const float by, const float kev,
    float **transr, float **transi, const long nx, const long ny,
    double *phirms, long *nbeams, const float k2max )
{
    int idx, idy, i, ixo, iyo, ix, iy, ixw, iyw, nx1, nx2, ny1, ny2;
    float k2;
    double r, rx2, rsq, vz, rmin, rmin2, sum, scale, scalex, scaley;

    const double rmax=3.0, rmax2=rmax*rmax; /* max atomic radius in Angstroms */

    scale = sigma( kev ) / 1000.0;  /* in 1/(volt-Angstroms) */

    scalex = ax/nx;
    scaley = by/ny;

    /* min radius to avoid  singularity */
    rmin = ax/((double)nx);
    r = by/((double)ny);
    rmin =  0.25 * sqrt( 0.5*(rmin*rmin + r*r) );
    rmin2 = rmin*rmin;

    idx = (int) ( nx*rmax/ax ) + 1;
    idy = (int) ( ny*rmax/by ) + 1;

    for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++)
        transr[ix][iy] = 0.0F;
    
    /*  force LUT init. to avoid redundant init in parallel form */ 
    rsq = 0.5;  /* arbitrary position */   
    for( i=0; i<natom; i++) vz =  vzatomLUT( Znum[i], rsq );
    
/*  paralleling this loop has little effect   */
/*#pragma omp parallel for private(ix,iy,ixo,iyo,nx1,nx2,ny1,ny2,rx2,r,ixw,iyw,vz,rsq)*/
    for( i=0; i<natom; i++) {
       ixo = (int) ( x[i]/scalex );
       iyo = (int) ( y[i]/scaley );
       nx1 = ixo - idx;
       nx2 = ixo + idx;
       ny1 = iyo - idy;
       ny2 = iyo + idy;

    /* add proj. atomic potential at a local region near its center
       taking advantage of small range of atomic potential */

       for( ix=nx1; ix<=nx2; ix++) {
            rx2 = x[i] - ((double)ix)*scalex;
            rx2 = rx2 * rx2;
            ixw = ix;
            while( ixw < 0 ) ixw = ixw + nx;
            ixw = ixw % nx;
            for( iy=ny1; iy<=ny2; iy++) {
               rsq = y[i] - ((double)iy)*scaley;
               rsq = rx2 + rsq*rsq;
               if( rsq <= rmax2 ) {
                iyw = iy;
                while( iyw < 0 ) iyw = iyw + ny;
                iyw = iyw % ny;
                if( rsq < rmin2 ) rsq = rmin2;
                /*r = sqrt( rx2 + r*r );
                  vz = occ[i] * vzatom( Znum[i], r ); slow */
                vz = occ[i] * vzatomLUT( Znum[i], rsq );
                transr[ixw][iyw] += (float) vz;
               }
            } /* end for(iy... */
       }  /* end for(ix... */

    } /* end for(i=0... */

    /* convert phase to a complex transmission function */
    sum = 0;
    for( ix=0; ix<nx; ix++) 
    for( iy=0; iy<ny; iy++) {
        vz = scale * transr[ix][iy];
        sum += vz;
        transr[ix][iy] = (float) cos( vz );
        transi[ix][iy] = (float) sin( vz );
    }

    *phirms = sum / ( ((double)nx)*((double)ny) );

    /* bandwidth limit the transmission function */
    *nbeams = 0;
    fft2d( transr, transi, nx, ny, +1);
    for( ix=0; ix<nx; ix++) for( iy=0; iy<ny; iy++) {
        k2 = ky2[iy] + kx2[ix];
        if (k2 < k2max) *nbeams += 1;
        else transr[ix][iy] = transi[ix][iy] = 0.0F;
    }
    fft2d( transr, transi, nx, ny, -1);
    
    return;

 }  /* end trlayer() */