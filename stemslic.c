/*          *** stemslic(e).c ***

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

   ANSI-C and TIFF version

  Calculate STEM image using repetitive multislice
  
  this file is formatted for a tab size of 8 characters

  rewritten in C, summer 1995 E. Kirkland
  removed declaration of unused var. 1-june-1995 ejk
  switch to TIFF file I/O 26-may-1996 ejk
  modified pointer arithmetic in STEMsignal() 22-jun-1996 ejk
  move layer expansion out of here and into atompot 4-aug-1996 ejk
  removed commas from input format and fixed minor bugs in
        image mode 21-jul-1997 ejk
  small cosmetic changes 17-jan-1998 ejk
  add astigmatism 4-feb-1998 ejk
  fix minor bug in astigmatism input, 1D comments
      and bad position error message    20-jul-1998 ejk
  updated memory allocator functions 14-nov-1999 ejk
  change void main() to int main() for better portability
             22-jan-2000 ejk
  update data type of nxl,nyl to be consistent with new tiffsubs
       18-jul-2007 ejk
  convert to GPL 3-jul-2008 ejk
  added function to read detector sensitivity 3-29-2016 cz

*/

#include <stdio.h>      /* ANSI C libraries used */
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "fft2dc.h"     /* FFT routines */
#include "slicelib.h"   /* misc. routines for multislice */
#include "tiffsubs.h"   /* file I/O routines in TIFF format */

#define BW (2.0F/3.0F)  /* antialiasing bandwidth limit factor */
#define ABERR 1.0e-5    /* max error for a,b */

#define NSMAX   1000    /* maximum number of total slices */
#define NCMAX   132     /* max number of characters per line */
#define NCINMAX 500     /* max number of characters in stacking spec */
#define NLMAX   52      /* maximum number of layers */
#define NPARAM  64      /* number of parameters */
#define TRUE 1

/* global specimen and probe parameters to share */

float **prober, **probei;       /* complex probe wavefunction */
float ***transr, ***transi;     /* complex transmission functions */
float **propxr, **propxi;       /* complex propagator vs x */
float **propyr, **propyi;       /* complex propagator vs y */
float **sensmap, ***sensmap_probe, ***sensmap_save; **sensmap_rotated; **sensmap_rotated_save;
float **mrad_sens;

int nx, ny, nxprobe, nyprobe, nslice, ndetect;
int *layer;
float *kx, *ky, *kx2, *ky2, *kxp, *kyp, *kxp2, *kyp2, *mradx, *mrady, *mradx_sens, *mrady_sens;
float  *rmin, *rmax;
double ax, by, wavlen, k2maxp, Cs,df,apert1, apert2, pi, keV;
double Df0, cDf_sigma;
float dfa2, dfa2phi, dfa3, dfa3phi;     /* astigmatism parameters */
double *almin, *almax, *k2max, *k2min, *detect;
float *ratio;
double **aber;

/* define functions at end of this file (i.e. so main can be 1st) */
void STEMsignal( double x, double y, double *detect, double* sum, float ***sensmap_probe);

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
    c += 0.25*pow(wl,3)*aber[5][0]*(pow(kxr, 4) - 6*pow(kxr,2)*pow(kyr, 2) + pow(kyr,4));   /* A3, 4-fold astigmasitm */
    kxr = kx*aber[6][3] + ky*aber[6][2];  kyr = ky*aber[6][3] - kx*aber[6][2];
    c += pow(wl,3)*aber[6][0]*(pow(kxr,4) - pow(kyr,4));                                        /* S3, star aberration */ 
    kxr = kx*aber[7][3] + ky*aber[7][2];  kyr = ky*aber[7][3] - kx*aber[7][2];
    c += 0.2*pow(wl,4)*aber[7][0]*(pow(kxr, 5) - 10*pow(kxr,3)*pow(kyr,2) + 5*kxr*pow(kyr,4));  /* A4, 5-fold astigmatism */
    kxr = kx*aber[8][3] + ky*aber[8][2];  kyr = ky*aber[8][3] - kx*aber[8][2];
    c += pow(wl,4)*aber[8][0]*(pow(kxr, 5) - 2*pow(kxr,3)*pow(kyr,2) - 3*kxr*pow(kyr, 4));  /* D4, 3-lobe aberration */
    kxr = kx*aber[9][3] + ky*aber[9][2];  kyr = ky*aber[9][3] - kx*aber[9][2];
    c += pow(wl,4)*aber[9][0]*(pow(kxr,5) + 2*pow(kxr,3)*pow(kyr,2) + kxr*pow(kyr,4));          /* B4, axial coma */

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
        char datetime[24];
        char *cin, *cin2, **filein, **fileout, **Sensitivity;
        const char version[] = "3-jul-2008 (ejk)";

        int ix, iy, nx2, ny2, i, islice, idetect,
                nout, nxout, nyout, ilayer, nlayer, l1d,
                nbits[3], samples, npix, senx, seny, npixx=1, j, minIdx, minIdy, centerx, centery;
        int px_in_1, px_out_1, px_in_2, px_out_2;

        long nbeamt, nbeamp, nxl, nyl;
        long32 nxl2, nyl2;   /* tiffsubs data type */

        float *x2, *y2, *param, *cz, ***pixr, *param_detector, reference, mindifference, difference;
        float temp_sens, px_in_1f, px_out_1f, px_in_2f, px_out_2f;
        float mradx_rot, mrady_rot, reference_x, reference_y, rot_angle, length, new_angle, ori_angle;

        double scale, v0, vz, mm0, x,y, ax2,by2, sum, k2, w,
           tctx, tcty, xi,xf, yi,yf, dx, dy, cztot, totmin, totmax,
           ctiltx, ctilty, time;
        double val, i_new, j_new;
        char **probe_file;
        int i_new_int,j_new_int;


        FILE *fp;

        /*  set up symbolic mapping
                this must be the same as in parlay  */

        const char cname[] =
                "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ";

/* start by announcing version etc */
        printf("stemslic(e) version dated %s\n", version );
        printf("Copyright (C) 1998, 2008 Earl J. Kirkland\n" );
        printf( "This program is provided AS-IS with ABSOLUTELY NO WARRANTY\n "
            " under the GNU general public license\n\n" );

        printf("perform STEM multislice\n\n");

        pi = 4.0 * atan( 1.0 );

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
           scanf("%s",  filein[i] );
        }


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
        printf("val = %lf, rot_angle = %lf\n",val, rot_angle );
        printf("cos = %lf, sin = %lf\n", cos(rot_angle*val), sin(rot_angle*val));
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
/*  get more parameter etc */
    aber = (double**)malloc2D(12, 4, sizeof(double), "Aberration array");

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
           printf("Size of probe wavefunction"
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

        printf("Number of detector geometries :\n");
        do scanf( "%d", &ndetect );
        while (ndetect <= 0);

        almin  = (double*) malloc1D( ndetect, sizeof(double), "almin" );
        almax  = (double*) malloc1D( ndetect, sizeof(double), "almax" );
        fileout = (char**) malloc2D( ndetect, NCMAX, sizeof(char), "fileout" );
        probe_file = (char**) malloc2D(ndetect, NCMAX, sizeof( char ),"probe_file" );
        ratio = (float*) malloc1D( ndetect,sizeof(float), "ratio");

        for( idetect=0; idetect<ndetect; idetect++) {
                if( (l1d == 1) && (idetect==0) ) {
                        printf("Name of file to get output"
                                " of 1D STEM multislice result :\n");
                        scanf("%s", fileout[idetect]);
                }
                printf("Detector %3d: Type, min,max angles(mrad)"
                        " of collector :\n", idetect+1);
                scanf("%lg %lg",
                        &almin[idetect], &almax[idetect] );
                if( l1d == 0 ) {
                        printf("Name of file to get output"
                                " of result for this detector:\n");
                        scanf("%s", fileout[idetect]);
                }       
                printf("Name of file to save adapted sensitivity map for detector %d: \n",idetect+1);
                scanf("%s", probe_file[idetect]);
                px_in_1 = 0;
                px_out_1 = 0;
                px_in_2 = 0;
                px_out_2 = 0;
                for (i = centery; i >0; i--)
                {
                    temp_sens = sensmap[centerx][i];
                    if ((temp_sens>0.5) && (px_in_1 == 0))
                    {
                        px_in_1 = i;
                    }
                    if ((temp_sens<0.5) && (px_in_1!=0) &&(px_out_1==0))
                    {
                        px_out_1 = i;
                    }
                }
                for (i = centery; i <seny; i++)
                {
                    temp_sens = sensmap[centerx][i];
                    if ((temp_sens>0.5) && (px_in_2 == 0))
                    {
                        px_in_2 = i;
                    }
                    if ((temp_sens<0.5) && (px_in_2!=0) &&(px_out_2==0))
                    {
                px_out_2 = i;
                    }
                }
                px_in_1f = (float)px_in_1;
                px_out_1f = (float)px_out_1;
                px_in_2f = (float)px_in_2;
                px_out_2f = (float)px_out_2;
                ratio[idetect] = 2*(almax[idetect]-almin[idetect])/(fabs(px_out_1f-px_in_1f)+fabs(px_out_2f-px_in_2f));
                printf("ratio of detector %d is %f mrad/px\n", idetect+1, ratio[idetect]);
                }  
                /* end for(idetect=.. */
        /*here xi and xf marks the border of final image, in Ang, nxout nyout is the px resolution*/

        if( l1d == 1 ) {
                printf("xi, xf, yi, yf, nout :\n");
                scanf("%lg %lg %lg %lg %d", &xi, &xf, &yi, &yf, &nout);
        }else {
                printf("xi,xf,yi,yf, nxout,nyout :\n");
                scanf("%lg %lg %lg %lg %d %d",
                         &xi, &xf, &yi, &yf, &nxout, &nyout);
        }

        /*  read in atomic potential and specimen parameters
                and calculate specimen transmission function
                for a single slice in transr,i
        */

        time = cputim();   /* get initial CPU time */

        param = (float*) malloc1D( NPARAM, sizeof(float), "param" );
        for( i=0; i<NPARAM; i++) param[i] = 0.0F;
        transr = (float***) malloc( nlayer* sizeof( float**) );
        transi = (float***) malloc( nlayer* sizeof( float**) );
        if( (transr==NULL) || (transi==NULL) ) {
                printf("Cannot allocate transmission function memory\n");
                exit( 0 );
        }
        cz   = (float*) malloc1D( nlayer, sizeof(float), "cz" );

        for( ilayer=0; ilayer<nlayer; ilayer++ ) {

                if( topenFloat( filein[ilayer] ) != 1 ) {
                        printf("Cannot open file %s.\n", filein[ilayer] );
                        exit( 0 );
                }
                tsize( &nxl2, &nyl2, nbits, &samples );
                nx2 = (int) nxl2;
                ny2 = (int) nyl2;
                transr[ilayer] = (float**) malloc2D( nx2,  ny2, 
                        sizeof(float), "transr" );
                transi[ilayer] = (float**) malloc2D( nx2,  ny2,
                        sizeof(float), "transi" );

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
                /*printf("param[pDX] = %f, pDX = %d\n\n",param[pDX], pDX );*/
                cz[ilayer] = param[pC];
                printf("layer %c, cz = %f\n\n", cname[ilayer], cz[ilayer]);
               /*printf("\nax = %f\n by = %f\n\n", ax, by);*/

                if (  ilayer != 0 ) {
                    /*To make sure each layer has same ax and by as previous layer*/
                        if ( (nx!=nx2) || (ny!=ny2) ) {
                                printf("pix size incompatible.\n");
                                printf("old size = %d, %d\n", nx, ny);
                                printf("new size = %d, %d\n", nx2, ny2 );
                                printf("layer = %1c\n", cname[ilayer]);
                                exit( 0 );
                        }
                        ax2 = param[pDX] * nx;
                        by2 = param[pDY] * ny;
                        if( ( fabs( ax-ax2 ) > fabs(ABERR*ax) ) ||
                            ( fabs( by-by2 ) > fabs(ABERR*by) ) ) {
                                printf("incompatible lattice constants\n");
                                printf("potential    a,b,c = %g, %g, %g",
                                       ax2,by2,cz[ilayer] );
                                printf("previous pix a,b = %g, %g\n",
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
                                         cname[ilayer], nx, ny  );
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
               /* printf("\nax = %f\n by = %f\n\n", ax, by); */

        }  /* end for(ilayer=... */

 /*       printf("Size in pixels Nx x Ny = %d x %d = %ld"
                " total pixels, \n lattice constants a,b = %f x %f\n param[pDX] = %f nx = %d\n\n",
                nx,ny, ((long)nx)*((long)ny), ax, by, param[pDX], nx);*/

        /*  check for valid scan coordinates  */

        if(     (xi < 0.0) || (xi > ax) ||
                (xf < 0.0) || (xf > ax) ||
                (yi < 0.0) || (yi > by) ||
                (yf < 0.0) || (yf > by) ) {
            printf("Coordinates out of range please try again.\n");
            exit( 0 );
        }

        /*  calculate the total specimen thickness and echo */

        cztot = 0.0F;
        for( islice=0; islice<nslice; islice++)
                cztot +=  cz[ layer[islice] ];
        printf("Total specimen thickness = %g Angstroms\n", cztot);

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

        /*  calculate spatial frequencies for future use
                (one set for transmission function and one for probe
                wavefunction)
        NOTE: zero freg is in the bottom left corner and
                expands into all other corners - not in the center
                this is required for FFT - don't waste time rearranging

        remember : the x values are the same for both sets
        
        x2, y2 are not used anywhere else - just for freqn()

        kxp, kxp2, x2, nxprobe is used for probe range while kx, kx2, x2, nx, is used for image pixel range
        
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
        printf("OK\n");
        
        /*kxp range only depends on ax, if nxprobe changes, ax*nxprobe/nx will also change and results remains same
          ax is the real space size of unit cell, in Angstrom*/

        freqn( kxp, kxp2, x2, nxprobe, ax*((double)nxprobe)/nx );
        freqn( kyp, kyp2, y2, nyprobe, by*((double)nyprobe)/ny );

        /*Initialize value in mrad list, both for detector sens and beam*/
        for (i = 0; i < nxprobe; i++)
        {
            mradx[i] = kxp[i]/0.001*wavlen;
            mrady[i] = kyp[i]/0.001*wavlen;
        }
        printf("OK\n");

        sensmap_probe = (float***)malloc(ndetect* sizeof(float**));
        sensmap_save = (float***)malloc(ndetect* sizeof(float**)); 
        printf("OK\n");

        if( sensmap_probe==NULL) {
            printf("Cannot allocate sensitivity map memory\n");
            exit( 0 );
        }
        for (idetect = 0; idetect < ndetect; idetect++)
        {
            sensmap_probe[idetect] = (float**)malloc2D(nxprobe, nyprobe, sizeof(float),"sensmap_probe");
            sensmap_save[idetect] = (float**)malloc2D(nxprobe, nyprobe, sizeof(float),"sensmap_save");
        }

        printf("OK\n");

        mrad_sens = (float**)malloc2D(senx, seny, sizeof(float),"mrad_sens");
        for (idetect = 0; idetect < ndetect; idetect++)
        {
                for (i = 0; i < senx; i++)
                {
                    mradx_sens[i] = (i-centerx)*ratio[idetect];
                    /*mradx_sens[i] = double(mradx_sens[i]);*/
                }
                for (i = 0; i < seny; i++)
                {
                    mrady_sens[i] = (i-centery)*ratio[idetect];
                    /*mrady_sens[i] = double(mrady_sens[i]);*/
                }
        }

            printf("OK\n");

        /* find value for each pixel in sensmap_probe, which is used to account for detector sensitivity */

            printf("senx = %d seny = %d\n",senx, seny );
            printf("nxprobe = %d nyprobe = %d\n",nxprobe, nyprobe );

      for (idetect = 0; idetect<ndetect; idetect++)
      {
        /*printf("ratio[%d] = %f\n",idetect+1, ratio[idetect] );*/
        for (i = 0; i < senx; i++)
        {
            mradx_sens[i] = (i-centerx)*ratio[idetect];
            mrady_sens[i] = (i-centery)*ratio[idetect];
            /*printf("mradx_sens[%d] = %f, mrady_sens[%d] = %f\n",i, mradx_sens[i], i, mrady_sens[i]);*/
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
            /*printf("sensmap_probe[%d][%d][%d] = sensmap[%d][%d] = %f\n",idetect, ix, iy, minIdx, minIdy, sensmap_probe[idetect][ix][iy]);*/
            /*printf("sensmap_probe[%d][%d][%d] = sensmap[%d][%d]\n", idetect, ix, iy, minIdx, minIdy);*/
            /*printf("mradx[%d] = %f\n",ix, mradx[ix] );*/
            }
        }
        param[1] = 1;
        param[3] = 0;

        tcreateFloatPixFile( probe_file[idetect+1], sensmap_save[idetect], (long)nxprobe, (long) nyprobe, 1, param);
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
        k2maxp = k2maxp * k2maxp;

        nxl = (long) nx;
        nyl = (long) ny;

        for( ilayer=0; ilayer<nlayer;ilayer++) {

           fft2d( transr[ilayer], transi[ilayer], nxl, nyl, +1);
           nbeamt = 0;

           for( ix=0; ix<nx; ix++ )
                for( iy=0; iy<ny; iy++) {
                        k2 = ky2[iy] + kx2[ix];
                        if( k2 < k2maxp ) nbeamt++;
                        else {
                                transr[ilayer][ix][iy] = 0.0F;
                                transi[ilayer][ix][iy] = 0.0F;
                        }
                } /* end for(iy... ) */
           fft2d( transr[ilayer], transi[ilayer], nxl, nyl, -1);

        } /* end for( ilayer */

        printf("Number of symmetrical anti-aliasing "
               "beams in trans. function = %ld\n", nbeamt);
        printf("  with a resolution of %g Angstroms.\n",
               1.0/sqrt(k2maxp));

        /* calculate propagator functions with probe sample size
                impose anti-aliasing bandwidth limit
        */
        propxr = (float**) malloc2D( nlayer, nxprobe,
                sizeof(float), "propxr" );
        propxi = (float**) malloc2D( nlayer, nxprobe,
                sizeof(float), "propxi" );
        propyr = (float**) malloc2D( nlayer, nyprobe,
                sizeof(float), "propyr" );
        propyi = (float**) malloc2D( nlayer, nyprobe,
                sizeof(float), "propyi" );

        tctx = 2.0 * tan(ctiltx);
        tcty = 2.0 * tan(ctilty);

        for( ilayer=0; ilayer<nlayer; ilayer++) {
                scale = cz[ilayer];
                for( ix=0; ix<nxprobe; ix++) {
                   w = pi * scale * 
                        ( kxp2[ix] * wavlen - kxp[ix]*tctx );
                   propxr[ilayer][ix]= (float)  cos(w);
                   propxi[ilayer][ix]= (float) -sin(w);
                }

                for( iy=0; iy<nyprobe; iy++) {
                   w = pi * scale * 
                        ( kyp2[iy] * wavlen - kyp[iy]*tcty );
                   propyr[ilayer][iy]= (float)  cos(w);
                   propyi[ilayer][iy]= (float) -sin(w);
                }

        }  /* end for(ilayer..) */

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
        printf("wavelength = %f\n\n", wavlen);

        /*  init the min/max record of total integrated intensity */

        totmin =  10.0;
        totmax = -10.0;
        detect = (double*) malloc1D( ndetect, sizeof(double), "detect" );
        rmin = (float*) malloc1D( ndetect, sizeof(float), "rmin" );
        rmax = (float*) malloc1D( ndetect, sizeof(float), "rmax" );

                /* probe wavefunction */
        prober = (float**) malloc2D( nxprobe, nyprobe, 
                sizeof(float), "prober" );
        probei = (float**) malloc2D( nxprobe, nyprobe,
                sizeof(float), "probei" );

        /*  start here for a full image output */

        if( l1d == 0 ) {
           printf("output file size in pixels is %d x %d\n",
                  nxout, nyout );
           pixr = (float***) malloc( ndetect * sizeof(float**) );
           if( pixr == NULL ) {
                printf("Cannot allocate pixr memory.\n");
                exit( 0 );
           }
           for( i=0; i<ndetect; i++) pixr[i] = 
                (float**) malloc2D( nxout, nyout, sizeof(float), "pixr" );

           /*  iterate the multislice algorithm proper for each position of
                the focused probe */

           dx = (xf-xi)/((double)(nxout-1));
           dy = (yf-yi)/((double)(nyout-1));

           for( iy=0; iy<nyout; iy++) {
                y = yi + dy * ((double)iy);
                for( ix=0; ix<nxout; ix++) {
                   x = xi + dx * ((double) ix);

                   STEMsignal( x, y, detect, &sum, sensmap_probe); 

  
                   if( sum < totmin ) totmin = sum;
                   if( sum > totmax ) totmax = sum;
                   for( i=0; i<ndetect; i++) {
                        pixr[i][ix][iy] = (float) detect[i];
                        if( (ix==0) && (iy==0) )
                           rmin[i] = rmax[i] = (float) detect[i];
                        else {
                           if( detect[i] < rmin[i] )
                                rmin[i] = (float) detect[i];
                           if( detect[i] > rmax[i] )
                                rmax[i] = (float) detect[i];
                        }
                   } /* end for(i..) */
                   if( sum < 0.9 ) printf("Warning: "
                       "integrated intensity too small, = %g, %d, %d\n",
                           sum, ix, iy );
                   if( sum > 1.1 ) printf("Warning: "
                       "integrated intensity too big, = %g, %d, %d\n",
                           sum, ix, iy );
                } /* end for(ix...) */

                printf("iy= %d, min[0]= %f, max[0]= %f\n", 
                       iy, rmin[0], rmax[0] );

        } /* end for(iy... ) */

        /*  store min and max */
        param[pIMAX]    = 0.0F;
        param[pIMIN]    = 0.0F;
        param[pXCTILT]  = (float) ctiltx;
        param[pYCTILT]  = (float) ctilty;
        param[pDEFOCUS] = (float) df;
        param[pDX]      = (float) dx;
        param[pDY]      = (float) dy;
        param[pENERGY]  = (float) v0;
        param[pOAPERT]  = (float) apert2;
        param[pCS]      = (float) Cs;
        param[pWAVEL]   = (float) wavlen;
        param[pNSLICES] = (float) nslice;

        for( i=0; i<ndetect; i++) {
           printf("output pix range for detector #%d : %g to %g\n", i, rmin[i], rmax[i]);
           param[pRMAX] = rmax[i];
           param[pRMIN] = rmin[i];
           param[pMINDET] = (float) ( almin[i] * 0.001 );
           param[pMAXDET] = (float) ( almax[i] * 0.001 );
           if( tcreateFloatPixFile( fileout[i], pixr[i], (long) nxout,
                         (long) nyout, 1, param ) != 1 ) {
                printf("Cannot write output file %s.\n", fileout[i] );
           }
        }

/*  start here for 1d line scan */

        }else if ( l1d == 1 ) {

           dx = (xf-xi)/((double)(nout-1));
           dy = (yf-yi)/((double)(nout-1));
           x = xi;
           y = yi;
           fp = fopen( fileout[0], "w+" );
           if( fp == NULL ) {
                printf("Cannot open output file %s.\n", fileout[0] );
                exit( 0 );
           }

           fprintf(fp, "C\n");
           fprintf(fp, "C   output of stemslice version %s\n", version);
           fprintf(fp, "C\n");
           fprintf(fp, "C   nslice= %d\n", nslice);
           for( i=0; i<nlayer; i++)
                fprintf(fp,"cz= %g, file in= %s\n", 
                        cz[i], filein[i] );
           fprintf(fp, "V0= %g, Cs= %g, df= %g\n", v0, Cs, df );
           fprintf(fp, "Apert= %g mrad to %g mrad\n", apert1, apert2 );
           fprintf(fp, "Crystal tilt x,y= %lg, %lg\n", ctiltx,ctilty);

           for(  idetect=0; idetect<ndetect; idetect++) {
                fprintf(fp, "Detector %d, Almin= %g mrad, Almax= %g mrad\n",
                        idetect, almin[idetect], almax[idetect] );
           }

           fprintf(fp, "ax= %g A, by= %g A\n", ax,by);
           fprintf(fp, "Number of symmetrical anti-aliasing "
                "beams in trans. function= %ld\n", nbeamt );
           fprintf(fp, "with a resolution (in Angstroms) = %g\n",
                1.0/sqrt(k2maxp) );
           fprintf(fp, "Number of symmetrical anti-aliasing "
                "beams in probe = %ld\n", nbeamp );
           fprintf(fp, "C         x              y         signal\n");

           for( ix=0; ix<nout; ix++) {

                x = xi + dx * ((double)ix);
                y = yi + dy * ((double)ix);
                STEMsignal( x, y, detect, &sum ,sensmap_probe);
                if( sum < totmin ) totmin = sum;
                if( sum > totmax ) totmax = sum;

                if( sum < 0.9) 
                   printf("Warning integrated intensity too small, = %g\n", 
                          sum );
                if( sum > 1.1) 
                   printf("Warning integrated intensity too large, = %g\n", 
                          sum );

                printf("%5d %14.7g %14.7g", ix+1, x, y);
                for(i=0; i<ndetect; i++) 
                        printf("%14.7g", detect[i] );
                printf("\n");
                
                fprintf(fp, "%14.7g %14.7g", x, y);
                for(i=0; i<ndetect; i++) 
                        fprintf(fp, "%14.7g", detect[i] );
                fprintf(fp, "\n");
           }

           fclose( fp );

        } /* end if( l1d.. ) */

        /*  echo min/max of total integrated intensity */
        printf("The total integrated intensity range was:\n");
        printf("   %g to %g\n\n",  totmin, totmax );

        printf("CPU time = %g sec.\n", cputim()-time);

        return 0;

}  /* end main() */

/*------------------------ STEMsignal() ---------------------*/
/*
  subroutine to calculate the stem signal arising from a given
  probe position

  iterate the multislice algorithm proper for each position of
  the focused probe

   note zero freg is in the bottom left corner and
     expands into all other corners - not in the center
     this is required for fft - don't waste time rearranging

  x,y       = real position of the incident probe
  detect[]  = real array to get signal into each detector
  sum       = real total integrated intensity
  
  the assumed global variables are:
  
  nxprobe,nyprobe       = int size of probe wavefunction in pixels
  nx,ny                 = int size of transmission function in pixels
  ndetect               = int number of detector geometries
  nslice                = int number of slices
  layer[]               = int array with slice layer indecies
  prober[][], probei[][] = float real,image probe wavefunction
  transr[][], transi[][] = float real,imag transmission function
  propxr[][], propxi[][] = float real,imag propagator vs x
  propyr[][], propyi[][] = float real,imag propagator vs y
  ax,by                 = float unit cell size in Angs
  kxp[], kyp[]          = float spatial frequencies vs x, y
  kxp2[], kyp2[]        = float square of kxp[], kyp[]
  apert1, apert2        = double min,max objective aperture (mrad)
  k2maxp                = double max spatial freq of probe squared
  pi                    = double constant PI
  wavlen                = double electron wavelength in Angs
  df                    = double defocus (in Ang)
  Cs                    = double spherical aberration (in mm)
  
*/
void STEMsignal( double x, double y, double *detect, double *sum, float ***sensmap_probe)
{
        int ix, iy, ixt, iyt, islice, ilayer,
                idetect, ixoff, iyoff, ixmid, iymid;

        long nxprobel, nyprobel;

        float scale, prr, pri, tr, ti, *trr, *tri;

        double xoff,yoff, chi1, chi2, k2maxa, k2maxb, chi,
                w, k2, phi;

        double sum0, delta;

/*   make sure x,y are ok */

        if( (x < 0.0) || (x > ax) ||
            (y < 0.0) || (y > by) ) {
                *sum = -1.2345;
                printf("bad x=%f,y=%f in STEMsignal()\n", x, y);
                return;
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
        ixoff = (int) floor( x*((double)nx) / ax ) - ixmid;
        xoff  = x - ax*((double)ixoff)/((double)nx);

        iymid = nyprobe/2;
        iyoff = (int) floor( y*((double)ny) / by ) - iymid;
        yoff  = y - by*((double)iyoff)/((double)ny);

        chi1 = pi * wavlen;
        chi2 = 0.5 * Cs * 1.e7*wavlen*wavlen;
        k2maxa = apert1 * 0.001/wavlen;
        k2maxa = k2maxa *k2maxa;
        k2maxb = apert2 * 0.001/wavlen;
        k2maxb = k2maxb * k2maxb;
        /* build wavefunction before hitting sample*/
        /* no chromatic aberration correction here*/
        sum0 = 0.0;
        for( ix=0; ix<nxprobe; ix++)
        for( iy=0; iy<nyprobe; iy++) {
           k2 = kxp2[ix] + kyp2[iy];
           if( (k2 >= k2maxa) && (k2 <= k2maxb) ) {
                chi = aphase(aber, wavlen, kxp[ix], kyp[iy], xoff, yoff);
                prober[ix][iy] = tr = (float) cos( chi );
                probei[ix][iy] = ti = (float) sin( chi );
                sum0 += (double) (tr*tr + ti*ti);
           } 
           else {
              prober[ix][iy] = 0.0F;
              probei[ix][iy] = 0.0F;
           }
        }

        scale = (float) ( 1.0/sqrt(sum0) );
        for( ix=0; ix<nxprobe; ix++)
        for( iy=0; iy<nyprobe; iy++) {
              prober[ix][iy] *= scale;
              probei[ix][iy] *= scale;
        }

        /*  transmit thru nslice layers
                beware ixoff,iyoff must be with one nx,ny
                        of the array bounds
        */

        nxprobel = (long) nxprobe;
        nyprobel = (long) nyprobe;
        
        for( islice=0; islice<nslice; islice++) {

           fft2d( prober, probei, nxprobel, nyprobel, -1);
           ilayer = layer[islice];

           for( ix=0; ix<nxprobe; ix++) {       
                ixt = ix + ixoff;
                if( ixt >= nx ) ixt = ixt - nx;
                        else if( ixt < 0 ) ixt = ixt + nx;
                trr = transr[ilayer][ixt];
                tri = transi[ilayer][ixt];
                for( iy=0; iy<nyprobe; iy++) {
                  iyt = iy + iyoff;
                  if( iyt >= ny ) iyt = iyt - ny;
                        else if( iyt < 0 ) iyt = iyt + ny;
                   prr = prober[ix][iy];
                   pri = probei[ix][iy];
                   prober[ix][iy] = prr*trr[iyt] - pri*tri[iyt];
                   probei[ix][iy] = prr*tri[iyt] + pri*trr[iyt];
                } /* end for(iy...) */
           }  /* end for(ix...) */

           fft2d( prober, probei, nxprobel, nyprobel, +1);

           /*  multiplied by the propagator function */
           propagate( prober, probei, propxr[ilayer], propxi[ilayer],
                propyr[ilayer], propyi[ilayer],
                kxp2, kyp2, (float)k2maxp, nxprobe, nyprobe );  

        }  /* end for( islice...) */

        /*  at this point the probe could be padded with zeros
                (in real space) and transformed back to get better
                sampling in reciprocal space but initial tests with
                this embedding technique did not indicate that it improved
                matters (26-jan-1995 ejk)

                this should be tested further sometime
        */

        /*  zero sum count */

        sum0 = 0.0;
        for(ix=0; ix<ndetect; ix++) detect[ix] = 0.0;

        /*  sum intensity incident on the detector
                and calculate total integrated intensity

                normal mode 
        */

        for( ix=0; ix<nxprobe; ix++)    /*loop over whole probe wavefunction, 512px*/
        for( iy=0; iy<nyprobe; iy++) {  /*all intensities are calculated in frequency domain*/
           prr = prober[ix][iy];
           pri = probei[ix][iy];
           delta = prr*prr + pri*pri; 
           sum0 += delta;
           k2 = kxp2[ix] + kyp2[iy];
           for( idetect=0; idetect<ndetect; idetect++) {
              if( (k2 >= k2min[idetect] ) &&    /*use k^2 to define the detect range as detector is a ring*/
                  (k2 <= k2max[idetect] ) )
                        detect[idetect] += delta*sensmap_probe[idetect][ix][iy];   /*within detect range*/
           }
        } /* end for(ix..) */

        *sum = sum0;

}/* end STEMsignal() */

