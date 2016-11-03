/*      *** fft2dc.h ***

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

        1D and 2D Fast Fourier Transforms ( FFT's )
        including a fast radix4,2 complex to complex FFT in C
        
        started 8-jan-1995 E. Kirkland
        added real to complex 1D and 2D FFT's 2-jun-1995 ejk
        reorganized rfft2d() 13-jun-1995 ejk
        added fft42t() and put in fft2d() 16-jul-1995 ejk
        added conditional compilation of fft42t() vs.
                fft42() because some CPUS are faster one
                way and some are faster the other 28-jul-1995 ejk
        add invert2D() and powerof2() 7-sep-1995 ejk
        add IBM/ESSL scft() 8-sep-1995 ejk
        removed possible free() of unallocated memory
                        in fft42t() 21-sep-1997 ejk
        removed some non-essential routines 2-mar-1998 ejk
        fixed small typo in "not a power of 2" error message
                in fft2d()  11-oct-1999 ejk
        fixed bug in fft2dc() when nx and ny are not the same
                   11-oct-1999 ejk

fft2d()   : 2D complex to complex FFT
rfft2d()  : 2D real to complex FFT
fft42()   : 1D complex to complex radix 4,2 FFT
fft42t()  : 1D complex to complex radix 4,2 FFT
                similar to fft42() but with look-up-tables
                faster for multiple calls of same length
invert2D(): move FFT center from corners to the center
powerof2(): return nearest power of 2 >= argument

*/

/*-----------------  fft2d() ---------------- */
/*
        2D complex to complex FFT

        pixr[ix][iy] = real part of 2D pix to Fourier transform
        pixi[ix][iy] = imag part of 2D pix to Fourier transform
        nx,ny = (long int) size of array
                ix=0,1,2...(nx-1)
                iy=0,1,2...(ny-1)
        inverse = if > 0 then forward transform
                  if < 0 then forward transform

        On exit pixr and pixi will have the Fourier transform 
        and the original data will be lost.      
*/

void fft2d( float **pixr, float **pixi, const long nx,
                   const long ny, int inverse );


/*--------------------  rfft2d() -------------------- */
/*
        2D real to complex FFT

        The fwd FFT produces the positive half plane with kx=0
        and the the ky frequencies are stored in the normal
        wrap around manner iy=0 is DC, iy=1 is the lowest positive
        frequency, iy=(ny-1) is the lowest neg. frequency.
        The complex values are stored at adjacent ix values.
        for example (pixr[0][iy],pixr[1][iy]) = (real,imag)

        pixr[ix][iy] = real 2D pix to Fourier transform
           NOTE: array size must be (nx+2) x ny
        nx,ny = (long int) size of original data (before FFT)
                ix=0,1,2...(nx-1)
                iy=0,1,2...(ny-1)
        inverse = if > 0 then forward transform
                  if < 0 then forward transform

        On exit pixr[][] will have the Fourier transform 
        and the original data will be lost.

        Although it is possible to pack the fwd FFT data into the
        same sized output array it is arranged in an awkward
        order so add 2 to the array size in the x direction for
        the high frequency component to make it easy to index
        the frequencies.
*/
void rfft2d( float **pixr, long nx, long ny, int inverse );


/*------------------------ fft42 --------------------------

        fft42( fr[], fi[], n )       radix-4,2 FFT in C

        fr[], fi[] = (float) real and imaginary array with input data
        n          = (long) size of array

  calculate the complex fast Fourier transform of (fr,fi)   
  input array fr,fi (real,imaginary) indexed from 0 to (n-1)
  on output fr,fi contains the transform 

  started 8-jan-1995 E. Kirkland
  
*/

void fft42 ( float *fr, float *fi, long n );

/*------------------------ fft42t --------------------------

        fft42t( fr[], fi[], n )       radix-4,2 FFT in C

        fr[], fi[] = (float) real and imaginary array with input data
        n          = (long) size of array

  calculate the complex fast Fourier transform of (fr,fi)   
  input array fr,fi (real,imaginary) indexed from 0 to (n-1)
  on output fr,fi contains the transform 

        this is similar to fft42() but uses a look-up-table that
        speed it up a little if there are many calls of the same
        length (i.e. in multidimensional FFT's)

  started from fft42() 16-jul-1995 E. Kirkland
  
*/

void fft42t ( float *fr, float *fi, long n  );

/*------------------------- invert2D() ----------------------*/
/*
        rearrange pix with corners moved to center (a la FFT's)

         pix[ix][iy] = real array with image
         nx,ny = range of image 0<ix<(nx-1) and 0<iy<(ny-1)

*/
void invert2D( float** pix, long nx, long ny );

/*---------------------------- powerof2() ---------------------------*/
/*
        return the nearest power of 2 greater than or equal to 
        the argument
*/
long powerof2( long n );
