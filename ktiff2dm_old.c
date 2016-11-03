/* tiff2dm: translate a Kirkland double tiff file into
   real 4-byte fixed format readable directly by Digital
   Micrograph.

   begun 1/8/01 pmv
   updated to ansi c, a new version of Kirkland libraries and for 64-bit machines  11-27-08 pmv
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tiffsubs.h"
#include "slicelib.h"

int WriteDMFixedFormat(FILE *fp, float **data, long32 size_x, long32 size_y);

/* fixed format data types enum from Gatan web site
   http://www.gatan.com/~software/ScriptingLanguage/ScriptingWebPages/importexport.html */
enum { NULL_DATA, SIGNED_INT16_DATA, REAL4_DATA, COMPLEX8_DATA, OBSELETE_DATA,
	 PACKED_DATA, UNSIGNED_INT8_DATA, SIGNED_INT32_DATA, RGB_DATA, SIGNED_INT8_DATA, 
	 UNSIGNED_INT16_DATA, UNSIGNED_INT32_DATA , REAL8_DATA, COMPLEX16_DATA, BINARY_DATA, 
	 RGBA_FLOAT32_DATA, RGB_UINT16_DATA , RGB_FLOAT64_DATA, RGBA_FLOAT64_DATA, 
	 RGBA_UINT16_DATA, RGB_UINT8_DATA , RGBA_UINT8_DATA, LAST_DATA, 
	 OS_RGBA_UINT8_DATA = RGB_DATA }; 

int main(int argc, char *argv[]) {

	FILE *out_f;
	long32 size_x, size_y;
	int nbits[3], nsamples, npix;
	float param[64];
	char datetime[20];
	float **data;

	char in_filename[255], out_filename[255];
	
	if(argc > 1) {
		if(!strcmp(argv[1], "-help") || !strcmp(argv[1], "-h")) {
			printf("tiff2dm translates a Kirkland 32-bit float tiff file into a fixed\n");
			printf("format (.gfx) DigitalMicrograph file.   The usage is:\n");
			printf("\nktiff2dm <tiff file> <DigitalMicrograph file>\n\n");
			exit(1);
		}
	}

	if(argc != 3) {
		printf("tiff2dm translates a Kirkland 32-bit float tiff file into a fixed\n");
		printf("format (.gfx) DigitalMicrograph file.  You've called it with the wrong\n");
		printf("number of arguments.  The usage is:\n");
		printf("\nktiff2dm <tiff file> <DigitalMicrograph file>\n\n");
		exit(1);
	}

	strcpy(in_filename, argv[1]);
	strcpy(out_filename, argv[2]);

	/* open and read the TIFF file. */
	if(topenFloat(in_filename) != 1) exit(0);

	tsize(&size_x, &size_y, nbits, &nsamples);
	data = (float**) malloc2D(size_x, size_y, sizeof(float), "tiff data");
	if(treadFloatPix(data, size_x, size_y, &npix, datetime, param) != 1) exit(0);
	tclose();

	/* open and write the DM file */
	if( !(out_f = fopen(out_filename, "w+b")) ) {
		printf("Cannot open the output file.  Exiting.\n");
		exit(1);
	}
	WriteDMFixedFormat(out_f, data, size_x, size_y);
	fclose(out_f);

	/*free2D(data);*/
}


int WriteDMFixedFormat(FILE *fp, float **in_data, long32 size_x, long32 size_y) {

  unsigned long32 endian = 0x0000FFFF;	/* indicates non-byte swapped ordering */
  long32 size[3];				/* x, y, z dimensions of the image */
  long32 depth = 4;			/* byte depth of the data (redundant in this case) */
  long32 type = REAL4_DATA;		/* data type */
	float *data;

	int i,j;


	/*  write the fixed format heaers */
	fwrite(&endian, sizeof(unsigned long32), 1, fp);
	size[0] = size_x;
	size[1] = size_y;
	size[2] = 1;
	fwrite(size, sizeof(long32), 3, fp);
	fwrite(&depth, sizeof(long32), 1, fp);
	fwrite(&type, sizeof(long32), 1, fp);


	/* get the data in one big block for an fwrite */
	data = (float*) malloc1D((size_x*size_y), sizeof(float), "output data block");
	for(j=0; j<size_y; j++) {
		for(i=0; i<size_x; i++) {
		  data[size_x*j + i] = in_data[i][(size_y-1) - j];  /*y goes from size_y-1 to 0.*/
		}
	}

	/* write the image data */
	fwrite(data, sizeof(float), (size_t) size_x*size_y, fp);

	free(data);
	return 1;
}

