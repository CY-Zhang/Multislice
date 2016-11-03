/* tiff2dm: translate a Kirkland double tiff file into
   real 4-byte fixed format readable directly by Digital
   Micrograph.

   begun 1/8/01 pmv
   updated to ansi c, a new version of Kirkland libraries and for 64-bit machines  11-27-08 pmv
   remove WriteDMFixedFormat from this file and include writegfx.h  07-04-13 pmv
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tiffsubs.h"
#include "slicelib.h"
#include "writegfx.h"

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


