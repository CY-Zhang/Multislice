/* 

WriteGFX.c: write Gatan Fixed Format data files for floating point data
from C

broken out from ktiff2dm.c  10/29/10, pmv


*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "tiffsubs.h"
#include "slicelib.h"
#include "writegfx.h"

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

