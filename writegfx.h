/*

WriteGFX.h: write Gatan Fixed Format file for 2D floating point data.

broken out of ktiff2dm by pmv, 10/30/10

*/

int WriteDMFixedFormat(FILE *fp, float **data, long32 size_x, long32 size_y);

/* fixed format data types enum from Gatan web site
   http://www.gatan.com/~software/ScriptingLanguage/ScriptingWebPages/importexport.html */
enum { NULL_DATA, SIGNED_INT16_DATA, REAL4_DATA, COMPLEX8_DATA, OBSELETE_DATA,
	 PACKED_DATA, UNSIGNED_INT8_DATA, SIGNED_INT32_DATA, RGB_DATA, SIGNED_INT8_DATA, 
	 UNSIGNED_INT16_DATA, UNSIGNED_INT32_DATA , REAL8_DATA, COMPLEX16_DATA, BINARY_DATA, 
	 RGBA_FLOAT32_DATA, RGB_UINT16_DATA , RGB_FLOAT64_DATA, RGBA_FLOAT64_DATA, 
	 RGBA_UINT16_DATA, RGB_UINT8_DATA , RGBA_UINT8_DATA, LAST_DATA, 
	 OS_RGBA_UINT8_DATA = RGB_DATA }; 

