// FITSread.h

enum valueType {None, Text, CharString, Logical, Integer, LongInteger, Float, Double, ComplexInt, Complex};

#include "DataStructures.h" // for struct Axis
#include <stddef.h> // for size_t

#define SINGLE -32 // precision of image data
#define DOUBLE -64

int FITSread_primaryHeaderInfo(char filename[],
  int *pDataSize, int *pNDim, struct Axis axis[], double *pBScale, double *pBZero, char bUnit[]);
int FITSread_primaryHeaderVerify(char filename[], int dataSize, int nDim, int dim[]);
int FITSread_startHeader(char filename[]);
int FITSfind_headerLine(char keyword[], enum valueType type, void *pValue, char comment[]);
int FITSfind_LastHeaderLine(char keyword[], char text[], char lastLine[], int *nTimes);
int FITSread_headerLine(char keyword[], enum valueType type, void *pValue, char comment[], int mustBeNext);
int FITSread_endHeader(void);
int FITSread_array2D(double *pArray, int datasize, size_t nx, size_t ny);
int FITSread_endFile(void);
int FITSclose(void);
