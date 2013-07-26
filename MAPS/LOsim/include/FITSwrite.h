// FITSwrite.h

#include "FITSread.h"
#include <stddef.h> // for size_t

int FITSwrite_startHeader(char filename[]);
int FITSwrite_headerLine(char keyword[], enum valueType type, void *pValue, char comment[]);
int FITSwrite_endHeader(void);
int FITSwrite_array2D(void *pArray, size_t nx, size_t ny, size_t dsize);
int FITSwrite_endFile(void);

