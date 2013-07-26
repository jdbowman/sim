#include "DataStructures.h"

int SourceListRead(int whatToDo, char sourceListFilename[], 
		   struct Source *pSource, int *pEOF);

// Values for the whatDoDo argument:
// Open the source list file and read all entries, 
// checking values (but not storing)
#define SLR_CHECK 1 

// Open the source list file and read, 
// check and return the first entry
#define SLR_FIRST 2 

#define SLR_NEXT  3
