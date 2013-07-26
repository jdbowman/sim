// RefreshRow.c: Read a row from the visibility file into the cache after writing any previous contents

// Author:	Peter Sherwood	iti@world.std.com	(617) 244-0836

// 7/4/01	Create

#include "RefreshRow.h" // includes LargeFiles.h
#include "Parameters.h"
#include "ArrayDimensions.h"
#include "GlobalReferences.h"
#include <math.h>

int RefreshRow(FileHandle visibilityFH, struct Complex rowCache[], long *pICached, long iWanted)
{
	long n;
	long long filePosition;
	// If there is an old row in the cache, write it
	if(*pICached>=0)
	{
		filePosition=((long long)(*pICached)*visibilityXN+0)*sizeof(struct Complex);
		fseekL(visibilityFH,filePosition,SEEK_SET);
		n=fwriteL(rowCache,sizeof(struct Complex),visibilityXN,visibilityFH);
		if(n!=visibilityXN) printf("** ERROR Wrote %ld (should be %ld **)",n,visibilityXN);
	}
	// If a new row is wanted, read it
	if(iWanted>=0)
	{
		filePosition=((long long)iWanted*visibilityXN+0)*sizeof(struct Complex);
		fseekL(visibilityFH,filePosition,SEEK_SET);
		n=freadL(rowCache,sizeof(struct Complex),visibilityXN,visibilityFH);
		if(n!=visibilityXN) printf("** ERROR Wrote %ld (should be %ld **)",n,visibilityXN);
		else *pICached=iWanted;
	}
	return 0;
}


