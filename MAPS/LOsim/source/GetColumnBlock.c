// GetColumnBlock.c: Get a single block of column of row-transformed brightness data

// Author:	Peter Sherwood	iti@world.std.com	(617) 244-0836

// 5/1/01	Create
// 7/10/01	v1.0

#include "Parameters.h"
#include "DataStructures.h"
#include "LargeFiles.h" // includes <stdio.h>, which must come before GetColumnBlock.h
#include "GetColumnBlock.h"
#include "GlobalReferences.h"
#include <string.h>

// Get a single block of column of row-transformed brightness data, representing columns j0-(j0+ELEMENTSPERBLOCK-1), all from row i
// The input matrix is conceptually visibilityXN x visibilityYN, but the first brightnessYNzero0 and last brightnessYNzero1 rows are all 0

int GetColumnBlock(FileHandle inFileL, long i, long j0, struct Complex s[])
{
	size_t nr; long long filePosition;
	if(i>=brightnessYNzero0 && i<(visibilityYN-brightnessYNzero1))
	{
		filePosition=((long long)(i-brightnessYNzero0)*visibilityXN+j0)*sizeof(struct Complex);
		fseekL(inFileL,filePosition,SEEK_SET);
		nr=freadL(s,sizeof(struct Complex),ELEMENTSPERBLOCK,inFileL);
	}
	else
	{
		nr=ELEMENTSPERBLOCK*sizeof(struct Complex); // break into 2 stmts to check w/debugger
		memset(s,0,nr); // in an all-zero row--return 0s
	}

	return 0;
}
