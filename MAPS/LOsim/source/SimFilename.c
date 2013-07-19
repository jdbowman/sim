//SimFilename.c: Generate a standard filename based on the simulation name

// Author:	Peter Sherwood	iti@world.std.com	(617) 244-0836

// 7/11/01	Create

#include "SimFilename.h"
#include "GlobalReferences.h"
#include <string.h>

// Creates filenames like <simulation name>_<category>.<extension> 
// (e g, bigTest_Description.txt)
int SimFilename(char filename[], 
		char category[], 
		char extension[])
{
  strcpy(filename,simName); strcat(filename,"_"); strcat(filename,category); 
  strcat(filename,"."); strcat(filename,extension);
  return 0;
}

// Creates filenames like <category>_<name>.<extension> 
// (e g, SourceList_3sources.txt)
int CategoryFilename(char filename[], 
		     char category[], 
		     char name[], 
		     char extension[])
{
  strcpy(filename,category); strcat(filename,"_"); strcat(filename,name); 
  strcat(filename,"."); strcat(filename,extension);
  return 0;
}
