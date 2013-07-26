// DescriptionFileRead.c: Read the description file for this simulation

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// 7/6/01	Create
// 8/22/01	Move ra, dec to observe file; and observation name
// 9/19/01	Allow small brightness grids, and give proper error messages
// 11/26/01	remove field size (now in the observation file)
// 12/27/01	v1.11: debug options

#include <stdio.h> // must appear before Utilites.h
#include "GlobalReferences.h"
#include "Utilities.h"
#include "Parameters.h"
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

// 0=description file read OK
// 1=file could not be opened
// 2=error in parsing file

#define MAX_FIELDNAME 50

/*---------------------------------------------------------------------------*/

int DescriptionFileRead(char simName[], char filename[])
{
#define MAX_LINE 150
  char lineIn[MAX_LINE], line[MAX_LINE], *p, fieldname[MAX_FIELDNAME+1];
  char value[L_EXTENDEDLINE], temp[100]; 
  int i, ln, nLine=0, error;
  FILE *descFile;
	
  /* To add a field
     1) Increase NFIELDNAMES
     2) Add the string the user puts in the description file at the 
        end of the list below
     3) Add code to process the value in the switch statement below
  */

#define NFIELDNAMES 22
  const char *name[NFIELDNAMES]={	// Files marked 'M' are mandatory; others may be omitted, or included with blank values
    "Title",						// M short phrase describing this simulation (e g, "Random 245 sources at full resolution")
    "Simulation name",				// M single word describing this simulation (e g, "Rand245")
    "When created",					// mm/dd/yy hh:mm:ss when this description file was first created
    "Who created",					// the name of the author of this simulation description file
    "Institution",					// the institutional affiliation of the author of this simulation description file
    "E-mail",						// the e-mail address of the author of this simulation description file
    "Phone",						// the telephone number of the author of this simulation description file
    "When last changed",			// similar information for the person who last modified this description file
    "Who last changed",
    "Description",					// a one-line description of the simulation, providing additional information from the title
    "Comments",						// any further remarks about this simulation
    "Brightness grid x-N",			// M number of points in the sky plane grid, x
    "Brightness grid y-N",			// M                                         y
    "Visibility grid x-log2N",		// M log2 of number of points in the visibility grid, x
    "Visibility grid y-log2N",		// M                                                  y
    "Source list name",				// M in-beam source list is in file SourceList_<name>.txt
    "Out-of-beam source list name", // M out-of-beam source list is in file OOBSourceList_<name>.txt
    "Station list name",			// M station list is in file StationList_<name>.txt
    "Ionosphere parameters name",	// M parameters for model of ionosphere are in file Ionosphere_<name>.txt
    "Observation name",				// M observing parameters are in file Observation_<name>.txt
    "Image name",					// external image parameters are in Image_<name>,txt
    "Debug option"					//   used for all debugging options; value is the name of the (boolean) option
  };
	
  // Open the description file
  strcpy(filename,simName); strcat(filename,"_Description.txt");
  if((descFile=fopen(filename,"r"))==NULL)
    {
      sprintf(temp,"Can't find the description file '%s'",filename);
      ErrorPrint(temp);
      return 1;
    }
	
  // Blank all the fields (except simName)
  title[0]=creator[0]=institution[0]=email[0]=phone[0]=changer[0]=description[0]=comments[0]=sourceListName[0]='\0';
  OOBSourceListName[0]=stationListName[0]=ionosphereParametersName[0]=imageParametersName[0]='\0';
  createDateTime=changeDateTime=0;
  brightnessXN = brightnessYN = visibilityXLn2 = visibilityYLn2 = -1;
  fieldXSize = fieldYSize = fieldRA = fieldDec = obsFreq = -1.0L;
  nErrors = 0;
	
  // Process the lines of the description file
  while(1)
    { // get a line from the file, stopping on error or EOF
      if((ln = GetLine(descFile, lineIn, line, MAX_LINE))<0) break; 
      nLine++;
      printf("\n%3d: %s",nLine,lineIn);
		
      // Remove leading and trailing blanks and tabs, and ignore blank lines or lines beginning with "//"
      if(!ln) continue; // skip a blank line
		
      /* Lines in the description file have the format
	 <fieldname> = <value>
	 The blanks around the "=" are optional
      */
      if((p=strchr(line,'='))==NULL)
	{
	  ErrorPrint("Line format should be <fieldname> = <value>");
	  continue; // ignore this line
	}
      strncpy(fieldname,line,p-line); 
      fieldname[p-line]=0; RemoveWhitespace(fieldname);
      strcpy(value,p+1);   RemoveWhitespace(value); // skip the '='
		
      // Identify the field
      for (i=0; i<NFIELDNAMES; i++)
	{
	  if(FieldMatches(fieldname,name[i])) break;
	}
		
      // Get the field value and store it
      if(i<NFIELDNAMES)
	{
	  switch(i)
	    {
	    case(0): // Title
	      if(StringValue(title,value,L_ONELINE)) ErrorPrint("Title is too long");
	      break;
	    case(1): // Simulation name
	      if(StringValue(temp,value,L_BRIEFNAME)) ErrorPrint("Simulation name is too long");
	      if(strcmp(simName,temp)) ErrorPrint("Simulation name field in description file doesn't match that in the name of the file");
	      break;
	    case(2): // When created
	      DateAndTime(&createDateTime,value); break;
	    case(3): // Who created
	      if(StringValue(creator,value,L_NORMAL)) ErrorPrint("Creator name is too long");
	      break;
	    case(4): // Institution
	      if(StringValue(institution,value,L_NORMAL)) ErrorPrint("Institution name is too long");
	      break;
	    case(5): // E-mail
	      if(StringValue(email,value,L_NORMAL)) ErrorPrint("E-mail address is too long");
	      break;
	    case(6): // Phone
	      if(StringValue(phone,value,L_BRIEFNAME)) ErrorPrint("Phone is too long");
	      break;
	    case(7): // When last changed
	      DateAndTime(&changeDateTime,value); break;
	    case(8): // Who last changed
	      if(StringValue(changer,value,L_NORMAL)) ErrorPrint("Who last changed name is too long");
	      break;
	    case(9): // Description
	      if(StringValue(description,value,L_EXTENDEDLINE)) ErrorPrint("Description is too long");
	      break;
	    case(10): // Comments
	      if(StringValue(comments,value,L_EXTENDEDLINE)) ErrorPrint("Comments string is too long");
	      break;
	    case(11): // Brightness grid x-size
	      if((error=IntegerValue(&brightnessXN,4L,30001,value)))
		{
		  if(error==1) ErrorPrint("No valid integer found");
		  if(error==2) ErrorPrint("Not within limits");
		}
	      break;
	    case(12): // Brightness grid y-size
	      if((error=IntegerValue(&brightnessYN,4L,30001,value)))
		{
		  if(error==1) ErrorPrint("No valid integer found");
		  if(error==2) ErrorPrint("Not within limits");
		}
	      break;
	    case(13): // Visibility grid x-size
	      if((error=IntegerValue(&visibilityXLn2,2L,15L,value)))
		{
		  if(error==1) ErrorPrint("No valid integer found");
		  if(error==2) ErrorPrint("Not within limits");
		}
	      break;
	    case(14): // Visibility grid y-size
	      if((error=IntegerValue(&visibilityYLn2,2L,15L,value)))
		{
		  if(error==1) ErrorPrint("No valid integer found");
		  if(error==2) ErrorPrint("Not within limits");
		}
	      break;
	    case(15): // Source list name
	      if(StringValue(sourceListName,value,L_BRIEFNAME)) ErrorPrint("Observation parameters name is too long");
	      break;
	    case(16): // Out-of-beam source list name
	      if(StringValue(OOBSourceListName,value,L_BRIEFNAME)) ErrorPrint("Observation parameters name is too long");
	      break;
	    case(17): // Station list name
	      if(StringValue(stationListName,value,L_BRIEFNAME)) ErrorPrint("Observation parameters name is too long");
	      break;
	    case(18): // Ionosphere parameters name
	      if(StringValue(ionosphereParametersName,value,L_BRIEFNAME)) ErrorPrint("Observation parameters name is too long");
	      break;
	    case(19): // observation parameters name
	      if(StringValue(observationName,value,L_BRIEFNAME)) ErrorPrint("Observation parameters name is too long");
	      break;
	    case(20): // external image parameters name
	      if(StringValue(imageParametersName,value,L_BRIEFNAME)) ErrorPrint("Image parameters name is too long");
	      break;
	    case(21): // debugging option
	      // Compare the debugging option name to, using lower case
	      if(StringValue(temp,value,L_BRIEFNAME)) ErrorPrint("Debugging option name is too long");
	      MultipleBlankRemove(temp); LowerCase(temp);
	      if(!strcmp(temp,"no source truncation")) {debugOption_noSourceTruncate=1; break;}
	      // ADD MORE DEBUGGING OPTIONS HERE
	      ErrorPrint("Debugging option name unknown");
	      break;
	    }
	}
      else
	{
	  sprintf(temp,"Field name '%s' is not known",fieldname);
	  ErrorPrint(temp);
	}
    }
	
  // Check that all the mandatory fields have been specified
  CheckStr(title,"Title"); CheckStr(sourceListName,"Source list name"); CheckStr(OOBSourceListName,"Out-of-beam source list name");
  CheckStr(observationName,"Observation name");
  CheckStr(stationListName,"Station list name"); CheckStr(ionosphereParametersName,"Ionosphere parameters name");
  CheckInt(brightnessXN,"Brightness grid x-N"); CheckInt(brightnessYN,"Brightness grid y-N");
  CheckInt(visibilityXLn2,"Visibility grid x-log2N"); CheckInt(visibilityYLn2,"Visibility grid y-log2N");

  // Check that the field values are compatible
  if(1<<visibilityXLn2 < brightnessXN) ErrorPrint("Visibility grid x must be at least as large as brightness grid x");	
  if(1<<visibilityYLn2 < brightnessYN) ErrorPrint("Visibility grid y must be at least as large as brightness grid y");
	
  // Do conversions to radians
  fieldXSize *= oneArcsec; 
  fieldYSize *= oneArcsec; 
  fieldDec   *= oneArcsec; 
  fieldRA    *= oneArcsec;
		
  if(nErrors) return 100+nErrors;
  else return 0;
}

