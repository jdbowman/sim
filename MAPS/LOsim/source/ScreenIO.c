#include "ScreenIO.h"
#include <stdio.h>
#include <string.h>

char textBuffer[MAX_LINES][LINE_LENGTH];
int textBufferLines;	// number of lines now in the text buffer
int textBufferSize;		// capacity of the text buffer, in lines
int lineSpacing, charHeight,charsPerLine,linesPerScreen;
int pos;				// line number currently at the top of the screen

// Initialize the text buffer
void ScreenInit(void)
{
	textBufferLines=0;	// number of lines now in the text buffer
	textBufferSize=MAX_LINES;		// capacity of the text buffer, in lines
	pos=0;
#if 0
	//DEBUG: make some text to play with
	textBufferLines=500;
	for (i=0; i<textBufferLines; i++) sprintf(textBuffer[i],"This is line %d",i);
#endif
	return;
	}
	
// Add a string to the screen
void OutputToScreen(char x[], int sameLine)
{
	extern FILE *outFile; char xLF[200]; int i;
	if(sameLine) xLF[0]=0; else strcpy(xLF,"\n"); // start with NL to make file readable
	strcat(xLF,x);
	if(outFile!=NULL) fwrite(xLF,sizeof(char),strlen(xLF),outFile);
	if(textBufferLines>=textBufferSize) // lose least-recent screenful (MAKE CIRCULAR BUFFER)
	{
		for (i=linesPerScreen; i<textBufferSize-1; i++) strcpy(textBuffer[i-linesPerScreen],textBuffer[i]);
		textBufferLines=textBufferSize-linesPerScreen;
	}
	strcpy(textBuffer[textBufferLines++],x);
	if(sameLine) printf("%s",x);
	else printf("\n%s",x);
	fflush(stdout);
	return;
}

// Skip a line
void SkipLine(void)
{
	char x[1]={0};
	OutputToScreen(x,0);
	return;
}
