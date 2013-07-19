// Utilities.c: Miscellaneous utility programs for LOsim

// Author:	Peter Sherwood	sherwood@computer.org	(617) 244-0836

// 8/4/01	split off from SourceListRead.c
// 8/22/01	add ErrorPrint

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <math.h>
#include "Utilities.h"
#include "Parameters.h"

// Functions dealing with time
// Validate year, month, etc, and store in a single long long result, 'time'
// Returns 0 if input is valid
int PackTime(int year, int month, int day, int hr, int min, int sec, int msec, long long *time)
{
	int daysInMo[12]={31,28,31,30,31,30,31,31,30,31,30,31}; // for a non-leap year
	
	// Validate the inputs
	if(year<1901 || year>2099) return 1;
	if(month<1 || month>12) return 2;
	if(year%4==0) daysInMo[1]++; // leap year => 29 days in February (valid from 1901-2099, since 2000 was a leap year)
	if(day<1 || day>daysInMo[month]) return 3;
	if(hr<0 || hr>23) return 4;
	if(min<0 || min>59) return 5;
	if(sec<0 || sec>59) return 6;
	if(msec<0 || msec>999) return 7;
	*time=(((((year*100LL+month)*100+day)*100+hr)*100+min)*100+sec)*1000+msec;
	return 0;
}

void UnpackTime(long long time, int *year, int *month, int *day, int *hr, int *min, int *sec, int *msec)
{
	*msec=time%1000; *sec=time/1000%100; *min=time/100000L%100; *hr=time/10000000L%100;
	*day=time/1000000000L%100; *month=time/100000000000L%100; *year=time/10000000000000L;
// CHECK FOR ERRONEOUS INPUT HERE
	return;
}


// Returns  time1-time0 in seconds
double DiffTime(long long time1, long long time0)
{
	int yr0,yr1,mo0,mo1,day0,day1,hr0,hr1,min0,min1,sec0,sec1,msec0,msec1,sec,msec,days0,days1;
	UnpackTime(time0,&yr0,&mo0,&day0,&hr0,&min0,&sec0,&msec0);
	UnpackTime(time1,&yr1,&mo1,&day1,&hr1,&min1,&sec1,&msec1);
	// Calculate day numbers for each date
	days0=DayNumber(yr0,mo0,day0); days1=DayNumber(yr1,mo1,day1);
	sec=(((days1-days0)*24+(hr1-hr0))*60+(min1-min0))*60+(sec1-sec0); msec=msec1-msec0;
	return sec+(msec/1000.);
}

int DayNumber(int year, int month, int day)
{
	int y,d;
	int cumDays[12]={0,31,59,90,120,151,181,212,243,273,304,334}; // number of days in non-leap-year which precede the first of the month
	y=year-1841; d=(y/4)*1461+((y%4)*365)+cumDays[month-1]+day; if(month>2 && y%4==3) d++;
	if(d>21608) --d; // 1900 was not a leap year
	return d;
}

// Add time in seconds to unpacked time
// Returns 0 if addition succeeded, 1 if incrementSec>=86400, 2 if year out of range 1901-2099
int AddTime(double incrementSec, int *year, int *month, int *day, int *hour, int *min, int *sec, int *msec)
{
	double fSec,iSec; int h,m,s,d;
	int daysInMo[12]={31,28,31,30,31,30,31,31,30,31,30,31}; // for a non-leap year
	
	fSec=modf(incrementSec,&iSec); if(iSec>=86400.) return 1; // input too large
	if(*year<1901 || *year>2099) return 2; // input year out of range
	(*msec)+=((int)(fSec*1000.0+0.5)); if((*msec)>=1000) {(*msec)-=1000; (*sec)++;}
	s=(int)iSec; m=s/60; s=s%60; h=m/60; m=m%60; d=h/24; h=h%24;
	(*sec)+=s; if((*sec)>=60) {(*sec)-=60; (*min)++;} // since input sec<60 and s<60, sec+carry+s<(59+60)
	(*min)+=m; if((*min)>=60) {(*min)-=60; (*hour)++;}
	(*hour)+=h; if((*hour)>=24) {(*hour)-=24; (*day)++;} // input hour<24 and h<24, so hour+carry+h<(23+24)
	else return 0; // no carry means we're done
	
	// Calculate number of days in February
	if((*year)%4==0) daysInMo[1]++; // leap year => 29 days (valid from 1901-2099, since 2000 was a leap year)
	if((*day)<=daysInMo[*month]) return 0;
	(*day)-=daysInMo[*month]; (*month)++; if((*month)<13) return 0;
	(*month)-=12; (*year)++; return 0;
}



extern int nErrors;

// Get a line of input, returning its length (-1 => EOF)
// lineIn=exact text   line=text w/comments, leading and trailing whitespace removed
int GetLine(FILE *inFile, char lineIn[], char line[], int maxLine)
{
	int i,n,commentSeen;
	fgets(lineIn,maxLine,inFile);
	if(feof(inFile)) {fclose(inFile); return -1;}
	n=strlen(lineIn);
	while((lineIn[n-1]=='\n')||(lineIn[n-1]=='\r')) {lineIn[n-1]=0; --n;} // remove a trailing return or newline
	strcpy(line,lineIn);
	// Remove a trailing comment
	commentSeen=0; i=0;
	while(i<n && !commentSeen)
	{
		switch(line[i++])
		{
			case('"'):	while(i<n) if(line[i++]=='"') break;
						break;
			case('/'):	commentSeen=line[i]=='/';
						break;
			default:	break;
		}
	}
	if(commentSeen) line[i-1]=0; // terminate the string at the start of the comment
	n=RemoveWhitespace(line); // remove the trailing and leading whitespace
	return n;
}

		
// Remove leading and trailing whitespace
int RemoveWhitespace(char line[])
{
	int i,j;
	i=0; while(line[i]==' ' || line[i]=='\t') i++; // i -> 1st non-blank character
	j=strlen(line)-1; while(j>=i && (line[j]==' ' || line[j]=='\t')) --j; // j -> last non-blank character or j<i
	if(j<i) return 0; // blank line
	line[j+1]=0; memmove(line,&line[i],j+1-i+1); // truncate the trailing and leading whitespace
	return j-i+1; // new length
}

/* See if a string matches a pattern
   Pattern elements are
	"<literal characters>", in which the character " is represented by \" (short for 1.1"<literal characters>")
	m.n"<literal characters>" from m>=0 through n>=m repetitions of the literal characters
	mD		exactly m>0 digits
	m.nD	from m>0 to n>=m digits
	m.D		m>=0 or more digits
*/

int PatternMatch(char string[], char pattern[])
{
	int i,j,k,c,m,n,lLit,dotSeen,seenM,seenN,patternLength,stringLength,iD; char literal[100]; char *p;
	int GetInteger(char **p); int GetString(char **p, char literal[]);
	void AdjustCounts(int dotSeen, int seenM, int seenN, int *pM, int *pN);
		
	if(!(patternLength=strlen(pattern))) return 0; // empty pattern matches nothing
	stringLength=strlen(string); // empty string may matches certain patterns, e g "" or D
	i=j=0; m=n=0; dotSeen=0; seenM=seenN=0;
	while(j<patternLength || i<stringLength)
	{
		if((c=pattern[j++])=='"')
		{ // advance in both pattern and string--chars must be identical
			p=&pattern[j];
			lLit=GetString(&p,literal);
			AdjustCounts(dotSeen,seenM,seenN,&m,&n);
			// Find the string m times (mandatory)
			for (iD=0,k=0; iD<m; iD++) while(k<lLit && i<stringLength) if(literal[k++]!=string[i++]) return 0; // mismatch-pattern fails
			// Find the string up to (n-m) more times (optional)
			for (iD=m,k=0; iD<n; iD++) while(k<lLit && i<stringLength) if(literal[k++]!=string[i++]) {i-=k; break;}; // mismatch-go back
			j=p-pattern;
			if(i>stringLength || j>patternLength) return 0; // string didn't contain entire literal || pattern had no close "
			m=n=0; dotSeen=0; seenM=seenN=0;
			continue;
		}
		else if(c=='D')
		{
			AdjustCounts(dotSeen,seenM,seenN,&m,&n);
			// Try to match from m to n digits
			for (iD=0; iD<m; iD++) if(i>=stringLength || !isdigit(string[i++])) return 0; // first match the m digits (mandatory)
			for (iD=m; iD<n; iD++) if(i>=stringLength || !isdigit(string[i++])) {i=i-1; break;} // skip up to (n-m) more digits (optional)
			m=n=0; dotSeen=0; seenM=seenN=0;
			continue;
		}
		else if(isdigit(c))
		{
			p=&pattern[j-1];
			if(!dotSeen) {m=GetInteger(&p); seenM=1;} else {n=GetInteger(&p); seenN=1;}
			j=p-pattern;
		}
		else if(c=='.')
		{
			if(dotSeen) return 0; // extra . in pattern
			dotSeen=1;	// next number is n
		}
		else return 0; // illegal pattern
	}

	
	return 1;
}
			
// Get a string literal starting at p
// On return, p->char following close "
int GetString(char **p, char literal[])
{
	int n=0; char c;
	while(**p && (c=*((*p)++))!='\"')
	{
		if(c=='\\')
		{
			if(**p) c=*((*p)++);
			// Escaped characters are \ and " only-others are copied, though
			else break; // trailing \ is dropped
		}
		literal[n++]=c;
	}
	literal[n]=0;
	return n;
}

// Get a non-negative integer starting at p
// On return, p->1st non-digit
int GetInteger(char **p)
{
	int k,n=0; char c;
	while(isdigit(c=**p))
	{
		(*p)++;
		k=c-'0'; // theoretically this is non-portable
		n=n*10+k;
	}
	return n;
}
			
void AdjustCounts(int dotSeen, int seenM, int seenN, int *pM, int *pN)
{
	if(!seenM) {*pM=*pN=1;} // special case for literal without counts
	else if(!dotSeen) *pN=*pM; // no '.' seen: mD
	else if(!seenN) *pN=INT_MAX; // .D or m.D
	return;
}

/*******************************************************************
      Functions used by the ASCII input file parsers
********************************************************************/
// Replace all multiple blanks and tabs with a single blank
void MultipleBlankRemove(char s[])
{
	int i=0;
	while(i<strlen(s))
	{
		if(s[i]==' ' || s[i]=='\t')
		{
			if((s[i+1]==' ' || s[i+1]=='\t'))
			{
				memmove(&s[i+1],&s[i+2],strlen(s)-(i+1)); --i;
			}
		}
		i++;
	}
	return;
}

// Convert a string to lower case
void LowerCase(char s[])
{
	int i,l=strlen(s);
	for (i=0; i<l; i++) s[i]=(char)tolower((int)s[i]);
	return;
}

// See if the field is the same as the name
int FieldMatches(char field[], const char name[])
{
	char fieldTest[MAX_FIELDNAME+1],nameTest[MAX_FIELDNAME+1]; int i,j,c,lField,lName;
	
	if((lField=strlen(field))>MAX_FIELDNAME || (lName=strlen(name))>MAX_FIELDNAME) return 0; // don't deal with bad input
	
	// Compare after removing all punctuation and making upper-case
	for (i=0,j=0; i<lField; i++) {c=field[i]; if(isalnum(c)) fieldTest[j++]=toupper(c);}
	fieldTest[j]=0;
	for (i=0,j=0; i<lName; i++) {c=name[i]; if(isalnum(c)) nameTest[j++]=toupper(c);}
	nameTest[j]=0;
	return strcmp(fieldTest,nameTest)==0;
}

// Validate and copy a string value
int StringValue(char var[], const char value[], int maxLen)
{
	if(strlen(value)>maxLen) return 1;
	strcpy(var,value);
	return 0;
}

// Validate and copy a date and time
// Returns 0=no error in input
int DateAndTime(long long *pDateTimeVar, const char value[])
{
	int mo=0,day=0,yr=0,hr=0,min=0,sec=0,ms=0;
	sscanf(value,"%d-%d-%d %d:%d:%d.%d",&yr,&mo,&day,&hr,&min,&sec,&ms);
	return PackTime(yr,mo,day,hr,min,sec,ms,pDateTimeVar);
}

// Validate and copy a long integer value
int IntegerValue(long *pVar, const long loLimit, const long hiLimit, const char value[])
{
	long val; int n;
	n=sscanf(value,"%ld",&val);
	if(n!=1) return 1; // no valid integer found
	if(val<loLimit || val>hiLimit) return 2; // not within limits
	*pVar=val;
	return 0;
}

// Validate and copy a double value
int DoubleValue(double *pVar, const double loLimit, const double hiLimit, const char value[])
{
	double val; int n;
	n=sscanf(value,"%le",&val);
	if(n!=1) return 1; // no valid real found
	if(val<loLimit || val>hiLimit) return 2; // not within limits
	*pVar=val;
	return 0;
}
	
void CheckStr(char field[], char fieldname[])
{
	char temp[100];
	if(!strlen(field)) {sprintf(temp,"The field '%s' must be defined in the description file",fieldname); ErrorPrint(temp);}
	return;
}
	
void CheckInt(long field, char fieldname[])
{
	char temp[100];
	if(field<0) {sprintf(temp,"The field '%s' must be defined in the description file",fieldname); ErrorPrint(temp);}
	return;
}
	
void CheckDbl(double field, char fieldname[])
{
	char temp[100];
	if(field<0.0L) {sprintf(temp,"The field '%s' must be defined in the description file",fieldname); ErrorPrint(temp);}
	return;
}


void ErrorPrint(char message[])
{
	fflush(stdout);
	fprintf(stderr,"\n*** ERROR: %s",message);
	nErrors++;
	return;
}
