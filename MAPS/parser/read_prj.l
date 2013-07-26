
%option yylineno noyywrap

%{
/* 
 * Parsing MAPS simulation project file 
 * Created 25 Jan 2011 by L. Benkevitch
 */

#include <stdio.h>
#include <string.h> 
#include "read_prj.tab.h"
#include "read_prj.h"
#include "brigen.h"
#include "visgen.h"

/* links to brigen and visgen parameter structures */  
extern struct bg_param *bgp;
extern struct vg_param *vgp; 

/* enum inptype {FNone = 0, FProject = 1, FSrclist, FObsspec, FArray}; */

/*
 * Linked list of the saved inpit files
 */
struct bufstack {
  struct bufstack *prev; /* previous entry */
  YY_BUFFER_STATE bs; /* saved buffer */
  int             lineno; /* saved line number */
  char            *filename; /* name of this file */
  FILE            *f; /* current file */
  enum inptype 	  ftype; /* type of include file: srclist, obsspec etc. */
} *curbs = 0;


char *curfilename = 0;  /* name of current input file */
//int newfile(enum inptype ftype, char *fn);
//enum inptype popfile(void);

/* NOTE: T_angle token is returned with the character value in yylval.s: 
 * angle and units. The first char, yylval.s[0], contains the info
 * about the angle format: 'a' - deg, ', "; 'x' - dd:mm:ss, hh:mm:ss 
 */


char *tl = NULL;

%}

D		[0-9]
L		[a-zA-Z_]
E		[Ee][+-]?{D}+
NAME		{L}({L}|{D})*
FLOAT0		[+-]?{D}+{E}	
FLOAT1		[+-]?{D}*"."{D}+{E}?
FLOAT2		[+-]?{D}+"."{D}*{E}?
FLT		{FLOAT0}|{FLOAT1}|{FLOAT2}
INT		[+-]?{D}+
I12		[0-9]{1,2}
I13		[0-9]{1,3}
YEAR		[0-9]{4}
NUM		{FLT}|{INT}
ANGLE		{NUM}[ \t]*("\""|"''"|"'"|"deg"|"rad")
FREQ		{NUM}[ \t]*("mHz"|"Hz"|"kHz"|"MHz"|"GHz")
DUR		{NUM}[ \t]*("us"|"usec"|"ms"|"msec"|"s"|"sec"|"min"|"hr")
XMS		[+-]?{I13}":"{I12}":"{I12}("\."{INT}?)?[ \t]*("dms"|"hms")
YDHMS		{YEAR}":"{I13}":"{I12}":"{I12}":"{I12}("\."{INT}?)?
ABSTIME		{YDHMS}[ \t]*"ydhms"	
FNAME		[0-9a-zA-Z\.\_\-\+\/]+

%x SRCLIST OBSSPEC ARRAY

%%

^[ \t]*"srclist"[ \t]+	   { BEGIN SRCLIST; }
<SRCLIST>{FNAME}  {
		     yylineno++;
		     { int c; while((c = input()) && c != '\n'); }
		     BEGIN INITIAL;		  
		     if (strcmp(yytext, "here") == 0) {
		       return T_srclist_here;
		     }
		     /* Only if called from brigen, dive into the srclist */
		     if (bgp && bgp->is_set.srclist != 1) {
		       bgp->srcfilename = yytext;
		       if (!newfile(FSrclist, yytext))
		         yyterminate(); /* no such file */
		       return T_srclist_here; /* imitate 'here' */
		     }
		   }

^[ \t]*"obsspec"[ \t]+	   { BEGIN OBSSPEC; }
<OBSSPEC>{FNAME}  {
		     yylineno++;
		     { int c; while((c = input()) && c != '\n'); }
		     BEGIN INITIAL;		  
		     if (strcmp(yytext, "here") == 0) {
		       return T_obsspec_here;
		     }
		     /* Only if called from visgen, dive into the obsspec */
		     if (vgp) {
		       if(!newfile(FObsspec, yytext))
		         yyterminate(); /* no such file */
		       return T_obsspec_here; /* imitate 'here' */
		     }
		   }

^[ \t]*"array"[ \t]+	   { BEGIN ARRAY; }
<ARRAY>{FNAME}  {
		     yylineno++;
		     { int c; while((c = input()) && c != '\n'); }
		     BEGIN INITIAL;		  
		     if (strcmp(yytext, "here") == 0) {
		       return T_array_here;
		     }
		     /* Only if called from visgen, dive into the array spec */
		     if (vgp) {
		       if(!newfile(FArray, yytext))
		         yyterminate(); /* no such file */
		       return T_array_here; /* imitate 'here' */
		     }
		   }

<<EOF>>		   { enum inptype ftype;
                     int nomorefiles = popfile(&ftype);
                     if (nomorefiles) {
		       curbs = 0; /* Must be NULL to run yylex() again! */
		       yyterminate(); //==============>>>
		     }
                     /* printf("EOF: ftype = %d\n", ftype); */
		     /* printf("EOF: curbs = %p\n", curbs); */
		     switch (ftype) {
		       case FNone:     break; 
                       case FProject:  break; 
		       case FSrclist:  return T_end_srclist;
		       case FObsspec:  return T_end_obsspec;		 
		       case FArray:    return T_end_array;
		     }		 
		   }

[ \t]+
"/*"		   { comment_c(); }
("//"[^\n]*[\n]+)|("#"[^\n]*[\n]+)  {unput('\n');
				    /*printf("COMMENT: %s", yytext);*/}

"simulation"       return T_simulation;
"section"    	   return T_section;
"end"[ \t]+"srclist"	  { return T_end_srclist; }
"end"[ \t]+"obsspec"	  { return T_end_obsspec; }
	  
"end"[ \t]+"array"	  return T_end_array;
"end"        	   { return T_end;}

"brigen"     	   return T_brigen;
"imgrid_pixels"    return T_imgrid; 
"imsize_radec" 	   return T_imsize;
"imcenter_radec"   return T_imcenter;
"imfitsout"        return T_imfitsout;
"addimage"         return T_addimage;
"visgen"     	   return T_visgen;
"array"      	   return T_array;
"visout"     	   return T_visout;
"obsspec"    	   return T_obsspec;
"textout"    	   return T_textout;
"FOV_center_RA"	   return T_FOV_center_RA;
"FOV_center_Dec"   return T_FOV_center_Dec;
"FOV_size_RA"	   return T_FOV_size_RA;
"FOV_size_Dec"	   return T_FOV_size_Dec;
"Corr_int_time"	   return T_Corr_int_time;
"Corr_chan_bw"	   return T_Corr_chan_bw;
"Time_cells"	   return T_Time_cells;
"Freq_cells"	   return T_Freq_cells;
"Scan_start"	   return T_Scan_start;
"Scan_duration"	   return T_Scan_duration;
"Central_freq"	   return T_Central_freq;
"Bandwidth"	   return T_Bandwidth;
"Endscan"	   return T_Endscan;
"PNT_center_RA"	   return T_PNT_center_RA;
"PNT_center_DEC"   return T_PNT_center_DEC;
"GHA"   	   return T_GHA;
"ydhms"   	   return T_ydhms;

","      	   return ',';
[\n]+	     	   return '\n';
"="	     	   return '=';


{NAME}		   { yylval.s = strdup(yytext); return T_name; }
{ANGLE}		   { yylval.s = strdup(yytext); return T_angle; }
{XMS}		   { yylval.s = strdup(yytext); return T_xms; }
{ABSTIME}	   { yylval.s = strdup(yytext); return T_abstime; }
{FREQ}		   { yylval.s = strdup(yytext); return T_freq; }
{DUR}		   { yylval.s = strdup(yytext); return T_dur; }
{FLT}		   { yylval.d = atof(yytext); return T_float; }
{INT}		   { yylval.n = atoi(yytext); return T_int; }
{FNAME}		   { yylval.s = strdup(yytext); return T_file; }

.           	   { printf("\nError: Illegal character '%s', line %d\n" \
			    "in file %s\n", yytext, yylineno, curfilename);
                     exit(1);
		   }

%%

void comment_c(void) { /* Skip C-style comment */
  char c = input(), c1 = 0;
  while (!(c1 == '*' && c == '/') && c != 0) {
    c1 = c;
    c = input();
  }
}


//void comment(void) { /* "//" or "#" -- up to eol */
//  char c;
//  do { 
//    c = input();
//  }
//  while (c != '\n' && c != 0);
//  unput(c);   /* put '\n' back to input stream */  
//}

/*********************************************
 * Skip leading whitespace characters,       *
 * '\n' and tabs '\t', in string str[]       *
 * IN PLACE.                                 *
 * Returns pointer to the first non-         *
 * whitespace character in str[].            *
 * Note: the routine works correctly ONLY if *
 * str[] is GUARANTEED to contain at least   *
 * one non-whitespace character.             *  
 * Created 25 Jan 2011 by L. Benkevitch      *
 *********************************************/
char *skiw(char* str) {
  while(*str == ' ' || *str == '\t') str++;
  return str; 
}

/****************************************************
 * newfile() prepares to read from the file named   *
 * fn, whose type can be one of the enum inptype.   *
 * It saves any previous input file and its         *
 * attributes, including the file type ftype. It    *
 * maintains a linked list of bufstack structures,  *
 * each of which has a link to the previous         *
 * bufstack along with the saved ftype, yylineno    *
 * and filename.                                    *
 * It opens the file; creates and switches to a     *
 * flex buffer; and saves the previous open file,   *
 * filename, file type, and buffer.                 *
 * Created 02 Feb 2011 by L. Benkevitch             * 
 ****************************************************/
 
/* enum inptype {FNone = 0, FProject = 1, FSrclist, FObsspec, FArray}; */

int newfile(enum inptype ftype, char *fname) {
  FILE *f = fopen(fname, "r");
  struct bufstack *bs = malloc(sizeof(struct bufstack));

  /* die if no file or no room */
  /* if(!f) { perror(fname); return 0; } // Expect input from stdin */
  if(!f)  { perror(fname); exit(1); }
  if(!bs) { perror("malloc"); exit(1); }

  /* remember state */
  if (curbs) curbs->lineno = yylineno;
  bs->prev = curbs;

  /* set up current entry */
  bs->bs = yy_create_buffer(f, YY_BUF_SIZE);
  bs->f = f;
  bs->filename = fname;
  bs->ftype = ftype;

  /* printf("NEWFILE: ftype = %d\n", ftype); */
  /* printf("NEWFILE: curbs = %p\n", curbs); */

  yy_switch_to_buffer(bs->bs);
  curbs = bs;
  yylineno = 1;
  curfilename = fname;
  return 1;
}

/****************************************************
 * popfile() undoes what newfile did. It closes the *
 * open file, deletes the current flex buffer, and  *
 * then restores the buffer, filename, and line     *
 * number from the prior stack entry. It returns 1  *
 * if that was the last file, and 0 if there still  *
 * are files on stack. Also, it sets the            *
 * enum inpfile *ftype parameter to the type of the *
 * previous file.                                   *
 *   Note that it doesn’t call yyrestart() when it  *
 * restores the prior buffer; if it did, it would   *
 * lose any input that had already been read into   *
 * the buffer.                                      *
 * Created 02 Feb 2011 by L. Benkevitch             * 
 ****************************************************/

/* enum inptype {FNone = 0, FProject = 1, FSrclist, FObsspec, FArray}; */

int popfile(enum inptype *ftype) {
  struct bufstack *bs = curbs;
  struct bufstack *prevbs;

  if(!bs) return 1; /* Empty buffer stack */
  /* get rid of current entry */
  fclose(bs->f);
  yy_delete_buffer(bs->bs);

  /* switch back to previous */
  prevbs = bs->prev;
  *ftype = bs->ftype;

  /* printf("POPFILE: ftype = %d\n", *ftype); */
  /* printf("POPFILE: curbs = %p\n", curbs); */
  /* printf("POPFILE: prevbs = %p\n", prevbs); */

  free(bs);
  if(!prevbs) return 1; /* No more files */
  yy_switch_to_buffer(prevbs->bs);
  curbs = prevbs;
  yylineno = curbs->lineno;
  curfilename = curbs->filename;
  return 0;
}
