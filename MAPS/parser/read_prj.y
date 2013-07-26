/*
 * read_prj.y
 * The grammar in BNF for the MAPS simulation project parser.
 * Created 25 Jan 2011 by L. Benkevitch
 * Modified 17 Feb 2011
 */

%{
#ifndef YYDEBUG
#  define YYDEBUG 0
#endif
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "read_prj.h"
#include "brigen.h"
#include "sky_source.h"

extern int yydebug;

int yylex(void);

/* enum inptype {FNone = 0, FProject = 1, FSrclist, FObsspec, FArray}; */

//extern double const pi; /* 3.1415926535897931; */
extern double const pi; /* 3.1415926535897931; */

extern char *yytext;
extern int yylineno;
extern char *curfilename;
extern int isrc;  /* Index into the list of celestial sources */
extern int iimgf; /* Index into the list of fits files with external images */

/* links to brigen and visgen parameter structures */  
extern struct bg_param *bgp;
extern struct vg_param *vgp; 


%}

/*
 * Parse stack element
 */
%union {
  char  *s;
  int    n;
  double d;
}

%token <s> T_name 
%token <s> T_file 
%token <d> T_float 
%token <n> T_int 
%token <s> T_angle T_xms T_abstime T_freq T_dur
%token T_simulation T_section T_end T_brigen 
%token T_imgrid T_imsize T_imcenter T_imfitsout T_addimage
%token T_srclist T_srclist_here  T_end_srclist
%token T_visgen T_visout T_textout 
%token T_obsspec T_obsspec_here T_end_obsspec
%token T_array  T_array_here  T_end_array
%token T_GHA T_ydhms 
%token T_FOV_center_RA T_FOV_center_Dec T_FOV_size_RA T_FOV_size_Dec
%token T_Corr_int_time  T_Corr_chan_bw T_Time_cells T_Freq_cells
%token T_Scan_start  T_Scan_duration T_Central_freq T_Bandwidth T_Endscan
%token T_PNT_center_RA T_PNT_center_DEC

%type <s> filename 
%type <s> angle
%type <d> num

%glr-parser
%expect-rr 3

%start combined

%%

combined
	: project
	| src_lines       { if (bgp) bgp->is_set.srclist = 1; }
        | obsspec         // { if (vgp) vgp->is_set.obsspec = 1; }  
        | array_lines     // { if (vgp) vgp->is_set.array = 1; }  
	;
 
project
	: /* empty */
	| project '\n'
        | project T_name T_simulation '\n' simulation T_end T_name '\n' 
            { if (strcmp($2, $7)) { 
                *yytext = 0; 
                printf("Error: '%s simulation' and 'end %s' do not match\n", 
		       $2, $7);
                printf("at line %d or earlier in file '%s'\n",
		       yylineno, curfilename);
		YYERROR;
              }
            }
	;

simulation
	: /* empty */
	| simulation '\n'
	| simulation brigen_section  
	| simulation visgen_section  
	;

brigen_section
	: T_brigen T_section '\n' brigen_statements T_end T_brigen '\n'
	;


brigen_statements
	: /* empty */
	| brigen_statements '\n'
        | brigen_statements T_imgrid T_int  T_int '\n' 
                            { if (bgp && set_imgrid($3, $4)) YYERROR; }
        | brigen_statements T_imgrid T_int  ',' T_int '\n' 
	                    { if (bgp && set_imgrid($3, $5)) YYERROR; }
        | brigen_statements T_imcenter angle angle '\n' 
	                    { if (bgp && set_imcenter($3, $4)) YYERROR; }
        | brigen_statements T_imcenter angle ',' angle '\n' 
	                    { if (bgp && set_imcenter($3, $5)) YYERROR; }
        | brigen_statements T_imsize angle angle '\n' 
	                    { if (bgp && set_imsize($3, $4)) YYERROR; }
        | brigen_statements T_imsize angle ',' angle '\n'
	                    { if (bgp && set_imsize($3, $5)) YYERROR; }
        | brigen_statements T_imfitsout filename '\n'
	                    { if (bgp && set_imfitsout($3)) YYERROR; }
        | brigen_statements T_addimage filename '\n'
	                    { if (bgp && set_extimage($3)) YYERROR; }
        | brigen_statements srclist_here
	;


srclist_here 
        : T_srclist_here  src_lines  T_end_srclist 
	        { if (bgp && !bgp->is_set.srclist) {
		      bgp->is_set.srclist = 2; 
		  }
		  else if (bgp->is_set.srclist == 2) {
		      *yytext = 0;
		      yyerror("multiple source lists"); YYERROR;
		  }
                }
	;

src_lines
	: /*empty */
        | src_lines '\n' 
	| src_lines src_data '\n'

src_data
	: T_name num num num num num num num num num num
	    { if (bgp && !bgp->is_set.srclist) {
	        struct Source *src = bgp->source;
  	        if (strlen($1) > MAX_SRCID_LEN) { 
	          *yytext = 0;
	          yyerror("sky source ID length esceeds MAX_SRCID_LEN"); 
	          YYERROR; }
	        else strcpy(src[isrc].id, $1);
	        src[isrc].fluxDensity = $2; 
		src[isrc].Q = $3; 
		src[isrc].U = $4; 
		src[isrc].V = $5; 
		src[isrc].specIndex = $6; 
		src[isrc].xOffset = $7; 
		src[isrc].yOffset = $8; 
		src[isrc].majorAxis = $9; 
		src[isrc].minorAxis = $10; 
		src[isrc].positionAngle = $11; 
		if (++isrc >= MAX_SOURCES) { 
		  *yytext = 0;
		  yyerror("Number of celestial sources exceeds MAX_SOURCES"); 
		  /*	    MAX_SOURCES); */ 
		  YYERROR; }
	      }
	    }
	; 
	
visgen_section
	: T_visgen T_section '\n' visgen_statements  T_end T_visgen '\n'
	;

visgen_statements 
	: /*empty */
	| visgen_statements '\n'
	| visgen_statements obsspec_here
	| visgen_statements T_obsspec filename '\n' 
	| visgen_statements T_array filename '\n' 
	| visgen_statements array_here
	| visgen_statements T_visout filename '\n'
        | visgen_statements T_textout filename '\n'
	| visgen_statements T_end T_array '\n'
	| visgen_statements T_end T_obsspec '\n'
	;

obsspec_here 
	: T_obsspec_here  obsspec  T_end T_obsspec '\n' 
	| T_obsspec_here  obsspec  T_end_obsspec 
	;

obsspec
	: 
	| obsspec '\n'
	| obsspec obsspec_glob
	| obsspec obsspec_scan
	;

obsspec_glob
	: T_FOV_center_RA  '='  angle '\n'
	| T_FOV_center_Dec '='  angle '\n' 
	| T_FOV_size_RA    '='  angle '\n'
	| T_FOV_size_Dec   '='  angle '\n'
	| T_Corr_int_time  '='  T_dur   '\n'
	| T_Corr_chan_bw   '='  T_freq  '\n'   
	;

obsspec_scan
	: T_Scan_start     '='  T_abstime '\n'
	| T_Scan_start     '='  num T_GHA '\n'
	| T_Scan_start     '='  T_GHA num '\n'
	| T_Scan_duration  '='  T_dur '\n'
	| T_Central_freq   '='  T_freq '\n'
	| T_Bandwidth      '='  T_freq '\n'
	| T_Endscan '\n'
	;

array_here 
	:  T_array_here  array_lines  T_end T_array '\n' 
	|  T_array_here  array_lines  T_end_array  
	;

array_lines
	: /*empty */
	| array_lines '\n'
	| array_lines array_data '\n'

array_data
	: T_name num num num T_name num num num 
	| num num T_name  

filename
	: T_file { $$ = $1; }
	| T_name { $$ = $1; }
	;

num	
	: T_int   { $$ = (double)$1; } 
	| T_float { $$ = $1; } 

angle
        : T_angle { $$ = $1; }
        | T_xms   { $$ = $1; }

%%


/* 
 * Convert T_angle value to degrees
 * There are two kinds of angle specs: rad, deg, ', "
 * and hh:mm:ss hms, dd:mm:ss dms.
 * The latter is recognized by units: dms or hms.
 */
double angle2deg(char *angle, int *err) {
  char *angle1 = strdup(angle);
  //char blank[] = " \t";
  char *units = NULL; 
  int au = 0;

  *err = 0; /* Assume no errors */
  /* Set au = 2 if the angle is 'hh:mm:ss hms' or 'dd:mm:ss dms' 
  * Otherwise au = 1 */
  if (strstr(angle1, "hms")) { au = 2; units = strdup("hms"); } else
  if (strstr(angle1, "dms")) { au = 2; units = strdup("dms"); } else
    au = 1;
  /* if (au == 2) printf("au = %d, units = '%s'\n", au, units); */


  if (au == 1) { /* rad, deg, ', or " */
    char *tail;
    double x = strtod(angle, &tail);
    units = skiw(tail);
    if      (!strcmp(units, "'"))    x =   x/60.0;
    else if (!strcmp(units, "\""))   x = x/3600.0;
    else if (!strcmp(units, "''"))   x = x/3600.0;
    else if (!strcmp(units, "rad"))  x = x*180.0/pi;
    else if (!strcmp(units, "deg")) /* Do nothing */;
    else { /* Impossible for good lexer: neither rad nor deg nor " etc. */
      *err = 2;
      return 0.0;
    }
    /* printf("units = '%s'\n", units); */
    return x;
  }

  else { /* if (au == 2) { // hh:mm:ss hms, dd:mm:ss dms */
    const char colon[] = ":"; /* delimiter */
    char *angle1 = strdup(angle); /* Once again, for angle1 is spoiled */
    double xx, mm, ss;
    xx = atof(strtok(angle1, colon));   /* hours or degrees */
    mm = atof(strtok(NULL, colon));  /* minutes */
    ss = atof(strtok(NULL, colon));  /* seconds */
    /* angle1 = strdup(angle); */
    /* units = strtok(angle1, blank); */
    /* units = strtok(NULL, blank); /\* Second after the blanks *\/ */

    if (mm > 60.0 || ss > 60.0) { /* Minutes and seconds are [0..60] */
      *err = 1;
      return 0.0;
    }
    if (xx < 0) { mm = - mm; ss = -ss; } /* Adjust if negative */
    if (!strcmp(units, "hms")) {    /* Hours:minutes:seconds */ 
      if (fabs(xx) > 24.0) { /* Hours must be 0 to 24 */
	*err = 1;
	return 0.0;
      }
      return 15.0*(xx + mm/60.0 + ss/3600.);
    }
    else if (!strcmp(units, "dms"))    /* Degrees:minutes:seconds */ 
      return xx + mm/60.0 + ss/3600.;
    else { /* Impossible: either dms or hms */
      *err = 2;
      return 0.0;

    }
  }
}

  
int set_imgrid(int ncol, int nrow) {
  if (!bgp->is_set.ncol) { 
    bgp->ncol = ncol; 
    bgp->is_set.ncol = 2; 
  }
  else if (bgp->is_set.ncol == 2) { 
    *yytext = 0;
    yyerror("multiple assignment to 'imgrid_radec'"); 
    return 1;
  }
  if (!bgp->is_set.nrow) {
    bgp->nrow = nrow; 
    bgp->is_set.nrow = 2; 
  }
  else if (bgp->is_set.nrow == 2) { 
    *yytext = 0;
    yyerror("multiple assignment to 'imgrid_radec'"); 
    return 2;
  }
  return 0;
}


/*
 * Set FOV RA and DEC center in degrees
 */
int set_imcenter(char *ra, char *dec) {
  char *ra1 = strdup(ra), *dec1 = strdup(dec);
  int err;

  if (!bgp->is_set.xcenter) {   
    bgp->xcenter = angle2deg(ra1, &err);
    if (err) {
      yytext = strdup(ra);
      yyerror("Wrong 'imcenter_radec' format");
      return 1;
    }
    bgp->is_set.xcenter = 2; 
  }
  else if (bgp->is_set.xcenter == 2) { 
    *yytext = 0;
    yyerror("multiple assignment to 'imcenter_radec'"); 
    return 1;
  }

  if (!bgp->is_set.ycenter) {
    bgp->ycenter = angle2deg(dec1, &err);
    if (err) {
      yytext = strdup(dec);
      yyerror("Wrong 'imcenter_radec' format");
      return 1;
    }
    bgp->is_set.ycenter = 2; 
  }
  else if (bgp->is_set.ycenter == 2) { 
    *yytext = 0;
    yyerror("multiple assignment to 'imcenter_radec'"); 
    return 2;
  }
  return 0;
}

/*
 * Set FOV RA and DEC eXtent in arcseconds
 */
int set_imsize(char *ra, char *dec) {
  char *ra1 = strdup(ra), *dec1 = strdup(dec);
  int err;

  if (!bgp->is_set.xsize) {
    bgp->xsize = angle2deg(ra1, &err)*3600.0;
    if (err) {
      yytext = strdup(ra);
      yyerror("Wrong 'imsize_radec' format");
      return 1;
    }
    bgp->is_set.xsize = 2; 
  }
  else if (bgp->is_set.xsize == 2) { 
    *yytext = 0;
    yyerror("multiple assignment to 'imsize_radec'"); 
    return 1;
  }

  if (!bgp->is_set.ysize) {
    bgp->ysize = angle2deg(dec1, &err)*3600.0;
    if (err) {
      yytext = strdup(dec);
      yyerror("Wrong 'imsize_radec' format");
      return 1;
    }
    bgp->is_set.ysize = 2; 
  }
  else if (bgp->is_set.ysize == 2) { 
    *yytext = 0;
    yyerror("multiple assignment to 'imsize_radec'"); 
    return 2;
  }
  return 0;
}


/*
 * Set the FITs file name for output 
 */
int set_imfitsout(char *fitsfilename) {
  if (!bgp->is_set.fitsfilename) { 
    bgp->fitsfilename = strdup(fitsfilename); 
    bgp->is_set.fitsfilename = 2; 
  }
  else if (bgp->is_set.fitsfilename == 2) {
    *yytext = 0;
    yyerror("multiple output FITs files"); 
    return 1;
  }
  return 0;
}


/*
 * Set the FITs file name with an external image 
 * (may be several such files) 
 */
int set_extimage(char *imgfilename) {
  if (iimgf < MAXEXTIMAGES)
    bgp->imgfilename[iimgf++] = imgfilename;
  else { /* maximum number of additioanl image files exceeded */
    printf("Number of additional image files exceeds %d\n", iimgf);
    return 1;
  }
  return 0;
}


void yyerror(const char *msg) {
  //printf("yytext = %d\n", (int)(*yytext));
  if (*yytext != 0)
    printf("Error: %s at '%s', line %d\nin file '%s'\n", 
	   msg, yytext, yylineno, curfilename);
  else
    printf("Error: %s at line %d or earlier, \nin file '%s'\n", 
	   msg, yylineno, curfilename);
}












