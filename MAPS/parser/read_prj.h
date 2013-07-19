/***********************************************************
 * read_prj.h                                              *
 * Used by the MAPS lexer read_prj.l and parser read_prj.y *
 * Created 02 Feb 2011 by Leonid Benkevitch                *
 ***********************************************************/

#ifndef INPTYPE
#define INPTYPE
/*
 * enum inptype is used by the lexer read_prj.l to determine
 * what kind of file it has finished scanning to return to the parser
 * a correct end-token: T_end_srclist, T_end_obsspec, or T_end_array.    
 */
enum inptype {FNone = 0, FProject = 1, FSrclist, FObsspec, FArray};
#endif

void yyerror(const char *msg);

char *skiw(char* str);  /* Skip leading whitespace characters */ 
void comment_c(void);   /* C-style comment */ 
int newfile(enum inptype ftype, char *fname);
int popfile(enum inptype *ftype);
int yyparse (void);
double angle2deg(char *angle, int *err); /* Convert T_angle value to degrees */
int set_imgrid(int ncol, int nrow);
int set_imcenter(char *ra, char *dec);
int set_imsize(char *ra, char *dec);
int set_imfitsout(char *fitsfilename);
int set_extimage(char *imgfilename); 
double xms2deg(char *xms, char *units, int *err);

