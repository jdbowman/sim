#define SET_INPUT(text,r,type,pat,ver0,ver1,ver2,er0,er1,er2,variable) {pText[nFieldNames]=text; rq[nFieldNames]=r; \
  dataType[nFieldNames]=type; pResult[nFieldNames]=&variable; lim[nFieldNames][0]=ver0; lim[nFieldNames][1]=ver1; lim[nFieldNames][2]=ver2; \
  pPat[nFieldNames]=pat; pErr[nFieldNames][0]=er0; pErr[nFieldNames][1]=er1; pErr[nFieldNames][2]=er2; nFieldNames++;}

// Data types
#define DT_STR 1	// character string
#define DT_INT 2	// integer
#define DT_HMS 3	// hour:min:sec
#define DT_DEG 4	// deg, arcmin, arcsec
#define DT_DBL 8	// real number, double precision
#define DT_TIME 10	// date and time
#define DT_BOOL 5	// boolean (0/1, Y/N, T/F, etc)

	int nFieldNames=0,fieldSeen[MAX_FIELDS],rq[MAX_FIELDS];
	char *pText[MAX_FIELDS]; char *pPat[MAX_FIELDS]; char *pErr[MAX_FIELDS][3]; double lim[MAX_FIELDS][3];
	void *pResult[MAX_FIELDS]; int dataType[MAX_FIELDS];
	#define MAX_LINE 150
	char filename[100+1],lineIn[MAX_LINE],line[MAX_LINE],*p,fieldname[MAX_FIELDNAME+1],value[L_EXTENDEDLINE];
