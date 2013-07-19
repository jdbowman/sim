/************************************
 ! MODULE: uvfits
 !         utilities for writing/reading UVFITS files in AIPS/MIRIAD format
 !         Author: Randall Wayth. Sep, 2006.
 !         Feb, 2007. Added uvfits reader.
 !         Nov, 2008. Added multi-channel, multi-time, multi IF.
 !         Dec, 2008. Fixed bugs in multi-IF functionality.
 *   RCS: $Id: uvfits.c,v 1.8 2007/04/20 16:00:46 rwayth Exp rwayth $
 ************************************/
 
/* a note about sign conventions etc:
  The de-facto standard for UVFITS files follows the AIPS convention:
  - a baseline is defined as location(ant1)-location(ant2) where the antenna indices are such that ant2 > ant1
  - using this definition of a baseline, a visibility (for a point source) is V(u,v) = I*exp(-2*pi*(ul+vm))
  - this means that if a baseline is defined as (ant2-ant1),
    then the exponent in the exponential must be +2pi, not -2pi
    
  This writer does NOT do any reformatting of data to/from the de-facto standard. It is up to the user of this
  code to encode/decode visibilities as appropriate before the reader/writer is called.
*/

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <float.h>
#include "uvfits.h"

// NAXIS represents: (a subset of these are N_CTYPES)
// 1 (random groups)
// 2 (set to 3 for real/imaginary/weight)
// 3 (polarisations)
// 5 (number of IFs)
// 4 (number of bands within an IF)
#define NAXIS 7
#define N_CTYPES 5

// 6 if there are multiple IFs, otherwise, the number is reduced to 5 in code
#define N_GRP_PARAMS 6
#define MAX_CSIZE 9

// For comparing floating point numbers
#ifndef CMP_FLT
#define CMP_FLT
#define cmpflt(f1,f2) (fabs(f1-f2) < 1e-5)
#endif

/* private function prototypes */
static int writeAntennaData(fitsfile *fptr, uvdata *data);
static int writeFQData(fitsfile *fptr, uvdata *data);
void printHeader(fitsfile *fptr, int header_number);
int readAntennaTable(fitsfile *fptr,uvdata *obj);
int readFQTable(fitsfile *fptr,uvdata *obj);
int writeSourceData(fitsfile *fptr, uvdata *data);

/* private global vars */
static int debug=0;

/******************************
 ! NAME:      uvfitsSetDebugLevel(
 ! PURPOSE:   set the module global debug level.
 ! ARGUMENTS: in_debug - integer global debug setting. 0= no debugging messages.
 ! RETURNS:   void
******************************/
void uvfitsSetDebugLevel(int in_debug) {
    debug = in_debug;
}


/******************************
 ! NAME:      readUVFITS
 ! PURPOSE:   read a UVFITS file and create/populate a uvdata data structure
 ! ARGUMENTS: filename: name of UVFITS file to read
 !            data: pointer to pointer data struct. The struct will be allocated by this function.
 ! RETURNS:   integer. 0 success, nonzero error. Returns CFITSIO error codes if a CFITSIO error
 !            occurred. returns -1 for malloc failure.
******************************/
int readUVFITS(char *filename, uvdata **data) {
    fitsfile *fptr;
    int i,j,k,l,m,status=0,retcode=0,num_hdu=0,hdu_type=0,index=0, grp_small_row_size;
    int bitpix,naxis,pcount,gcount,num_IF,num_pol,num_freq,num_bl=0,time_index=0,bl_index=0;
    int date_pindex=4,uscale_pindex=0,vscale_pindex=1,wscale_pindex=2,fscale_pindex=5;
    int complex_cindex=1, pol_cindex=2, if_cindex=3, freq_cindex=4, ctype_index; // Default axes index, into naxes
    long len_bl, len_if, len_freq, group_num; // line lengths for array indexing
    long naxes[NAXIS];
    char temp[84];
    char *ptype[MAX_CSIZE],*ctype[MAX_CSIZE];
    float crval[MAX_CSIZE],crpix[MAX_CSIZE],cdelt[MAX_CSIZE],*grp_par=NULL, *grp_small_row=NULL;
    float currtime=FLT_MAX,mintime=FLT_MAX, curr_if=FLT_MAX;
    double pscal[MAX_CSIZE],pzero[MAX_CSIZE];
    uvdata *obj;
    
    /* init */
    for (i=0; i<MAX_CSIZE; i++ ) {
        ptype[i]=NULL;
        ctype[i]=NULL;
    }
    
    /* open the file */
    fits_open_file(&fptr,filename, READONLY, &status);
    if(status!=0) {
        fprintf(stderr,"readUVFITS: failed to open file <%s>\n",filename);
        return status;
    }
    
    /* get data size, number of extensions etc */
    fits_get_num_hdus(fptr,&num_hdu,&status);
    if(debug) fprintf(stdout,"readUVFITS: there are %d HDUs\n",num_hdu);
    if (num_hdu > 3) {
        fprintf(stderr,"readUVFITS ERROR: %s appears to be a multi-source file. This is not supported yet.\n",filename);
        //    return 1;
    }
    fits_get_hdu_type(fptr,&hdu_type,&status);
    if (hdu_type != IMAGE_HDU) {
        fprintf(stderr,"readUVFITS: first HDU in file is not an image. type: %d\n",hdu_type);
    }
    
    /* print out all the headers for debugging*/
    if (debug) {
        for(i=1;i<=num_hdu; i++) {
            printHeader(fptr,i);
        }
        fits_movabs_hdu(fptr, 1, NULL, &status);
    }
    
    /* get data size */
    fits_get_img_param(fptr,NAXIS,&bitpix,&naxis,naxes,&status);
    if (debug) {
        fprintf(stdout,"readUVFITS: NAXIS is %d\n",naxis);
        for(i=0;i<naxis;i++) fprintf(stdout,"\tNAXIS%d is %d\n",i+1,(int)naxes[i]);
    }
    if (status !=0 ) {
        fprintf(stderr,"readUVFITS: status %d after get img params\n",status);
    }
    
    /* basic stuff seems to be in order, allocate space for data structure (not visibilites yet) */
    *data = calloc(1,sizeof(uvdata));
    if(*data==NULL) {
        fprintf(stderr,"readUVFITS: no malloc for main struct\n");
        status = -1;
        goto EXIT;
    }
    obj = *data;
    obj->visdata=NULL;
    obj->weightdata=NULL;
    obj->baseline=NULL;
    obj->u=NULL;
    obj->v=NULL;
    obj->w=NULL;
    obj->n_baselines=NULL;
    obj->source = calloc(1,sizeof(source_table)); /* no source table yet... */
    obj->array = NULL; /* array table created in readAntennaTable */
    obj->date = calloc(1,sizeof(double));
    if (obj->source==NULL || obj->date==NULL) {
        fprintf(stderr,"readUVFITS: no malloc for source table\n");
        status=-1; goto EXIT;
    }

    /* find out how many groups there are. There is one group for each source,time,baseline combination.
     we don't know how many baselines there are yet- have to work this out while reading the data.
     The rows must be time ordered or everything will break.  */
    fits_read_key(fptr,TSTRING,"GROUPS",temp,NULL,&status);
    if (temp[0] != 'T') fprintf(stderr,"readUVFITS: GROUPS fits keyword not set\n");
    fits_read_key(fptr,TINT,"PCOUNT",&pcount,NULL,&status);
    fits_read_key(fptr,TINT,"GCOUNT",&gcount,NULL,&status);
    if (debug) fprintf(stdout,"PCOUNT is %d, GCOUNT is %d\n",pcount,gcount);
    fits_read_key(fptr,TSTRING,"OBJECT",obj->source->name,NULL,&status);
    fits_read_key(fptr,TDOUBLE,"OBSRA",&(obj->source->ra),NULL,&status);
    fits_read_key(fptr,TDOUBLE,"OBSDEC",&(obj->source->dec),NULL,&status);
    if (status !=0) {
        fprintf(stderr,"readUVFITS WARNING: status %d while getting basic keywords.\n",status);
        status=0;
    }

    /* get the details of the observations PTYPE/CTYPE keyword from header*/
    for(i=0; i<pcount; i++) {
        sprintf(temp,"PTYPE%d",i+1);
        ptype[i] = malloc(sizeof(char)*84);
        if (ptype[i]==NULL) {
            fprintf(stderr,"readUVFITS: no malloc for ptype char array\n");
            status = -1; goto EXIT;
        }
        ptype[i][0]='\0';
        fits_read_key(fptr,TSTRING,temp,ptype[i],NULL,&status);
        sprintf(temp,"PSCAL%d",i+1);
        fits_read_key(fptr,TDOUBLE,temp,pscal+i,NULL,&status);
        sprintf(temp,"PZERO%d",i+1);
        fits_read_key(fptr,TDOUBLE,temp,pzero+i,NULL,&status);
        
        /* remember the stuff we care about */
        if (strncmp("DATE",ptype[i],4)==0) {
            date_pindex = i;
        }
        if (strncmp("UU",ptype[i],2)==0) {
            uscale_pindex = i;
        }
        if (strncmp("VV",ptype[i],2)==0) {
            vscale_pindex = i;
        }
        if (strncmp("WW",ptype[i],2)==0) {
            wscale_pindex = i;
        }
        if (strncmp("FREQSEL",ptype[i],7)==0) {
            fscale_pindex = i;
        }
        
        if(debug) fprintf(stdout,"PZERO%d: %f, PSCAL%d: %g, PTYPE%d: %s\n",i+1,pzero[i],i+1,pscal[i],i+1,ptype[i]);
    }
    
    // get the CRs, starting with CTYPE2 (the complex vis axis, and following the number of ctypes)
    for(i=0; i<=N_CTYPES; i++) {
        ctype_index = i+2; // values start at CTYPE2
        sprintf(temp,"CTYPE%d",ctype_index);
        ctype[i] = malloc(sizeof(char)*84); /* strings can be 80 chars not including terminator */
        if (ctype[i]==NULL) {
            fprintf(stderr,"readUVFITS: no malloc for ctype char array\n");
            status = -1; goto EXIT;
        }
        ctype[i][0] = '\0';
        fits_read_key(fptr,TSTRING,temp,ctype[i],NULL,&status);
        if (status != 0) {
            if(debug) fprintf(stdout,"no more C keys. status: %d\n",status);
            status =0;
            break;
        }
        sprintf(temp,"CRVAL%d",ctype_index);
        fits_read_key(fptr,TFLOAT,temp,crval+i,NULL,&status);
        sprintf(temp,"CRPIX%d",ctype_index);
        fits_read_key(fptr,TFLOAT,temp,crpix+i,NULL,&status);
        sprintf(temp,"CDELT%d",ctype_index);
        fits_read_key(fptr,TFLOAT,temp,cdelt+i,NULL,&status);
        if(debug) fprintf(stdout,"CRVAL%d: %g,\tCRPIX%d: %g,\tCDELT%d: %g,\tCTYPE%d: %s\n",
                          ctype_index,crval[i],ctype_index,crpix[i],ctype_index,cdelt[i],ctype_index,ctype[i]);
        
        /* check for the stuff we're interested in */
        if(strncmp("COMPLEX", ctype[i], 7)==0) {
            complex_cindex = i+1;
        }
        // TODO: need to check convention for other polarization types
        else if(strncmp("STOKES", ctype[i], 7)==0 || strncmp("LINEAR POLARIZATION", ctype[i], 19)==0 || 
                strncmp("CIRCULAR POLARIZATION", ctype[i], 21)==0) {
            pol_cindex = i+1;
        }
        else if(strncmp("IF", ctype[i], 2)==0) {
            if_cindex = i+1;
        }
        else if (strncmp("FREQ",ctype[i],4)==0) {
            /* centre freq is CRVAL, BW is CDELT */
            freq_cindex = i+1;
            obj->cent_freq = crval[i];
            obj->freq_delta = cdelt[i];
        }
    }
  
    if(debug) fprintf(stdout,"NAXES: Complex(%d)=%ld, Pol(%d)=%ld, IF(%d)=%ld, Freq(%d)=%ld\n",complex_cindex,
                    naxes[complex_cindex], pol_cindex, naxes[pol_cindex], if_cindex, naxes[if_cindex],
                    freq_cindex, naxes[freq_cindex]);

    // Takes care of all axes options (complex is always 3: real, imag, weight)
    obj->n_pol = num_pol = naxes[pol_cindex];
    obj->n_IF = num_IF = naxes[if_cindex];
    obj->n_freq = num_freq = naxes[freq_cindex];
  
    // Calculate line lengths
    len_freq = num_pol;
    len_if = num_freq*len_freq;
    len_bl = num_IF*len_if;
    
    // The total row size, for all axes, and for all baselines, times (IFs included in the groups)
    grp_small_row_size = naxes[complex_cindex]*num_pol*num_freq*num_IF;
    grp_small_row = malloc(sizeof(float)*grp_small_row_size);
    
    /* each row in the group table has this size */
    grp_par = malloc(sizeof(float)*pcount);
    if (grp_par==NULL) {
        fprintf(stderr,"readUVFITS: no malloc for big arrays\n");
        status = -1; goto EXIT;
    }

    /* for reasons I don't understand, it doesn't work to only read a section
       of this file at a time. So read the whole lot... */
    //  fits_read_img(fptr,TFLOAT,1,grp_row_size,NULL,grp_row,NULL,&status);
    //  if (status != 0) {
    //    fprintf(stderr,"readUVFITS: error reading data %d. status: %d\n",i,status);
    //    status = -1; goto EXIT;
    //  }
    
    /* copy data into data structure. We know how many groups there are, but we don't know
     how many different baselines or time units there are. The data must be ordered by increasing
     time, so we keep track of the time until it changes. This then tells us how many baselines
     there are. If there are missing baselines in the first time step, it will be bad.... */
    for (i=0; i< gcount; i++) {
        
        /* read a row from the parameters. contains U,V,W,baseline,time, and sometimes freqsel for IF */
        fits_read_grppar_flt(fptr,i+1,1,pcount, grp_par, &status);
        if (status != 0) {
            fprintf(stderr,"readUVFITS: error reading group %d. status: %d\n",i+1,status);
            status = -1; goto EXIT;
        }
        /* if the time changes, make some space for more data */
        if (currtime != grp_par[4]) {
            
            /* new set of baselines */
            currtime = grp_par[4];
            if (debug) printf("new time: %f\n",currtime);
            /* check that time is increasing only */
            if (currtime < mintime) {
                if (mintime != FLT_MAX) {
                    fprintf(stderr,"readUVFITS ERROR: data is not time ordered. current: %f, min: %f\n",
                            currtime,mintime);
                }
                else {
                    /* first time around, set mintime */
                    mintime = currtime;
                }
            }
            
            /* allocate space for visibility data. Increase
             the size of these arrays to have one for each time sample */
            obj->visdata = realloc(obj->visdata,sizeof(float *)*(time_index+1));
            obj->weightdata = realloc(obj->weightdata,sizeof(float *)*(time_index+1));
            obj->baseline = realloc(obj->baseline,sizeof(float *)*(time_index+1));
            obj->n_baselines = realloc(obj->n_baselines,sizeof(int)*(time_index+1));
            obj->u = realloc(obj->u,sizeof(double *)*(time_index+1));
            obj->v = realloc(obj->v,sizeof(double *)*(time_index+1));
            obj->w = realloc(obj->w,sizeof(double *)*(time_index+1));
            obj->date = realloc(obj->date,sizeof(double)*(time_index+1));

            if (obj->weightdata==NULL || obj->visdata==NULL || obj->baseline==NULL || 
                obj->u==NULL || obj->v==NULL || obj->w==NULL || obj->date==NULL ) {
                fprintf(stderr,"readUVFITS: no malloc in param loop\n");
                status = -1; goto EXIT;
            }
            
            obj->visdata[time_index]=NULL;
            obj->weightdata[time_index]=NULL;
            obj->baseline[time_index]=NULL;
            obj->u[time_index]=NULL;
            obj->v[time_index]=NULL;
            obj->w[time_index]=NULL;
            obj->n_baselines[time_index]=0;
            
            /* there may be less than the full number of baselines in the first and last scan from VLA data.
             so keep track of this so we don't run over the array end below */
            bl_index=0;  /* reset for next scan */
            time_index++;
            
            // Reset the current IF
            if(pcount > 5)
                curr_if = grp_par[5];
        } // end if time changes
    
        // Only do if the IF has changed (as baseline is the same for multiple IFs)
        if(pcount <= 5 || grp_par[5] == curr_if) {
            /* put the UVW data into the data structure now */
            /* first, resize the arrays for the next baseline */
            obj->baseline[time_index-1] = realloc(obj->baseline[time_index-1],sizeof(float)*(bl_index+1));
            obj->u[time_index-1] = realloc(obj->u[time_index-1],sizeof(double)*(bl_index+1));
            obj->v[time_index-1] = realloc(obj->v[time_index-1],sizeof(double)*(bl_index+1));
            obj->w[time_index-1] = realloc(obj->w[time_index-1],sizeof(double)*(bl_index+1));
            if (obj->u[time_index-1]==NULL || obj->v[time_index-1]==NULL || obj->w[time_index-1]==NULL) {
                fprintf(stderr,"readUVFITS: realloc failed for UVW arrays\n");
                status = -1; goto EXIT;
            }
            /* copy the UVW, baseline, time */
            obj->u[time_index-1][bl_index] = grp_par[0]*pscal[uscale_pindex] + pzero[uscale_pindex];
            obj->v[time_index-1][bl_index] = grp_par[1]*pscal[vscale_pindex] + pzero[vscale_pindex];
            obj->w[time_index-1][bl_index] = grp_par[2]*pscal[wscale_pindex] + pzero[wscale_pindex];
            obj->baseline[time_index-1][bl_index] = grp_par[3];
            obj->date[time_index-1] = grp_par[4]+pzero[date_pindex];
            obj->n_baselines[time_index-1] = bl_index+1;

            bl_index++;
        }
        
        /* keep track of the maximum number of baselines that are listed */
        if (num_bl < bl_index) num_bl = bl_index;
        if(debug) {
            printf("PARAMS %d: %g,\t%g,\t%g,\t%g,\t%g",
                   i+1,grp_par[0],grp_par[1],grp_par[2],grp_par[3],grp_par[4]);
            if(gcount > 5)
                printf("\t%g",grp_par[5]);
            printf("\n");
        }
    } // end for each gcount
    
    /* now we know how many time steps there are and how many groups there are.
     The number of baselines has been tracked in the loop. */
    obj->n_vis = time_index;
    if (debug) printf("There are %d time steps. Counted maximum %d baselines\n",time_index,num_bl);
    
    /* now alloc space for the visibilities, now includes multiple IF channels */
    for (i=0; i< time_index; i++) {
        obj->visdata[i] = malloc(sizeof(float)*num_IF*num_freq*num_pol*num_bl*2);
        obj->weightdata[i] = malloc(sizeof(float)*num_IF*num_freq*num_pol*num_bl);
        if (obj->visdata[i]==NULL || obj->weightdata[i]==NULL) {
            fprintf(stderr,"No malloc in time loop\n");
            status = -1; goto EXIT;
        }
    }
    
    /* pack the visibility data structure */
    for (i=0; i< time_index; i++) {
        for (j=0; j<obj->n_baselines[i]; j++) {
            
            group_num = i*obj->n_baselines[i] + j + 1;
            fits_read_img_flt(fptr, group_num, 1, grp_small_row_size, 0, grp_small_row, NULL, &status);
            index = 0;
            
            /* vis data and weights. U,V,W already done above. now includes multiple IF channels */
            for(m=0; m<num_IF; m++) {
                for(k=0; k<num_freq; k++) {
                    for(l=0; l< num_pol; l++) {
                        if(debug) printf("VisR/I/W/A/P [%d]: %f, %f, %f, %f, P\n",
                                         index, 
                                         grp_small_row[index], 
                                         grp_small_row[index+1],
                                         grp_small_row[index+2],
                                         sqrt(grp_small_row[index]*grp_small_row[index]+
                                              grp_small_row[index+1]*grp_small_row[index+1]));
                        obj->visdata[i][(j*len_bl + m*len_if + k*len_freq + l)*2+0] = grp_small_row[index++];
                        obj->visdata[i][(j*len_bl + m*len_if + k*len_freq + l)*2+1] = grp_small_row[index++];
                        obj->weightdata[i][j*len_bl + m*len_if + k*len_freq + l] = grp_small_row[index++];
                    } // n_pol
                } // n_freq
            } // n_IF
        } // n_baselines
    } // n_times

    /* read antenna table. Search thru HDUs to find it. */
    for(i=2;i<=num_hdu; i++) {
        fits_movabs_hdu(fptr, i, NULL, &status);
        temp[0] = '\0';
        fits_read_key(fptr,TSTRING,"EXTNAME",temp,NULL,&status);
        if (strncmp(temp,"AIPS AN",7) == 0) {
            if(debug) {
                fprintf(stdout,"AN table is in HDU %d\n",i);
                fflush(stdout);
            }
            status = readAntennaTable(fptr,obj);
        }
        else if(strncmp(temp,"AIPS FQ",7) == 0) {
            if(debug) {
                fprintf(stdout,"FQ table is in HDU %d\n",i);
                fflush(stdout);
            }
            status = readFQTable(fptr,obj);
        }
        else {
            if (debug) fprintf(stdout,"readUVFITS: Ignoring extension type: %s\n",temp);
            continue;
        }
    }
    
EXIT:
    retcode = status;
    status=0;
    /* all done, close file */
    fits_close_file(fptr, &status);
    if (status!=0) fprintf(stderr,"readUVFITS: problem closing file. status: %d\n",status);
    /* done with raw data, so free everything we have allocated */
    if (grp_par!=NULL) free(grp_par);
    if(grp_small_row!=NULL) free(grp_small_row);
    for (i=0; i< MAX_CSIZE; i++) {
        if (ptype[i]!=NULL) free(ptype[i]);
        if (ctype[i]!=NULL) free(ctype[i]);
    }
    return retcode;
}


/********************************
 ! NAME:      writeUVFITS
 ! PURPOSE:   write a UVFITS file
 ! ARGUMENTS: filename: name of output file (will be created)
 !            data: pointer to uvdata struct that contains the data
 ! RETURNS:   integer: 0 success, nonzero failure. returns CFITSIO errors in the case of a CFITSIO error.
*********************************/
int writeUVFITS(char *filename, uvdata *data) {
    
    fitsfile *fptr;
    FILE *fp;
    int status=0;
    long naxes[NAXIS],i,j,k,l,m,rec_in,rec_out;
    long len_baseline, len_if, len_freq, len_pol; // line lengths for array indexing
    float *array=NULL,temp;
    long nelements=0;
    char *ctypes[N_CTYPES] = {"STOKES","FREQ","IF","RA","DEC"}; // order changed to fit with Miriad
    char tempstr[40];
    
    // The FREQSEL parameter is only used if there are multiple IFs, otherwise it isn't necessary
    // (we only allow for one set of IFs, so the configuration is always "1" but I think we can't
    // leave out the parameter because it triggers Miriad to load the FQ table)
    char *params[N_GRP_PARAMS] = {"UU","VV","WW","BASELINE","DATE","FREQSEL"}; // FQID
    int nGroupParams, nCTypes=N_CTYPES; // axes and parameter lengths to adjust for multiple IFs
    int mon=0,day=0,year=0;
    double jd_day,jd_day_trunc,dtemp=0;
    double *u_ptr, *v_ptr, *w_ptr;
    float  *vis_ptr, *wt_ptr;
    
    // The number of group parameters depending on whether there are multiple IFs or not
    if(data->n_IF > 1) {
        nGroupParams = N_GRP_PARAMS;
    }
    else {
        nGroupParams = N_GRP_PARAMS - 1;
    }
    
    /* set up for cfitsio interface */
    // The IF/freq parameters have been set according to Miriad convention
    // RA/DEC axes need to be there to be read into Miriad
    naxes[0] = 0;
    naxes[1] = 3;
    naxes[2] = data->n_pol;
    naxes[3] = data->n_freq; // swapped IF and freq so major axis second (for Miriad compatibility)
    naxes[4] = data->n_IF;
    naxes[5] = 1; // ra
    naxes[6] = 1; // dec
    
    // number of elements in a group (separate for each baseline and time)
    nelements = naxes[1]*naxes[2]*naxes[3]*naxes[4];
    
    if (data->n_baselines[0] < 1) {
        fprintf(stderr,"There are no baselines.\n");
        return 1;
    }
    
    array = calloc(nelements+nGroupParams,sizeof(float));
    if(array==NULL) {
        fprintf(stderr,"writeUVFITS: no malloc\n");
        exit(1);
    }
    
    /* open the file after checking to see if it already exists and clobbering */
    if ((fp=fopen(filename,"r"))!=NULL) {
        fclose(fp);
        remove(filename);
    }
    fits_create_file(&fptr,filename,&status);
    if (status !=0) {
        fits_report_error(stderr, status);
        return 1;
    }
    
    /* set up for writing fits "random groups". There is one group for each baseline
     of visibilities per time sample including all pols and channels. The group consists of:
     1- the (N_GRP_PARAMS) preamble of: U,V,W (nanoseconds),baseline (see below), time offset (days)
     2- the vis data as real,imag,weight for pol1, real,imag,weight for pol2 etc
     for each pol, then repeated for each freq channel */
    /* using this function creates the keywords PCOUNT=5,GROUPS=T and GCOUNT=n_vis*data->n_baselines*data->n_IF */
    fits_write_grphdr(fptr,TRUE,FLOAT_IMG,NAXIS,naxes,nGroupParams,data->n_vis*data->n_baselines[0],TRUE,&status);
    
    if(status != 0) {
        fits_report_error(stderr, status);
        return 1;
    }
    
    /* Write a keyword; must pass pointers to values */
    fits_update_key(fptr,TSTRING, "OBJECT", data->source->name, NULL, &status);
    dtemp = data->source->ra*15.0;
    fits_update_key(fptr,TDOUBLE, "OBSRA", &dtemp, NULL, &status);
    fits_update_key(fptr,TDOUBLE, "OBSDEC", &(data->source->dec), NULL, &status);
    fits_update_key(fptr,TSTRING, "TELESCOP", "MWA-LFD" , NULL, &status);
    fits_update_key(fptr,TSTRING, "INSTRUME", "MWA-LFD" , NULL, &status);
    temp=2000.0;
    fits_update_key(fptr,TFLOAT, "EPOCH", &temp, NULL, &status);
    /* the following is a fix for an AIPS bug. Technically, shouldn't need it */
    temp=1.0;
    fits_update_key(fptr,TFLOAT, "BSCALE", &temp, NULL, &status);
    
    /* get the JD of the first obs. Then truncate the time within the day.
     The times associated with the visibilities are then increments to
     this truncated start-of-jd time. JDs start at 12:00h on a day. */
    JD_to_Cal(data->date[0],&year,&mon,&day); /* get year, month, day for this JD */
    Cal_to_JD(year,mon,day,&jd_day);          /* gets the JD for the start of this day */
    jd_day_trunc = ((int)jd_day)+0.5;         /* truncate */
    sprintf(tempstr,"%d-%02d-%02dT00:00:00.0",year,mon,day);
    fits_update_key(fptr,TSTRING, "DATE-OBS", tempstr , NULL, &status);
    
    /* UV parameters */
    /* this sets the UU,VV,WW,BASELINE and DATE header params, and possible FREQSEL */
    for (i=0; i<nGroupParams; i++) {
        sprintf(tempstr,"PTYPE%d",(int)i+1);
        fits_update_key(fptr,TSTRING, tempstr, params[i] , NULL, &status);
        sprintf(tempstr,"PSCAL%d",(int)i+1);
        temp=1.0;
        fits_update_key(fptr,TFLOAT, tempstr,&temp , NULL, &status);
        sprintf(tempstr,"PZERO%d",(int)i+1);
        temp=0.0;
        fits_update_key(fptr,TFLOAT, tempstr, &temp , NULL, &status);
    }
    /* except the PZERO for the date, which must match the OBS-DATE in Julian Days
     the following relies on the tempstr from previous loop. */
    sprintf(tempstr, "PZERO%d", 5);
    fits_update_key(fptr,TDOUBLE,tempstr,&jd_day_trunc, NULL, &status);
    
    /* write the coord system keywords CRVAL, CRPIX, CDELT, CTYPE.*/
    temp=1.0;
    fits_update_key(fptr,TFLOAT, "CRVAL2", &temp , NULL, &status);
    fits_update_key(fptr,TFLOAT, "CRPIX2", &temp , NULL, &status);
    fits_update_key(fptr,TFLOAT, "CDELT2", &temp , NULL, &status);
    fits_update_key(fptr,TSTRING,"CTYPE2", "COMPLEX" , NULL, &status);
    
    for (i=0; i<nCTypes; i++) {
        /* this fills in the CRVAL etc arrays with some generic values, then
         specific cases are updated below. */
        sprintf(tempstr,"CRVAL%d",(int)i+3);
        fits_update_key(fptr,TFLOAT, tempstr, &temp , NULL, &status);
        sprintf(tempstr,"CRPIX%d",(int)i+3);
        fits_update_key(fptr,TFLOAT, tempstr, &temp , NULL, &status);
        sprintf(tempstr,"CDELT%d",(int)i+3);
        fits_update_key(fptr,TFLOAT, tempstr, &temp , NULL, &status);
        sprintf(tempstr,"CTYPE%d",(int)i+3);
        fits_update_key(fptr,TSTRING, tempstr, ctypes[i] , NULL, &status);
    }
    
    /* set the input pol type. for circular, CRVAL3 should be -1, CDELT3 -1
     for linear CRVAL3 should be -5, CDELT3 -1
     for Stokes CRVAL3 should be  1, CDELT3  1 */
    temp = data->pol_type;
    fits_update_key(fptr,TFLOAT, "CRVAL3", &temp , NULL, &status);
    temp = (data->pol_type < 0 ? -1.0 : 1.0);
    fits_update_key(fptr,TFLOAT, "CDELT3", &temp , NULL, &status);
    
    // IF axes settings in CR*5 but no scaling or offset is required
    
    /* specific updates for freq channels. FIXME: check this is correct!!
     need to check that delta might want to be negative. */
    /* set the reference pixel for the central freq to be the center of the channels */
    /* FIXME: is there supposed to be a 1/2 channel offset here for an even number of channels? */
    //temp = data->n_freq/2 + 1; /* pixel indexes start at 1... */
    
    // Setting pixel offset to always be the first pixel
    // (leave cent_freq to be the centre of the first band)
    temp = 1;
    
    fits_update_key(fptr,TFLOAT,"CRPIX4",&temp , NULL, &status);
    fits_update_key(fptr,TFLOAT,"CRVAL4",&(data->cent_freq),NULL,&status);
    fits_update_key(fptr,TFLOAT,"CDELT4",&(data->freq_delta) , NULL, &status);
    
    /* specific updates for RA/DEC */
    temp = data->source->ra*15; fits_update_key(fptr,TFLOAT,"CRVAL6",&temp,NULL,&status);
    temp = data->source->dec;   fits_update_key(fptr,TFLOAT,"CRVAL7",&temp,NULL,&status);
    
    /* write the HISTORY AIPS WTSCAL item which seems to be required */
    fits_write_history(fptr,"AIPS WTSCAL =  1.0",&status);
    fits_write_comment(fptr,"written by the UV FITS writer of RBW.",&status);
    
    
    /* write the acutal UV data to the main array.
     this is done in "random groups" with one group for each
     instant of visibilities. Each instant contains all freqs,
     pols and baselines */
    
    // Calculate the line lengths for each level to make array indexing simpler
    // Each line length is basically the number of elements within that loop level
    // Nothing within polarisation line, just the 3 data items (the factor of 3 in output only)
    len_pol = 1;
    len_freq = len_pol*data->n_pol;
    len_if = len_freq*data->n_freq;
    len_baseline = len_if*data->n_IF;
    
    /* pack the array for writing to the file */
    for (i=0; i<data->n_vis; i++) {
        /*    fprintf(stderr,"doing vis set %d of %d. There are %d baselines\n",i+1,data->n_vis,data->n_baselines);*/
        u_ptr = data->u[i];
        v_ptr = data->v[i];
        w_ptr = data->w[i];
        vis_ptr = data->visdata[i];
        wt_ptr = data->weightdata[i];
        
        for(l=0; l<data->n_baselines[i]; l++) {
            array[0] = u_ptr[l];
            array[1] = v_ptr[l];
            array[2] = w_ptr[l];
            /*      EncodeBaseline(data->array->baseline_ant1[l],data->array->baseline_ant2[l],array+3); */
            array[3] = data->baseline[i][l];
            /* last, the time offset */
            array[4] = data->date[i]-jd_day_trunc;
            
            // If multiple IFs, need to add the FREQSEL parameter here (so separate group for each IF)
            // Changed to start at 1, as that's probably what's expected in Miriad
            // hardcode this value for now, only one frequency configuration ie. IFs don't change)
            if(data->n_IF > 1)
                array[5] = 1;
            
            // Loop through each of the NAXIS dimensions of the data (3,n_pol,n_IF,n_freq)
            for(m=0; m<data->n_IF; m++) {
                for (j=0; j<data->n_freq; j++){
                    for (k=0; k<data->n_pol; k++) {
                        // Got rid of m*len_if in the output, as now a group is written for every IF
                        rec_in = l*len_baseline + m*len_if + j*len_freq + k; // no len_pol, index 3 items separately
                        rec_out = m*len_if + j*len_freq + k*len_pol; // array per baseline so no l index
                        
                        // naxes[1] represents the 3 component representation
                        array[rec_out*naxes[1] + nGroupParams] = vis_ptr[rec_in*2]; // real
                        array[rec_out*naxes[1] + nGroupParams + 1] = vis_ptr[rec_in*2 + 1]; // imaginary
                        array[rec_out*naxes[1] + nGroupParams + 2] = wt_ptr[rec_in]; // weight
                    } // n_pol
                } // n_freq
            } // n_IF
            
            // Try writing group just for each baseline (removed all reference to m's)
            fits_write_grppar_flt(fptr,i*data->n_baselines[i]+l+1,1,naxes[1]*len_baseline+nGroupParams, array, &status);
            
        } // n_baselines
    } // n_vis
    
    /* create and write the antenna table */
    status = writeAntennaData(fptr, data);
    if (status !=0) {
        fprintf(stderr,"ERROR: writeAntennaData failed with status %d\n",status);
    }
    
    /* write out the multiple frequency table, if there are multiple IFs */
    if(data->n_IF > 1) {
        status = writeFQData(fptr, data);
        
        if (status !=0) {
            fprintf(stderr,"ERROR: writeFQData failed with status %d\n",status);
        }
    }
    
    /* tidy up */
    fits_close_file(fptr, &status);            /* close the file */
    fits_report_error(stderr, status);  /* print out any error messages */
    free(array);
    return 0;
}

/**************************
 ! NAME:      writeSourceData
 ! PURPOSE:   write the source ("SU") table in the UVFITS file.
 ! ARGUMENTS: fptr: an open fits file.
 !            data: uvdata struct with existing data.
 ! RETURNS:   integer. 0 success.
 **************************/
int writeSourceData(fitsfile *fptr, uvdata *data) {
    fprintf(stderr, "ERROR: writeSourceData not yet implemented.\n");
    return -1;
}


/**************************
 ! NAME:      writeAntennaData
 ! PURPOSE:   write the source ("AN") table in the UVFITS file.
 ! ARGUMENTS: fptr: an open fits file.
 !            data: uvdata struct with existing data.
 ! RETURNS:   integer. 0 success.
 **************************/
int writeAntennaData(fitsfile *fptr, uvdata *data) {
    
    int status=0,i,itemp=0,year,mon,day;
    double temp=0;
    char *ptemp;
    char *col_names[] = {"ANNAME", "STABXYZ", "ORBPARM","NOSTA", "MNTSTA", "STAXOF",
    "POLTYA", "POLAA", "POLCALA", "POLTYB", "POLAB", "POLCALB"};
    char *col_format[] = {"8A","3D","0D","1J","1J","1E",
    "1A","1E","3E","1A","1E","3E"};
    char *col_units[] = {"","METERS","","","","METERS","","DEGREES","","","DEGREES",""};
    char tempstr[40];
    
    /* convert JD date to Calendar date. I don't think
     sub-day accuracy is needed here. */
    JD_to_Cal(data->date[0],&year,&mon,&day);
    
    /* create a new binary table */
    fits_create_tbl(fptr,BINARY_TBL,0, 12 ,col_names, col_format,col_units,"AIPS AN", &status);
    
    /* write table header info */
    fits_update_key(fptr,TDOUBLE,"ARRAYX", &(data->array->xyz_pos[0]), NULL, &status);
    fits_update_key(fptr,TDOUBLE,"ARRAYY", &(data->array->xyz_pos[1]), NULL, &status);
    fits_update_key(fptr,TDOUBLE,"ARRAYZ", &(data->array->xyz_pos[2]), NULL, &status);
    fits_update_key(fptr,TFLOAT,"FREQ", &(data->cent_freq) , NULL, &status);
    
    /* TODO FIXME: finish this */
    temp = 1.32584e2;
    fits_update_key(fptr,TDOUBLE,"GSTIA0",&temp , NULL, &status);
    temp = 3.60985e2;
    fits_update_key(fptr,TDOUBLE,"DEGPDY",&temp, NULL, &status);
    
    sprintf(tempstr,"%d-%02d-%02dT00:00:00.0",year,mon,day);
    fits_update_key(fptr,TSTRING,"RDATE",tempstr , NULL, &status);
    
    temp=0.0;
    fits_update_key(fptr,TDOUBLE,"POLARX", &temp, NULL, &status);
    fits_update_key(fptr,TDOUBLE,"POLARY", &temp, NULL, &status);
    fits_update_key(fptr,TDOUBLE,"UT1UTC", &temp, NULL, &status);
    fits_update_key(fptr,TDOUBLE,"DATUTC", &temp, NULL, &status);
    fits_update_key(fptr,TSTRING,"TIMSYS","UTC" , NULL, &status);
    fits_update_key(fptr,TSTRING,"ARRNAM", data->array->name, NULL, &status);
    itemp=0;
    fits_update_key(fptr,TINT,"NUMORB",&itemp , NULL, &status);
    itemp =3;
    fits_update_key(fptr,TINT,"NOPCAL",&itemp , NULL, &status);
    itemp =-1;
    fits_update_key(fptr,TINT,"FREQID",&itemp, NULL, &status);
    temp=33.0;
    fits_update_key(fptr,TDOUBLE,"IATUTC",&temp , NULL, &status);
    
    /* write data row by row. CFITSIO automatically adjusts the size of the table */
    for(i=0; i<data->array->n_ant; i++) {
        ptemp = data->antennas[i].name; /* need this to get a **char in the next line */
        fits_write_col(fptr,TSTRING,1,i+1,1,1,&ptemp, &status);  /* ANNAME */
        fits_write_col_dbl(fptr,2,i+1,1,3,data->antennas[i].xyz_pos, &status); /* STABXYZ */
        itemp = i+1;
        fits_write_col_int(fptr,4,i+1,1,1,&itemp,&status);  /* NOSTA */
        itemp = data->antennas[i].mount_type;
        fits_write_col_int(fptr,5,i+1,1,1,&itemp,&status);  /* MNTSTA */
        ptemp = data->antennas[i].pol_typeA;
        fits_write_col_str(fptr,7,i+1,1,1,&ptemp,&status);  /* POLTYA */
        fits_write_col_flt(fptr,8,i+1,1,1,&(data->antennas[i].pol_angleA),&status); /* POLAA */
        fits_write_col_flt(fptr,9,i+1,1,1,&(data->antennas[i].pol_calA),&status); /* POL calA */
        ptemp = data->antennas[i].pol_typeB;
        fits_write_col_str(fptr,10,i+1,1,1,&ptemp,&status);
        fits_write_col_flt(fptr,11,i+1,1,1,&(data->antennas[i].pol_angleB),&status); /*POLAB */
        fits_write_col_flt(fptr,12,i+1,1,1,&(data->antennas[i].pol_calB),&status); /*POL calB */
    }
    return status;
}

/**************************
 ! NAME:      writeFQData
 ! PURPOSE:   write the frequency ("FQ") table in the UVFITS file.
 ! ARGUMENTS: fptr: an open fits file.
 !            data: uvdata struct with existing data.
 ! RETURNS:   integer. 0 success.
 **************************/
int writeFQData(fitsfile *fptr, uvdata *data) {
    double f;
    int status=0,i, id;
    int NUM_COLS = 5;
    char *col_names[] = {"FRQSEL", "IF FREQ", "CH WIDTH", "TOTAL BANDWIDTH", "SIDEBAND"};
    char *col_format[] = {"1I","1D","1E","1E","1I"};
    char *col_units[] = {"","","","",""};
    
    /* create a new binary table */
    fits_create_tbl(fptr,BINARY_TBL,0, NUM_COLS ,col_names, col_format,col_units,"AIPS FQ", &status);
    
    /* write table header info: just the number of IFs */
    fits_update_key(fptr,TINT,"NO_IF", &(data->n_IF), NULL, &status);
    
    // write data row by row. CFITSIO automatically adjusts the size of the table
    // Array index is used to create the identifier (index+1 so starts at 1 not 0)
    // Arg order: fptr, col#, first row, first elem, num elems, value, status
    // Value is passed as an address, as uses arrays for multiple row writing
    for(i=0; i<data->n_IF; i++) {
        // Changed from centre of entire IF to just centre of first channel
        f = data->fq[i].freq + data->fq[i].chbw/2 - data->cent_freq;
        //printf("if #%d, f=%f\n", (i+1), f);
        id = 1;//i+1;
        
        fits_write_col_int(fptr,1,i+1,1,1,&(id),&status); /* FRQSEL */
        fits_write_col_dbl(fptr,2,i+1,1,1,&(f), &status); /* IF FREQ */
        fits_write_col_flt(fptr,3,i+1,1,1,&(data->fq[i].chbw), &status); /* CH WIDTH */
        fits_write_col_flt(fptr,4,i+1,1,1,&(data->fq[i].bandwidth), &status); /* TOTAL BANDWIDTH */
        fits_write_col_int(fptr,5,i+1,1,1,&(data->fq[i].sideband), &status); /* SIDEBAND */
    }
    
    return status;
}

/**************************
 ! NAME:      JD_to_Cal
 ! PURPOSE:   convert a Julian date to a Gregorian calendar date
 ! ARGUMENTS: jd: input julian date in double format
 !            year, month, day: output calendar date. day/months start at 1.
 ! RETURNS:   void
 **************************/
/* follows from: http://quasar.as.utexas.edu/BillInfo/JulianDatesG.html */
void JD_to_Cal(double jd, int *year, int *month, int *day) {
    int z,w,x,a,b,c,d,e,f;
    
    z=jd+0.5;
    w = (z-1867216.25)/36524.25;
    x=w/4;
    a=z+1+w-x;
    b=a+1524;
    c= (b-122.1)/365.25;
    d=365.25*c;
    e=(b-d)/30.6001;
    f=30.6001*e;
    *day = b-d-f;
    while (e-1 > 12) e-=12;
    *month = e-1;
    *year = c-4715-((e-1)>2?1:0);
    /*  printf("z:%d, w:%d, x: %d, a: %d, b: %d, c: %d, e: %d, f: %d\n",z,w,x,a,b,c,e,f);*/
}


/**************************
 ! NAME:      Cal_to_JD
 ! PURPOSE:   convert a Gregorian calendar date to a Julian date
 ! ARGUMENTS: jd: output julian date in double format
 !            year, month, day: input calendar date. day/months start at 1.
 ! RETURNS:   void
 **************************/
/* follows from: http://quasar.as.utexas.edu/BillInfo/JulianDatesG.html */
void Cal_to_JD(int year, int month, int day, double *jd) {
    int L;
    
    if (month < 1 || day < 1) {
        fprintf(stderr,"Cal_to_JD WARNING: month (%d) or day (%d) in conversion is less than 1\n",month, day);
    }
    
    /* this one ripped off from IDL */
    L = (month-14)/12;        //   In leap years, -1 for Jan, Feb, else 0
    *jd = day - 32075 + 1461*(year+4800+L)/4 +  367*(month - 2-L*12)/12 - 3*((year+4900+L)/100)/4 - 0.5;
    
}


/**************************
 ! NAME:      printUVData
 ! PURPOSE:   debugging utility to print the uvdata structure.
 ! ARGUMENTS: data: input data struct
 !            fp: open file pointer. can also be "stderr" or "stdout".
 ! RETURNS:   void
 **************************/
void printUVData(uvdata *data,FILE *fp) {
    int i,j,k,l,m;
    int len_bl, len_if, len_freq;
    
    // Calcualte line lengths for indexing into array
    len_freq = data->n_pol;
    len_if = data->n_freq*len_freq;
    len_bl = data->n_IF*len_if;
    
    fprintf(fp,"Source:  \t<%s>\n",data->source->name);
    fprintf(fp,"OBSRA:   \t%f\n",data->source->ra);
    fprintf(fp,"OBSDEC:  \t%f\n",data->source->dec);
    fprintf(fp,"Num Pols:\t%d\n",data->n_pol);
    fprintf(fp,"Num IF:  \t%d\n",data->n_IF);
    fprintf(fp,"Num Freq:\t%d\n",data->n_freq);
    fprintf(fp,"Cen Freq:\t%g\n",data->cent_freq);
    fprintf(fp,"Del Freq:\t%g\n",data->freq_delta);
    fprintf(fp,"Num times:\t%d\n",data->n_vis);
    
    // Print the IFs
    if(data->n_IF > 1) {
        for(i=0; i<data->n_IF; i++) {
            fprintf(fp,"IF%d: freq=%lf, chbw=%f, bandwidth=%f, sideband=%d\n",i+1, data->fq[i].freq,data->fq[i].chbw,
                    data->fq[i].bandwidth, data->fq[i].sideband);
        }
    }
    
    // Print all the visibilities (for all times, baselines, IFs, freqs, polarizations)
    for (i=0; i< data->n_vis; i++) {
        fprintf(fp,"\nScan %d. Julian date:\t%g, Num BL:\t\t%d ",i,data->date[i],data->n_baselines[i]);
        for (j=0; j<data->n_baselines[i]; j++) {
            fprintf(fp,"\nUVW: %g,%g,%g. Baseline: %g\n",data->u[i][j],data->v[i][j],
                    data->w[i][j],data->baseline[i][j]);
            for(m=0; m<data->n_IF; m++) {
                fprintf(fp,"\tIF: %d\n",(m+1));
                for(k=0; k<data->n_freq; k++) { 
                    fprintf(fp,"\t\tFreq: %d\n",k);
                    for(l=0;l<data->n_pol; l++) {
                        fprintf(fp,"\t\tPol %d Val,wt: (%g,%g),%g\n",l,
                                data->visdata[i][(j*len_bl+m*len_if+k*len_freq+l)*2  ],
                                data->visdata[i][(j*len_bl+m*len_if+k*len_freq+l)*2+1],
                                data->weightdata[i][j*len_bl+m*len_if+k*len_freq+l]);
                    } // n_pol
                } // n_freq
            } // n_IF
        } // n_baselines
    } // n_vis
}


/**************************
 ! NAME:      printAntennaData
 ! PURPOSE:   debugging utility to print the uvdata antenna structure.
 ! ARGUMENTS: data: input data struct
 !            fp: open file pointer. can be "stderr" or "stdout".
 ! RETURNS:   void
 **************************/
void printAntennaData(uvdata *data,FILE *fp) {
    int i;
    
    fprintf(fp,"Array X,Y,Z: %.10g,%.10g,%.10g\n",
            data->array->xyz_pos[0],data->array->xyz_pos[1],data->array->xyz_pos[2]);
    for (i=0; i< data->array->n_ant; i++) {
        fprintf(fp,"Index: %d, Name: %s, id: %d, XYZ: (%g,%g,%g), Mount: %d\n",
                i,data->antennas[i].name, data->antennas[i].station_num, data->antennas[i].xyz_pos[0],
                data->antennas[i].xyz_pos[1],data->antennas[i].xyz_pos[2],
                data->antennas[i].mount_type);
    }
}


/**************************
 ! NAME:      printHeader
 ! PURPOSE:   debugging utility to print a header from an open FITS file.
 ! ARGUMENTS: fptr: pointer to open FITS file.
 !            header_number: HDU number in the file. starts at 1.
 ! RETURNS:   void
 **************************/
void printHeader(fitsfile *fptr, int header_number) {
    int i,status=0,nkeys=0;
    char *str=NULL,line[81];
    
    fits_movabs_hdu(fptr, header_number, NULL, &status);
    fits_hdr2str(fptr,TRUE,NULL,0,&str,&nkeys,&status);
    printf("Header %d. There are %d keys\n",header_number,nkeys);
    for(i=0; i<nkeys; i++) {
        memset(line,'\0',81);
        strncpy(line,str+i*80,80);
        printf("%s\n",line);
    }
    if (str!=NULL) free(str);
}


/**************************
 ! NAME:      readAntennaTable
 ! PURPOSE:   extract data we care about from the antenna table in the FITS file
 ! ARGUMENTS: fptr: pointer to open FITS file.
 !            obj: pointer to existing uvdata struct. Assumes no existing antenna/array data in the struct.
 ! RETURNS:   integer: 0 for success. returns CFITSIO error codes or -1 for memory allocation problems.
 **************************/
int readAntennaTable(fitsfile *fptr,uvdata *obj) {
    int status=0,ncols,i,curr_col;
    long nrows;
    char **str_names=NULL;
    
    /* get the properties of the table */
    fits_get_num_rows(fptr, &nrows, &status);
    fits_get_num_cols(fptr, &ncols, &status);
    if (status !=0) {
        fprintf(stderr,"readAntennaTable: status %d reading number of rows/cols\n",status);
        goto EXIT;
    }
    if (debug) printf("AN table has %d rows, %d columns\n",(int)nrows,ncols);
    
    /* make space for 'array' struct and an array of antenna structs */
    if (obj->array==NULL) obj->array = calloc(1,sizeof(array_data));
    if(obj->antennas != NULL) {
        fprintf(stderr,
            "readAntennaTable WARNING: existing antennas array will be overwritten and lost. This is a memory leak.\n");
    }
    obj->antennas = calloc(nrows,sizeof(ant_table));
    str_names = calloc(nrows+1,sizeof(char *));
    if (obj->array==NULL || obj->antennas==NULL || str_names==NULL) {
        fprintf(stderr,"readAntennaTable: no malloc\n");
        status=-1; goto EXIT;
    }
    obj->array->n_ant = nrows;
    
    /* load the array X,Y,Z location */
    fits_read_key(fptr,TDOUBLE,"ARRAYX",&(obj->array->xyz_pos[0]),NULL,&status);
    fits_read_key(fptr,TDOUBLE,"ARRAYY",&(obj->array->xyz_pos[1]),NULL,&status);
    fits_read_key(fptr,TDOUBLE,"ARRAYZ",&(obj->array->xyz_pos[2]),NULL,&status);
    if (status != 0) {
        fprintf(stderr,"readAntennaTable: status %d reading array XYZ\n",status);
        goto EXIT;
    }
    
    /* get ANNAME */
    /* create pointers to the start of the antenna name arrays as required by CFITSIO */
    for(i=0;i<nrows; i++) {
        str_names[i] = obj->antennas[i].name;
    }
    curr_col = 0;
    fits_get_colnum(fptr,CASEINSEN, "ANNAME", &curr_col, &status);
    if (curr_col==0 || status !=0) {
        fprintf(stderr,"readAntennaTable: no ANNAME column in antenna table\n");
        goto EXIT;
    }
    
    fits_read_col_str(fptr, curr_col ,1, 1, nrows, NULL, str_names, NULL, &status);
    if (debug)  for (i=0; i<nrows; i++) printf("Antenna %d name: %s\n",i+1,obj->antennas[i].name);
    if (status !=0) {
        fprintf(stderr,"readAntennaTable: status %d reading column ANNAME\n",status);
    }
    
    /* get STABXYZ */
    curr_col = 0;
    fits_get_colnum(fptr,CASEINSEN, "STABXYZ", &curr_col, &status);
    if (curr_col==0 || status !=0) {
        fprintf(stderr,"readAntennaTable: no STABXYZ column in antenna table\n");
        goto EXIT;
    }
    
    for (i=0; i<nrows; i++) {
        fits_read_col_dbl(fptr, curr_col ,i+1, 1, 3, 0.0, obj->antennas[i].xyz_pos, NULL, &status);
        if (debug) printf("Antenna %d pos: %f,%f,%f\n",i+1,obj->antennas[i].xyz_pos[0],
                          obj->antennas[i].xyz_pos[1],obj->antennas[i].xyz_pos[2]);
        if (status !=0) {
            fprintf(stderr,"readAntennaTable: status %d reading column for STABXYZ\n",status);
        }
    }
    
    /* get NOSTA */
    curr_col = 0;
    fits_get_colnum(fptr,CASEINSEN, "NOSTA", &curr_col, &status);
    if (curr_col==0 || status !=0) {
        fprintf(stderr,"readAntennaTable: no NOSTA column in antenna table\n");
        goto EXIT;
    }
    
    for (i=0; i<nrows; i++) {
        fits_read_col_int(fptr, curr_col ,i+1, 1, 1, 0, &(obj->antennas[i].station_num), NULL, &status);
        if (debug) printf("Antenna %d num: %d\n",i+1,obj->antennas[i].station_num);
        if (status !=0) {
            fprintf(stderr,"readAntennaTable: status %d reading column for NOSTA\n",status);
        }
    }
    
    /* get MNTSTA */
    curr_col = 0;
    fits_get_colnum(fptr,CASEINSEN, "MNTSTA", &curr_col, &status);
    if (curr_col==0 || status !=0) {
        fprintf(stderr,"readAntennaTable: no MNTSTA column in antenna table\n");
        goto EXIT;
    }
    
    for (i=0; i<nrows; i++) {
        fits_read_col_int(fptr, curr_col ,i+1, 1, 1, 0, &(obj->antennas[i].mount_type), NULL, &status);
        if (debug) printf("Antenna %d mount type: %d\n",i+1,obj->antennas[i].mount_type);
        if (status !=0) {
            fprintf(stderr,"readAntennaTable: status %d reading column for MNTSTA\n",status);
        }
    }
    
    /* get POLTYA */
    /* create pointers to the start of the antenna name arrays as required by CFITSIO */
    for(i=0;i<nrows; i++) {
        str_names[i] = obj->antennas[i].pol_typeA;
    }
    curr_col = 0;
    fits_get_colnum(fptr,CASEINSEN, "POLTYA", &curr_col, &status);
    if (curr_col==0 || status !=0) {
        fprintf(stderr,"readAntennaTable: no POLTYA column in antenna table\n");
        goto EXIT;
    }
    
    fits_read_col_str(fptr, curr_col ,1, 1, nrows, NULL, str_names, NULL, &status);
    if (debug)  for (i=0; i<nrows; i++) printf("Antenna %d POLTYA: %s\n",i+1,obj->antennas[i].pol_typeA);
    if (status !=0) {
        fprintf(stderr,"readAntennaTable: status %d reading column POLTYA\n",status);
    }
    
    /* get POLTYB */
    /* create pointers to the start of the antenna name arrays as required by CFITSIO */
    for(i=0;i<nrows; i++) {
        str_names[i] = obj->antennas[i].pol_typeB;
    }
    curr_col = 0;
    fits_get_colnum(fptr,CASEINSEN, "POLTYB", &curr_col, &status);
    if (curr_col==0 || status !=0) {
        fprintf(stderr,"readAntennaTable: no POLTYB column in antenna table\n");
        goto EXIT;
    }
    
    fits_read_col_str(fptr, curr_col ,1, 1, nrows, NULL, str_names, NULL, &status);
    if (debug)  for (i=0; i<nrows; i++) printf("Antenna %d POLTYB: %s\n",i+1,obj->antennas[i].pol_typeB);
    if (status !=0) {
        fprintf(stderr,"readAntennaTable: status %d reading column POLTYB\n",status);
    }
    
    /* get POLAA */
    curr_col = 0;
    fits_get_colnum(fptr,CASEINSEN, "POLAA", &curr_col, &status);
    if (curr_col==0 || status !=0) {
        fprintf(stderr,"readAntennaTable: no POLAA column in antenna table\n");
        goto EXIT;
    }
    
    for (i=0; i<nrows; i++) {
        fits_read_col_flt(fptr, curr_col ,i+1, 1, 1, 0, &(obj->antennas[i].pol_angleA), NULL, &status);
        if (debug) printf("Antenna %d pol angle A: %f\n",i+1,obj->antennas[i].pol_angleA);
        if (status !=0) {
            fprintf(stderr,"readAntennaTable: status %d reading column for POLAA\n",status);
        }
    }
    
    /* get POLAB */
    curr_col = 0;
    fits_get_colnum(fptr,CASEINSEN, "POLAB", &curr_col, &status);
    if (curr_col==0 || status !=0) {
        fprintf(stderr,"readAntennaTable: no POLAB column in antenna table\n");
        goto EXIT;
    }
    
    for (i=0; i<nrows; i++) {
        fits_read_col_flt(fptr, curr_col ,i+1, 1, 1, 0, &(obj->antennas[i].pol_angleB), NULL, &status);
        if (debug) printf("Antenna %d pol angle B: %f\n",i+1,obj->antennas[i].pol_angleB);
        if (status !=0) {
            fprintf(stderr,"readAntennaTable: status %d reading column for POLAB\n",status);
        }
    }
    
EXIT:
    if(str_names!=NULL) free(str_names);
    return status;
}

/**************************
 ! NAME:      readFQTable
 ! PURPOSE:   read the IF FQ data from the table "AIPS GQ"
 ! ARGUMENTS: fptr: pointer to open FITS file.
 !            obj: pointer to existing uvdata struct. Assumes no existing fq array.
 ! RETURNS:   integer: 0 for success. returns CFITSIO error codes or -1 for memory allocation problems.
 **************************/
int readFQTable(fitsfile *fptr,uvdata *obj) {
    int status=0,ncols,i,j,curr_col, num_IF=0;
    long nrows;
    int NUM_COLS = 5;
    char *col_names[] = {"FRQSEL", "IF FREQ", "CH WIDTH", "TOTAL BANDWIDTH", "SIDEBAND"};
    
    /* get the properties of the table */
    fits_get_num_rows(fptr, &nrows, &status);
    fits_get_num_cols(fptr, &ncols, &status);
    if (status !=0) {
        fprintf(stderr,"readFQTable: status %d reading number of rows/cols\n",status);
        goto EXIT;
    }
    if (debug) printf("FQ table has %d rows, %d columns\n",(int)nrows,ncols);
    
    /* make space for fqtable struct and an array of antenna structs */
    if(obj->fq != NULL)
        fprintf(stderr,"readFQTable WARNING: existing fq array will be overwritten and lost. This is a memory leak.\n");
    
    obj->fq = calloc(nrows,sizeof(fq_table));
    if (obj->fq==NULL) {
        fprintf(stderr,"readFQTable: no malloc for fqtable\n");
        status=-1; goto EXIT;
    }
    
    /* load the NO_IF parameter */
    fits_read_key(fptr,TINT,"NO_IF",&num_IF,NULL,&status);
    if (status != 0) {
        fprintf(stderr,"readFQTable: status %d reading NO_IF\n",status);
        goto EXIT;
    }
    
    // Check the number of rows matches the IFs already read, and the IF parameter matches also
    if(nrows != obj->n_IF || num_IF != obj->n_IF) {
        fprintf(stderr, 
                "readFQTable: ERROR: n_ows in IF table or NO_IF value does not match number of IFs in main data\n");
        status=-1; goto EXIT;
    }
    
    // Check for the existence of the necessary columns
    for(i=0; i<NUM_COLS; i++) {
        curr_col = 0;
        fits_get_colnum(fptr,CASEINSEN, col_names[i], &curr_col, &status);
        if (curr_col==0 || status !=0) {
            fprintf(stderr,"readFQTable: no %s column in FQ table\n", col_names[i]);
            goto EXIT;
        }
        
        // Read the data out of each row
        for (j=0; j<nrows; j++) {
            if(i==1) {
                fits_read_col_dbl(fptr, curr_col,j+1, 1, 1, 0.0, &(obj->fq[j].freq), NULL, &status);
                if (debug) printf("FQ%d, %s: %lf\n", j+1, col_names[i], obj->fq[j].freq);
            }
            else if(i==2) {
                fits_read_col_flt(fptr, curr_col,j+1, 1, 1, 0.0, &(obj->fq[j].chbw), NULL, &status);
                if (debug) printf("FQ%d, %s: %f\n", j+1, col_names[i], obj->fq[j].chbw);
            }
            else if(i==3) {
                fits_read_col_flt(fptr, curr_col,j+1, 1, 1, 0.0, &(obj->fq[j].bandwidth), NULL, &status);
                if (debug) printf("FQ%d, %s: %f\n", j+1, col_names[i], obj->fq[j].bandwidth);
            }
            else if(i==4) {
                fits_read_col_int(fptr, curr_col,j+1, 1, 1, 0.0, &(obj->fq[j].sideband), NULL, &status);
                if (debug) printf("FQ%d, %s: %d\n", j+1, col_names[i], obj->fq[j].sideband);
            }
            
            if (status !=0) {
                fprintf(stderr,"readFQTable: status %d reading column for %s\n",status, col_names[i]);
            }
        }
    }
    
    // Convert back to absolute frequencies (in file as relative)
    for(i=0; i<obj->n_IF; i++) {
        obj->fq[i].freq = obj->fq[i].freq - obj->fq[i].bandwidth/2 + obj->cent_freq;
    }
    
EXIT:
    return status;
}


/**************************
 ! NAME:      freeUVFITSdata
 ! PURPOSE:   free space associated with a uvfits data struct
 ! ARGUMENTS: data: pointer to data struct
 ! RETURNS:   void
 **************************/
void freeUVFITSdata(uvdata *data) {
    int i;
    
    if (data->date !=NULL) free(data->date);
    if (data->antennas!=NULL) free(data->antennas);
    if (data->array !=NULL) free(data->array);
    if (data->source != NULL) free(data->source);
    if (data->n_baselines !=NULL) free(data->n_baselines);
    for(i=0; i<data->n_vis; i++) {
        if (data->visdata!=NULL && data->visdata[i] !=NULL) free(data->visdata[i]);
        if (data->weightdata!=NULL && data->weightdata[i] !=NULL) free(data->weightdata[i]);
        if (data->u!=NULL && data->u[i] !=NULL) free(data->u[i]);
        if (data->v!=NULL && data->v[i] !=NULL) free(data->v[i]);
        if (data->w!=NULL && data->w[i] !=NULL) free(data->w[i]);
        if (data->baseline!=NULL && data->baseline[i]!=NULL) free(data->baseline[i]);
    }
    if (data->visdata!=NULL) free(data->visdata);
    if (data->weightdata!=NULL) free(data->weightdata);
    if (data->u != NULL) free(data->u);
    if (data->v != NULL) free(data->v);
    if (data->w != NULL) free(data->w);
    if (data->baseline !=NULL) free(data->baseline);
    //free(data);
}


/* encode the baseline. Use the miriad convention to handle more than 255 antennas (up to 2048).
   this is backwards compatible with the standard UVFITS convention */
void EncodeBaseline(int b1, int b2, float *result) {
    if (b2 > 255) {
        *result = b1*2048 + b2 + 65536;
    } else {
        *result = b1*256 + b2;
    }
}


/* decode a baseline using same (extended) miriad format as above */
void DecodeBaseline(float blcode, int *b1, int *b2) {
    if (blcode > 65536) {
        blcode -= 65536;
        *b2 = (int) blcode % 2048;
        *b1 = (blcode - *b2)/2048;
    }
    else {
        *b2 = (int) blcode % 256;
        *b1 = (blcode - *b2)/256;
    }
}

/* Compares all the data within two uvdata structures. Takes into account the last significant digit
 * may change due to floating-point rounding errors.
 */
int compareUVFITS(uvdata *data1, uvdata *data2) {
    int i, j, k, l, m;
    int n_vis, n_baselines, n_IF, n_freq, n_pol;
    int index, len_bl, len_if, len_freq, len_pol;
    
    // Set basic dimensions
    n_vis = data1->n_vis;
    n_baselines = data1->n_baselines[0];
    n_IF = data1->n_IF;
    n_freq = data1->n_freq;
    n_pol = data1->n_pol;
    
    // Calculate line lengths for array indexing
    len_pol = 1;
    len_freq = len_pol*data1->n_pol;
    len_if = len_freq*data1->n_freq;
    len_bl = len_if*data1->n_IF;
    
    // Check basic dimensions
    if(n_vis != data2->n_vis) {
        printf("Number of time steps different.\n");
        return 1;
    }
    if(n_baselines != data2->n_baselines[0]) {
        printf("Number of baselines different.\n");
        return 1;
    }
    if(n_IF != data2->n_IF) {
        printf("Number of IFs different.\n");
        return 1;
    }
    if(n_freq != data2->n_freq) {
        printf("Number of frequencies different.\n");
        return 1;
    }
    if(n_pol != data2->n_pol) {
        printf("Number of polarizations different.\n");
        return 1;
    }
    
    // Compare frequencies
    if(!cmpflt(data1->cent_freq, data2->cent_freq)  || !cmpflt(data1->freq_delta, data2->freq_delta)) {
        printf("Center frequency or delta differ.\n");
        return 1;
    }
    
    // Check IF table if present
    if(n_IF > 1) {
        for(i=0; i<n_IF; i++) {
            if(!cmpflt(data1->fq[i].freq, data2->fq[i].freq)!=0 || 
               !cmpflt(data1->fq[i].chbw, data2->fq[i].chbw) ||
               !cmpflt(data1->fq[i].bandwidth, data2->fq[i].bandwidth) || 
               data1->fq[i].sideband != data2->fq[i].sideband) {
                printf("Freq/bandiwidth in IF table different.\n");
            }
        }
    }
    
    // Check all visibility/time/baseline data
    for(i=0; i<n_vis; i++) {
        if(!cmpflt(data1->date[i],data2->date[i])) {
            printf("n_vis dates differ\n");
            return 1;
        }
        if(data1->n_baselines[i] != data2->n_baselines[i]) {
            printf("n_baselines differ for n_vis=%d\n", i+1);
            return 1;
        }
        
        for(j=0; j<n_baselines; j++) {
            // Compare baseline,u,v,w
            if(!cmpflt(data1->baseline[i][j], data2->baseline[i][j])) {
                printf("Baselines differ at n_vis=%d, n_baseline=%d\n", i+1, j+1);
                return 1;
            }
            if(!cmpflt(data1->u[i][j], data2->u[i][j]) || !cmpflt(data1->v[i][j], data2->v[i][j]) ||
               !cmpflt(data1->u[i][j], data2->u[i][j])) {
                printf("UVW different at n_vis=%d, n_baseline=%d\n", i+1, j+1);
                return 1;
            }
            
            // Check the visibility and weight data
            for(k=0; k<n_IF; k++) {
                for(l=0; l<n_freq; l++) {
                    for(m=0; m<n_pol; m++) {
                        index = j*len_bl + k*len_if + l*len_freq + m;
                        if(!cmpflt(data1->visdata[i][2*index],data2->visdata[i][2*index]) || 
                           !cmpflt(data1->visdata[i][2*index+1],data2->visdata[i][2*index+1])) {
                            printf("Visdata different at n_vis=%d, n_baseline=%d, n_IF=%d, n_freq=%d, n_pol=%d\n",
                                   i+1, j+1, k+1, l+1, m+1);
                            return 1;
                        }
                        if(!cmpflt(data1->weightdata[i][index],data2->weightdata[i][index])) {
                            printf("Weightdata different at n_vis=%d, n_baseline=%d, n_IF=%d, n_freq=%d, n_pol=%d\n",
                                   i+1, j+1, k+1, l+1, m+1);
                            return 1;
                        }
                    } // n_pol
                } // n_freq
            } // n_IF
            
        } // n_baselines
    } // n_vis
    
    // TODO: also check source and antenna table
    
    
    // If get to here, everything equal
    return 0;
}

/* convert Geodetic lat/lon/height to XYZ coords */
/* uses constants from WGS84 system, defined in header */
void Geodetic2XYZ(double lat_rad, double lon_rad, double height_meters, double *X, double *Y, double *Z) {
    double s_lat,s_lon,c_lat,c_lon;
    double chi;

    s_lat = sin(lat_rad); c_lat = cos(lat_rad);
    s_lon = sin(lon_rad); c_lon = cos(lon_rad);
    chi = sqrt(1.0 - E_SQUARED*s_lat*s_lat);

    *X = (EARTH_RAD_WGS84/chi + height_meters)*c_lat*c_lon;
    *Y = (EARTH_RAD_WGS84/chi + height_meters)*c_lat*s_lon;
    *Z = (EARTH_RAD_WGS84*(1.0-E_SQUARED)/chi + height_meters)*s_lat;
}

/* convert local topocentric East, North, Up units to Earth-centred earth-fixed cartesian
   coords relative to a specified center. Latitude is geodetic, not geocentric */
/* lat,lon in radian. E,N,U in same units as X,Y,Z */
void ENH2XYZ_absolute(double E,double N, double H, double lat_rad, double lon_rad, double *X, double *Y, double *Z) {
    double c_lat,c_lon,s_lat,s_lon;
 
    s_lat = sin(lat_rad); c_lat = cos(lat_rad);
    s_lon = sin(lon_rad); c_lon = cos(lon_rad);

    *X = -s_lon*E - s_lat*c_lon*N + c_lat*c_lon*H;
    *Y =  c_lon*E - s_lat*s_lon*N + c_lat*s_lon*H;
    *Z =            c_lat*N       + s_lat*H;
}

/*********************************
  convert coords in local topocentric East, North, Height units to
  'local' XYZ units. Local means Z point north, X point to the equator through the local
  meridian and Y is East. This is like the absolute system except that zero lon is now
  the local meridian rather than prime meridian.
  Latitude is geodetic.
  This is what you want for constructing the local antenna positions in a UVFITS antenna table.
**********************************/
void ENH2XYZ_local(double E,double N, double H, double lat, double *X, double *Y, double *Z) {
  double sl,cl;

  sl = sin(lat);
  cl = cos(lat);
  *X = -N*sl + H*cl;
  *Y = E;
  *Z = N*cl + H*sl;
}
