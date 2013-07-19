#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <fitsio.h>
#include <wcslib.h>

#define TRUE 1
#define FALSE 0

#define MAX_DIM 3
#define MAX_POL 4
enum pol_indices {POL_I=0,POL_Q,POL_U,POL_V};
enum celestial_coords {EQUATORIAL, GALACTIC, ECLIPTIC, HELIOECLIPTIC, SUPERGALACTIC };

/* private function prototypes */
static void usage(const char *progname);
static int read_input_image(char *filename, float **data, long imsize[3], struct wcsprm **wcs, int *nwcs);
int make_output_wcs(struct wcsprm **wcs, long *naxes, double lst, double lat, double pix_ang_size,float flipx);
static int add_text_source_list(long naxes[],double pix_ang_size,float *data,char *sourcelistfilename);
int reproject(struct wcsprm *wcs_in, struct wcsprm *wcs_out, long in_naxes[3], long out_naxes[3],
                float *in_data, float *outdata[MAX_POL], int npol);
void altaz2hadec(double alt, double az, double lat, double *ha, double *dec);

/* private global variables */
static int debug=0;
static const char *pol_names[] = {"I","Q","U","V"};
static FILE *fpd;   // debugging output to this file handle.
  /* for multiple polarisations */
static int npol=1;
static float percent = 0.0,Q_val=1.0,U_val=1.0;

/* magic constants for use with sphs2x() for transforming from equatorial to Galactic coords or vice-versa */
const double equ2gal[] = {192.85948125, 62.87174874, 122.93191814,0.455983795747, 0.889988077457};
const double gal2equ[] = {122.93191814, -62.87174874, 192.85948125,0.455983795747, -0.889988077457};

int main(int argc, char *argv[]) {
  char *infilename=NULL,*outfilename=NULL,*sourcelistfilename=NULL;
  char optstring[] = "i:o:w:t:l:s:a:xdP:QU",*tempname=NULL,*header=NULL;
  int  imsize=0,result;
  long out_naxes[MAX_DIM],in_naxes[MAX_DIM]={0,0,0},fpixel[MAX_DIM]={1,1,1};
  double lst=99,lat=99,pix_ang_size=0,flipx=1.0;
  struct wcsprm *wcs_in=NULL,*wcs_out=NULL;
  float  *in_image=NULL;
  int nwcs_in=0,nwcs_out=0, nkeys=0;

  fitsfile *outfp[MAX_POL]={NULL,NULL,NULL,NULL};
  float *outdata[MAX_POL]={NULL,NULL,NULL,NULL}; 
  int poln = 0;
  size_t namelen = 0;
  /* end pol defs */

  fpd = stdout;

  /* parse command line args */
  if (argc==1) usage(argv[0]);
  while ((result = getopt(argc, argv, optstring)) != -1){
    switch(result) {
      case 'i': infilename = optarg;
        break;
      case 'o': outfilename = optarg;
        break;
      case 's': sourcelistfilename = optarg;
        break;
      case 'd': debug=1;
        break;
      case 'w': imsize = atoi(optarg);
        break;
      case 't': lst = atof(optarg);
        break;
      case 'l': lat = atof(optarg);
        break;
      case 'a': pix_ang_size = atof(optarg);
        break;
      case 'x': flipx=-1.0;
        break;
      case 'P': 
        percent = atof(optarg);
        npol = 4;
        break;
      case 'Q':
        Q_val = 0.0;
        break;
      case 'U':
        U_val = 0.0;
        break;
      default:    
        break;

    }
  }
  /* sanity checks */
  if (imsize <= 0 || imsize > 65535) {
    fprintf(stderr,"invalid image size: %d. exiting.\n",imsize);
    exit(1);
  }
  if (infilename==NULL) {
    fprintf(stderr,"ERROR: input image file name not specified\n");
    exit(1);
  }
  if (outfilename==NULL) {
    fprintf(stderr,"ERROR: output image file name not specified\n");
    exit(1);
  }
  if (lst < 0 || lst > 24.0) {
    fprintf(stderr,"ERROR: LST is invalid or not set\n");
    exit(1);
  }
  if (lat < -90.0 || lat > 90.0) {
    fprintf(stderr,"ERROR: lat is invalid or unset\n");
    exit(1);
  }
  if (pix_ang_size <= 0.0 || pix_ang_size > 60000.0) {
    fprintf(stderr,"ERROR: pixel angular size is invalid or unset: %g\n",pix_ang_size);
    exit(1);
  }
  if (pix_ang_size/3600.0*imsize > 180) {
    fprintf(stderr,"WARNING: angular width of image (%g) is > 180 degrees.\n",pix_ang_size/3600.0*imsize);
  }
  if (npol > 1) {
    if (percent < 0.0 || percent > 100.00) {
      fprintf(stderr,"ERROR: %f percent polarisation non-physical\n",percent);
      exit(1);
    }
    else if (percent < 1.0) {
      fprintf(stderr,"WARNING: %f percent polarisation unlikely?\n",percent);
    }
  }

  if (debug) fprintf(fpd,"output image size: %dx%d\n",imsize,imsize);
  if (debug) fprintf(fpd,"source list file name: %s\n",sourcelistfilename);
  if (debug) fprintf(fpd,"output image file name: %s\n",outfilename);
  if (debug) fprintf(fpd,"input image file name: %s\n",infilename);
  if (debug) fprintf(fpd,"LST: %f hours\n",lst);
  if (debug) fprintf(fpd,"Lat: %f degrees\n",lat);
  if (debug) fprintf(fpd,"pixel angular size: %g (arcsec)\n",pix_ang_size);
  if (debug && flipx < 0) fprintf(fpd,"Flipping direction of RA axis in FITS header (RA increases to +X)\n");
  if (debug && npol > 1) fprintf(fpd,"Generating I,Q,U,V, PA angle points north, %f percent polarisation \n", percent);

  /* convert stuff to radians */
  pix_ang_size *= M_PI/(180.0*3600.0);
  lat *= M_PI/(180.0);
  lst *= M_PI*(15.0/180.0);

  out_naxes[0] = out_naxes[1] = imsize;
  
  /* see if we open the output files now before we do any work */
  /* clobber existing file by adding '!' to front of file name which
     is recognised by cfitsio */
 
  namelen = strlen(outfilename);
  assert(namelen);

  tempname = malloc(sizeof(char)*(namelen+9));
  for(poln=0; poln < npol; poln++) {

      // add a polarisation suffix if we're writing polarisations.
      if (npol > 1) {
          sprintf(tempname,"!%s_%s.fits",outfilename,pol_names[poln]);
      }
      else {
          sprintf(tempname,"!%s.fits",outfilename);
      }
      fits_create_file(&(outfp[poln]),tempname,&result);
      if (result !=0) {
        fprintf(stderr,"could not open output file %s\n",tempname);
        exit(1);
      }
      if (debug) fprintf(fpd,"Opened %s for polarisation %d\n",tempname,poln);
  }
  free(tempname);
  
  /* get the input sky image */
  result = read_input_image(infilename,&in_image,in_naxes, &wcs_in, &nwcs_in);
  if (result) {
    fprintf(stderr,"Failed to read input file <%s>\n",infilename);
    exit(1);
  }
  if (wcs_in->cdelt[0] ==0 || wcs_in->cdelt[1]==0) {
    fprintf(stderr,"ERROR: unknown WCS in input image.\n");
    exit(1);
  }

  for (poln=0; poln < npol; poln++) {
    outdata[poln] = calloc(imsize*imsize,sizeof(float));
    assert(outdata[poln] != NULL);
  }

  /*construct WCS for output image. */
  result = make_output_wcs(&wcs_out,out_naxes,lst,lat,pix_ang_size,flipx);
  nwcs_out =1;

  result = reproject(wcs_in,wcs_out, in_naxes,out_naxes,in_image, outdata,npol);

  /* add point sources -- only works for single pol output at the moment */
  if (npol == 1) {
    if (sourcelistfilename != NULL) result = add_text_source_list(out_naxes,pix_ang_size,outdata[0],sourcelistfilename);
  }

  /* create WCS for FITS header in output files */
  result = wcshdo(WCSHDO_safe,wcs_out,&nkeys,&header);

  /* create an image in the FITS output */
  for (poln=0; poln < npol; poln++) {
    int i;

    if (debug) fprintf(fpd,"Writing POL:%d\n",poln);
    
    fits_create_img(outfp[poln], FLOAT_IMG, 2, out_naxes, &result);

    if (result !=0) {
      fprintf(stderr,"fits_create_img failed with error code: %d\n",result);
      exit(1);
    }

    /* first write empty image */
    //fits_write_null_img(outfp[poln], 1, imsize*imsize, &result);

    /* write out the image */
    fits_write_pix(outfp[poln], TFLOAT, fpixel,imsize*imsize,outdata[poln],&result);
    if (result !=0) {
      fprintf(stderr,"error writing data to output file. code: %d\n",result);
    }

    /* update the header with the WCS info */
    for (i=0; i<nkeys; i++) {
        char keyname[10];
        strncpy(keyname,header+80*i,8);
        keyname[8]='\0';
        fits_update_card(outfp[poln],keyname,header+80*i, &result);
    }
    fits_close_file(outfp[poln],&result);
  }
  if (wcs_in !=NULL) wcsvfree(&nwcs_in, &wcs_in);
  if (wcs_out !=NULL) wcsvfree(&nwcs_out, &wcs_out);
  if (header !=NULL) free(header);
  return 0;
}

/*********************
  construct WCS for output image. This follows from the twcs.c test program distributed with wcslib
**********************/
int make_output_wcs(struct wcsprm **wcs, long *naxes, double lst, double lat, double pix_ang_size,float flipx) {
  int status=0,naxis=2;
  struct wcsprm *wcs_out = NULL;

  wcs_out = malloc(sizeof(struct wcsprm));
  assert(wcs_out != NULL);
  *wcs = wcs_out;
  wcs_out->flag = -1;
  if ((status = wcsini(1,naxis,wcs_out))) {
    printf("make_output_wcs: wcsini ERROR%3d\n", status);
  }

  // make CRPIX
  wcs_out->crpix[0] = naxes[0]/2;
  wcs_out->crpix[1] = naxes[1]/2;
  // make CDELT
  wcs_out->cdelt[0] = -flipx*pix_ang_size*180.0/(M_PI);
  wcs_out->cdelt[1] = pix_ang_size*180.0/(M_PI);
  // make CTYPE
  strcpy(wcs_out->ctype[0],"RA---SIN");
  strcpy(wcs_out->ctype[1],"DEC--SIN");
  // make CRVAL
  wcs_out->crval[0] = lst*180.0/(M_PI);
  wcs_out->crval[1] = lat*180.0/(M_PI);

  // do the initialisation to set up all the required internal wcs stuff
  if ((status = wcsset(wcs_out))) {
    printf("make_output_wcs: wcsset ERROR%3d\n", status);
  }

  return status;
}

/************************
*************************/
int decode_coord_type(struct wcsprm *wcs) {
    if ( wcs->ctype[0][0]=='R') return EQUATORIAL; // looking for RA---
    if ( wcs->ctype[0][0]=='G') return GALACTIC  ; // looking for GLON
    if ( wcs->ctype[0][0]=='E') return ECLIPTIC  ; // looking for ELON
    if ( wcs->ctype[0][0]=='H') return HELIOECLIPTIC; // looking for HLON
    if ( wcs->ctype[0][0]=='S') return SUPERGALACTIC; // looking for SLON

    return -1;
}


/************************
*************************/
int translate_coords(struct wcsprm *wcs_in, struct wcsprm *wcs_out, double celestial_coords[MAX_DIM]) {
    int type_in, type_out;

    type_in = decode_coord_type(wcs_in);
    type_out = decode_coord_type(wcs_out);

    if (type_in ==-1 || type_out==-1) {
        fprintf(stderr,"ERROR translate_coords: unknown celestial coord type %d, %d\n",type_in,type_out);
        return -1;
    }

    if (type_in == type_out) return 0;  // no translation necessary.

    if (type_out == ECLIPTIC || type_out==SUPERGALACTIC || type_out==HELIOECLIPTIC) {
        fprintf(stderr,"ERROR translate_coords: type %d output coords not implemented yet... sorry.\n",type_out);
        return 1;
    }

    /* implement the coord transform */
    if (type_out==EQUATORIAL) {
        switch (type_in) {
            case GALACTIC: {
                double ra,dec,glon,glat;

                ra = celestial_coords[0];
                dec= celestial_coords[1];
                sphs2x(equ2gal, 1, 1, 1, 1, &ra, &dec, &glon, &glat);
                celestial_coords[0] = glon;
                celestial_coords[1] = glat;
                break;
            }
            default:
                fprintf(stderr,"translate_coords: transform %d to %d not implemented yet. How about doing it?\n",
                            type_in, type_out);
                return 1;
                break;
        }
    }
    if (type_out == GALACTIC) {
        switch (type_in) {
            case EQUATORIAL: {
                double ra,dec,glon,glat;

                glon = celestial_coords[0];
                glat = celestial_coords[1];
                sphs2x(gal2equ, 1, 1, 1, 1, &glon, &glat, &ra, &dec);
                celestial_coords[0] = ra;
                celestial_coords[1] = dec;
                break;
            }
            default:
                fprintf(stderr,"translate_coords: transform %d to %d not implemented yet. How about doing it?\n",
                            type_in, type_out);
                return 1;
                break;
        }
    }


    return 0;
}


/************************
*************************/
int reproject(struct wcsprm *wcs_in, struct wcsprm *wcs_out, long in_naxes[3], long out_naxes[3],
                float *in_data, float *outdata[MAX_POL], int npol) {
  int i,j,pol,status=0,stat[1];
  int x1,x2,y1,y2;
  double f11,f12,f21,f22; /* for interpolation function */
  double pixcrd[MAX_DIM],imgcrd[MAX_DIM],phi,theta,world[MAX_DIM],xpos,ypos,ra,dec;

  /* work through the image, calculating the angular location of
     each pixel and sampling the input image */
  for(j=0; j < out_naxes[1]; j++) {
    for(i=0; i < out_naxes[0]; i++) {

        /* calculate the sky position of this pixel in the output image */
        pixcrd[0] = i;
        pixcrd[1] = j;
        pixcrd[2] = 0;
        status = wcsp2s(wcs_out,1,1,pixcrd,imgcrd,&phi,&theta,world,stat);
        if(status != 0 && status != 8) {
            fprintf(stderr,"reproject: status %d from wcsp2s. Pixel %d, %d\n",status,i,j);
            return status;
        }
        if(status == 8) {
            if (debug) fprintf(fpd,"reproject: pixel %d, %d in output does not map to sky\n",i,j);
            for(pol=0; pol < npol; pol++) outdata[pol][i+j*out_naxes[0]] = 0.0;
            continue;
        }
        ra = world[0];
        dec= world[1];
        if(debug) {
            fprintf(fpd,"i,j: %d,%d. RA,DEC: %g,%g\n",i,j,ra,dec);
        }

        /* translate between different sky coord systems if necessary. This uses the current values in world and
           puts the results there too. If no translation is necessary,
           then the results in world are unchanged. */
        status = translate_coords(wcs_in, wcs_out, world);
        if (status) return status;

        /* calculate the pixel position for this sky location in the input image */
        status = wcss2p(wcs_in,1,1,world,&phi,&theta,imgcrd,pixcrd,stat);
        xpos = pixcrd[0];
        ypos = pixcrd[1];

        if (xpos < 0 || xpos > in_naxes[0] || ypos < 0 || ypos > in_naxes[1]) {
          fprintf(fpd,"WARNING: somehow calculated input image pixel that is outside image. i,j: %d,%d ra: %g, dec: %g, x: %.1f, y:%.1f\n",i,j,ra,dec,xpos,ypos);
        }
        else {
          /* sample the sky */
          /* use bilinear interpolation */
          x1 = floor(xpos); x2 = x1+1;
          y1 = floor(ypos); y2 = y1+1;
          /* check to see if we need to wrap around */
          f11 = in_data[x1+in_naxes[0]*y1];
          f12 = in_data[x1+in_naxes[0]*(y2%in_naxes[1])];
          f21 = in_data[(x2%in_naxes[0])+in_naxes[0]*y1];
          f22 = in_data[(x2%in_naxes[0])+in_naxes[0]*(y2%in_naxes[1])];

          outdata[0][i+j*out_naxes[0]] = f11*(x2-xpos)*(y2-ypos) + f21*(xpos-x1)*(y2-ypos) + f12*(x2-xpos)*(ypos-y1) + f22*(xpos-x1)*(ypos-y1);
          /* data[i+j*naxes[0]]= in_image[(int)xpos+in_naxes[0]*(int)ypos]; */ /* this is direct sampling */
          if (debug) printf("i,j: %d,%d, x: %g, y: %g, x1,y1: %d,%d, x2,y2: %d,%d ,RA: %g, DEC: %g, val: %g\n",
              i,j,xpos,ypos,x1,y1,x2,y2,ra,dec,outdata[0][i+j*out_naxes[0]]);
          
        
          if (npol > 1) {
            /* make some fake polarisation data. Stokes I already done above */
            outdata[POL_Q][i+j*out_naxes[0]] = Q_val*outdata[POL_I][i+j*out_naxes[0]]*(0.01*percent);
            outdata[POL_U][i+j*out_naxes[0]] = U_val*outdata[POL_Q][i+j*out_naxes[0]];
            outdata[POL_V][i+j*out_naxes[0]] = 0.0;
          }
        }
    }
  }
  return status;
}

/************************
*************************/
static int add_text_source_list(long naxes[], double pix_ang_size,float *data,char *filename) {
  FILE *fp;
  double norm=0,flux=0,az,el,x,y;
  char name[20];

  /* open the file */
  if ((fp=fopen(filename,"r")) ==NULL) {
    fprintf(stderr,"add_text_source_list: could not open file %s\n",filename);
    return 1;
  }

  /* calculate the normalisation factor */
  /* need Jansky unit point sources to be spread over pixels
     with units Jansky/steradian. */
  /* represent each point as a Gaussian spread over a few
     pixels, so also need to divide by 2pi */
  norm = 1./(pix_ang_size*pix_ang_size*2.0*M_PI);

  /* process each line in the file, adding sources */
  while (fscanf(fp,"%lf %lf %lf %s",&flux,&az,&el,name) != EOF) {
    if(debug) fprintf(fpd,"adding source %s at pos (%f,%f) with flux density %f Jy\n",name,az,el,flux);
    /* find the pixel location in the image */
    /* just undo direction cosines, so x = sin(az)cos(el), y = cos(az)cos(el) */
    x = cos(az)*cos(el);
    y = sin(az)*cos(el);

    /* calc the total in each pixel */
  }
  return 0;
}

/****************************
convert an alt/az coord to HA/DEC given a latitude.
   All quantities are RADIANS
*****************************/
void altaz2hadec(double alt, double az, double lat, double *ha, double *dec) {
  double sindec=0;
  /* find local hour angle */
  *ha = atan2(-sin(az)*cos(alt),-cos(az)*sin(lat)*cos(alt)+sin(alt)*cos(lat));
  if (*ha < -M_PI) *ha += 2.0*M_PI;
  if (*ha > M_PI) *ha -= 2.0*M_PI;

  /* Find declination (positive if north of Celestial Equator, negative if south) */
  sindec = sin(lat)*sin(alt)+cos(lat)*cos(alt)*cos(az);
  *dec = asin(sindec);
}

/*****************************
read the FITS input sky image into an array and collect extra info about the 
   coord system if it is present 
*******************************/
static int read_input_image(char *filename,float **data, long imsize[3], struct wcsprm **wcs_in, int *nwcs) {
  fitsfile *infp=NULL;
  long fpixel[MAX_DIM] = {1,1,1};
  int i,result=0,nkeys=0,nreject=0,fix_status[NWCSFIX],anynul=0;
  char *hdr=NULL;

  fits_open_image(&infp,filename,READONLY,&result);
  if (result !=0) {
    fprintf(stderr,"could not open input file %s\n",filename);
    exit(1);
  }
  fits_get_img_size(infp, 3, imsize, &result);
  if (result !=0) {
    fprintf(stderr,"read_input_image: fits_get_img_size result code: %d\n",result);
    return result;
  }
  if (debug) fprintf(fpd,"input image is %ldx%ldx%ld pixels.\n",imsize[0],imsize[1],imsize[2]);
  if (imsize[2] > 1) {
    fprintf(stderr,"read_input_image ERROR: cannot handle input images with more than 2 dimensions.\n");
    return 1;
  }
  *data = malloc(sizeof(float)*imsize[0]*imsize[1]);
  assert(*data !=NULL);

  /* load the data */
  fits_read_pix(infp, TFLOAT,fpixel, imsize[0]*imsize[1],NULL, *data, &anynul, &result);
  if (result !=0) {
    fprintf(stderr,"read_input_image: fits_read_pix result code: %d\n",result);
    return result;
  }

  /* extract the WCS */
  // get the FITS header as one big string
  fits_hdr2str(infp,1, NULL, 0,&hdr, &nkeys, &result);
  if(result) {
    fprintf(stderr,"Status %d getting header string\n",result);
    return result;
  }
  //extract wcs from header string
  result = wcspih(hdr,nkeys,WCSHDR_all,2,&nreject,nwcs,wcs_in);
  if(result) {
    fprintf(stderr,"Status %d getting WCS\n",result);
    exit(1);
  }
  if (debug) fprintf(fpd,"Got %d WCSes\n",*nwcs);
  // fix the WCS if necessary. print diagnostics
  if( (result = wcsfix(7,NULL,*wcs_in,fix_status)) ) {
    for (i = 0; i < NWCSFIX; i++) {
      if (fix_status[i] > 0) {
        fprintf(stderr, "wcsfix ERROR %d: %s.\n", result,wcsfix_errmsg[fix_status[i]]);
      }
    }
  }
  // set internal WCS data
  if ((result = wcsset(*wcs_in))) {
    fprintf(stderr, "wcsset ERROR %d: %s.\n", result, wcs_errmsg[result]);
    return 1;
  }

  if (debug) {
    for (i=0; i<(*wcs_in)->naxis; i++) fprintf(fpd,"CRVAL%d: %g, CRPIX%d: %g, CDELT%d: %g (degrees)\n",
                                i,(*wcs_in)->crval[i],i,(*wcs_in)->crpix[i],i,(*wcs_in)->cdelt[i]);
  }

  if(debug) wcsprt(*wcs_in);

  fits_close_file(infp,&result);
  if (hdr!=NULL) free(hdr);
  return result;
}

/***********************************
  trusty usage message 
*********************************/
static void usage(const char *progname) {
  fprintf(stderr,"Usage:\n%s\n",progname);
  fprintf(stderr,"\t-i input_sky_fits_file\n");
  fprintf(stderr,"\t-o output_image_name\n");
  fprintf(stderr,"\t-w width(and height)_output_image\n");
  fprintf(stderr,"\t-t LST (hours)\n");
  fprintf(stderr,"\t-l latitude (degrees)\n");
  fprintf(stderr,"\t-s text_source_file_name\n");
  fprintf(stderr,"\t-a angluar_size_of_pixels (arcsec)\n");
  fprintf(stderr,"\t-x (flip X/RA axis)\n");
  fprintf(stderr,"\t-d (turns debugging on)\n");
  exit(0);
}
