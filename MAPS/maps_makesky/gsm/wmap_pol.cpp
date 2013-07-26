#include <iostream>
#include <vector>

#include <string>
#include <fstream>
#include <sstream>
#include <math.h>

#include "fitsio.h"
#include "utils.h"

using namespace std;

#define NFILES 4


void usage () {
  cout << "Simple program to convert the GSM HealPix files to Plan-carre all sky maps" << endl;
  cout << " Usage: " << endl;
  cout << "wmap_pol -I <input stokes I> -i <root> -o <outfile> [-p PA file]" << endl;
  cout << "\nWe assume filename is of the form <root>_I.fits,<root>_Q.fits, <root>_U.fits" << endl;
  cout << "The output format is the same as the input but the 2 binary tables are percentage polarisation\
 and polarisation position angle\n"; 
  cout << "PA file is optional and can be constructed from the input Q and U. This allows a smoothed PA to be used if required" << endl;
}

static int read_input_image(const char *filename,float **data,long *imsize);

static int unload_file(const char *inputfile, const char *outputfile, float *data, long *imsize);

int main (int argc, char **argv) {    

 

  int gotc; 
  bool got_in = false, got_out = false;
  float sig_0=1.5; 

  string *input=NULL,*output=NULL,*stokes_I=NULL,*pa_file=NULL;
  bool got_stokes_I = false;
  bool got_noise = true;
  bool got_pa = false;
  
  while ((gotc = getopt( argc,argv,"hI:i:No:p:s:")) != -1) {

    switch (gotc) {
      case 'h':
	usage();
	exit(0);
	break;
      case 'i':
	input = new string(optarg);
	got_in = true;
	break;
	case 'N':
	got_noise = false;
	break;
      case 'o':
	output= new string(optarg);
	got_out = true;
	break;
      case 'p':
	got_pa= true;
	pa_file = new string (optarg);
	break;
      case 's':
	sig_0 = atof(optarg);
	break;
      case 'I':
	stokes_I = new string(optarg);
        got_stokes_I = true;
	break;
      default:
	break;
    }
  }
  
  if (!(got_in || got_out)) {
    usage();
    exit(1);
  }    
  if (!got_stokes_I) {
    usage();
    exit(1);
  }
  /* have to allocate storage for and open I,Q,U input files */
  string stokes[NFILES];
  string derived[2];
  string output_files[4];
  
  float *data[NFILES]; // array of NFILES pointers to floats
  float *angelica;
  
  long imsize[2];
  float *pa=NULL;
  float *degree=NULL;
  
  char *tempname = new char [64];
  sprintf(tempname,"%s_I.fits",input->c_str()); 
  stokes[0] = tempname;

  sprintf(tempname,"%s_Q.fits",input->c_str()); 
  stokes[1] = tempname;

  sprintf(tempname,"%s_U.fits",input->c_str()); 
  stokes[2] = tempname;
  
  if (got_noise == true) {
    sprintf(tempname,"%s_NOBS.fits",input->c_str()); 
    stokes[3] = tempname;
  }
  
  sprintf(tempname,"%s_DEG.fits",output->c_str());
  derived[0] = tempname;

 
  sprintf(tempname,"%s_PA.fits",output->c_str());
  derived[1] = tempname;

  sprintf(tempname,"%s_I.fits",output->c_str()); 
  output_files[0] = tempname;

  sprintf(tempname,"%s_Q.fits",output->c_str()); 
  output_files[1] = tempname;

  sprintf(tempname,"%s_U.fits",output->c_str()); 
  output_files[2] = tempname;
 
  sprintf(tempname,"%s_V.fits",output->c_str()); 
  output_files[3] = tempname;



  delete(tempname);

  /* begin loop to open all files */
  int fcount = 0;
  cout << "Going to open ..." << endl;
  int nfiles = NFILES;
  if (!got_noise) {
    nfiles--;
  }
  while (fcount < nfiles) {
    
    cout << stokes[fcount] << endl;  
    
    // the images are all the same so I can avoid the WCS overhead
    // by just opening up and doing this pixel by pixel
    read_input_image(stokes[fcount].c_str(),&data[fcount],imsize);
    
    fcount++;
  }
  
  read_input_image(stokes_I->c_str(),&angelica,imsize);
  if (got_pa) {
    read_input_image(pa_file->c_str(),&pa,imsize);
  }
  else {
    pa = new float [(imsize[0]*imsize[1])];
  }
  degree = new float [(imsize[0]*imsize[1])];
  
  float * stokes_q = new float [(imsize[0]*imsize[1])];
  float * stokes_u = new float [(imsize[0]*imsize[1])];
  float * stokes_v = new float [(imsize[0]*imsize[1])];
  /* the images are square .... */

  long pcount = 0;

  while (pcount < imsize[0]) {
    for (int i=0l ; i < imsize[1]; i++) {
      size_t offset = (pcount*imsize[1] + i);
      float noise = 0;
      
      if (got_noise) {
	noise = sig_0/sqrtf(data[3][offset]);
      }
      else {
	noise = data[0][offset];
      }
      
      float snr=0; 
      float pol_ang = 0;
      if (isnan(data[0][offset])) {
	degree[offset] = nanf("0");
	if (!got_pa) {
	  pa[offset] = nanf("0");
	}
      }
      else {
        snr = data[0][offset]/noise;
        if (snr >= 1.0) {	
	  if (!got_pa) {
	    pol_ang =  0.5*atan2(data[2][offset],data[1][offset]);
	    pa[offset] = pol_ang;
	  }
	  degree[offset]  = (sqrtf(data[1][offset]*data[1][offset] + data[2][offset]*data[2][offset])/(0.9*fabs(data[0][offset])));
	}
	else {
	    degree[offset] = 0;
	    pa[offset] = 0;
	}
      }
    }
    pcount++;
  }  
  pcount = 0; 
  while (pcount < imsize[0]) {
    for (int i=0l ; i < imsize[1]; i++) {
      size_t offset = (pcount*imsize[1] + i);
      if (isnan(angelica[offset])) {
	stokes_q[offset] = nanf("0");
	stokes_u[offset] = nanf("0");
	stokes_v[offset] = nanf("0");
      }
      else {
	float k= sqrtf(1+(tanf(2*pa[offset])*tanf(2*pa[offset])));
	float q = 0;
	if (fabs(k)>1e-6) {
	  q = angelica[offset]*degree[offset]/k;
        
	  stokes_q[offset] = q;

	  stokes_u[offset] = stokes_q[offset]*tanf(2*pa[offset]);

	  stokes_v[offset] = 0.0;
	}
	else {
	  stokes_q[offset] = 0;

	  stokes_u[offset] = 0.0;

	  stokes_v[offset] = 0.0;
        }

      }
  
    }
    pcount++;
  }
      
  unload_file(stokes[0].c_str(),derived[0].c_str(),degree,imsize); 
  unload_file(stokes[0].c_str(),derived[1].c_str(),pa,imsize); 
  
  
  unload_file(stokes[0].c_str(),output_files[0].c_str(),angelica,imsize); 
  unload_file(stokes[0].c_str(),output_files[1].c_str(),stokes_q,imsize); 
  unload_file(stokes[0].c_str(),output_files[2].c_str(),stokes_u,imsize); 
  unload_file(stokes[0].c_str(),output_files[3].c_str(),stokes_v,imsize); 

  return EXIT_SUCCESS;

}

static int read_input_image(const char *filename,float **data,long *imsize) {
  fitsfile *infp=NULL;
  int result=0,n_nulls=0;
  long fpixel[3]={1,1,1}; /* use 3 just in case we have a wavelength dimension */
  float nullval=0;
  
  int debug = 1;

  fits_open_image(&infp,filename,READONLY,&result);
  if (result !=0) {
    fprintf(stderr,"could not open input file %s\n",filename);
    exit(1);
  }
  fits_get_img_size(infp, 2, imsize, &result);
  if (result !=0) {
    fprintf(stderr,"error getting size of image. err code: %d\n",result);
    goto EXIT;
  }
  if (debug) fprintf(stderr,"input image is %ldx%ld pixels.\n",imsize[0],imsize[1]);
  *data = (float *) malloc(sizeof(float)*imsize[0]*imsize[1]);
  if (*data ==NULL) {
    fprintf(stderr,"read_input_image: no malloc\n");
    exit(1);
  }
  /* extract some info about the coord system */
  
  fits_read_pix(infp,TFLOAT,fpixel,imsize[0]*imsize[1],&nullval,*data,&n_nulls,&result);
  
  if (result !=0) {
    fprintf(stderr,"read_input_image: error reading data. err code: %d\n",result);
  }

EXIT:
  fits_close_file(infp,&result);
  return result;
}


static int unload_file(const char *inputfile, const char *outputfile, float *data, long *imsize){

  fitsfile *fptr1, *fptr2;
  long fpixel[3]={1,1,1};
  long npix=0;
  int status=0;

  if (fits_open_file(&fptr1,inputfile,0,&status)) {
    fprintf(stderr,"fits_open_file has failed\n");
    fits_report_error(stderr,status);
    return EXIT_FAILURE;
  } 

  if (fits_create_file(&fptr2,outputfile,&status)) {
    fprintf(stderr,"fits_create_file has failed\n");
    fits_report_error(stderr,status);
    return EXIT_FAILURE;
  }

  if (fits_copy_header(fptr1,fptr2,&status)) {
    fprintf(stderr,"fits_copy_header has failed\n");
    fits_report_error(stderr,status);
    return EXIT_FAILURE;
  }

  npix = imsize[0]*imsize[1];

  if (fits_write_pix(fptr2,TFLOAT,fpixel,npix,data,&status)) {
    fprintf(stderr,"fits_write_pix has failed\n");
    fits_report_error(stderr,status);
    return EXIT_FAILURE;
  }

  if (fits_close_file(fptr2,&status)) {
    fprintf(stderr,"fits_close_file has failed on the new file\n");
    fits_report_error(stderr,status);
    return EXIT_FAILURE;
  }

  if (fits_close_file(fptr1,&status)) {
    fprintf(stderr,"fits_close_file has failed on the old file\n");
    fits_report_error(stderr,status);
    return EXIT_FAILURE;
  }

  return EXIT_SUCCESS;

}
