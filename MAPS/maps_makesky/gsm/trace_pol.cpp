#include <iostream>
#include <string>
#include <math.h>
extern "C" {
#include "cblas.h"
}
#include "fitsio.h"

using namespace std;

#define SOL 299792458

static int read_input_image(const char *filename,float **data,long *imsize);

static int unload_file(const char *inputfile, const char *outputfile, float *data, long *imsize);


void usage() {
  cout << "trace_pol <config_file>" << endl;
  cout << "-q <file> input Stokes Q" << endl;
  cout << "-u <file> input Stokes U" << endl;
  cout << "-r <file> RM screen\n" << endl;
  cout << "----------------------------" << endl;
  cout << "-     Processing options   -" << endl;
  cout << "----------------------------" << endl;
  cout << endl;
  cout << "-f <freq> Centre Frequency in *MHz* [100]" << endl;
  cout << "-b <bw> channel BW in *kHz* [40]" << endl;
  cout << "-c <>   number of subchannels" << endl;
  cout << "-s <>   number of Faraday screens" << endl;
}

int main (int argc, char **argv) {
   
  enum CBLAS_ORDER order;
  enum CBLAS_TRANSPOSE transa;

 string *q=NULL,*u=NULL,*r=NULL;
 float freq=10.0,bw=40.0,cbw=1.0;
 int nsub=20;
 char *optstring = "hq:u:r:f:b:c:s:";
 int c;
 int screens=1;
 
 while ((c = getopt(argc,argv,optstring)) != -1) {
   switch (c) {
     case 'q':
       q = new string(optarg);
       break;
     case 'u':
       u = new string(optarg);
       break;
     case 'r':
       r = new string(optarg);
       break;
     case 'f':
       freq = atof(optarg);
       break;
     case 'b':
       bw = atof(optarg);
       break;
     case 'c':
       nsub = atoi(optarg);
       break;
     case 's':
       screens = atoi(optarg);
       break;
     default:
       usage();
       exit(0);
   }
 }
 
 if (!q || !u || !r) {
   cerr << endl;
   cerr << "Need a Q and U and Faraday screen" << endl;
   usage();
   exit(0);
 }
   

 cout << "Q file is " << q->c_str() << endl;
 cout << "U file is " << u->c_str() << endl;
 cout << "Faraday rotation screen file is " << r->c_str() << endl;
 cout << "Channel centre (MHz) = " << freq << endl;
 cout << "Channel width (kHz) = " << bw << endl;
 cout << "Artificial subchannels (num) = " << nsub << endl;
 cbw = bw/nsub;
 cout << "Artificial subchannels (width) = " << cbw << endl;
 // Lets open up the files

 float *qpixels;
 float *upixels;
 float *rmpixels;
 long imsize[2];
 
 // the Q file
 read_input_image(q->c_str(),&qpixels,imsize);
  // the U file
 read_input_image(u->c_str(),&upixels,imsize);
 // the screen
 read_input_image(r->c_str(),&rmpixels,imsize);

 long npix = imsize[0]*imsize[1];

 // a lot of  these are nan
 
 float input_stokes[4];
 float output_stokes[4];
 float answer[4];
 float scale[4];
 float a[16];
 double lambda=0;
 double delta_x=0;
 int debug = 0; 
 int m,n,lda,incx,incy;
 float alpha,beta;
 
 m=4;
 n=4;
 lda=4;
 incx=1;
 incy=1;
 alpha=1;
 beta=0;
 

   
 order = CblasColMajor;
 transa = CblasNoTrans;

 float stokes_total[4];

 // The output arrays

 float *out_q;
 float *out_u;

 out_q = (float *) malloc(npix*sizeof(float));
 out_u = (float *) malloc(npix*sizeof(float));


 for (long pix=0;pix<npix;pix++) {
   
   stokes_total[0]=0.0;
   stokes_total[1]=0.0;
   stokes_total[2]=0.0;
   stokes_total[3]=0.0;


   if (isnan(rmpixels[pix])) {
     out_u[pix] = nan(NULL);
     out_q[pix] = nan(NULL);
     continue;
   }
   else {
     for (int nscreen = 1; nscreen<=screens; nscreen=nscreen+1) {

       input_stokes[0] = (1.0/screens + stokes_total[0])/nsub;
       input_stokes[1] = (qpixels[pix]/screens + stokes_total[1])/nsub;
       input_stokes[2] = (upixels[pix]/screens + stokes_total[2])/nsub;
       input_stokes[3] = (0.0 + stokes_total[3])/nsub;

       output_stokes[0] = 0.0;
       output_stokes[1] = 0.0;
       output_stokes[2] = 0.0;
       output_stokes[3] = 0.0;

       int m = 4;   
       for (int subc = 0 ; subc < nsub ; subc++) {
	 float subfreq = (freq*1000.0) - ((cbw*nsub/2.0)) + (subc*cbw); // kHz
	 lambda = SOL/(subfreq*1000.0);
	 // lambda = lambda/1000.0; // units ....

	 delta_x = lambda*lambda*rmpixels[pix]/screens;  

	 if (debug) {
	   cout << "Sub channel no." << subc << endl;
	   cout << "Frequency: " << subfreq << " kHz" << endl;
	   cout << "Wavelength: " << lambda << " m" << endl;
	   cout << "RM: " << rmpixels[pix]/screens << " radm-2 " << endl;
	   cout << "Rotation angle: " << delta_x << " rad " << endl;
	   
	   cout << "Rotation angle: " << delta_x*180/M_PI << " deg" << endl;

	   double turns = floor(delta_x/(2*M_PI));
	   double resid = delta_x - turns*2*M_PI;
           delta_x=resid;
	   cout << "Rotation angle: " << delta_x << " rad " << endl;
	   
	   cout << "Rotation angle: " << delta_x*180/M_PI << " deg" << endl;



	 }

	 // ok now we have to build the rotation matrix
	 a[0] = 1;
	 a[1] = 0;
	 a[2] = 0;
	 a[3] = 0;
	 /* The elements of the second column */
	 a[m] = 0;
	 a[m+1] = cos(2*delta_x);
	 a[m+2] = -1*sin(2*delta_x);
	 a[m+3] = 0;
	 /* The elements of the third column */
	 a[m*2] = 0;
	 a[m*2+1] = sin(2*delta_x);  
	 a[m*2+2] = cos(2*delta_x);
	 a[m*2+3] = 0;
	 /* The elements of the fourth column */
	 a[m*3] = 0;
	 a[m*3+1] = 0;
	 a[m*3+2] = 0;
	 a[m*3+3] = 1;
	 /* the answer */
	 answer[0] = 0;
	 answer[1] = 0;
	 answer[2] = 0;
	 answer[3] = 0;

	 cblas_sgemv(order,transa,m,n,alpha,a,lda,input_stokes,incx,beta,answer,incy);

	 output_stokes[0]=output_stokes[0] + answer[0];
	 output_stokes[1]=output_stokes[1] + answer[1];
	 output_stokes[2]=output_stokes[2] + answer[2];
	 output_stokes[3]=output_stokes[3] + answer[3];

       }
       stokes_total[0] = output_stokes[0];
       stokes_total[1] = output_stokes[1];
       stokes_total[2] = output_stokes[2];
       stokes_total[3] = output_stokes[3];
     }

     if (debug) {
       for (int i=0;i<4;i++) {
	 cout << "S[" << i << "] = " << stokes_total[i] << endl;
       }
      
       float inpol =  sqrtf(qpixels[pix]*qpixels[pix] + upixels[pix]*upixels[pix]);
       float outpol = sqrtf(stokes_total[1]*stokes_total[1] + stokes_total[2]*stokes_total[2]);

       if ((1 - outpol/inpol) > 0.1) {  
	 cout << "OLD " << sqrtf(qpixels[pix]*qpixels[pix] + upixels[pix]*upixels[pix]) << " NEW " << sqrtf(stokes_total[1]*stokes_total[1] + stokes_total[2]*stokes_total[2]) << " DEPOL: " << (1-outpol/inpol) * 100.0 << "%" << endl;;
       }
     }
     
     out_q[pix] = stokes_total[1]/stokes_total[0];
     out_u[pix] = stokes_total[2]/stokes_total[0];


   }
 }
	 	
 
 unload_file(q->c_str(),"!test_q.fits", out_q, imsize);
 unload_file(u->c_str(),"!test_u.fits", out_u, imsize);
	   
 exit(1);
}

void rotate(float *stokes, float RM, float *answer) {
  cout << "NOT YET ROTATING" << endl;
}
static int read_input_image(const char *filename,float **data,long *imsize) {
  fitsfile *infp=NULL;
  int result=0,n_nulls=0;
  long fpixel[3]={1,1,1}; /* use 3 just in case we have a wavelength dimension */
  float nullval=0;
  
  int debug = 0;

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
