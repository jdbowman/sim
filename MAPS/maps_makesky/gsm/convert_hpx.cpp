#include <iostream>
#include <vector>

#include <cstdlib>

#include <string>
#include <fstream>
#include <sstream>

#include "utils.h"

using namespace std;
void usage () {
  cout << "Simple program to convert the GSM HealPix files to Plan-carre all sky maps or FITs files for later conversion" << endl;
  cout << " Usage: " << endl;
  cout << "convert_hpx -f <freq in MHz> -i <infile> -o <outfile>" << endl;
  cout << "other options:"<< endl;
  cout << "-P\toutput plan carre" << endl;
  cout << "-F\tjust convert to FITs format for later conversion to HPX. The conversion is NOT done by this program but by HPXcvt." << endl;
  cout << "-H\tinput is in HPX projection." << endl;
}
int main (int argc, char **argv) {    



  int gotc; 
  bool got_in = false, got_out = false,got_freq=false;
  bool cvt = true,justFITS = false, isHPX = false;
  float freq;

  string *input=NULL,*output=NULL;

  while ((gotc = getopt( argc,argv,"hFf:Hi:o:P")) != -1) {

    switch (gotc) {
      case 'h':
	usage();
	exit(0);
	break;
      case 'H':
	cvt = false;
	isHPX = true;
	break;
      case 'i':
	input = new string(optarg);
	got_in = true;
	break;
      case 'o':
	output= new string(optarg);
	got_out = true;
	break;
      case 'f':
	freq = atof(optarg);
	got_freq= true;
	break;
      case 'F':
	cvt = false;
	justFITS = true;
	break;
      case 'P':
	cvt = true;
	break;
      default:
	break;
    }
  }

  if (!(got_in || got_out)) {
    usage();
    exit(1);
  }    
  if (!got_freq) {
    fprintf(stderr,"A frequency (in MHz) must be specified\n");
    usage();
    exit(EXIT_FAILURE);
  }
  long nside=0;
  long npix=0;
  vector<float> storage;
  if (cvt | justFITS) {
    ifstream infile;
    infile.open(input->c_str());
    if (!infile.is_open()) {
      cerr << "Cannot open " << input << endl;
    }  

    string s; 
    while (!infile.eof()) {

      getline(infile,s);
      stringstream ss (stringstream::in | stringstream::out);
      ss << s;

      float number;

      ss >> number; // cute c++ method to output the string to a double 
      storage.push_back(number);
    }

    if (infile.is_open()) {
      infile.close();
    }

    npix = long (storage.size());
    nside = npix2nside(npix);
    cout << "Read in " <<  npix << " pixels" << endl;
    cout << "Which implies nside of " <<  nside << " pixels" << endl;
    cout << "We are assuming ring format GSM data" << endl;
  }
  else {
    cout << "Must read in HPX" << endl;
  }
  if (cvt) {
    healpix2carre (&storage[0], nside, output->c_str(), 0,freq) ;
  }
  else if (justFITS) {
    hpxwrap (&storage[0], nside, output->c_str(), 0,freq);
  }
  else if (isHPX) {
    hpx2carre (input->c_str(), output->c_str(), freq) ;
  }
}


