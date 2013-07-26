#include <iostream>
#include <vector>
#include <math.h>

#include <cstdlib>


#include <string>
#include <fstream>
#include <sstream>

#include "utils.h"

using namespace std;
void usage () {
  cout << "Simple program to convert the 360 degree poln map to a fits file" << endl;
  cout << " Usage: " << endl;
  cout << "convert_dat -i <infile> -o <outfile>" << endl;
}
int main (int argc, char **argv) {    

 

  int gotc; 
  bool got_in = false, got_out = false;
  float freq;

  string *input=NULL,*output=NULL;
  
  while ((gotc = getopt( argc,argv,"hf:i:o:")) != -1) {

    switch (gotc) {
      case 'h':
	usage();
	exit(0);
	break;
      case 'f':
	freq = atof(optarg);
	break;
      case 'i':
	input = new string(optarg);
	got_in = true;
	break;
      case 'o':
	output= new string(optarg);
	got_out = true;
	break;
      default:
	break;
    }
  }
  
  if (!(got_in || got_out)) {
    usage();
    exit(1);
  }    

  vector<float> storage;

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
    
  long npix = long (storage.size());
  
  double nside = sqrt(double(npix));
  
  cout << "Read in " <<  npix << " pixels" << endl;
  cout << "Which implies nside of " <<  nside << " pixels" << endl;

  dat2carre (&storage[0], (long int) nside, output->c_str(),freq );
}


