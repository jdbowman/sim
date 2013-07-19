#include <iostream>
#include <vector>

#include <string>
#include <fstream>
#include <sstream>

#include "utils.h"

using namespace std;
void usage () {
  cout << "Simple program to convert the GSM HealPix files to Plan-carre all sky maps" << endl;
  cout << " Usage: " << endl;
  cout << "convert_hpx -i <infile> -o <outfile>" << endl;
}
int main (int argc, char **argv) {    

 

  int gotc; 
  bool got_in = false, got_out = false;


  string *input=NULL,*output=NULL;
  
  while ((gotc = getopt( argc,argv,"hi:o:")) != -1) {

    switch (gotc) {
      case 'h':
	usage();
	exit(0);
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
  long nside = npix2nside(npix);
  cout << "Read in " <<  npix << " pixels" << endl;
  cout << "Which implies nside of " <<  nside << " pixels" << endl;
  cout << "We are assuming ring format GSM data" << endl;


  hpx2carre (&storage[0], nside, output->c_str(), 0) ;

}


