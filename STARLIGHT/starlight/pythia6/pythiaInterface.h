//#include "Pythia.h"
//using namespace Pythia8; 

//==========================================================================

// Declare the Fortran subroutines that may be used.
// This code section is generic.

#include <string>

extern "C" {
  extern void pyinit_(const char*, const char*, const char*, double&, int, int, int);
  extern void pyevnt_();
  extern void pygive_(const char*, int);
  extern void pyfram_(int&);
  extern void pylist_(int&);
  extern void pystat_(int&);
  extern int pycomp_(int&);
  
  extern struct
        {
	  int n;
	  int npad;
	  int k[5][4000];
	  double p[5][4000];
	  double v[5][4000];
        } pyjets_;
    extern struct
    {
      int mdcy[3][500];
      int mdme[2][8000];
      double brat[8000];
      int kfpd[5][8000];
    } pydat3_;
}


class pythiaInterface {

public:

  // Give in a command to change a setting.
  static void pygive(const std::string cmnd) { 
    const char* cstring = cmnd.c_str(); int len = cmnd.length(); 
    pygive_(cstring, len);
  }

  // Initialize the generation for the given beam confiuration.
  static void pyinit(const std::string frame, const std::string beam, 
    const std::string target, double wIn) { 
    const char* cframe = frame.c_str(); int lenframe = frame.length();
    const char* cbeam = beam.c_str(); int lenbeam = beam.length();
    const char* ctarget = target.c_str(); int lentarget = target.length();
    pyinit_(cframe, cbeam, ctarget, wIn, lenframe, lenbeam, lentarget); 
  }
  
  static void pyevnt() {pyevnt_();}
  
  static void pyfram(int frame) { pyfram_(frame); }
  
  // List the event at the process level.
  static void pylist(int mode) {pylist_(mode);}

  // Print statistics on the event generation process.
  static void pystat(int mode) {pystat_(mode);}

  // Get compressed code (KC) from PDG code
  static int pycomp(int pdg) { return pycomp_(pdg);}

  
};
