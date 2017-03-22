// main91.cc is a part of the PYTHIA event generator.
// Copyright (C) 2013 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This is a simple test program illustrating how to link in Pythia 6.4 
// for the generation of hard processes. That part is now considered
// obsolete, but for some debug work this example has been collected
// from various code pieces that were separately available up until 
// Pythia 8.125. In addition to modifying the below code to fit your
// needs, you also have to modify examples/Makefile to link properly
// for main91, including the Pythia 6.4xx library to be used (where
// xx is the current subversion number). 

// All hard PYTHIA 6.4 processes should be available for full generation
// in PYTHIA 8, at least to the extent that they are defined for p p,
// pbar p or e+ e-collisions. Soft processes, i.e. elastic and diffractive 
// scattering, as well as minimum-bias events, require a different 
// kinematics machinery, and can only be generated with the internal 
// PYTHIA 8 processes.

// PYTHIA 6.4 does its own rejection of events internally, according to
// the strategy option 3. However, the current runtime interface does not 
// take cross-section information from PYTHIA 6.4. This means that both
// the initial maxima and the final cross sections printed by the PYTHIA 8
// routines are irrelevant in this case. Instead you have to study the
// standard PYTHIA 6.4 initialization printout and call on pystat(...)
// at the end of the run to obtain this information. It also means that 
// you cannot mix with internally generated PYTHIA 8 events.

//==========================================================================

#include "Pythia.h"
#include "LHAFortran.h"

using namespace Pythia8; 

//==========================================================================

// Declare the Fortran subroutines that may be used.
// This code section is generic.

#ifdef _WIN32
  #define pygive_ PYGIVE
  #define pyinit_ PYINIT
  #define pyupin_ PYUPIN
  #define pyupev_ PYUPEV
  #define pylist_ PYLIST
  #define pystat_ PYSTAT
#endif

extern "C" {
#ifdef _WIN32
  extern void pyinit_(const char*, int, const char*, int, const char*, 
    int, double&);
#else
  extern void pyinit_(const char*, const char*, const char*, double&, 
    int, int, int);
#endif
}

extern "C" {
  extern void pygive_(const char*, int);
  extern void pyupin_();
  extern void pyupev_();
  extern void pylist_(int&);
  extern void pystat_(int&);
}

//==========================================================================

// Interfaces to the above routines, to make the C++ calls similar to Fortran.
// This code section is generic.

class Pythia6Interface {

public:

  // Give in a command to change a setting.
  static void pygive(const string cmnd) { 
    const char* cstring = cmnd.c_str(); int len = cmnd.length(); 
    pygive_(cstring, len);
  }

  // Initialize the generation for the given beam confiuration.
  static void pyinit(const string frame, const string beam, 
    const string target, double wIn) { 
    const char* cframe = frame.c_str(); int lenframe = frame.length();
    const char* cbeam = beam.c_str(); int lenbeam = beam.length();
    const char* ctarget = target.c_str(); int lentarget = target.length();
#ifdef _WIN32
    pyinit_(cframe, lenframe, cbeam, lenbeam, ctarget, lentarget, wIn);
#else
    pyinit_(cframe, cbeam, ctarget, wIn, lenframe, lenbeam, lentarget); 
#endif
  }
  
  // Fill the initialization information in the HEPRUP commonblock.
  static void pyupin() {pyupin_();}

  // Generate the next hard process and 
  // fill the event information in the HEPEUP commonblock
  static void pyupev() {pyupev_();}

  // List the event at the process level.
  static void pylist(int mode) {pylist_(mode);}

  // Print statistics on the event generation process.
  static void pystat(int mode) {pystat_(mode);}

};

//==========================================================================

// Implement initialization fillHepRup method for Pythia6 example.
// This code section is specific to the kind of precesses you generate.

// Of all parameters that could be set with pygive, only those that 
// influence the generation of the hard processes have any impact, 
// since this is the only part of the Fortran code that is used. 

bool LHAupFortran::fillHepRup() { 

  // Set process to generate.
  // Example 1: QCD production; must set pTmin.  
  //Pythia6Interface::pygive("msel = 1"); 
  //Pythia6Interface::pygive("ckin(3) = 20.");
  // Example 2: t-tbar production.  
  //Pythia6Interface::pygive("msel = 6"); 
  // Example 3: Z0 production; must set mMin.
  Pythia6Interface::pygive("msel = 11"); 
  Pythia6Interface::pygive("ckin(1) = 50."); 

  // Speed up initialization: multiparton interactions only in C++ code.
  Pythia6Interface::pygive("mstp(81)=0");
    
  // Initialize for 14 TeV pp collider.
  Pythia6Interface::pyinit("cms","p","p",14000.);   

  // Fill initialization information in HEPRUP.
  Pythia6Interface::pyupin();

  // Done.
  return true;

}

//==========================================================================

// Implement event generation fillHepEup method for Pythia 6 example.
// This code section is generic.

bool LHAupFortran::fillHepEup() { 

  // Generate and fill the next Pythia6 event in HEPEUP.
  Pythia6Interface::pyupev();

  // Done.
  return true;

}

//==========================================================================

// The main program.
// This code section is specific to the physics study you want to do.

int main() {

  // Generator. Shorthand for the event and for settings.
  Pythia pythia;
  Event& event = pythia.event;
  Settings& settings = pythia.settings;

  // Set Pythia8 generation aspects. Examples only. 
  pythia.readString("BeamRemnants:primordialKThard = 2.");    
  pythia.readString("MultipartonInteractions:bProfile = 3");    
  pythia.readString("Next:numberShowInfo = 0"); 
  pythia.readString("Next:numberShowProcess = 0"); 
  pythia.readString("Next:numberShowEvent = 0"); 

  // Initialize to access Pythia6 generator by Les Houches interface.
  pythia.readString("Beams:frameType = 5"); 
  LHAupFortran pythia6;
  pythia.setLHAupPtr( &pythia6);
  pythia.init();    

  // Set some generation values.
  int nEvent = 100;
  int nList  = 1;
  int nAbort = 10;

  // List changed settings data.
  settings.listChanged();

  // Histograms.
  double eCM = 14000.;
  double epTol = 1e-7 * eCM;
  Hist epCons("deviation from energy-momentum conservation",100,0.,epTol);
  Hist nFinal("final particle multiplicity",100,-0.5,1599.5);
  Hist nChg("final charged multiplicity",100,-0.5,799.5);

  // Begin event loop.
  int iAbort = 0; 
  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

    // Generate events. Quit if too many failures.
    if (!pythia.next()) {
      if (++iAbort < nAbort) continue;
      cout << " Event generation aborted prematurely, owing to error!\n"; 
      break;
    }
 
    // List first few events, both hard process and complete events.
    if (iEvent < nList) { 
      pythia.info.list();
      // This call to Pythia6 is superfluous, but shows it can be done.
      Pythia6Interface::pylist(1);
      pythia.process.list();
      event.list();
    }

    // Reset quantities to be summed over event.
    int nfin = 0;
    int nch = 0;
    Vec4 pSum = - (event[1].p() + event[2].p());

    // Loop over final particles in the event. 
    for (int i = 0; i < event.size(); ++i) 
    if (event[i].isFinal()) {
      ++nfin;
      if (event[i].isCharged()) ++nch;
      pSum += event[i].p();
    }

    // Fill summed quantities. 
    double epDev = abs(pSum.e()) + abs(pSum.px()) + abs(pSum.py())
      + abs(pSum.pz());
    epCons.fill(epDev);
    nFinal.fill(nfin);
    nChg.fill(nch);

  // End of event loop.
  }

  // Final statistics. Must do call to Pythia6 explicitly.
  pythia.stat();
  Pythia6Interface::pystat(1);  

  // Histogram output.
  cout << epCons << nFinal<< nChg; 

  // Done.
  return 0;
}
