// LHAPDF5.h is a part of the PYTHIA event generator.
// Copyright (C) 2015 Torbjorn Sjostrand.
// PYTHIA is licenced under the GNU GPL version 2, see COPYING for details.
// Please respect the MCnet Guidelines, see GUIDELINES for details.

// This file contains the LHAPDF5 PDF plugin class.

#ifndef Pythia8_LHAPDF5_H
#define Pythia8_LHAPDF5_H

#include "Pythia8/PartonDistributions.h"

namespace Pythia8 {

//==========================================================================

// Interfaces to to make the C++ calls similar to f77.

//--------------------------------------------------------------------------

// Declare the LHAPDF5 f77 subroutines that are needed.

extern "C" {

  extern void initpdfsetm_(int&, const char*, int);

  extern void initpdfsetbynamem_(int&, const char*, int);

  extern void initpdfm_(int&, int&);

  extern void evolvepdfm_(int&, double&, double&, double*);

  extern void evolvepdfphotonm_(int&, double&, double&, double*, double&);

  extern void setlhaparm_(const char*, int);

}

//--------------------------------------------------------------------------

// Map the f77 routines to C++.

namespace LHAPDF5Interface {

  // Initialize set with full pathname, allowing multiple sets.
  void initPDFsetM( int& nSet, string name) {
    const char* cName = name.c_str(); int lenName = name.length();
    initpdfsetm_( nSet, cName, lenName);
  }

  // Initialize set with simple name, allowing multiple sets.
  void initPDFsetByNameM( int& nSet, string name) {
    const char* cName = name.c_str(); int lenName = name.length();
    initpdfsetbynamem_( nSet, cName, lenName);
  }

  // Initialize member of set.
  void initPDFM(int& nSet, int member) {
    initpdfm_(nSet, member);
  }

  // Evaluate x f_i(x, Q).
  void evolvePDFM( int& nSet, double x, double Q, double* xfArray) {
    evolvepdfm_( nSet, x, Q, xfArray);
  }

  // Evaluate x f_i(x, Q) including photon
  void evolvePDFPHOTONM(int& nSet, double x, double Q, double* xfArray,
    double& xPhoton) {
    evolvepdfphotonm_( nSet, x, Q, xfArray, xPhoton);
  }

  // Extrapolate PDF set beyond boundaries, or freeze them there.
  void setPDFparm(string name) {
    const char* cName = name.c_str(); int lenName = name.length();
    setlhaparm_( cName, lenName);
  }

  // Global tracking of opened PDF sets.
  map< int, pair<string, int> > initializedSets;

  // Method to find the nSet number corresponding to a name and member.
  // Returns -1 if no such LHAPDF5 set has been initialized.
  int findNSet(string setName, int member) {
    for (map<int, pair<string, int> >::const_iterator
         i = initializedSets.begin();
         i != initializedSets.end(); ++i) {
      int    iSet    = i->first;
      string iName   = i->second.first;
      int    iMember = i->second.second;
      if (iName == setName && iMember == member) return iSet;
    }
    return -1;
  }

  // Method to return the lowest non-occupied nSet number.
  int freeNSet() {
    for (int iSet = 1; iSet <= int(initializedSets.size()); ++iSet) {
      if (initializedSets.find(iSet) == initializedSets.end())
        return iSet;
    }
    return initializedSets.size() + 1;
  }

}

//==========================================================================

// Plugin interface to the LHAPDF5 library.

//==========================================================================

// Provide plugin interface to the LHAPDF5 library of parton densities.

class LHAPDF5 : public PDF {

public:

  // Constructor.
  LHAPDF5(int idBeamIn, string setName, int member,  int nSetIn = -1,
    Info* infoPtr = 0) : PDF(idBeamIn), nSet(nSetIn)
    {init(setName, member, infoPtr);}

  // Allow extrapolation beyond boundaries. This is optional.
  void setExtrapolate(bool extrapol);

private:

  // Initialization of PDF set.
  void init(string setName, int member, Info* infoPtr);

  // Update all PDF values.
  void xfUpdate(int , double x, double Q2);

  // Current set and pdf values.
  int    nSet;
  double xfArray[13];
  bool   hasPhoton;
  double xPhoton;

};

//--------------------------------------------------------------------------

// Initialize a parton density function from LHAPDF5.

void LHAPDF5::init(string setName, int member, Info*) {

  // Determine whether the pdf set contains the photon or not.
  // So far only MRST2004QED and  NNPDF2.3QED.
  if ( setName == "MRST2004qed.LHgrid"
    || setName == "NNPDF23_lo_as_0130_qed.LHgrid"
    || setName == "NNPDF23_lo_as_0130_qed_mem0.LHgrid"
    || setName == "NNPDF23_lo_as_0119_qed_mem0.LHgrid"
    || setName == "NNPDF23_lo_as_0130_qed_mem0.LHgrid"
    || setName == "NNPDF23_nlo_as_0119_qed_mc.LHgrid"
    || setName == "NNPDF23_nlo_as_0119_qed_mc_mem0.LHgrid"
    || setName == "NNPDF23_nnlo_as_0119_qed_mc.LHgrid"
    || setName == "NNPDF23_nnlo_as_0119_qed_mc_mem0.LHgrid" ) hasPhoton = true;
  else hasPhoton = false;

  // If already initialized then need not do anything further.
  pair<string, int> initializedNameMember =
    LHAPDF5Interface::initializedSets[nSet];
  string initializedSetName   = initializedNameMember.first;
  int    initializedMember    = initializedNameMember.second;
  if (setName == initializedSetName && member == initializedMember) return;

  // Initialize set. If first character is '/' then assume that name
  // is given with path, else not.
  if (setName[0] == '/') LHAPDF5Interface::initPDFsetM( nSet, setName);
  else LHAPDF5Interface::initPDFsetByNameM( nSet, setName);
  isSet = (nSet >= 0);

  // Initialize member.
  LHAPDF5Interface::initPDFM(nSet, member);

  // Do not collect statistics on under/overflow to save time and space.
  LHAPDF5Interface::setPDFparm( "NOSTAT" );
  LHAPDF5Interface::setPDFparm( "LOWKEY" );

  // Save values to avoid unnecessary reinitializations.
  if (nSet > 0) LHAPDF5Interface::initializedSets[nSet] =
                  make_pair(setName, member);

}

//--------------------------------------------------------------------------

// Allow optional extrapolation beyond boundaries.

void LHAPDF5::setExtrapolate(bool extrapol) {

   LHAPDF5Interface::setPDFparm( (extrapol) ? "EXTRAPOLATE" : "18" );

}

//--------------------------------------------------------------------------

// Give the parton distribution function set from LHAPDF5.

void LHAPDF5::xfUpdate(int, double x, double Q2) {

  // Let LHAPDF5 do the evaluation of parton densities.
  double Q = sqrt( max( 0., Q2));

  // Use special call if photon included in proton.
  if (hasPhoton) {
    LHAPDF5Interface::evolvePDFPHOTONM( nSet, x, Q, xfArray, xPhoton);
  }
  // Else use default LHAPDF5 call.
  else {
    LHAPDF5Interface::evolvePDFM( nSet, x, Q, xfArray);
    xPhoton=0.0;
  }

  // Update values.
  xg     = xfArray[6];
  xu     = xfArray[8];
  xd     = xfArray[7];
  xs     = xfArray[9];
  xubar  = xfArray[4];
  xdbar  = xfArray[5];
  xsbar  = xfArray[3];
  xc     = xfArray[10];
  xb     = xfArray[11];
  xgamma = xPhoton;

  // Subdivision of valence and sea.
  xuVal  = xu - xubar;
  xuSea  = xubar;
  xdVal  = xd - xdbar;
  xdSea  = xdbar;

  // idSav = 9 to indicate that all flavours reset.
  idSav = 9;

}

//--------------------------------------------------------------------------

// Define external handles to the plugin for dynamic loading.

extern "C" {

  LHAPDF5* newLHAPDF(int idBeamIn, string setName, int member,
                    Info* infoPtr) {
    int nSet = LHAPDF5Interface::findNSet(setName, member);
    if (nSet == -1) nSet = LHAPDF5Interface::freeNSet();
    return new LHAPDF5(idBeamIn, setName, member, nSet, infoPtr);
  }

  void deleteLHAPDF(LHAPDF5* pdf) {
    delete pdf;
  }

}

//==========================================================================

} // end namespace Pythia8

#endif // end Pythia8_LHAPDF5_H
