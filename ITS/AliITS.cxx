/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/*
$Log$
Revision 1.10  2000/01/19 17:16:51  fca
Introducing a list of lists of hits -- more hits allowed for detector now

Revision 1.9  1999/11/14 14:33:25  fca
Correct problems with distructors and pointers, thanks to I.Hrivnacova

Revision 1.8  1999/09/29 09:24:19  fca
Introduction of the Copyright and cvs Log

*/

///////////////////////////////////////////////////////////////////////////////
//
//      An overview of the basic philosophy of the ITS code development
// and analysis is show in the figure below.
//Begin_Html
/*
<img src="picts/ITS/ITS_Analysis_schema.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>Roberto Barbera is in charge of the ITS Offline code (1999).
<a href="mailto:roberto.barbera@ct.infn.it">Roberto Barbera</a>.
</font>
<pre>
*/
//End_Html
//
//  AliITS. Inner Traking System base class.
//  This class contains the base procedures for the Inner Tracking System
//
//Begin_Html
/*
<img src="picts/ITS/AliITS_Class_Diagram.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>This show the class diagram of the different elements that are part of
the AliITS class.
</font>
<pre>
*/
//End_Html
//
// Version: 0
// Written by Rene Brun, Federico Carminati, and Roberto Barbera
//
// Version: 1
// Modified and documented by Bjorn S. Nilsen
// July 11 1999
//
// AliITS is the general base class for the ITS. Also see AliDetector for
// futher information.
//
///////////////////////////////////////////////////////////////////////////////
 
#include <TMath.h>
#include <TRandom.h>
#include <TVector.h>
#include <TGeometry.h>
#include <TNode.h>
#include <TTUBE.h>

#include "AliITSmodule.h"
#include "AliDetector.h"
#include "AliITS.h"
#include "TClonesArray.h"
#include "TObjArray.h"
#include "AliITShit.h"
#include "AliITSdigit.h"
#include "AliRun.h"

ClassImp(AliITS)

////////////////////////////////////////////////////////////////////////
//
//      An overview of the basic philosophy of the ITS code development
// and analysis is show in the figure below.
//Begin_Html
/*
<img src="picts/ITS/ITS_Analysis_schema.gif">
</pre>
<br clear=left>
<font size=+2 color=red>
<p>Roberto Barbera is in charge of the ITS Offline code (1999).
<a href="mailto:roberto.barbera@ct.infn.it">Roberto Barbera</a>.
</font>
<pre>
*/
//End_Html
//
// Version: 0
// Written by Rene Brun, Federico Carminati, and Roberto Barbera
//
// Version: 1
// Modified and documented by Bjorn S. Nilsen
// July 11 1999
//
// AliITS is the general base class for the ITS. Also see AliDetector for
// futher information.
//
// Data members:
//
// AliITSgeom *fITSgeom
//     All of the geometry information about the active volumes that
// make up the ITS are described in the AliITSgeom class. This includes
// the transformation functions between the local and global coordinate
// systems and the like. See the full description found in the AliITSgeom
// class. Here in the AliITS class is kept the pointer to the geometry
// used in the simulations or that thought to be the correct one for the
// data. Until a more general class is define and a more appropriate
// place is found to keep this geometry information, it is kept here in
// AliITS.
//
// TObjArray *fITSpoints
//     This is a pointer to the points, to be used by the tracking algorithms
// for example, found in the detectors of the ITS. To allow for the most
// general points structure it is defined to be a pointer to a TObjArray where
// each array element would be one point found in the ITS detectors. An
// Addpoints function is defined below. By default an array of 16 TObjects are
// defined during the initialization of AliITS. This is automatically expanded
// when necessary by the Addpoints function.
//
// Bool_t fEuclidOut
//     This is a flag used to indicate that an Euclid compatible CAD
// file will be created upon the creation of the ITS Monte Carlo
// geometry definition, in the function CreatGeometry. If fEuclidOut is
// true, then a file called ITSgeometry.euc will be created.
//
// Int_t fIdN
//     This variable contains the number of layers defined in the ITS
// geometry. It is primarily used as a size indicator for fIdSens and
// fIdName described below. In general the number of layers, ladders, or
// detectors should be gotten from the AliITSgeom functions. Upon
// creating the AliITS object it is set to zero.
//
// Int_t *fIdSens
//     This is a pointer to an array containing the Monte Carlo volume
// numbers for the different layers of the ITS. These numbers are needed
// by the StepManager function to determine what layer a hit was on. It
// is sized and initialized in the Init function and the AliITSv? Init
// function, called after a call to CreateGeometry. Upon creating the
// AliITS object it points to zero. This variable is made a pointer
// in order to keep the maximum flexibility at this level of the code.
//
// char **fIdName
//     This is a pointer to an array of characters containing the names of
// the different ITS layers as defined in the Monte Carlo geometry data-
// base. It is sized and filled in the AliITSv? Init function, called
// after a call to CreatGeometry. Upon creating the AliITS object it
// points to zero. This variable is make a pointer in order to keep the
// maximum flexibility at this level of the code.
//
// Member Functions:
//
// AliITS()
//     The default constructor of the AliITS class. In addition to
// creating the AliITS class it zeros the variables fIshunt (a member
// of AliDetector class), fEuclidOut, and fIdN, and zeros the pointers
// fITSpoints, fIdSens, and fIdName. The AliDetector default constructor
// is also called.
//
// AliITS(const char *name, const char *title)
//     The constructor of the AliITS class. In addition to creating the
// AliITS class, it allocates memory for the TClonesArrays fHits and
// fDigits, and for the TObjArray fITSpoints. It also zeros the variables
// fIshunt (a member of AliDetector class), fEuclidOut, and fIdN, and zeros
// the pointers fIdSens and fIdName. To help in displaying hits via the ROOT
// macro display.C AliITS also sets the marker color to red. The variables
// passes with this constructor, const char *name and *title, are used by
// the constructor of AliDetector class. See AliDetector class for a
// description of these parameters and its constructor functions.
//
// ~AliITS()
//     The default destructor of the AliITS class. In addition to deleting
// the AliITS class it deletes the memory pointed to by the fHits, fDigits,
// fIdSens, fIdName, and fITSpoints.
//
// AddHit(Int_t track, Int_t *vol, Float_t *hits)
//     The function to add information to the AliITShit class. See the
// AliITShit class for a full description. This function allocates the
// necessary new space for the hit information and passes the variable
// track, and the pointers *vol and *hits to the AliITShit constructor
// function.
//
// AddDigit(Int_t *track, Int_t *digits)
//     The function to add information to the AliITSdigits class. See the
// AliITSdigits class for a full description. This function allocates the
// necessary new space for the digits information and passes the pointers
// *track and *digits to the AliITSdigits constructor function.
//
// BuildGeometry()
//     This function builds a simple ITS geometry used by the ROOT macro
// display.C. In general the geometry as coded is wrong.
//
// CreateGeometry()
//     This function builds the detailed geometry used by the Geant
// Monte Carlo. As defined here it is a dummy routine to be replaced
// by the version coded up in AliITSv? where the specific geometry to
// be used by the simulation is defined. See the definition of AliITSv5
// or the other routines for a complete definition.
//
// CreateMaterials()
//     This function defines the default materials used in the Geant
// Monte Carlo simulations. In general it is automatically replaced by
// the CreatMaterials routine defined in AliITSv?. Should the function
// CreateMaterials not exist for the geometry version you are using this
// one is used. See the definition found in AliITSv5 or the other routine
// for a complete definition.
//
// IsVersion()
//     Returns the version number of the AliITS class. At present it is
// version 1.
//
// DistancetoPrimitive(Int_t x, Int_t y)
//     A dummy routine used by the ROOT macro display.C to allow for the
// use of the mouse (pointing device) in the macro. In general this should
// never be called. If it is it returns the number 9999 for any value of
// x and y.
//
// Init()
//     This routine initializes the AliITS class. It is intended to be called
// from the Init function in AliITSv?. Besides displaying a banner
// indicating that it has been called it initializes the array fIdSens.
// Therefore it should be called after a call to CreateGeometry.
//
// MakeBranch(Option_t *Opt=" ")
//     Creates the TTree branch where the class AliITS is kept.
//
// SetEUCLID(bool_t euclid=1)
//     Sets the flag fEuclidOut to true (default) of false (euclid=0).
// By setting or clearing the fEuclidOut flag you can controls whether
// or not a euclid formatted output file of the ITS geometry is written.
// If fEuclidOut is set true then a file called ITSgeometry.euc will be
// written after the ITS geometry is defined in the Monte Carlo. If
// fEuclidOut is set false then no file is created.
//
// StepManager()
//     Dummy routine which is replaced by the routine StepManager() defined
// in AliITSv?. If no such routine exist then this routine returns zero.
// See AliITSv? for a detailed description of the step manager routines.
//
// GetITSgeom()
//     Returns the value of the pointer fITSgeom. This is used to get
// access to the ITS geometry stored in the file. See AliITSgeom for a
// full description of the geometry package.
//
// GetITSpoints()
//     Returns the value of the pointer fITSpoints. This is used to get
// access to the ITS cluster objects, if filled, stored in the file. See
// AliITSCluster for a full description of the cluster data.
////////////////////////////////////////////////////////////////////////
//_____________________________________________________________________________
AliITS::AliITS() {
  //
  // Default initialiser for ITS
  //     The default constructor of the AliITS class. In addition to
  // creating the AliITS class it zeros the variables fIshunt (a member
  // of AliDetector class), fEuclidOut, and fIdN, and zeros the pointers
  // fITSpoints, fIdSens, and fIdName.
  //
  fITSpoints  = 0;
  fIshunt     = 0;
  fEuclidOut  = 0;
  fIdN        = 0;
  fIdName     = 0;
  fIdSens     = 0;
  fITSmodules = 0;

}
//_____________________________________________________________________________
AliITS::AliITS(const char *name, const char *title):AliDetector(name,title){
  //
  // Default initialiser for ITS
  //     The constructor of the AliITS class. In addition to creating the
  // AliITS class, it allocates memory for the TClonesArrays fHits and
  // fDigits, and for the TObjArray fITSpoints. It also zeros the variables
  // fIshunt (a member of AliDetector class), fEuclidOut, and fIdN, and zeros
  // the pointers fIdSens and fIdName. To help in displaying hits via the ROOT
  // macro display.C AliITS also sets the marker color to red. The variables
  // passes with this constructor, const char *name and *title, are used by
  // the constructor of AliDetector class. See AliDetector class for a
  // description of these parameters and its constructor functions.
  //

  fHits       = new TClonesArray("AliITShit", 1560);
  gAlice->AddHitList(fHits);
  fDigits     = new TClonesArray("AliITSdigit",1000);
  fITSpoints  = new TObjArray();
  fITSmodules = 0; //new AliITSmodules();

  fIshunt     = 0;
  fEuclidOut  = 0;
  fIdN        = 0;
  fIdName     = 0;
  fIdSens     = 0;

  SetMarkerColor(kRed);

}

//_____________________________________________________________________________
AliITS::~AliITS(){
  //
  // Default distructor for ITS
  //     The default destructor of the AliITS class. In addition to deleting
  // the AliITS class it deletes the memory pointed to by the fHits, fDigits,
  // fIdSens, fIdName, and fITSpoints.
  //
  delete fHits;
  delete fDigits;
  if(fIdName!=0) delete[] fIdName;
  if(fIdSens!=0) delete[] fIdSens;
  delete fITSmodules;
  if(fITSpoints!=0) delete fITSpoints;
}

//_____________________________________________________________________________
void AliITS::AddDigit(Int_t *tracks, Int_t *digits){
  //
  // Add an ITS Digit
  //     The function to add information to the AliITSdigits class. See the
  // AliITSdigits class for a full description. This function allocates the
  // necessary new space for the digits information and passes the pointers
  // *track and *digits to the AliITSdigits constructor function.
  //
  TClonesArray &ldigits = *fDigits;
  new(ldigits[fNdigits++]) AliITSdigit(tracks,digits);
}

Int_t AliITS::AddDigit(AliITSdigit* d) {

   fDigits->Add(d);
   fNdigits = fDigits->GetEntriesFast();
   return fNdigits;
}

//_____________________________________________________________________________
void AliITS::AddHit(Int_t track, Int_t *vol, Float_t *hits){
  //
  // Add an ITS hit
  //     The function to add information to the AliITShit class. See the
  // AliITShit class for a full description. This function allocates the
  // necessary new space for the hit information and passes the variable
  // track, and the pointers *vol and *hits to the AliITShit constructor
  // function.
  //
  TClonesArray &lhits = *fHits;
  new(lhits[fNhits++]) AliITShit(fIshunt,track,vol,hits);
}
//_____________________________________________________________________________
Int_t AliITS::DistancetoPrimitive(Int_t , Int_t ){
  //
  // Distance from mouse to ITS on the screen. Dummy routine
  //     A dummy routine used by the ROOT macro display.C to allow for the
  // use of the mouse (pointing device) in the macro. In general this should
  // never be called. If it is it returns the number 9999 for any value of
  // x and y.
  //
  return 9999;
}

//_____________________________________________________________________________
void AliITS::Init(){
  //
  // Initialise ITS after it has been built
  //     This routine initializes the AliITS class. It is intended to be called
  // from the Init function in AliITSv?. Besides displaying a banner
  // indicating that it has been called it initializes the array fIdSens.
  // Therefore it should be called after a call to CreateGeometry.
  //
  Int_t i;
  //
  printf("\n");
  for(i=0;i<35;i++) printf("*");
  printf(" ITS_INIT ");
  for(i=0;i<35;i++) printf("*");
  printf("\n");
  //
  //
  for(i=0;i<fIdN;i++) fIdSens[i] = gMC->VolId(fIdName[i]);
  //
  for(i=0;i<80;i++) printf("*");
  printf("\n");
}

//_____________________________________________________________________________
void AliITS::MakeBranch(Option_t* option){
  //
  // Create Tree branches for the ITS.
  // Creates the TTree branch where the class AliITS is kept.
  //
  Int_t buffersize = 4000;
  char branchname[10];
  sprintf(branchname,"%s",GetName());

  AliDetector::MakeBranch(option);

  char *D = strstr(option,"D");

  if (fDigits   && gAlice->TreeD() && D) {
    gAlice->TreeD()->Branch(branchname,&fDigits, buffersize);
    printf("Making Branch %s for digits\n",branchname);
  } // end if
}

//____________________________________________________________________________
void AliITS::Streamer(TBuffer &R__b){
   // Stream an object of class AliITS.
    Int_t i,j,l;

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); 
      if (R__v == 1) {
	  AliDetector::Streamer(R__b);
	  R__b >> fITSgeom;
//        R__b >> fITSmodules; //We do not write out modules so don't read them
	  R__b >> fITSpoints;
	  R__b >> fEuclidOut;
	  R__b >> fIdN;
	  if(fIdSens!=0) delete[] fIdSens;
	  if(fIdName!=0) delete[] fIdName;
	  fIdSens = new Int_t[fIdN];
	  fIdName = new char*[fIdN];
	  for(i=0;i<fIdN;i++) R__b >> fIdSens[i];
	  for(i=0;i<fIdN;i++){
	      R__b >> l;
	      fIdName[i] = new char[l+1]; // add room for null character.
	      for(j=0;j<l;j++) R__b >> fIdName[i][j];
	      fIdName[i][l] = '\0'; // Null terminate this string.
	  } // end for i
	  R__b >> fMajorVersion;
	  R__b >> fMinorVersion;
      } // end if (R__v)
   } else {
      R__b.WriteVersion(AliITS::IsA());
      AliDetector::Streamer(R__b);
      R__b << fITSgeom;
//    R__b << fITSmodules; //We don't want to write out the modules class.
      R__b << fITSpoints;
      R__b << fEuclidOut;
      R__b << fIdN;
      for(i=0;i<fIdN;i++) R__b <<fIdSens[i];
      for(i=0;i<fIdN;i++){
	  l = strlen(fIdName[i]);
	  R__b << l;
	  for(j=0;j<l;j++) R__b << fIdName[i][j];
      } // end for i
      R__b << fMajorVersion;
      R__b << fMinorVersion;
   }
}
