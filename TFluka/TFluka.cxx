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
Revision 1.4  2002/11/04 16:00:46  iglez2
The conversion between ID and PDG now uses Fluka routines and arrays which is more consistent.

Revision 1.3  2002/10/22 15:12:14  alibrary
Introducing Riostream.h

Revision 1.2  2002/10/14 14:57:40  hristov
Merging the VirtualMC branch to the main development branch (HEAD)

Revision 1.1.2.8  2002/10/08 16:33:17  iglez2
LSOUIT is set to true before the second call to flukam.

Revision 1.1.2.7  2002/10/08 09:30:37  iglez2
Solved stupid missing ;

Revision 1.1.2.6  2002/10/07 13:40:22  iglez2
First implementations of the PDG <--> Fluka Id conversion routines

Revision 1.1.2.5  2002/09/26 16:26:03  iglez2
Added verbosity
Call to gAlice->Generator()->Generate()

Revision 1.1.2.4  2002/09/26 13:22:23  iglez2
Naive implementation of ProcessRun and ProcessEvent
Opening/Closing of input file (fInputFileName) with FORTRAN unit 5 before/after the first call to flukam inside Init()

Revision 1.1.2.3  2002/09/20 15:35:51  iglez2
Modification of LFDRTR. Value is passed to FLUKA !!!

Revision 1.1.2.2  2002/09/18 14:34:44  iglez2
Revised version with all pure virtual methods implemented

Revision 1.1.2.1  2002/07/24 08:49:41  alibrary
Adding TFluka to VirtualMC

Revision 1.1  2002/07/05 13:10:07  morsch
First commit of Fluka interface.

*/

#include <Riostream.h>

#include "TFluka.h"
#include "TCallf77.h"      //For the fortran calls
#include "Fdblprc.h"       //(DBLPRC) fluka common
#include "Fiounit.h"       //(IOUNIT) fluka common
#include "Fepisor.h"       //(EPISOR) fluka common
#include "Fpart.h"         //(PART)   fluka common
#include "TVirtualMC.h"

#include "TG4GeometryManager.h" //For the geometry management
#include "TG4DetConstruction.h" //For the detector construction

#include "FGeometryInit.hh"

// Fluka methods that may be needed.
#ifndef WIN32 
# define flukam  flukam_
# define fluka_openinp fluka_openinp_
# define fluka_closeinp fluka_closeinp_
# define mcihad mcihad_
# define mpdgha mpdgha_
#else 
# define flukam  FLUKAM
# define fluka_openinp FLUKA_OPENINP
# define fluka_closeinp FLUKA_CLOSEINP
# define mcihad MCIHAD
# define mpdgha MPDGHA
#endif

extern "C" 
{
  //
  // Prototypes for FLUKA functions
  //
  void type_of_call flukam(const int&);
  void type_of_call fluka_openinp(const int&, DEFCHARA);
  void type_of_call fluka_closeinp(const int&);
  int  type_of_call mcihad(const int&);
  int  type_of_call mpdgha(const int&);
}

//
// Class implementation for ROOT
//
ClassImp(TFluka)

//
//----------------------------------------------------------------------------
// TFluka constructors and destructors.
//____________________________________________________________________________ 
TFluka::TFluka()
  :TVirtualMC(),
   fVerbosityLevel(0),
   fInputFileName(""),
   fDetector(0)
{ 
  //
  // Default constructor
  //
} 
 
TFluka::TFluka(const char *title, Int_t verbosity)
  :TVirtualMC("TFluka",title),
   fVerbosityLevel(verbosity),
   fInputFileName(""),
   fDetector(0)
{
  if (fVerbosityLevel >=3)
    cout << "==> TFluka::TFluka(" << title << ") constructor called." << endl;

  
  // create geometry manager
  if (fVerbosityLevel >=2)
    cout << "\t* Creating G4 Geometry manager..." << endl;
  fGeometryManager = new TG4GeometryManager();
  if (fVerbosityLevel >=2)
    cout << "\t* Creating G4 Detector..." << endl;
  fDetector = new TG4DetConstruction();
  FGeometryInit* geominit = FGeometryInit::GetInstance();
  if (geominit)
    geominit->setDetConstruction(fDetector);
  else {
    cerr << "ERROR: Could not create FGeometryInit!" << endl;
    cerr << "       Exiting!!!" << endl;
    abort();
  }

  if (fVerbosityLevel >=3)
    cout << "<== TFluka::TFluka(" << title << ") constructor called." << endl;
}

TFluka::~TFluka() {
  if (fVerbosityLevel >=3)
    cout << "==> TFluka::~TFluka() destructor called." << endl;

  delete fGeometryManager;

  if (fVerbosityLevel >=3)
    cout << "<== TFluka::~TFluka() destructor called." << endl;
}

//
//_____________________________________________________________________________
// TFluka control methods
//____________________________________________________________________________ 
void TFluka::Init() {
  if (fVerbosityLevel >=3)
    cout << "==> TFluka::Init() called." << endl;

  if (fVerbosityLevel >=2)
    cout << "\t* Changing lfdrtr = (" << (GLOBAL.lfdrtr?'T':'F')
	 << ") in fluka..." << endl;
  GLOBAL.lfdrtr = true;

  if (fVerbosityLevel >=2)
    cout << "\t* Opening file " << fInputFileName << endl;
  const char* fname = fInputFileName;
  fluka_openinp(lunin, PASSCHARA(fname));

  if (fVerbosityLevel >=2)
    cout << "\t* Calling flukam..." << endl;
  flukam(1);

  if (fVerbosityLevel >=2)
    cout << "\t* Closing file " << fInputFileName << endl;
  fluka_closeinp(lunin);

  if (fVerbosityLevel >=3)
    cout << "<== TFluka::Init() called." << endl;
}

void TFluka::FinishGeometry() {
  if (fVerbosityLevel >=3)
    cout << "==> TFluka::FinishGeometry() called." << endl;

  fGeometryManager->Ggclos();

  if (fVerbosityLevel >=3)
    cout << "<== TFluka::FinishGeometry() called." << endl;
} 

void TFluka::BuildPhysics() {
  if (fVerbosityLevel >=3)
    cout << "==> TFluka::BuildPhysics() called." << endl;


  if (fVerbosityLevel >=3)
    cout << "<== TFluka::BuildPhysics() called." << endl;
}  

void TFluka::ProcessEvent() {
  if (fVerbosityLevel >=3)
    cout << "==> TFluka::ProcessEvent() called." << endl;

  if (fVerbosityLevel >=3)
    cout << "<== TFluka::ProcessEvent() called." << endl;
}


void TFluka::ProcessRun(Int_t nevent) {
  if (fVerbosityLevel >=3)
    cout << "==> TFluka::ProcessRun(" << nevent << ") called." 
	 << endl;

  if (fVerbosityLevel >=2) {
    cout << "\t* GLOBAL.fdrtr = " << (GLOBAL.lfdrtr?'T':'F') << endl;
    cout << "\t* Calling flukam again..." << endl;
  }
  fApplication->GeneratePrimaries();
  EPISOR.lsouit = true;
  flukam(1);

  if (fVerbosityLevel >=3)
    cout << "<== TFluka::ProcessRun(" << nevent << ") called." 
	 << endl;
}

//_____________________________________________________________________________
// methods for building/management of geometry
//____________________________________________________________________________ 
// functions from GCONS 
void TFluka::Gfmate(Int_t imat, char *name, Float_t &a, Float_t &z,  
		    Float_t &dens, Float_t &radl, Float_t &absl,
		    Float_t* ubuf, Int_t& nbuf) {
//
  fGeometryManager->Gfmate(imat, name, a, z, dens, radl, absl, ubuf, nbuf);
} 

void TFluka::Gfmate(Int_t imat, char *name, Double_t &a, Double_t &z,  
		    Double_t &dens, Double_t &radl, Double_t &absl,
		    Double_t* ubuf, Int_t& nbuf) {
//
  fGeometryManager->Gfmate(imat, name, a, z, dens, radl, absl, ubuf, nbuf);
} 

// detector composition
void TFluka::Material(Int_t& kmat, const char* name, Double_t a, 
		      Double_t z, Double_t dens, Double_t radl, Double_t absl,
		      Float_t* buf, Int_t nwbuf) {
//
  fGeometryManager
    ->Material(kmat, name, a, z, dens, radl, absl, buf, nwbuf); 
} 
void TFluka::Material(Int_t& kmat, const char* name, Double_t a, 
		      Double_t z, Double_t dens, Double_t radl, Double_t absl,
		      Double_t* buf, Int_t nwbuf) {
//
  fGeometryManager
    ->Material(kmat, name, a, z, dens, radl, absl, buf, nwbuf); 
} 

void TFluka::Mixture(Int_t& kmat, const char *name, Float_t *a, 
		     Float_t *z, Double_t dens, Int_t nlmat, Float_t *wmat) {
//
   fGeometryManager
     ->Mixture(kmat, name, a, z, dens, nlmat, wmat); 
} 
void TFluka::Mixture(Int_t& kmat, const char *name, Double_t *a, 
		     Double_t *z, Double_t dens, Int_t nlmat, Double_t *wmat) {
//
   fGeometryManager
     ->Mixture(kmat, name, a, z, dens, nlmat, wmat); 
} 

void TFluka::Medium(Int_t& kmed, const char *name, Int_t nmat, 
		    Int_t isvol, Int_t ifield, Double_t fieldm, Double_t tmaxfd, 
		    Double_t stemax, Double_t deemax, Double_t epsil, 
		    Double_t stmin, Float_t* ubuf, Int_t nbuf) { 
  //
  fGeometryManager
    ->Medium(kmed, name, nmat, isvol, ifield, fieldm, tmaxfd, stemax, deemax, 
	     epsil, stmin, ubuf, nbuf);
} 
void TFluka::Medium(Int_t& kmed, const char *name, Int_t nmat, 
		    Int_t isvol, Int_t ifield, Double_t fieldm, Double_t tmaxfd, 
		    Double_t stemax, Double_t deemax, Double_t epsil, 
		    Double_t stmin, Double_t* ubuf, Int_t nbuf) { 
  //
  fGeometryManager
    ->Medium(kmed, name, nmat, isvol, ifield, fieldm, tmaxfd, stemax, deemax, 
	     epsil, stmin, ubuf, nbuf);
} 

void TFluka::Matrix(Int_t& krot, Double_t thetaX, Double_t phiX, 
		    Double_t thetaY, Double_t phiY, Double_t thetaZ, 
		    Double_t phiZ) {
//		     
  fGeometryManager
    ->Matrix(krot, thetaX, phiX, thetaY, phiY, thetaZ, phiZ); 
} 

void TFluka::Gstpar(Int_t itmed, const char *param, Double_t parval) {
//
  fGeometryManager->Gstpar(itmed, param, parval); 
}    

// functions from GGEOM 
Int_t TFluka::Gsvolu(const char *name, const char *shape, Int_t nmed,  
		     Float_t *upar, Int_t np)  {
//
  return fGeometryManager->Gsvolu(name, shape, nmed, upar, np); 
}
Int_t TFluka::Gsvolu(const char *name, const char *shape, Int_t nmed,  
		     Double_t *upar, Int_t np)  {
//
  return fGeometryManager->Gsvolu(name, shape, nmed, upar, np); 
}
 
void TFluka::Gsdvn(const char *name, const char *mother, Int_t ndiv, 
		   Int_t iaxis) {
//
  fGeometryManager->Gsdvn(name, mother, ndiv, iaxis); 
} 

void TFluka::Gsdvn2(const char *name, const char *mother, Int_t ndiv, 
		    Int_t iaxis, Double_t c0i, Int_t numed) {
//
  fGeometryManager->Gsdvn2(name, mother, ndiv, iaxis, c0i, numed); 
} 

void TFluka::Gsdvt(const char *name, const char *mother, Double_t step, 
		   Int_t iaxis, Int_t numed, Int_t ndvmx) {
//			
  fGeometryManager->Gsdvt(name, mother, step, iaxis, numed, ndvmx); 
} 

void TFluka::Gsdvt2(const char *name, const char *mother, Double_t step, 
		    Int_t iaxis, Double_t c0, Int_t numed, Int_t ndvmx) { 
//
  fGeometryManager->Gsdvt2(name, mother, step, iaxis, c0, numed, ndvmx); 
} 

void TFluka::Gsord(const char *name, Int_t iax) {
//
  fGeometryManager->Gsord(name, iax); 
} 

void TFluka::Gspos(const char *name, Int_t nr, const char *mother,  
		   Double_t x, Double_t y, Double_t z, Int_t irot, 
		   const char *konly) {
//
  fGeometryManager->Gspos(name, nr, mother, x, y, z, irot, konly); 
} 

void TFluka::Gsposp(const char *name, Int_t nr, const char *mother,  
		    Double_t x, Double_t y, Double_t z, Int_t irot,
		    const char *konly, Float_t *upar, Int_t np)  {
  //
  fGeometryManager->Gsposp(name, nr, mother, x, y, z, irot, konly, upar, np); 
} 
void TFluka::Gsposp(const char *name, Int_t nr, const char *mother,  
		    Double_t x, Double_t y, Double_t z, Int_t irot,
		    const char *konly, Double_t *upar, Int_t np)  {
  //
  fGeometryManager->Gsposp(name, nr, mother, x, y, z, irot, konly, upar, np); 
} 

void TFluka::Gsbool(const char* onlyVolName, const char* manyVolName) {
//
  fGeometryManager->Gsbool(onlyVolName, manyVolName);
}

void TFluka::SetCerenkov(Int_t itmed, Int_t npckov, Float_t *ppckov,
			 Float_t *absco, Float_t *effic, Float_t *rindex) {
//
  fGeometryManager->SetCerenkov(itmed, npckov, ppckov, absco, effic, rindex);
}  
void TFluka::SetCerenkov(Int_t itmed, Int_t npckov, Double_t *ppckov,
			 Double_t *absco, Double_t *effic, Double_t *rindex) {
//
  fGeometryManager->SetCerenkov(itmed, npckov, ppckov, absco, effic, rindex);
}  
    
// Euclid
void TFluka::WriteEuclid(const char* fileName, const char* topVol, 
                          Int_t number, Int_t nlevel) {
//
  fGeometryManager->WriteEuclid(fileName, topVol, number, nlevel); 
} 






//____________________________________________________________________________ 
// ID <--> PDG transformations
//_____________________________________________________________________________
Int_t TFluka::IdFromPDG(Int_t pdg) const 
{
  //
  // Return Fluka code from PDG and pseudo ENDF code

  // MCIHAD() goes from pdg to fluka internal.
  Int_t intfluka = mcihad(pdg);
  // KPTOIP array goes from internal to official
  return GetFlukaKPTOIP(intfluka);
}

Int_t TFluka::PDGFromId(Int_t id) const 
{
  //
  // Return PDG code and pseudo ENDF code from Fluka code

  //IPTOKP array goes from official to internal
  Int_t intfluka = GetFlukaIPTOKP(id);
  //MPKDHA() goes from internal to PDG
  return mpdgha(intfluka);
  
}
