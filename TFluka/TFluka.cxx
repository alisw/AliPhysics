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
Revision 1.6  2002/11/21 18:40:06  iglez2
Media handling added

Revision 1.5  2002/11/07 17:59:10  iglez2
Included the geometry through geant4_vmc/FLUGG

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
#include "Fepisor.h"       //(EPISOR) fluka common
#include "Ffinuc.h"        //(FINUC) fluka common
#include "Fiounit.h"       //(IOUNIT) fluka common
#include "Fpaprop.h"       //(PAPROP) fluka common
#include "Fpart.h"         //(PART)   fluka common
#include "Ftrackr.h"       //(TRACKR) fluka common
#include "Ffheavy.h"       //(FHEAVY) fluka common

#include "TVirtualMC.h"
#include "TG4GeometryManager.h" //For the geometry management
#include "TG4DetConstruction.h" //For the detector construction

#include "FGeometryInit.hh"
#include "TLorentzVector.h"

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
   fDetector(0),
   fCurrentFlukaRegion(-1)
{ 
  //
  // Default constructor
  //
} 
 
TFluka::TFluka(const char *title, Int_t verbosity)
  :TVirtualMC("TFluka",title),
   fVerbosityLevel(verbosity),
   fInputFileName(""),
   fDetector(0),
   fCurrentFlukaRegion(-1)
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

  FGeometryInit* flugg = FGeometryInit::GetInstance();
  map<TString, Int_t, less<TString> >::iterator i;
  for (fVolumeMediaMap.begin(); i != fVolumeMediaMap.end(); i++) {
    TString volName = (*i).first;
    Int_t   media   = (*i).second;
    Int_t   region  = flugg->GetRegionFromName(volName);
    fMediaByRegion[region] = media;
  }
  
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
  fVolumeMediaMap[TString(name)] = nmed;
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



//_____________________________________________________________________________
// methods needed by the stepping
//____________________________________________________________________________ 
Int_t TFluka::GetMedium() const {
  return fMediaByRegion[fCurrentFlukaRegion];
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



//_____________________________________________________________________________
// methods for step management
//____________________________________________________________________________ 
//
// dynamic properties
//
void TFluka::TrackPosition(TLorentzVector& position) const
{
// Return the current position in the master reference frame of the
// track being transported
// TRACKR.atrack = age of the particle
// TRACKR.xtrack = x-position of the last point
// TRACKR.ytrack = y-position of the last point
// TRACKR.ztrack = z-position of the last point
  position.SetX(TRACKR.xtrack[TRACKR.ntrack]);
  position.SetY(TRACKR.ytrack[TRACKR.ntrack]);
  position.SetZ(TRACKR.ztrack[TRACKR.ntrack]);
  position.SetT(TRACKR.atrack);
}

void TFluka::TrackMomentum(TLorentzVector& momentum) const
{
// Return the direction and the momentum (GeV/c) of the track
// currently being transported
// TRACKR.ptrack = momentum of the particle (not always defined, if
//               < 0 must be obtained from etrack) 
// TRACKR.cx,y,ztrck = direction cosines of the current particle
// TRACKR.etrack = total energy of the particle
// TRACKR.jtrack = identity number of the particle
// PAPROP.am[TRACKR.jtrack] = particle mass in gev
  if (TRACKR.ptrack >= 0) {
    momentum.SetPx(TRACKR.ptrack*TRACKR.cxtrck);
    momentum.SetPy(TRACKR.ptrack*TRACKR.cytrck);
    momentum.SetPz(TRACKR.ptrack*TRACKR.cztrck);
    momentum.SetE(TRACKR.etrack);
    return;
  }
  else {
    Double_t p = sqrt(TRACKR.etrack*TRACKR.etrack - PAPROP.am[TRACKR.jtrack]*PAPROP.am[TRACKR.jtrack]);
    momentum.SetPx(p*TRACKR.cxtrck);
    momentum.SetPy(p*TRACKR.cytrck);
    momentum.SetPz(p*TRACKR.cztrck);
    momentum.SetE(TRACKR.etrack);
    return;
  }
}

Double_t TFluka::TrackStep() const
{
// Return the length in centimeters of the current step
// TRACKR.ctrack = total curved path
    return TRACKR.ctrack;
}

Double_t TFluka::TrackLength() const
{
// It is wrong
// should be the sum of all steps starting from the beginning of the track
// for the time being returns only the length in centimeters of the current step
    return TRACKR.ctrack;
}

Double_t TFluka::TrackTime() const
{
// Return the current time of flight of the track being transported
// TRACKR.atrack = age of the particle
  return TRACKR.atrack;
}

Double_t TFluka::Edep() const
{
// Energy deposition
// if TRACKR.ntrack = 0, TRACKR.mtrack = 0:
// -->local energy deposition (the value and the point are not recorded in TRACKR)
//    but in the variable "rull" of the procedure "endraw.cxx"
// if TRACKR.ntrack > 0, TRACKR.mtrack = 0:
// -->no energy loss along the track
// if TRACKR.ntrack > 0, TRACKR.mtrack > 0:
// -->energy loss distributed along the track
// TRACKR.dtrack = energy deposition of the jth deposition even
  if (TRACKR.ntrack == 0 && TRACKR.mtrack == 0)
    return fRull;
  else {
    Double_t sum = 0;
    for ( Int_t j=0;j<TRACKR.mtrack;j++) {
      sum +=TRACKR.dtrack[j];  
    }
    return sum;
  }
}

Int_t TFluka::TrackPid() const
{
// Return the id of the particle transported
// TRACKR.jtrack = identity number of the particle
  return PDGFromId(TRACKR.jtrack);
}

Double_t TFluka::TrackCharge() const
{
// Return charge of the track currently transported
// PAPROP.ichrge = electric charge of the particle
  return PAPROP.ichrge[TRACKR.jtrack+6];
}

Double_t TFluka::TrackMass() const
{
// PAPROP.am = particle mass in GeV
  return PAPROP.am[TRACKR.jtrack+6];
}

Double_t TFluka::Etot() const
{
// TRACKR.etrack = total energy of the particle
  return TRACKR.etrack;
}

//
// track status
//
Bool_t   TFluka::IsNewTrack() const
{
// ???????????????,
// True if the track is not at the boundary of the current volume
// Not true in some cases in bxdraw - to be solved
  return 1;
}

Bool_t   TFluka::IsTrackInside() const
{
// True if the track is not at the boundary of the current volume
// In Fluka a step is always inside one kind of material
// If the step would go behind the region of one material,
// it will be shortened to reach only the boundary.
// Therefore IsTrackInside() is always true.
// Not true in some cases in bxdraw - to be solved
  return 1;
}

Bool_t   TFluka::IsTrackEntering() const
{
// True if this is the first step of the track in the current volume
// Boundary- (X) crossing
// Icode = 19: boundary crossing - call from Kaskad
// Icode = 29: boundary crossing - call from Emfsco
// Icode = 39: boundary crossing - call from Kasneu
// Icode = 49: boundary crossing - call from Kashea
// Icode = 59: boundary crossing - call from Kasoph
  if (fIcode == 19 ||
      fIcode == 29 ||
      fIcode == 39 ||
      fIcode == 49 ||
      fIcode == 59) return 1;
  else return 0;
}

Bool_t   TFluka::IsTrackExiting() const
{
// True if this is the last step of the track in the current volume
// Boundary- (X) crossing
// Icode = 19: boundary crossing - call from Kaskad
// Icode = 29: boundary crossing - call from Emfsco
// Icode = 39: boundary crossing - call from Kasneu
// Icode = 49: boundary crossing - call from Kashea
// Icode = 59: boundary crossing - call from Kasoph
  if (fIcode == 19 ||
      fIcode == 29 ||
      fIcode == 39 ||
      fIcode == 49 ||
      fIcode == 59) return 1;
  else return 0;
}

Bool_t   TFluka::IsTrackOut() const
{
// True if the track is out of the setup
// means escape
// Icode = 14: escape - call from Kaskad
// Icode = 23: escape - call from Emfsco
// Icode = 32: escape - call from Kasneu
// Icode = 40: escape - call from Kashea
// Icode = 51: escape - call from Kasoph
  if (fIcode == 14 ||
      fIcode == 23 ||
      fIcode == 32 ||
      fIcode == 40 ||
      fIcode == 51) return 1;
  else return 0;
}

Bool_t   TFluka::IsTrackDisappeared() const
{
// means all inelastic interactions and decays
// fIcode from usdraw
  if (fIcode == 101 || // inelastic interaction
      fIcode == 102 || // particle decay
      fIcode == 214 || // in-flight annihilation
      fIcode == 215 || // annihilation at rest
      fIcode == 217 || // pair production
      fIcode == 221) return 1;
  else return 0;
}

Bool_t   TFluka::IsTrackStop() const
{
// True if the track energy has fallen below the threshold
// means stopped by signal or below energy threshold
// Icode = 12: stopping particle       - call from Kaskad
// Icode = 15: time kill               - call from Kaskad
// Icode = 21: below threshold, iarg=1 - call from Emfsco
// Icode = 22: below threshold, iarg=2 - call from Emfsco
// Icode = 24: time kill               - call from Emfsco
// Icode = 31: below threshold         - call from Kasneu
// Icode = 33: time kill               - call from Kasneu
// Icode = 41: time kill               - call from Kashea
// Icode = 52: time kill               - call from Kasoph
  if (fIcode == 12 ||
      fIcode == 15 ||
      fIcode == 21 ||
      fIcode == 22 ||
      fIcode == 24 ||
      fIcode == 31 ||
      fIcode == 33 ||
      fIcode == 41 ||
      fIcode == 52) return 1;
  else return 0;
}

Bool_t   TFluka::IsTrackAlive() const
{
// means not disappeared or not out
  if (IsTrackDisappeared() || IsTrackOut() ) return 0;
  else return 1;
}

//
// secondaries
//

Int_t TFluka::NSecondaries() const
// Number of secondary particles generated in the current step
// fIcode = 100 = elastic interaction
// fIcode = 101 = inelastic interaction
// fIcode = 102 = particle decay
// fIcode = 103 = delta ray
// fIcode = 104 = pair production
// fIcode = 105 = bremsstrahlung
{
  if (fIcode >= 100 && fIcode <= 105)
    return FINUC.np + FHEAVY.npheav;
  else 
    return -1;
}

void     TFluka::GetSecondary(Int_t isec, Int_t& particleId,
		TLorentzVector& position, TLorentzVector& momentum)
// 
// fIcode = 100 = elastic interaction
// fIcode = 101 = inelastic interaction
// fIcode = 102 = particle decay
// fIcode = 103 = delta ray
// fIcode = 104 = pair production
// fIcode = 105 = bremsstrahlung
{
  if (fIcode >= 100 && fIcode <= 105) {
    if (isec >= 0 && isec < FINUC.np) {
      particleId = PDGFromId(FINUC.kpart[isec]);
      position.SetX(fXsco);
      position.SetY(fYsco);
      position.SetZ(fZsco);
      position.SetT(FINUC.agesec[isec]);
      momentum.SetPx(FINUC.plr[isec]*FINUC.cxr[isec]);
      momentum.SetPy(FINUC.plr[isec]*FINUC.cyr[isec]);
      momentum.SetPz(FINUC.plr[isec]*FINUC.czr[isec]);
      momentum.SetE(FINUC.tki[isec] + PAPROP.am[FINUC.kpart[isec]+6]);
    }
    if (isec >= FINUC.np && isec < FINUC.np + FHEAVY.npheav) {
      Int_t jsec = isec - FINUC.np;
      particleId = FHEAVY.kheavy[jsec]; // this is Fluka id !!!
      position.SetX(fXsco);
      position.SetY(fYsco);
      position.SetZ(fZsco);
      position.SetT(FHEAVY.agheav[jsec]);
      momentum.SetPx(FHEAVY.pheavy[jsec]*FHEAVY.cxheav[jsec]);
      momentum.SetPy(FHEAVY.pheavy[jsec]*FHEAVY.cyheav[jsec]);
      momentum.SetPz(FHEAVY.pheavy[jsec]*FHEAVY.czheav[jsec]);
      if (FHEAVY.tkheav[jsec] >= 3 && FHEAVY.tkheav[jsec] <= 6) 
        momentum.SetE(FHEAVY.tkheav[jsec] + PAPROP.am[jsec+6]);
      else if (FHEAVY.tkheav[jsec] > 6)
        momentum.SetE(FHEAVY.tkheav[jsec] + FHEAVY.amnhea[jsec]); // to be checked !!!
    }
  }
}

//TMCProcess ProdProcess(Int_t isec) const
// Name of the process that has produced the secondary particles
// in the current step
//{
// will come from FINUC when called from USDRAW
//}

//Int_t StepProcesses(TArrayI &proc) const
// Return processes active in the current step
//{
//ck = total energy of the particl ???????????????? 
//}


// ===============================================================
void TFluka::FutoTest() 
{
  Int_t icode, mreg, newreg, particleId;
//  Int_t medium;
  Double_t rull, xsco, ysco, zsco;
  TLorentzVector position, momentum;
  icode = GetIcode();
  if (icode == 0) {
    cout << " icode=" << icode << endl;
    /*
    cout << "TLorentzVector positionX=" << position.X()
       << "positionY=" << position.Y()
       << "positionZ=" << position.Z()
       << "timeT=" << position.T() << endl;
    cout << "TLorentzVector momentumX=" << momentum.X()
       << "momentumY=" << momentum.Y()
       << "momentumZ=" << momentum.Z()
       << "energyE=" << momentum.E() << endl;
    cout << "TrackPid=" << TrackPid() << endl;
    */
  }

  else if (icode > 0 && icode <= 5) {
// mgdraw
    mreg = GetMreg();
//    medium = GetMedium();
    cout << " icode=" << icode
	 << " mreg=" << mreg
//	 << " medium=" << medium
	 << endl;
  TrackPosition(position);
  TrackMomentum(momentum);
  cout << "TLorentzVector positionX=" << position.X()
       << "positionY=" << position.Y()
       << "positionZ=" << position.Z()
       << "timeT=" << position.T() << endl;
  cout << "TLorentzVector momentumX=" << momentum.X()
       << "momentumY=" << momentum.Y()
       << "momentumZ=" << momentum.Z()
       << "energyE=" << momentum.E() << endl;
  cout << "TrackStep=" << TrackStep() << endl;
  cout << "TrackLength=" << TrackLength() << endl;
  cout << "TrackTime=" << TrackTime() << endl;
  cout << "Edep=" << Edep() << endl;
  cout << "TrackPid=" << TrackPid() << endl;
  cout << "TrackCharge=" << TrackCharge() << endl;
  cout << "TrackMass=" << TrackMass() << endl;
  cout << "Etot=" << Etot() << endl;
  cout << "IsNewTrack=" << IsNewTrack() << endl;
  cout << "IsTrackInside=" << IsTrackInside() << endl;
  cout << "IsTrackEntering=" << IsTrackEntering() << endl;
  cout << "IsTrackExiting=" << IsTrackExiting() << endl;
  cout << "IsTrackOut=" << IsTrackOut() << endl;
  cout << "IsTrackDisappeared=" << IsTrackDisappeared() << endl;
  cout << "IsTrackAlive=" << IsTrackAlive() << endl;
  }

  else if((icode >= 10 && icode <= 15) ||
          (icode >= 20 && icode <= 24) ||
          (icode >= 30 && icode <= 33) ||
          (icode >= 40 && icode <= 41) ||
          (icode >= 50 && icode <= 52)) {
// endraw
    mreg = GetMreg();
//    medium = GetMedium();
    rull = GetRull();
    xsco = GetXsco();
    ysco = GetYsco();
    zsco = GetZsco();
    cout << " icode=" << icode
         << " mreg=" << mreg
//	 << " medium=" << medium
	 << " rull=" << rull
	 << " xsco=" << xsco
	 << " ysco=" << ysco
	 << " zsco=" << zsco << endl;
  TrackPosition(position);
  TrackMomentum(momentum);
  cout << "Edep=" << Edep() << endl;
  cout << "Etot=" << Etot() << endl;
  cout << "TrackPid=" << TrackPid() << endl;
  cout << "TrackCharge=" << TrackCharge() << endl;
  cout << "TrackMass=" << TrackMass() << endl;
  cout << "IsTrackOut=" << IsTrackOut() << endl;
  cout << "IsTrackDisappeared=" << IsTrackDisappeared() << endl;
  cout << "IsTrackStop=" << IsTrackStop() << endl;
  cout << "IsTrackAlive=" << IsTrackAlive() << endl;
  }

  else if((icode >= 100 && icode <= 105) ||
           (icode == 208) ||
           (icode == 210) ||
           (icode == 212) ||
           (icode >= 214 && icode <= 215) ||
           (icode == 217) ||
           (icode == 219) ||
           (icode == 221) ||
           (icode == 225) ||
           (icode == 300) ||
           (icode == 400)) {
// usdraw
    mreg = GetMreg();
//    medium = GetMedium();
    xsco = GetXsco();
    ysco = GetYsco();
    zsco = GetZsco();
    cout << " icode=" << icode
         << " mreg=" << mreg
//	 << " medium=" << medium
	 << " xsco=" << xsco
	 << " ysco=" << ysco
	 << " zsco=" << zsco << endl;
    cout << "TrackPid=" << TrackPid() << endl;
    cout << "NSecondaries=" << NSecondaries() << endl;
    for (Int_t isec=0; isec< NSecondaries(); isec++) {
//void     TFluka::GetSecondary(Int_t isec, Int_t& particleId,
//                 TLorentzVector& position, TLorentzVector& momentum)
      TFluka::GetSecondary(isec, particleId, position, momentum);
      cout << "TLorentzVector positionX=" << position.X()
           << "positionY=" << position.Y()
           << "positionZ=" << position.Z()
           << "timeT=" << position.T() << endl;
      cout << "TLorentzVector momentumX=" << momentum.X()
           << "momentumY=" << momentum.Y()
           << "momentumZ=" << momentum.Z()
           << "energyE=" << momentum.E() << endl;
      cout << "TrackPid=" << particleId << endl;

    }
  }

  else if((icode == 19) ||
          (icode == 29) ||
          (icode == 39) ||
          (icode == 49) ||
          (icode == 59)) {
    mreg = GetMreg();
//    medium = GetMedium();
    newreg = GetNewreg();
    xsco = GetXsco();
    ysco = GetYsco();
    zsco = GetZsco();
    cout << " icode=" << icode
         << " mreg=" << mreg
//	 << " medium=" << medium
	 << " newreg=" << newreg
	 << " xsco=" << xsco
	 << " ysco=" << ysco
	 << " zsco=" << zsco << endl;
  }
//
// ====================================================================
//

  

} // end of FutoTest

