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

/* $Id$ */

#include <Riostream.h>

#include "TClonesArray.h"
#include "TFlukaGeo.h"
#include "TCallf77.h"      //For the fortran calls
#include "Fdblprc.h"       //(DBLPRC) fluka common
#include "Fepisor.h"       //(EPISOR) fluka common
#include "Ffinuc.h"        //(FINUC) fluka common
#include "Fiounit.h"       //(IOUNIT) fluka common
#include "Fpaprop.h"       //(PAPROP) fluka common
#include "Fpart.h"         //(PART)   fluka common
#include "Ftrackr.h"       //(TRACKR) fluka common
#include "Fpaprop.h"       //(PAPROP) fluka common
#include "Ffheavy.h"       //(FHEAVY) fluka common

#include "TVirtualMC.h"
#include "TGeoManager.h"
#include "TFlukaMCGeometry.h"

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
//______________________________________________________________________________
TFluka::TFluka()
  :TVirtualMC(),
   fVerbosityLevel(0),
   sInputFileName("")
{ 
  //
  // Default constructor
  //
   fNVolumes = 0;
   fCurrentFlukaRegion = -1;
   fGeom = 0;
} 
 
//______________________________________________________________________________ 
TFluka::TFluka(const char *title, Int_t verbosity, Bool_t isRootGeometrySupported)
  :TVirtualMC("TFluka",title, isRootGeometrySupported),
   fVerbosityLevel(verbosity),
   sInputFileName(""),
   fTrackIsEntering(0),
   fTrackIsExiting(0)
{
  // create geometry interface
  if (fVerbosityLevel >=3)
    cout << "<== TFluka::TFluka(" << title << ") constructor called." << endl;

   fNVolumes      = 0;
   fCurrentFlukaRegion = -1;
   fGeom = new TFlukaMCGeometry("geom", "ALICE geometry");
}

//______________________________________________________________________________ 
TFluka::~TFluka() {
  if (fVerbosityLevel >=3)
    cout << "==> TFluka::~TFluka() destructor called." << endl;

  delete fGeom;

  if (fVerbosityLevel >=3)
    cout << "<== TFluka::~TFluka() destructor called." << endl;
}

//
//______________________________________________________________________________
// TFluka control methods
//______________________________________________________________________________ 
void TFluka::Init() {

  if (fVerbosityLevel >=3)
    cout << "==> TFluka::Init() called." << endl;

   if (!gGeoManager) new TGeoManager("geom", "FLUKA geometry");
   fApplication->ConstructGeometry();
   TGeoVolume *top = (TGeoVolume*)gGeoManager->GetListOfVolumes()->First();
   gGeoManager->SetTopVolume(top);
   gGeoManager->CloseGeometry("di");
   gGeoManager->DefaultColors();  // to be removed
   fNVolumes = fGeom->NofVolumes();
   printf("== Number of volumes: %i\n ==", fNVolumes);
   fGeom->CreateFlukaMatFile("flukaMat.inp");   
  cout << "\t* InitPhysics() - Prepare input file to be called" << endl; 
  // now we have TGeo geometry created and we have to patch alice.inp
  // with the material mapping file FlukaMat.inp
  InitPhysics(); // prepare input file with the current physics settings
    cout << "\t* InitPhysics() - Prepare input file was called" << endl; 

  if (fVerbosityLevel >=2)
    cout << "\t* Changing lfdrtr = (" << (GLOBAL.lfdrtr?'T':'F')
	 << ") in fluka..." << endl;
  GLOBAL.lfdrtr = true;

  if (fVerbosityLevel >=2)
    cout << "\t* Opening file " << sInputFileName << endl;
  const char* fname = sInputFileName;
  fluka_openinp(lunin, PASSCHARA(fname));

  if (fVerbosityLevel >=2)
    cout << "\t* Calling flukam..." << endl;
  flukam(1);

  if (fVerbosityLevel >=2)
    cout << "\t* Closing file " << sInputFileName << endl;
  fluka_closeinp(lunin);

  FinishGeometry();

  if (fVerbosityLevel >=3)
    cout << "<== TFluka::Init() called." << endl;
}

//______________________________________________________________________________ 
void TFluka::FinishGeometry() {
//
// Build-up table with region to medium correspondance
//
  if (fVerbosityLevel >=3)
    cout << "==> TFluka::FinishGeometry() called." << endl;

   printf("----FinishGeometry - nothing to do with TGeo\n");
  
  if (fVerbosityLevel >=3)
    cout << "<== TFluka::FinishGeometry() called." << endl;
} 

//______________________________________________________________________________ 
void TFluka::BuildPhysics() {
  if (fVerbosityLevel >=3)
    cout << "==> TFluka::BuildPhysics() called." << endl;


  if (fVerbosityLevel >=3)
    cout << "<== TFluka::BuildPhysics() called." << endl;
}  

//______________________________________________________________________________ 
void TFluka::ProcessEvent() {
  if (fVerbosityLevel >=3)
    cout << "==> TFluka::ProcessEvent() called." << endl;
  fApplication->GeneratePrimaries();
  EPISOR.lsouit = true;
  flukam(1);
  if (fVerbosityLevel >=3)
    cout << "<== TFluka::ProcessEvent() called." << endl;
}

//______________________________________________________________________________ 
void TFluka::ProcessRun(Int_t nevent) {
  if (fVerbosityLevel >=3)
    cout << "==> TFluka::ProcessRun(" << nevent << ") called." 
	 << endl;

  if (fVerbosityLevel >=2) {
    cout << "\t* GLOBAL.fdrtr = " << (GLOBAL.lfdrtr?'T':'F') << endl;
    cout << "\t* Calling flukam again..." << endl;
  }
  fApplication->InitGeometry();
  fApplication->BeginEvent();
  ProcessEvent();
  fApplication->FinishEvent();
  if (fVerbosityLevel >=3)
    cout << "<== TFluka::ProcessRun(" << nevent << ") called." 
	 << endl;

}

//_____________________________________________________________________________
// methods for building/management of geometry

// functions from GCONS 
//____________________________________________________________________________ 
void TFluka::Gfmate(Int_t imat, char *name, Float_t &a, Float_t &z,  
		    Float_t &dens, Float_t &radl, Float_t &absl,
		    Float_t* ubuf, Int_t& nbuf) {
//
  fGeom->Gfmate(imat, name, a, z, dens, radl, absl, ubuf, nbuf);
} 

//______________________________________________________________________________ 
void TFluka::Gfmate(Int_t imat, char *name, Double_t &a, Double_t &z,  
		    Double_t &dens, Double_t &radl, Double_t &absl,
		    Double_t* ubuf, Int_t& nbuf) {
//
  fGeom->Gfmate(imat, name, a, z, dens, radl, absl, ubuf, nbuf);
} 

// detector composition
//______________________________________________________________________________ 
void TFluka::Material(Int_t& kmat, const char* name, Double_t a, 
		      Double_t z, Double_t dens, Double_t radl, Double_t absl,
		      Float_t* buf, Int_t nwbuf) {
//
   fGeom->Material(kmat, name, a, z, dens, radl, absl, buf, nwbuf);
} 

//______________________________________________________________________________ 
void TFluka::Material(Int_t& kmat, const char* name, Double_t a, 
		      Double_t z, Double_t dens, Double_t radl, Double_t absl,
		      Double_t* buf, Int_t nwbuf) {
//
   fGeom->Material(kmat, name, a, z, dens, radl, absl, buf, nwbuf);
} 

//______________________________________________________________________________ 
void TFluka::Mixture(Int_t& kmat, const char *name, Float_t *a, 
		     Float_t *z, Double_t dens, Int_t nlmat, Float_t *wmat) {
//
   fGeom->Mixture(kmat, name, a, z, dens, nlmat, wmat); 
} 

//______________________________________________________________________________ 
void TFluka::Mixture(Int_t& kmat, const char *name, Double_t *a, 
		     Double_t *z, Double_t dens, Int_t nlmat, Double_t *wmat) {
//
   fGeom->Mixture(kmat, name, a, z, dens, nlmat, wmat); 
} 

//______________________________________________________________________________ 
void TFluka::Medium(Int_t& kmed, const char *name, Int_t nmat, 
		    Int_t isvol, Int_t ifield, Double_t fieldm, Double_t tmaxfd, 
		    Double_t stemax, Double_t deemax, Double_t epsil, 
		    Double_t stmin, Float_t* ubuf, Int_t nbuf) { 
  //
  fGeom->Medium(kmed, name, nmat, isvol, ifield, fieldm, tmaxfd, stemax, deemax, 
	     epsil, stmin, ubuf, nbuf);
} 

//______________________________________________________________________________ 
void TFluka::Medium(Int_t& kmed, const char *name, Int_t nmat, 
		    Int_t isvol, Int_t ifield, Double_t fieldm, Double_t tmaxfd, 
		    Double_t stemax, Double_t deemax, Double_t epsil, 
		    Double_t stmin, Double_t* ubuf, Int_t nbuf) { 
  //
  fGeom->Medium(kmed, name, nmat, isvol, ifield, fieldm, tmaxfd, stemax, deemax, 
	     epsil, stmin, ubuf, nbuf);
} 

//______________________________________________________________________________ 
void TFluka::Matrix(Int_t& krot, Double_t thetaX, Double_t phiX, 
		    Double_t thetaY, Double_t phiY, Double_t thetaZ, 
		    Double_t phiZ) {
//		     
  fGeom->Matrix(krot, thetaX, phiX, thetaY, phiY, thetaZ, phiZ); 
} 

//______________________________________________________________________________ 
void TFluka::Gstpar(Int_t /*itmed*/, const char */*param*/, Double_t /*parval*/) {
//
// Is it needed with TGeo ??? - to clear-up
   Warning("Gstpar", "Not implemented with TGeo");
}    

// functions from GGEOM 
//_____________________________________________________________________________
void TFluka::Gsatt(const char *name, const char *att, Int_t val)
{ 
   fGeom->Gsatt(name,att, val);
}

//______________________________________________________________________________ 
Int_t TFluka::Gsvolu(const char *name, const char *shape, Int_t nmed,  
		     Float_t *upar, Int_t np)  {
//
    return fGeom->Gsvolu(name, shape, nmed, upar, np); 
}

//______________________________________________________________________________ 
Int_t TFluka::Gsvolu(const char *name, const char *shape, Int_t nmed,  
		     Double_t *upar, Int_t np)  {
//
    return fGeom->Gsvolu(name, shape, nmed, upar, np); 
}
 
//______________________________________________________________________________ 
void TFluka::Gsdvn(const char *name, const char *mother, Int_t ndiv, 
		   Int_t iaxis) {
//
    fGeom->Gsdvn(name, mother, ndiv, iaxis); 
} 

//______________________________________________________________________________ 
void TFluka::Gsdvn2(const char *name, const char *mother, Int_t ndiv, 
		    Int_t iaxis, Double_t c0i, Int_t numed) {
//
    fGeom->Gsdvn2(name, mother, ndiv, iaxis, c0i, numed); 
} 

//______________________________________________________________________________ 
void TFluka::Gsdvt(const char *name, const char *mother, Double_t step, 
		   Int_t iaxis, Int_t numed, Int_t ndvmx) {
//	
    fGeom->Gsdvt(name, mother, step, iaxis, numed, ndvmx); 
} 

//______________________________________________________________________________ 
void TFluka::Gsdvt2(const char *name, const char *mother, Double_t step, 
		    Int_t iaxis, Double_t c0, Int_t numed, Int_t ndvmx) { 
//
    fGeom->Gsdvt2(name, mother, step, iaxis, c0, numed, ndvmx); 
} 

//______________________________________________________________________________ 
void TFluka::Gsord(const char * /*name*/, Int_t /*iax*/) {
//
// Nothing to do with TGeo
} 

//______________________________________________________________________________ 
void TFluka::Gspos(const char *name, Int_t nr, const char *mother,  
		   Double_t x, Double_t y, Double_t z, Int_t irot, 
		   const char *konly) {
//
  fGeom->Gspos(name, nr, mother, x, y, z, irot, konly); 
} 

//______________________________________________________________________________ 
void TFluka::Gsposp(const char *name, Int_t nr, const char *mother,  
		    Double_t x, Double_t y, Double_t z, Int_t irot,
		    const char *konly, Float_t *upar, Int_t np)  {
  //
  fGeom->Gsposp(name, nr, mother, x, y, z, irot, konly, upar, np); 
} 

//______________________________________________________________________________ 
void TFluka::Gsposp(const char *name, Int_t nr, const char *mother,  
		    Double_t x, Double_t y, Double_t z, Int_t irot,
		    const char *konly, Double_t *upar, Int_t np)  {
  //
  fGeom->Gsposp(name, nr, mother, x, y, z, irot, konly, upar, np); 
} 

//______________________________________________________________________________ 
void TFluka::Gsbool(const char* /*onlyVolName*/, const char* /*manyVolName*/) {
//
// Nothing to do with TGeo
   Warning("Gsbool", "Not implemented with TGeo");
}

//______________________________________________________________________________ 
void TFluka::SetCerenkov(Int_t /*itmed*/, Int_t /*npckov*/, Float_t */*ppckov*/,
			 Float_t * /*absco*/, Float_t * /*effic*/, Float_t * /*rindex*/) {
//
// Not implemented with TGeo - what G4 did ? Any FLUKA card generated?
   Warning("SetCerenkov", "Not implemented with TGeo");
}  

//______________________________________________________________________________ 
void TFluka::SetCerenkov(Int_t /*itmed*/, Int_t /*npckov*/, Double_t * /*ppckov*/,
			 Double_t * /*absco*/, Double_t * /*effic*/, Double_t * /*rindex*/) {
//
// Not implemented with TGeo - what G4 did ? Any FLUKA card generated?
   Warning("SetCerenkov", "Not implemented with TGeo");
}  
    
// Euclid
//______________________________________________________________________________ 
void TFluka::WriteEuclid(const char* /*fileName*/, const char* /*topVol*/, 
                          Int_t /*number*/, Int_t /*nlevel*/) {
//
// Not with TGeo
   Warning("WriteEuclid", "Not implemented with TGeo");
} 



//_____________________________________________________________________________
// methods needed by the stepping
//____________________________________________________________________________ 

Int_t TFluka::GetMedium() const {
//
//  Get the medium number for the current fluka region
//
    return fGeom->GetMedium(); // this I need to check due to remapping !!!
}



//____________________________________________________________________________ 
// particle table usage
// ID <--> PDG transformations
//_____________________________________________________________________________
Int_t TFluka::IdFromPDG(Int_t pdg) const 
{
    //
    // Return Fluka code from PDG and pseudo ENDF code
    
    // Catch the feedback photons
    if (pdg == 50000051) return (-1);
    // MCIHAD() goes from pdg to fluka internal.
    Int_t intfluka = mcihad(pdg);
    // KPTOIP array goes from internal to official
    return GetFlukaKPTOIP(intfluka);
}

//______________________________________________________________________________ 
Int_t TFluka::PDGFromId(Int_t id) const 
{
  //
  // Return PDG code and pseudo ENDF code from Fluka code

  // IPTOKP array goes from official to internal

    if (id == -1) {
// Cerenkov photon
	if (fVerbosityLevel >= 1)
	    printf("\n PDGFromId: Cerenkov Photon \n");
	return  50000050;
    }
// Error id    
    if (id == 0) {
	if (fVerbosityLevel >= 1)
	    printf("PDGFromId: Error id = 0\n");
	return -1;
    }
// Good id    
    Int_t intfluka = GetFlukaIPTOKP(id);
    if (intfluka == 0) {
	if (fVerbosityLevel >= 1)
	    printf("PDGFromId: Error intfluka = 0: %d\n", id);
	return -1;
    } else if (intfluka < 0) {
	if (fVerbosityLevel >= 1)
	    printf("PDGFromId: Error intfluka < 0: %d\n", id);
	return -1;
    }
    if (fVerbosityLevel >= 3)
	printf("mpdgha called with %d %d \n", id, intfluka);
    // MPDGHA() goes from fluka internal to pdg.
    return mpdgha(intfluka);
}

//_____________________________________________________________________________
// methods for physics management
//____________________________________________________________________________ 
//
// set methods
//

//______________________________________________________________________________ 
void TFluka::SetProcess(const char* flagName, Int_t flagValue)
{
  Int_t i;
  if (iNbOfProc < 100) {
    for (i=0; i<iNbOfProc; i++) {
      if (strcmp(&sProcessFlag[i][0],flagName) == 0) {
        iProcessValue[iNbOfProc] = flagValue;
	goto fin;
      }
    }
    strcpy(&sProcessFlag[iNbOfProc][0],flagName);
    iProcessValue[iNbOfProc++] = flagValue;
  }
  else
    cout << "Nb of SetProcess calls exceeds 100 - ignored" << endl;
fin:
  iNbOfProc = iNbOfProc;
}

//______________________________________________________________________________ 
void TFluka::SetCut(const char* cutName, Double_t cutValue)
{
  Int_t i;
  if (iNbOfCut < 100) {
    for (i=0; i<iNbOfCut; i++) {
      if (strcmp(&sCutFlag[i][0],cutName) == 0) {
        fCutValue[iNbOfCut] = cutValue;
	goto fin;
      }
    }
    strcpy(&sCutFlag[iNbOfCut][0],cutName);
    fCutValue[iNbOfCut++] = cutValue;
  }
  else
    cout << "Nb of SetCut calls exceeds 100 - ignored" << endl;
fin:
  iNbOfCut = iNbOfCut;
}

//______________________________________________________________________________ 
Double_t TFluka::Xsec(char*, Double_t, Int_t, Int_t)
{
  printf("WARNING: Xsec not yet implemented !\n"); return -1.;
}


//______________________________________________________________________________ 
void TFluka::InitPhysics()
{
// Last material number taken from the "corealice.inp" file, presently 31
// !!! it should be available from Flugg !!!
  Int_t i, j, k;
  Double_t fCut;
  Float_t fLastMaterial = 31.0;
 
// construct file names
  TString sAliceInp = getenv("ALICE_ROOT");
  TString sAliceTmp = sAliceInp;
  sAliceInp +="/TFluka/input/";
  sAliceTmp +="/tmp/flukaMat.inp";
  TString sAliceCoreInp = sAliceInp;
  sAliceInp += GetInputFileName();
  sAliceCoreInp += GetCoreInputFileName();
  ifstream AliceCoreInp(sAliceCoreInp.Data());
  ifstream AliceFlukaMat(sAliceTmp.Data());
  ofstream AliceInp(sAliceInp.Data());

// copy core input file 
  Char_t sLine[255];
  Float_t fEventsPerRun;

  while (AliceCoreInp.getline(sLine,255)) {
    if (strncmp(sLine,"GEOEND",6) != 0)
      AliceInp << sLine << endl; // copy until GEOEND card
    else {
      AliceInp << "GEOEND" << endl; // add GEOEND card
      goto flukamat;
    }
  } // end of while until GEOEND card

flukamat:
  while (AliceFlukaMat.getline(sLine,255)) { // copy flukaMat.inp file
        AliceInp << sLine << endl;
  }

  while (AliceCoreInp.getline(sLine,255)) {
    if (strncmp(sLine,"START",5) != 0)
      AliceInp << sLine << endl;
    else {
      sscanf(sLine+10,"%10f",&fEventsPerRun);
      goto fin;
    }
  } //end of while until START card

fin:
// in G3 the process control values meaning can be different for
// different processes, but for most of them is:
//   0  process is not activated
//   1  process is activated WITH generation of secondaries
//   2  process is activated WITHOUT generation of secondaries
// if process does not generate secondaries => 1 same as 2
//
// Exceptions:
//   MULS:  also 3
//   LOSS:  also 3, 4
//   RAYL:  only 0,1
//   HADR:  may be > 2
//
 
// Loop over number of SetProcess calls  
  AliceInp << "*----------------------------------------------------------------------------- "; 
  AliceInp << endl;
  AliceInp << "*----- The following data are generated from SetProcess and SetCut calls ----- "; 
  AliceInp << endl;
  AliceInp << "*----------------------------------------------------------------------------- "; 
    AliceInp << endl;
  for (i=0; i<iNbOfProc; i++) {

    // annihilation
    // G3 default value: 1
    // G4 processes: G4eplusAnnihilation/G4IeplusAnnihilation
    // Particles: e+
    // Physics:   EM
    // flag = 0 no annihilation
    // flag = 1 annihilation, decays processed
    // flag = 2 annihilation, no decay product stored
    // gMC ->SetProcess("ANNI",1); // EMFCUT   -1.   0.  0. 3. lastmat 0. ANNH-THR
    if (strncmp(&sProcessFlag[i][0],"ANNI",4) == 0) {
      if (iProcessValue[i] == 1 || iProcessValue[i] == 2) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Kinetic energy threshold (GeV) for e+ annihilation - resets to default=0."; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('ANNI',1) or SetProcess('ANNI',2)"; 
        AliceInp << endl;
        AliceInp << setw(10) << "EMFCUT    "; 
        AliceInp << setiosflags(ios::scientific) << setprecision(5);
        AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
        AliceInp << setw(10) << -1.0; // kinetic energy threshold (GeV) for e+ annihilation (resets to default=0)
        AliceInp << setw(10) << 0.0;  // not used
        AliceInp << setw(10) << 0.0;  // not used
        AliceInp << setw(10) << 3.0;  // lower bound of the material indices in which the respective thresholds apply
        AliceInp << setw(10) << setprecision(2);
        AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(1);
        AliceInp << setw(10) << 1.0;  // step length in assigning indices
        AliceInp << setw(8)  << "ANNH-THR"; 
        AliceInp << endl;
      }
      else if (iProcessValue[i] == 0) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*No annihilation - no FLUKA card generated"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('ANNI',0)"; 
        AliceInp << endl;
      }
      else  {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Illegal flag value in SetProcess('ANNI',?) call."; 
        AliceInp << endl;
        AliceInp << "*No FLUKA card generated"; 
        AliceInp << endl;
      }
    }
    
    // bremsstrahlung and pair production are both activated
    // G3 default value: 1
    // G4 processes: G4eBremsstrahlung/G4IeBremsstrahlung,
    //               G4MuBremsstrahlung/G4IMuBremsstrahlung,
    //               G4LowEnergyBremstrahlung
    // Particles: e-/e+; mu+/mu-
    // Physics:   EM
    // flag = 0 no bremsstrahlung
    // flag = 1 bremsstrahlung, photon processed
    // flag = 2 bremsstrahlung, no photon stored
    // gMC ->SetProcess("BREM",1); // PAIRBREM  2.   0.  0. 3. lastmat
                                 // EMFCUT   -1.   0.  0. 3. lastmat 0. ELPO-THR
    // G3 default value: 1
    // G4 processes: G4GammaConversion,
    //               G4MuPairProduction/G4IMuPairProduction
    //               G4LowEnergyGammaConversion
    // Particles: gamma, mu
    // Physics:   EM
    // flag = 0 no delta rays
    // flag = 1 delta rays, secondaries processed
    // flag = 2 delta rays, no secondaries stored
    // gMC ->SetProcess("PAIR",1); // PAIRBREM  1.   0.  0. 3. lastmat
                                 // EMFCUT    0.   0. -1. 3. lastmat 0. PHOT-THR
    else if ((strncmp(&sProcessFlag[i][0],"PAIR",4) == 0) && (iProcessValue[i] == 1 || iProcessValue[i] == 2)) {
      for (j=0; j<iNbOfProc; j++) {
        if ((strncmp(&sProcessFlag[j][0],"BREM",4) == 0) && (iProcessValue[j] == 1 || iProcessValue[j] == 2)) {
          AliceInp << "*"; 
          AliceInp << endl;
          AliceInp << "*Bremsstrahlung and pair production by muons and charged hadrons both activated"; 
          AliceInp << endl;
          AliceInp << "*Generated from call: SetProcess('BREM',1) and SetProcess('PAIR',1)"; 
          AliceInp << endl;
          AliceInp << "*Energy threshold set by call SetCut('BCUTM',cut) or set to 0."; 
          AliceInp << endl;
          AliceInp << "*Energy threshold set by call SetCut('PPCUTM',cut) or set to 0."; 
          AliceInp << endl;
          AliceInp << setw(10) << "PAIRBREM  "; 
          AliceInp << setiosflags(ios::scientific) << setprecision(5);
          AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
          AliceInp << setw(10) << 3.0; // bremsstrahlung and pair production by muons and charged hadrons both are activated
          // direct pair production by muons
          // G4 particles: "e-", "e+"
          // G3 default value: 0.01 GeV
          //gMC ->SetCut("PPCUTM",cut); // total energy cut for direct pair prod. by muons
          fCut = 0.0;
          for (k=0; k<iNbOfCut; k++) {
            if (strncmp(&sCutFlag[k][0],"PPCUTM",6) == 0) fCut = fCutValue[k];
          }
          AliceInp << setiosflags(ios::scientific) << setprecision(5);
          AliceInp << setw(10) << fCut; // e+, e- kinetic energy threshold (in GeV) for explicit pair production.
          // muon and hadron bremsstrahlung
          // G4 particles: "gamma"
          // G3 default value: CUTGAM=0.001 GeV
          //gMC ->SetCut("BCUTM",cut);  // cut for muon and hadron bremsstrahlung
          fCut = 0.0;
          for (k=0; k<iNbOfCut; k++) {
            if (strncmp(&sCutFlag[k][0],"BCUTM",5) == 0) fCut = fCutValue[k];
          }
          AliceInp << setiosflags(ios::scientific) << setprecision(5);
          AliceInp << setw(10) << fCut; // photon energy threshold (GeV) for explicit bremsstrahlung production
          AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
          AliceInp << setw(10) << 3.0; // lower bound of the material indices in which the respective thresholds apply
          AliceInp << setw(10) << setprecision(2);
          AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
          AliceInp << endl;
	  
          // for e+ and e-
          AliceInp << "*"; 
          AliceInp << endl;
          AliceInp << "*Kinetic energy threshold (GeV) for e+/e- bremsstrahlung - resets to default=0."; 
          AliceInp << endl;
          AliceInp << "*Generated from call: SetProcess('BREM',1);"; 
          AliceInp << endl;
          AliceInp << setw(10) << "EMFCUT    "; 
          fCut = -1.0;
          for (k=0; k<iNbOfCut; k++) {
            if (strncmp(&sCutFlag[k][0],"BCUTE",5) == 0) fCut = fCutValue[k];
          }
          AliceInp << setiosflags(ios::scientific) << setprecision(5);
          AliceInp << setw(10) << fCut; // kinetic energy threshold (GeV) for e+/e- bremsstrahlung (resets to default=0)
          AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
          AliceInp << setw(10) << 0.0;  // not used
          AliceInp << setw(10) << 0.0;  // not used
          AliceInp << setw(10) << 3.0;  // lower bound of the material indices in which the respective thresholds apply
          AliceInp << setw(10) << setprecision(2);
          AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
          AliceInp << setprecision(1);
          AliceInp << setw(10) << 1.0; // step length in assigning indices
          AliceInp << setw(8)  << "ELPO-THR"; 
          AliceInp << endl;
      
          // for e+ and e-
          AliceInp << "*"; 
          AliceInp << endl;
          AliceInp << "*Pair production by electrons is activated"; 
          AliceInp << endl;
          AliceInp << "*Generated from call: SetProcess('PAIR',1);"; 
          AliceInp << endl;
          AliceInp << setw(10) << "EMFCUT    "; 
          AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
          AliceInp << setw(10) << 0.0;  // energy threshold (GeV) for Compton scattering (= 0.0 : ignored)
          AliceInp << setw(10) << 0.0;  // energy threshold (GeV) for Photoelectric (= 0.0 : ignored)
          fCut = -1.0;
          for (j=0; j<iNbOfCut; j++) {
            if (strncmp(&sCutFlag[j][0],"CUTGAM",6) == 0) fCut = fCutValue[j];
          }
          AliceInp << setiosflags(ios::scientific) << setprecision(5);
          AliceInp << setw(10) << fCut; // energy threshold (GeV) for gamma pair production (< 0.0 : resets to default, = 0.0 : ignored)
          AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
          AliceInp << setw(10) << 3.0;  // lower bound of the material indices in which the respective thresholds apply
          AliceInp << setprecision(2);
          AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
          AliceInp << setprecision(1);
          AliceInp << setw(10) << 1.0;  // step length in assigning indices
          AliceInp << setw(8) << "PHOT-THR"; 
          AliceInp << endl;
	  goto BOTH;
        } // end of if for BREM
      } // end of loop for BREM

      // only pair production by muons and charged hadrons is activated
      AliceInp << "*"; 
      AliceInp << endl;
      AliceInp << "*Pair production by muons and charged hadrons is activated"; 
      AliceInp << endl;
      AliceInp << "*Generated from call: SetProcess('PAIR',1) or SetProcess('PAIR',2)"; 
      AliceInp << endl;
      AliceInp << "*Energy threshold set by call SetCut('PPCUTM',cut) or set to 0."; 
      AliceInp << endl;
      AliceInp << setw(10) << "PAIRBREM  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
      AliceInp << setw(10) << 1.0; // pair production by muons and charged hadrons is activated
      // direct pair production by muons
      // G4 particles: "e-", "e+"
      // G3 default value: 0.01 GeV
      //gMC ->SetCut("PPCUTM",cut); // total energy cut for direct pair prod. by muons
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
      AliceInp << setw(10) << 0.0; // e+, e- kinetic energy threshold (in GeV) for explicit pair production.
      AliceInp << setw(10) << 0.0; // no explicit bremsstrahlung production is simulated
      AliceInp << setw(10) << 3.0; // lower bound of the material indices in which the respective thresholds apply
      AliceInp << setprecision(2);
      AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
      AliceInp << endl;
      
      // for e+ and e-
      AliceInp << "*"; 
      AliceInp << endl;
      AliceInp << "*Pair production by electrons is activated"; 
      AliceInp << endl;
      AliceInp << "*Generated from call: SetProcess('PAIR',1) or SetProcess('PAIR',2)"; 
      AliceInp << endl;
      AliceInp << setw(10) << "EMFCUT    "; 
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
      AliceInp << setw(10) << 0.0;  // energy threshold (GeV) for Compton scattering (= 0.0 : ignored)
      AliceInp << setw(10) << 0.0;  // energy threshold (GeV) for Photoelectric (= 0.0 : ignored)

      fCut = -1.0;
      for (j=0; j<iNbOfCut; j++) {
        if (strncmp(&sCutFlag[j][0],"CUTGAM",6) == 0) fCut = fCutValue[j];
      }
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << fCut; // energy threshold (GeV) for gamma pair production (< 0.0 : resets to default, = 0.0 : ignored)
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
      AliceInp << setw(10) << 3.0;  // lower bound of the material indices in which the respective thresholds apply
      AliceInp << setprecision(2);
      AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
      AliceInp << setprecision(1);
      AliceInp << setw(10) << 1.0;  // step length in assigning indices
      AliceInp << setw(8) << "PHOT-THR"; 
      AliceInp << endl;

BOTH:
    k = 0;
    } // end of if for PAIR



    // bremsstrahlung
    // G3 default value: 1
    // G4 processes: G4eBremsstrahlung/G4IeBremsstrahlung,
    //               G4MuBremsstrahlung/G4IMuBremsstrahlung,
    //               G4LowEnergyBremstrahlung
    // Particles: e-/e+; mu+/mu-
    // Physics:   EM
    // flag = 0 no bremsstrahlung
    // flag = 1 bremsstrahlung, photon processed
    // flag = 2 bremsstrahlung, no photon stored
    // gMC ->SetProcess("BREM",1); // PAIRBREM  2.   0.  0. 3. lastmat
                                 // EMFCUT   -1.   0.  0. 3. lastmat 0. ELPO-THR
    else if (strncmp(&sProcessFlag[i][0],"BREM",4) == 0) {
      for (j=0; j<iNbOfProc; j++) {
        if ((strncmp(&sProcessFlag[j][0],"PAIR",4) == 0) && iProcessValue[j] == 1) goto NOBREM;
      }
      if (iProcessValue[i] == 1 || iProcessValue[i] == 2) { 
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Bremsstrahlung by muons and charged hadrons is activated"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('BREM',1) or SetProcess('BREM',2)"; 
        AliceInp << endl;
        AliceInp << "*Energy threshold set by call SetCut('BCUTM',cut) or set to 0."; 
        AliceInp << endl;
        AliceInp << setw(10) << "PAIRBREM  "; 
        AliceInp << setiosflags(ios::scientific) << setprecision(5);
        AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
        AliceInp << setw(10) << 2.0; // bremsstrahlung by muons and charged hadrons is activated
        AliceInp << setw(10) << 0.0; // no meaning
        // muon and hadron bremsstrahlung
        // G4 particles: "gamma"
        // G3 default value: CUTGAM=0.001 GeV
        //gMC ->SetCut("BCUTM",cut);  // cut for muon and hadron bremsstrahlung
        fCut = 0.0;
        for (j=0; j<iNbOfCut; j++) {
          if (strncmp(&sCutFlag[j][0],"BCUTM",5) == 0) fCut = fCutValue[j];
        }
        AliceInp << setw(10) << fCut; // photon energy threshold (GeV) for explicit bremsstrahlung production
        AliceInp << setw(10) << 3.0; // lower bound of the material indices in which the respective thresholds apply
        AliceInp << setw(10) << setprecision(2);
        AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
        AliceInp << endl;
 
        // for e+ and e-
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Kinetic energy threshold (GeV) for e+/e- bremsstrahlung - resets to default=0."; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('BREM',1);"; 
        AliceInp << endl;
        AliceInp << setw(10) << "EMFCUT    "; 
        AliceInp << setiosflags(ios::scientific) << setprecision(5);
        AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
        AliceInp << setw(10) << -1.0; // kinetic energy threshold (GeV) for e+/e- bremsstrahlung (resets to default=0)
        AliceInp << setw(10) << 0.0;  // not used
        AliceInp << setw(10) << 0.0;  // not used
        AliceInp << setw(10) << 3.0;  // lower bound of the material indices in which the respective thresholds apply
        AliceInp << setw(10) << setprecision(2);
        AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(1);
        AliceInp << setw(10) << 1.0; // step length in assigning indices
        AliceInp << setw(8)  << "ELPO-THR"; 
        AliceInp << endl;
      }
      else if (iProcessValue[i] == 0) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*No bremsstrahlung - no FLUKA card generated"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('BREM',0)"; 
        AliceInp << endl;
      }
      else  {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Illegal flag value in SetProcess('BREM',?) call."; 
        AliceInp << endl;
        AliceInp << "*No FLUKA card generated"; 
        AliceInp << endl;
      }
NOBREM:
      j = 0;
    } // end of else if (strncmp(&sProcessFlag[i][0],"BREM",4) == 0)

    
    // Cerenkov photon generation
    // G3 default value: 0
    // G4 process: G4Cerenkov
    // 
    // Particles: charged
    // Physics:   Optical
    // flag = 0 no Cerenkov photon generation
    // flag = 1 Cerenkov photon generation
    // flag = 2 Cerenkov photon generation with primary stopped at each step
    //xx gMC ->SetProcess("CKOV",1); // ??? Cerenkov photon generation
    else if (strncmp(&sProcessFlag[i][0],"CKOV",4) == 0) {
      if (iProcessValue[i] == 1 || iProcessValue[i] == 2) { 
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Cerenkov photon generation"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('CKOV',1) or SetProcess('CKOV',2)"; 
        AliceInp << endl;
        AliceInp << setw(10) << "OPT-PROD  "; 
        AliceInp << setiosflags(ios::scientific) << setprecision(5);
        AliceInp << setw(10) <<  2.07e-9 ; //  minimum Cerenkov photon emission energy (in GeV!). Default: 2.07E-9 GeV (corresponding to 600 nm)
        AliceInp << setw(10) << 4.96e-9;  // maximum Cerenkov photon emission energy (in GeV!). Default: 4.96E-9 GeV (corresponding to 250 nm)
        AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
        AliceInp << setw(10) << 0.0;  // not used
        AliceInp << setw(10) << 3.0;  // lower bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(2);
        AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(1);
        AliceInp << setw(10) << 1.0; // step length in assigning indices
        AliceInp << setw(8) << "CERENKOV"; 
        AliceInp << endl;
      }
      else if (iProcessValue[i] == 0) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*No Cerenkov photon generation"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('CKOV',0)"; 
        AliceInp << endl;
        AliceInp << setw(10) << "OPT-PROD  "; 
        AliceInp << setiosflags(ios::scientific) << setprecision(5);
        AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
        AliceInp << setw(10) << 0.0;  // not used
        AliceInp << setw(10) << 0.0;  // not used
        AliceInp << setw(10) << 0.0;  // not used
        AliceInp << setw(10) << 3.0;  // lower bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(2);
        AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(1);
        AliceInp << setw(10) << 1.0; // step length in assigning indices
        AliceInp << setw(8) << "CERE-OFF"; 
        AliceInp << endl;
      }
      else  {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Illegal flag value in SetProcess('CKOV',?) call."; 
        AliceInp << endl;
        AliceInp << "*No FLUKA card generated"; 
        AliceInp << endl;
      }
    } // end of else if (strncmp(&sProcessFlag[i][0],"CKOV",4) == 0)
	    

    // Compton scattering
    // G3 default value: 1
    // G4 processes: G4ComptonScattering,
    //               G4LowEnergyCompton,
    //               G4PolarizedComptonScattering
    // Particles: gamma
    // Physics:   EM
    // flag = 0 no Compton scattering
    // flag = 1 Compton scattering, electron processed
    // flag = 2 Compton scattering, no electron stored
    // gMC ->SetProcess("COMP",1); // EMFCUT   -1.   0.  0. 3. lastmat 0. PHOT-THR
    else if (strncmp(&sProcessFlag[i][0],"COMP",4) == 0) {
      if (iProcessValue[i] == 1 || iProcessValue[i] == 2) { 
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Energy threshold (GeV) for Compton scattering - resets to default=0."; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('COMP',1);"; 
        AliceInp << endl;
        AliceInp << setw(10) << "EMFCUT    "; 
        AliceInp << setiosflags(ios::scientific) << setprecision(5);
        AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
        AliceInp << setw(10) << -1.0; // energy threshold (GeV) for Compton scattering - resets to default=0.
        AliceInp << setw(10) << 0.0;  // not used
        AliceInp << setw(10) << 0.0;  // not used
        AliceInp << setw(10) << 3.0;  // lower bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(2);
        AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(1);
        AliceInp << setw(10) << 1.0; // step length in assigning indices
        AliceInp << setw(8) << "PHOT-THR"; 
        AliceInp << endl;
      }
      else if (iProcessValue[i] == 0) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*No Compton scattering - no FLUKA card generated"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('COMP',0)"; 
        AliceInp << endl;
      }
      else  {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Illegal flag value in SetProcess('COMP',?) call."; 
        AliceInp << endl;
        AliceInp << "*No FLUKA card generated"; 
        AliceInp << endl;
      }
    } // end of else if (strncmp(&sProcessFlag[i][0],"COMP",4) == 0)

    // decay
    // G3 default value: 1
    // G4 process: G4Decay
    // 
    // Particles: all which decay is applicable for
    // Physics:   General
    // flag = 0 no decays
    // flag = 1 decays, secondaries processed
    // flag = 2 decays, no secondaries stored
    //gMC ->SetProcess("DCAY",1); // not available
    else if ((strncmp(&sProcessFlag[i][0],"DCAY",4) == 0) && iProcessValue[i] == 1) 
      cout << "SetProcess for flag=" << &sProcessFlag[i][0] << " value=" << iProcessValue[i] << " not avaliable!" << endl;
      
    // delta-ray
    // G3 default value: 2
    // !! G4 treats delta rays in different way
    // G4 processes: G4eIonisation/G4IeIonization,
    //               G4MuIonisation/G4IMuIonization,
    //               G4hIonisation/G4IhIonisation
    // Particles: charged
    // Physics:   EM
    // flag = 0 no energy loss
    // flag = 1 restricted energy loss fluctuations
    // flag = 2 complete energy loss fluctuations
    // flag = 3 same as 1
    // flag = 4 no energy loss fluctuations
    // gMC ->SetProcess("DRAY",0); // DELTARAY 1.E+6 0.  0. 3. lastmat 0.
    else if (strncmp(&sProcessFlag[i][0],"DRAY",4) == 0) {
      if (iProcessValue[i] == 0 || iProcessValue[i] == 4) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Kinetic energy threshold (GeV) for delta ray production"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('DRAY',0) or SetProcess('DRAY',4)"; 
        AliceInp << endl;
        AliceInp << "*No delta ray production by muons - threshold set artificially high"; 
        AliceInp << endl;
        AliceInp << setw(10) << "DELTARAY  "; 
        AliceInp << setiosflags(ios::scientific) << setprecision(5);
        AliceInp << setw(10) << 1.0e+6; // kinetic energy threshold (GeV) for delta ray production (discrete energy transfer)
        AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
        AliceInp << setw(10) << 0.0; // ignored
        AliceInp << setw(10) << 0.0; // ignored
        AliceInp << setw(10) << 3.0; // lower bound of the material indices in which the respective thresholds apply
        AliceInp << setw(10) << setprecision(2);
        AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(1);
        AliceInp << setw(10) << 1.0; // step length in assigning indices
        AliceInp << endl;
      }
      else if (iProcessValue[i] == 1 || iProcessValue[i] == 2 || iProcessValue[i] == 3) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Kinetic energy threshold (GeV) for delta ray production"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('DRAY',flag), flag=1,2,3"; 
        AliceInp << endl;
        AliceInp << "*Delta ray production by muons switched on"; 
        AliceInp << endl;
        AliceInp << "*Energy threshold set by call SetCut('DCUTM',cut) or set to 0."; 
        AliceInp << endl;
        AliceInp << setw(10) << "DELTARAY  "; 
        AliceInp << setiosflags(ios::scientific) << setprecision(5);
        fCut = 1.0e+6;
        for (j=0; j<iNbOfCut; j++) {
          if (strncmp(&sCutFlag[j][0],"DCUTM",5) == 0) fCut = fCutValue[j];
        }
        AliceInp << setw(10) << fCut; // kinetic energy threshold (GeV) for delta ray production (discrete energy transfer)
        AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
        AliceInp << setw(10) << 0.0; // ignored
        AliceInp << setw(10) << 0.0; // ignored
        AliceInp << setw(10) << 3.0; // lower bound of the material indices in which the respective thresholds apply
        AliceInp << setw(10) << setprecision(2);
        AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(1);
        AliceInp << setw(10) << 1.0; // step length in assigning indices
        AliceInp << endl;
      }
      else  {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Illegal flag value in SetProcess('DRAY',?) call."; 
        AliceInp << endl;
        AliceInp << "*No FLUKA card generated"; 
        AliceInp << endl;
      }
    } // end of else if (strncmp(&sProcessFlag[i][0],"DRAY",4) == 0)
    
    // hadronic process
    // G3 default value: 1
    // G4 processes: all defined by TG4PhysicsConstructorHadron
    //  
    // Particles: hadrons
    // Physics:   Hadron
    // flag = 0 no multiple scattering
    // flag = 1 hadronic interactions, secondaries processed
    // flag = 2 hadronic interactions, no secondaries stored
    // gMC ->SetProcess("HADR",1); // ??? hadronic process
    //Select pure GEANH (HADR 1) or GEANH/NUCRIN (HADR 3) ?????
    else if (strncmp(&sProcessFlag[i][0],"HADR",4) == 0) {
      if (iProcessValue[i] == 1 || iProcessValue[i] == 2) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Hadronic interaction is ON by default in FLUKA"; 
        AliceInp << endl;
        AliceInp << "*No FLUKA card generated"; 
        AliceInp << endl;
      }
      else if (iProcessValue[i] == 0) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Hadronic interaction is set OFF"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('HADR',0);"; 
        AliceInp << endl;
        AliceInp << setw(10) << "MULSOPT  "; 
        AliceInp << setiosflags(ios::scientific) << setprecision(5);
        AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
        AliceInp << setw(10) << 0.0;  // ignored
        AliceInp << setw(10) << 3.0;  // multiple scattering for hadrons and muons is completely suppressed
        AliceInp << setw(10) << 0.0;  // no spin-relativistic corrections
        AliceInp << setw(10) << 3.0;  // lower bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(2);
        AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
        AliceInp << endl;

      }
      else  {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Illegal flag value in SetProcess('HADR',?) call."; 
        AliceInp << endl;
        AliceInp << "*No FLUKA card generated"; 
        AliceInp << endl;
      }
    } // end of else if (strncmp(&sProcessFlag[i][0],"HADR",4) == 0)


    // energy loss
    // G3 default value: 2
    // G4 processes: G4eIonisation/G4IeIonization,
    //               G4MuIonisation/G4IMuIonization,
    //               G4hIonisation/G4IhIonisation
    // 
    // Particles: charged
    // Physics:   EM
    // flag=0 no energy loss
    // flag=1 restricted energy loss fluctuations
    // flag=2 complete energy loss fluctuations
    // flag=3 same as 1
    // flag=4 no energy loss fluctuations
    // If the value ILOSS is changed, then (in G3) cross-sections and energy
    // loss tables must be recomputed via the command 'PHYSI'
    // gMC ->SetProcess("LOSS",2); // ??? IONFLUCT ? energy loss
    else if (strncmp(&sProcessFlag[i][0],"LOSS",4) == 0) {
      if (iProcessValue[i] == 2) { // complete energy loss fluctuations
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Complete energy loss fluctuations do not exist in FLUKA"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('LOSS',2);"; 
        AliceInp << endl;
        AliceInp << "*flag=2=complete energy loss fluctuations"; 
        AliceInp << endl;
        AliceInp << "*No input card generated"; 
        AliceInp << endl;
      }
      else if (iProcessValue[i] == 1 || iProcessValue[i] == 3) { // restricted energy loss fluctuations
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Restricted energy loss fluctuations"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('LOSS',1) or SetProcess('LOSS',3)"; 
        AliceInp << endl;
        AliceInp << setw(10) << "IONFLUCT  "; 
        AliceInp << setiosflags(ios::scientific) << setprecision(5);
        AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
        AliceInp << setw(10) << 1.0;  // restricted energy loss fluctuations (for hadrons and muons) switched on
        AliceInp << setw(10) << 1.0;  // restricted energy loss fluctuations (for e+ and e-) switched on
        AliceInp << setw(10) << 1.0;  // minimal accuracy
        AliceInp << setw(10) << 3.0;  // lower bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(2);
        AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
        AliceInp << endl;
      }
      else if (iProcessValue[i] == 4) { // no energy loss fluctuations
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*No energy loss fluctuations"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('LOSS',4)"; 
        AliceInp << endl;
        AliceInp << setw(10) << -1.0;  // restricted energy loss fluctuations (for hadrons and muons) switched off
        AliceInp << setw(10) << -1.0;  // restricted energy loss fluctuations (for e+ and e-) switched off
        AliceInp << setw(10) << 1.0;  // minimal accuracy
        AliceInp << setw(10) << 3.0;  // lower bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(2);
        AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
        AliceInp << endl;
      }
      else  {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Illegal flag value in SetProcess('LOSS',?) call."; 
        AliceInp << endl;
        AliceInp << "*No FLUKA card generated"; 
        AliceInp << endl;
      }
    } // end of else if (strncmp(&sProcessFlag[i][0],"LOSS",4) == 0)
	
 
    // multiple scattering
    // G3 default value: 1
    // G4 process: G4MultipleScattering/G4IMultipleScattering
    // 
    // Particles: charged
    // Physics:   EM
    // flag = 0 no multiple scattering
    // flag = 1 Moliere or Coulomb scattering
    // flag = 2 Moliere or Coulomb scattering
    // flag = 3 Gaussian scattering
    // gMC ->SetProcess("MULS",1); // MULSOPT multiple scattering
    else if (strncmp(&sProcessFlag[i][0],"MULS",4) == 0) {
      if (iProcessValue[i] == 1 || iProcessValue[i] == 2 || iProcessValue[i] == 3) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Multiple scattering is ON by default for e+e- and for hadrons/muons"; 
        AliceInp << endl;
        AliceInp << "*No FLUKA card generated"; 
        AliceInp << endl;
      }
      else if (iProcessValue[i] == 0) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Multiple scattering is set OFF"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('MULS',0);"; 
        AliceInp << endl;
        AliceInp << setw(10) << "MULSOPT  "; 
        AliceInp << setiosflags(ios::scientific) << setprecision(5);
        AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
        AliceInp << setw(10) << 0.0;  // ignored
        AliceInp << setw(10) << 3.0;  // multiple scattering for hadrons and muons is completely suppressed
        AliceInp << setw(10) << 3.0;  // multiple scattering for e+ and e- is completely suppressed
        AliceInp << setw(10) << 3.0;  // lower bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(2);
        AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
        AliceInp << endl;
      }
      else  {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Illegal flag value in SetProcess('MULS',?) call."; 
        AliceInp << endl;
        AliceInp << "*No FLUKA card generated"; 
        AliceInp << endl;
      }
    } // end of else if (strncmp(&sProcessFlag[i][0],"MULS",4) == 0)


    // muon nuclear interaction
    // G3 default value: 0
    // G4 processes: G4MuNuclearInteraction,
    // G4MuonMinusCaptureAtRest
    // 
    // Particles: mu
    // Physics:   Not set
    // flag = 0 no muon-nuclear interaction
    // flag = 1 nuclear interaction, secondaries processed
    // flag = 2 nuclear interaction, secondaries not processed
    // gMC ->SetProcess("MUNU",1); // MUPHOTON  1.   0.  0. 3. lastmat
    else if (strncmp(&sProcessFlag[i][0],"MUNU",4) == 0) {
      if (iProcessValue[i] == 1) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Muon nuclear interactions with production of secondary hadrons"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('MUNU',1);"; 
        AliceInp << endl;
        AliceInp << setw(10) << "MUPHOTON  "; 
        AliceInp << setiosflags(ios::scientific) << setprecision(5);
        AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
        AliceInp << setw(10) << 1.0;  // full simulation of muon nuclear interactions and production of secondary hadrons
        AliceInp << setw(10) << 0.0; // ratio of longitudinal to transverse virtual photon cross-section - Default = 0.25.
        AliceInp << setw(10) << 0.0; // fraction of rho-like interactions ( must be < 1) - Default = 0.75.
        AliceInp << setprecision(1);
        AliceInp << setw(10) << 3.0;  // lower bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(2);
        AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
        AliceInp << endl;
      }
      else if (iProcessValue[i] == 2) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Muon nuclear interactions without production of secondary hadrons"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('MUNU',2);"; 
        AliceInp << endl;
        AliceInp << setw(10) << "MUPHOTON  "; 
        AliceInp << setiosflags(ios::scientific) << setprecision(5);
        AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
        AliceInp << setw(10) << 2.0; // full simulation of muon nuclear interactions and production of secondary hadrons
        AliceInp << setw(10) << 0.0; // ratio of longitudinal to transverse virtual photon cross-section - Default = 0.25.
        AliceInp << setw(10) << 0.0; // fraction of rho-like interactions ( must be < 1) - Default = 0.75.
        AliceInp << setprecision(1);
        AliceInp << setw(10) << 3.0;  // lower bound of the material indices in which the respective thresholds apply
        AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
        AliceInp << endl;
      }
      else if (iProcessValue[i] == 0) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*No muon nuclear interaction - no FLUKA card generated"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('MUNU',0)"; 
        AliceInp << endl;
      }
      else  {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Illegal flag value in SetProcess('MUNU',?) call."; 
        AliceInp << endl;
        AliceInp << "*No FLUKA card generated"; 
        AliceInp << endl;
      }
    } // end of else if (strncmp(&sProcessFlag[i][0],"MUNU",4) == 0)


    // photofission
    // G3 default value: 0
    // G4 process: ??
    //
    // Particles: gamma
    // Physics:   ??
    // gMC ->SetProcess("PFIS",0); // PHOTONUC -1.   0.  0. 3. lastmat 0.
    // flag = 0 no photon fission
    // flag = 1 photon fission, secondaries processed
    // flag = 2 photon fission, no secondaries stored
    else if (strncmp(&sProcessFlag[i][0],"PFIS",4) == 0) {
      if (iProcessValue[i] == 0) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*No photonuclear interactions";
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('PFIS',0);"; 
        AliceInp << endl;
        AliceInp << setw(10) << "PHOTONUC  "; 
        AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
        AliceInp << setw(10) << -1.0; // no photonuclear interactions
        AliceInp << setw(10) << 0.0;  // not used
        AliceInp << setw(10) << 0.0;  // not used
        AliceInp << setw(10) << 3.0;  // upper bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(2); 
        AliceInp << setw(10) << fLastMaterial;
        AliceInp << setprecision(1);  // upper bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(1);
        AliceInp << setw(10) << 1.0;  // step length in assigning indices
        AliceInp << endl;
      }
      else if (iProcessValue[i] == 1) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Photon nuclear interactions are activated at all energies";
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('PFIS',1);"; 
        AliceInp << endl;
        AliceInp << setw(10) << "PHOTONUC  "; 
        AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
        AliceInp << setw(10) << 1.0; // photonuclear interactions are activated at all energies
        AliceInp << setw(10) << 0.0; // not used
        AliceInp << setw(10) << 0.0; // not used
        AliceInp << setprecision(2); 
        AliceInp << setw(10) << 3.0; // upper bound of the material indices in which the respective thresholds apply
        AliceInp << setw(10) << fLastMaterial;
        AliceInp << setprecision(1); // upper bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(1);
        AliceInp << setw(10) << 1.0; // step length in assigning indices
        AliceInp << endl;
      }
      else if (iProcessValue[i] == 0) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*No photofission - no FLUKA card generated"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('PFIS',0)"; 
        AliceInp << endl;
      }
      else {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Illegal flag value in SetProcess('PFIS',?) call."; 
        AliceInp << endl;
        AliceInp << "*No FLUKA card generated"; 
        AliceInp << endl;
      }
    }

 
    // photo electric effect
    // G3 default value: 1
    // G4 processes: G4PhotoElectricEffect
    //               G4LowEnergyPhotoElectric
    // Particles: gamma
    // Physics:   EM
    // flag = 0 no photo electric effect
    // flag = 1 photo electric effect, electron processed
    // flag = 2 photo electric effect, no electron stored
    // gMC ->SetProcess("PHOT",1); // EMFCUT    0.  -1.  0. 3. lastmat 0. PHOT-THR
    else if (strncmp(&sProcessFlag[i][0],"PHOT",4) == 0) {
      if (iProcessValue[i] == 1 || iProcessValue[i] == 2) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Photo electric effect is activated"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('PHOT',1);"; 
        AliceInp << endl;
        AliceInp << setw(10) << "EMFCUT    "; 
        AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
        AliceInp << setw(10) << 0.0;  // ignored
        AliceInp << setw(10) << -1.0; // resets to default=0.
        AliceInp << setw(10) << 0.0;  // ignored
        AliceInp << setw(10) << 3.0;  // upper bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(2);
        AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(1);
        AliceInp << setw(10) << 1.0;  // step length in assigning indices
        AliceInp << setw(8) << "PHOT-THR"; 
        AliceInp << endl;
      }
      else if (iProcessValue[i] == 0) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*No photo electric effect - no FLUKA card generated"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('PHOT',0)"; 
        AliceInp << endl;
      }
      else {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Illegal flag value in SetProcess('PHOT',?) call."; 
        AliceInp << endl;
        AliceInp << "*No FLUKA card generated"; 
        AliceInp << endl;
      }
    } // else if (strncmp(&sProcessFlag[i][0],"PHOT",4) == 0)


    // Rayleigh scattering
    // G3 default value: 0
    // G4 process: G4OpRayleigh
    // 
    // Particles: optical photon
    // Physics:   Optical
    // flag = 0 Rayleigh scattering off
    // flag = 1 Rayleigh scattering on
    //xx gMC ->SetProcess("RAYL",1);
    else if (strncmp(&sProcessFlag[i][0],"RAYL",4) == 0) {
      if (iProcessValue[i] == 1) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Rayleigh scattering is ON by default in FLUKA"; 
        AliceInp << endl;
        AliceInp << "*No FLUKA card generated"; 
        AliceInp << endl;
      }
      else if (iProcessValue[i] == 0) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Rayleigh scattering is set OFF"; 
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('RAYL',0);"; 
        AliceInp << endl;
        AliceInp << setw(10) << "EMFRAY    "; 
        AliceInp << setiosflags(ios::scientific) << setprecision(5);
        AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
        AliceInp << setw(10) << -1.0;  // no Rayleigh scattering and no binding corrections for Compton
        AliceInp << setw(10) << 3.0;  // lower bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(2);
        AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
        AliceInp << endl;
      }
      else  {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Illegal flag value in SetProcess('RAYL',?) call."; 
        AliceInp << endl;
        AliceInp << "*No FLUKA card generated"; 
        AliceInp << endl;
      }
    } // end of else if (strncmp(&sProcessFlag[i][0],"RAYL",4) == 0)


    // synchrotron radiation in magnetic field
    // G3 default value: 0
    // G4 process: G4SynchrotronRadiation
    // 
    // Particles: ??
    // Physics:   Not set
    // flag = 0 no synchrotron radiation
    // flag = 1 synchrotron radiation
    //xx gMC ->SetProcess("SYNC",1); // synchrotron radiation generation
    else if (strncmp(&sProcessFlag[i][0],"SYNC",4) == 0) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Synchrotron radiation generation is NOT implemented in FLUKA"; 
        AliceInp << endl;
        AliceInp << "*No FLUKA card generated"; 
        AliceInp << endl;
    }
	

    // Automatic calculation of tracking medium parameters
    // flag = 0 no automatic calculation
    // flag = 1 automatic calculation
    //xx gMC ->SetProcess("AUTO",1); // ??? automatic computation of the tracking medium parameters
    else if (strncmp(&sProcessFlag[i][0],"AUTO",4) == 0) {
        AliceInp << "*"; 
        AliceInp << endl;
        AliceInp << "*Automatic calculation of tracking medium parameters is always ON in FLUKA"; 
        AliceInp << endl;
        AliceInp << "*No FLUKA card generated"; 
        AliceInp << endl;
    }
	

    // To control energy loss fluctuation model
    // flag = 0 Urban model
    // flag = 1 PAI model
    // flag = 2 PAI+ASHO model (not active at the moment)
    //xx gMC ->SetProcess("STRA",1); // ??? energy fluctuation model
    else if (strncmp(&sProcessFlag[i][0],"STRA",4) == 0) {
      if (iProcessValue[i] == 0 || iProcessValue[i] == 2 || iProcessValue[i] == 3) {
        AliceInp << "*";
        AliceInp << endl;
        AliceInp << "*Ionization energy losses calculation is activated";
        AliceInp << endl;
        AliceInp << "*Generated from call: SetProcess('STRA',n);, n=0,1,2";
        AliceInp << endl;
        AliceInp << setw(10) << "IONFLUCT  ";
        AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
        AliceInp << setw(10) << 1.0;  // restricted energy loss fluctuations
                                      // (for hadrons and muons) switched on
        AliceInp << setw(10) << 1.0;  // restricted energy loss fluctuations
                                      // (for e+ and e-) switched on
        AliceInp << setw(10) << 1.0;  // minimal accuracy
        AliceInp << setw(10) << 3.0;  // upper bound of the material indices in
                                      // which the respective thresholds apply
        AliceInp << setprecision(2);
        AliceInp << setw(10) << fLastMaterial; // upper bound of the material indices in which the respective thresholds apply
        AliceInp << setprecision(1);
        AliceInp << setw(10) << 1.0;  // step length in assigning indices
        AliceInp << endl;
      }
      else {
        AliceInp << "*";
        AliceInp << endl;
        AliceInp << "*Illegal flag value in SetProcess('STRA',?) call.";
        AliceInp << endl;
        AliceInp << "*No FLUKA card generated";
        AliceInp << endl;
      }
    } // else if (strncmp(&sProcessFlag[i][0],"STRA",4) == 0)




    else { // processes not yet treated

    // light photon absorption (Cerenkov photons)
    // it is turned on when Cerenkov process is turned on
    // G3 default value: 0
    // G4 process: G4OpAbsorption, G4OpBoundaryProcess
    // 
    // Particles: optical photon
    // Physics:   Optical
    // flag = 0 no absorption of Cerenkov photons
    // flag = 1 absorption of Cerenkov photons
    // gMC ->SetProcess("LABS",2); // ??? Cerenkov light absorption



      cout << "SetProcess for flag=" << &sProcessFlag[i][0] << " value=" << iProcessValue[i] << " not yet implemented!" << endl;
    }
  } //end of loop number of SetProcess calls

 
// Loop over number of SetCut calls  
  for (Int_t i=0; i<iNbOfCut; i++) {

    // cuts used in SetProcess calls
    if (strncmp(&sCutFlag[i][0],"BCUTM",5) == 0) continue;
    else if (strncmp(&sCutFlag[i][0],"BCUTE",5) == 0) continue;
    else if (strncmp(&sCutFlag[i][0],"DCUTM",5) == 0) continue;
    else if (strncmp(&sCutFlag[i][0],"PPCUTM",6) == 0) continue;

    // gammas
    // G4 particles: "gamma"
    // G3 default value: 0.001 GeV
    //gMC ->SetCut("CUTGAM",cut); // cut for gammas
    else if (strncmp(&sCutFlag[i][0],"CUTGAM",6) == 0) {
      AliceInp << "*"; 
      AliceInp << endl;
      AliceInp << "*Cut for gamma"; 
      AliceInp << endl;
      AliceInp << "*Generated from call: SetCut('CUTGAM',cut);"; 
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
      AliceInp << setw(10) << 7.0;
      AliceInp << endl;
    }

    // electrons
    // G4 particles: "e-"
    // ?? positrons
    // G3 default value: 0.001 GeV
    //gMC ->SetCut("CUTELE",cut); // cut for e+,e-
    else if (strncmp(&sCutFlag[i][0],"CUTELE",6) == 0) {
      AliceInp << "*"; 
      AliceInp << endl;
      AliceInp << "*Cut for electrons"; 
      AliceInp << endl;
      AliceInp << "*Generated from call: SetCut('CUTELE',cut);"; 
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
      AliceInp << setw(10) << 3.0;
      AliceInp << setw(10) << 4.0;
      AliceInp << setw(10) << 1.0;
      AliceInp << endl;
    }

    // neutral hadrons
    // G4 particles: of type "baryon", "meson", "nucleus" with zero charge
    // G3 default value: 0.01 GeV
    //gMC ->SetCut("CUTNEU",cut); // cut for neutral hadrons
    else if (strncmp(&sCutFlag[i][0],"CUTNEU",6) == 0) {
      AliceInp << "*"; 
      AliceInp << endl;
      AliceInp << "*Cut for neutral hadrons"; 
      AliceInp << endl;
      AliceInp << "*Generated from call: SetCut('CUTNEU',cut);"; 
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
      AliceInp << setw(10) << 8.0; // Neutron
      AliceInp << setw(10) << 9.0; // Antineutron
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(2);
      AliceInp << setw(10) << 12.0; // Kaon zero long
      AliceInp << setw(10) << 12.0; // Kaon zero long
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(2);
      AliceInp << setw(10) << 17.0; // Lambda, 18=Antilambda
      AliceInp << setw(10) << 19.0; // Kaon zero short
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(2);
      AliceInp << setw(10) << 22.0; // Sigma zero, Pion zero, Kaon zero
      AliceInp << setw(10) << 25.0; // Antikaon zero
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(2);
      AliceInp << setw(10) << 32.0; // Antisigma zero
      AliceInp << setw(10) << 32.0; // Antisigma zero
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(2);
      AliceInp << setw(10) << 34.0; // Xi zero
      AliceInp << setw(10) << 35.0; // AntiXi zero
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(2);
      AliceInp << setw(10) << 47.0; // D zero
      AliceInp << setw(10) << 48.0; // AntiD zero
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(2);
      AliceInp << setw(10) << 53.0; // Xi_c zero
      AliceInp << setw(10) << 53.0; // Xi_c zero
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(2);
      AliceInp << setw(10) << 55.0; // Xi'_c zero
      AliceInp << setw(10) << 56.0; // Omega_c zero
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(2);
      AliceInp << setw(10) << 59.0; // AntiXi_c zero
      AliceInp << setw(10) << 59.0; // AntiXi_c zero
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(2);
      AliceInp << setw(10) << 61.0; // AntiXi'_c zero
      AliceInp << setw(10) << 62.0; // AntiOmega_c zero
      AliceInp << endl;
    }

    // charged hadrons
    // G4 particles: of type "baryon", "meson", "nucleus" with non-zero charge
    // G3 default value: 0.01 GeV
    //gMC ->SetCut("CUTHAD",cut); // cut for charged hadrons
    else if (strncmp(&sCutFlag[i][0],"CUTHAD",6) == 0) {
      AliceInp << "*"; 
      AliceInp << endl;
      AliceInp << "*Cut for charged hadrons"; 
      AliceInp << endl;
      AliceInp << "*Generated from call: SetCut('CUTHAD',cut);"; 
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
      AliceInp << setw(10) << 1.0; // Proton
      AliceInp << setw(10) << 2.0; // Antiproton
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(2);
      AliceInp << setw(10) << 13.0; // Positive Pion, Negative Pion, Positive Kaon
      AliceInp << setw(10) << 16.0; // Negative Kaon
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(2);
      AliceInp << setw(10) << 20.0; // Negative Sigma
      AliceInp << setw(10) << 16.0; // Positive Sigma
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(2);
      AliceInp << setw(10) << 31.0; // Antisigma minus
      AliceInp << setw(10) << 33.0; // Antisigma plus
      AliceInp << setprecision(1);
      AliceInp << setw(10) << 2.0;  // step length
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(2);
      AliceInp << setw(10) << 36.0; // Negative Xi, Positive Xi, Omega minus
      AliceInp << setw(10) << 39.0; // Antiomega
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(2);
      AliceInp << setw(10) << 45.0; // D plus
      AliceInp << setw(10) << 46.0; // D minus
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(2);
      AliceInp << setw(10) << 49.0; // D_s plus, D_s minus, Lambda_c plus
      AliceInp << setw(10) << 52.0; // Xi_c plus
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(2);
      AliceInp << setw(10) << 54.0; // Xi'_c plus
      AliceInp << setw(10) << 60.0; // AntiXi'_c minus
      AliceInp << setprecision(1);
      AliceInp << setw(10) << 6.0;  // step length
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(2);
      AliceInp << setw(10) << 57.0; // Antilambda_c minus
      AliceInp << setw(10) << 58.0; // AntiXi_c minus
      AliceInp << endl;
    }

    // muons
    // G4 particles: "mu+", "mu-"
    // G3 default value: 0.01 GeV
    //gMC ->SetCut("CUTMUO",cut); // cut for mu+, mu-
    else if (strncmp(&sCutFlag[i][0],"CUTMUO",6) == 0) {
      AliceInp << "*"; 
      AliceInp << endl;
      AliceInp << "*Cut for muons"; 
      AliceInp << endl;
      AliceInp << "*Generated from call: SetCut('CUTMUO',cut);"; 
      AliceInp << endl;
      AliceInp << setw(10) << "PART-THR  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
      AliceInp << setprecision(2);
      AliceInp << setw(10) << 10.0;
      AliceInp << setw(10) << 11.0;
      AliceInp << endl;
    }
    // delta-rays by electrons
    // G4 particles: "e-"
    // G3 default value: 10**4 GeV
    // gMC ->SetCut("DCUTE",cut);  // cut for deltarays by electrons ???????????????
    else if (strncmp(&sCutFlag[i][0],"DCUTE",5) == 0) {
      AliceInp << "*"; 
      AliceInp << endl;
      AliceInp << "*Cut for delta rays by electrons ????????????";
      AliceInp << endl;
      AliceInp << "*Generated from call: SetCut('DCUTE',cut);";
      AliceInp << endl;
      AliceInp << setw(10) << "EMFCUT    ";
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << -fCutValue[i];
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
      AliceInp << setw(10) << 0.0;
      AliceInp << setw(10) << 0.0;
      AliceInp << setw(10) << 3.0;
      AliceInp << setprecision(2);
      AliceInp << setw(10) << fLastMaterial;
      AliceInp << setprecision(1);
      AliceInp << setw(10) << 1.0;
      AliceInp << endl;
    }
    
    //
    // time of flight cut in seconds
    // G4 particles: all
    // G3 default value: 0.01 GeV
    //gMC ->SetCut("TOFMAX",tofmax); // time of flight cuts in seconds
    else if (strncmp(&sCutFlag[i][0],"TOFMAX",6) == 0) {
      AliceInp << "*"; 
      AliceInp << endl;
      AliceInp << "*Time of flight cuts in seconds"; 
      AliceInp << endl;
      AliceInp << "*Generated from call: SetCut('TOFMAX',tofmax);"; 
      AliceInp << endl;
      AliceInp << setw(10) << "TIME-CUT  "; 
      AliceInp << setiosflags(ios::scientific) << setprecision(5);
      AliceInp << setw(10) << fCutValue[i]*1.e9;
      AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint) << setprecision(1);
      AliceInp << setw(10) << 0.0;
      AliceInp << setw(10) << 0.0;
      AliceInp << setw(10) << -6.0; // lower bound of the particle numbers for which the transport time cut-off and/or the start signal is to be applied
      AliceInp << setprecision(2);
      AliceInp << setw(10) << 64.0; // upper bound of the particle numbers for which the transport time cut-off and/or the start signal is to be applied
      AliceInp << setprecision(1);
      AliceInp << setw(10) << 1.0; // step length in assigning numbers
      AliceInp << endl;
    }

    else {
      cout << "SetCut for flag=" << &sCutFlag[i][0] << " value=" << fCutValue[i] << " not yet implemented!" << endl;
    }
  } //end of loop over SeCut calls
    
// Add START and STOP card
   AliceInp << setw(10) << "START     "; 
   AliceInp << setiosflags(ios::fixed) << setiosflags(ios::showpoint);
   AliceInp << setw(10) << fEventsPerRun;
   AliceInp << endl;
   AliceInp << setw(10) << "STOP      "; 
   AliceInp << endl;

} // end of InitPhysics


//______________________________________________________________________________ 
void TFluka::SetMaxStep(Double_t)
{
// SetMaxStep is dummy procedure in TFluka !
  if (fVerbosityLevel >=3)
  cout << "SetMaxStep is dummy procedure in TFluka !" << endl;
}

//______________________________________________________________________________ 
void TFluka::SetMaxNStep(Int_t)
{
// SetMaxNStep is dummy procedure in TFluka !
  if (fVerbosityLevel >=3)
  cout << "SetMaxNStep is dummy procedure in TFluka !" << endl;
}

//______________________________________________________________________________ 
void TFluka::SetUserDecay(Int_t)
{
// SetUserDecay is dummy procedure in TFluka !
  if (fVerbosityLevel >=3)
  cout << "SetUserDecay is dummy procedure in TFluka !" << endl;
}

//
// dynamic properties
//
//______________________________________________________________________________ 
void TFluka::TrackPosition(TLorentzVector& position) const
{
// Return the current position in the master reference frame of the
// track being transported
// TRACKR.atrack = age of the particle
// TRACKR.xtrack = x-position of the last point
// TRACKR.ytrack = y-position of the last point
// TRACKR.ztrack = z-position of the last point
  Int_t caller = GetCaller();
  if (caller == 1 || caller == 3 || caller == 6 || caller == 11 || caller == 12) { //bxdraw,endraw,usdraw
    position.SetX(GetXsco());
    position.SetY(GetYsco());
    position.SetZ(GetZsco());
    position.SetT(TRACKR.atrack);
  }
  else if (caller == 4) { // mgdraw
    position.SetX(TRACKR.xtrack[TRACKR.ntrack]);
    position.SetY(TRACKR.ytrack[TRACKR.ntrack]);
    position.SetZ(TRACKR.ztrack[TRACKR.ntrack]);
    position.SetT(TRACKR.atrack);
  }
  else if (caller == 5) { // sodraw
    position.SetX(TRACKR.xtrack[TRACKR.ntrack]);
    position.SetY(TRACKR.ytrack[TRACKR.ntrack]);
    position.SetZ(TRACKR.ztrack[TRACKR.ntrack]);
    position.SetT(0);
  }
  else
    Warning("TrackPosition","position not available");
}

//______________________________________________________________________________ 
void TFluka::TrackPosition(Double_t& x, Double_t& y, Double_t& z) const
{
// Return the current position in the master reference frame of the
// track being transported
// TRACKR.atrack = age of the particle
// TRACKR.xtrack = x-position of the last point
// TRACKR.ytrack = y-position of the last point
// TRACKR.ztrack = z-position of the last point
  Int_t caller = GetCaller();
  if (caller == 1 || caller == 3 || caller == 6 || caller == 11 || caller == 12) { //bxdraw,endraw,usdraw
    x = GetXsco();
    y = GetYsco();
    z = GetZsco();
  }
  else if (caller == 4) { // mgdraw
    x = TRACKR.xtrack[TRACKR.ntrack];
    y = TRACKR.ytrack[TRACKR.ntrack];
    z = TRACKR.ztrack[TRACKR.ntrack];
  }
  else if (caller == 5) { // sodraw
    x = TRACKR.xtrack[TRACKR.ntrack];
    y = TRACKR.ytrack[TRACKR.ntrack];
    z = TRACKR.ztrack[TRACKR.ntrack];
  }
  else
    Warning("TrackPosition","position not available");
}

//______________________________________________________________________________ 
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
  Int_t caller = GetCaller();
  if (caller != 2) { // not eedraw 
    if (TRACKR.ptrack >= 0) {
      momentum.SetPx(TRACKR.ptrack*TRACKR.cxtrck);
      momentum.SetPy(TRACKR.ptrack*TRACKR.cytrck);
      momentum.SetPz(TRACKR.ptrack*TRACKR.cztrck);
      momentum.SetE(TRACKR.etrack);
      return;
    }
    else {
      Double_t p = sqrt(TRACKR.etrack*TRACKR.etrack - PAPROP.am[TRACKR.jtrack+6]*PAPROP.am[TRACKR.jtrack+6]);
      momentum.SetPx(p*TRACKR.cxtrck);
      momentum.SetPy(p*TRACKR.cytrck);
      momentum.SetPz(p*TRACKR.cztrck);
      momentum.SetE(TRACKR.etrack);
      return;
    }
  }
  else
    Warning("TrackMomentum","momentum not available");
}

//______________________________________________________________________________ 
void TFluka::TrackMomentum(Double_t& px, Double_t& py, Double_t& pz, Double_t& e) const
{
// Return the direction and the momentum (GeV/c) of the track
// currently being transported
// TRACKR.ptrack = momentum of the particle (not always defined, if
//               < 0 must be obtained from etrack) 
// TRACKR.cx,y,ztrck = direction cosines of the current particle
// TRACKR.etrack = total energy of the particle
// TRACKR.jtrack = identity number of the particle
// PAPROP.am[TRACKR.jtrack] = particle mass in gev
  Int_t caller = GetCaller();
  if (caller != 2) { // not eedraw 
    if (TRACKR.ptrack >= 0) {
      px = TRACKR.ptrack*TRACKR.cxtrck;
      py = TRACKR.ptrack*TRACKR.cytrck;
      pz = TRACKR.ptrack*TRACKR.cztrck;
      e = TRACKR.etrack;
      return;
    }
    else {
      Double_t p = sqrt(TRACKR.etrack*TRACKR.etrack - PAPROP.am[TRACKR.jtrack+6]*PAPROP.am[TRACKR.jtrack+6]);
      px = p*TRACKR.cxtrck;
      py = p*TRACKR.cytrck;
      pz = p*TRACKR.cztrck;
      e = TRACKR.etrack;
      return;
    }
  }
  else
    Warning("TrackMomentum","momentum not available");
}

//______________________________________________________________________________ 
Double_t TFluka::TrackStep() const
{
// Return the length in centimeters of the current step
// TRACKR.ctrack = total curved path
  Int_t caller = GetCaller();
  if (caller == 1 || caller == 3 || caller == 6) //bxdraw,endraw,usdraw
    return 0.0;
  else if (caller == 4) //mgdraw
    return TRACKR.ctrack;
  else
    return -1.0;
}

//______________________________________________________________________________ 
Double_t TFluka::TrackLength() const
{
// TRACKR.cmtrck = cumulative curved path since particle birth
  Int_t caller = GetCaller();
  if (caller == 1 || caller == 3 || caller == 4 || caller == 6) //bxdraw,endraw,mgdraw,usdraw
    return TRACKR.cmtrck;
  else 
    return -1.0;
}

//______________________________________________________________________________ 
Double_t TFluka::TrackTime() const
{
// Return the current time of flight of the track being transported
// TRACKR.atrack = age of the particle
  Int_t caller = GetCaller();
  if (caller == 1 || caller == 3 || caller == 4 || caller == 6) //bxdraw,endraw,mgdraw,usdraw
    return TRACKR.atrack;
  else 
    return -1;
}

//______________________________________________________________________________ 
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
  Double_t sum = 0;
  for ( Int_t j=0;j<TRACKR.mtrack;j++) {
    sum +=TRACKR.dtrack[j];  
  }
  if (TRACKR.ntrack == 0 && TRACKR.mtrack == 0)
    return fRull + sum;
  else {
    return sum;
  }
}

//______________________________________________________________________________ 
Int_t TFluka::TrackPid() const
{
// Return the id of the particle transported
// TRACKR.jtrack = identity number of the particle
  Int_t caller = GetCaller();
  if (caller != 2)  // not eedraw 
    return PDGFromId(TRACKR.jtrack);
  else
    return -1000;
}

//______________________________________________________________________________ 
Double_t TFluka::TrackCharge() const
{
// Return charge of the track currently transported
// PAPROP.ichrge = electric charge of the particle
// TRACKR.jtrack = identity number of the particle
  Int_t caller = GetCaller();
  if (caller != 2)  // not eedraw 
    return PAPROP.ichrge[TRACKR.jtrack+6];
  else
    return -1000.0;
}

//______________________________________________________________________________ 
Double_t TFluka::TrackMass() const
{
// PAPROP.am = particle mass in GeV
// TRACKR.jtrack = identity number of the particle
  Int_t caller = GetCaller();
  if (caller != 2)  // not eedraw 
    return PAPROP.am[TRACKR.jtrack+6];
  else
    return -1000.0;
}

//______________________________________________________________________________ 
Double_t TFluka::Etot() const
{
// TRACKR.etrack = total energy of the particle
  Int_t caller = GetCaller();
  if (caller != 2)  // not eedraw
    return TRACKR.etrack;
  else
    return -1000.0;
}

//
// track status
//
//______________________________________________________________________________ 
Bool_t   TFluka::IsNewTrack() const
{
// True if the track has positive cummulative length
  Int_t caller = GetCaller();
  if (caller != 2) { // not eedraw
    if (TRACKR.cmtrck > 0.0)
      return 1; 
    else
      return 0; 
  }
  else
    return 0;
}

//______________________________________________________________________________ 
Bool_t   TFluka::IsTrackInside() const
{
// True if the track is not at the boundary of the current volume
// In Fluka a step is always inside one kind of material
// If the step would go behind the region of one material,
// it will be shortened to reach only the boundary.
// Therefore IsTrackInside() is always true.
  Int_t caller = GetCaller();
  if (caller == 1)  // bxdraw
    return 0;
  else
    return 1;
}

//______________________________________________________________________________ 
Bool_t   TFluka::IsTrackEntering() const
{
// True if this is the first step of the track in the current volume

  Int_t caller = GetCaller();
  if (caller == 11)  // bxdraw entering
    return 1;
  else return 0;
}

//______________________________________________________________________________ 
Bool_t   TFluka::IsTrackExiting() const
{
  Int_t caller = GetCaller();
  if (caller == 12)  // bxdraw exiting
    return 1;
  else return 0;
}

//______________________________________________________________________________ 
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

//______________________________________________________________________________ 
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

//______________________________________________________________________________ 
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

//______________________________________________________________________________ 
Bool_t   TFluka::IsTrackAlive() const
{
// means not disappeared or not out
  if (IsTrackDisappeared() || IsTrackOut() ) return 0;
  else return 1;
}

//
// secondaries
//

//______________________________________________________________________________ 
Int_t TFluka::NSecondaries() const
// Number of secondary particles generated in the current step
// FINUC.np = number of secondaries except light and heavy ions
// FHEAVY.npheav = number of secondaries for light and heavy secondary ions
{
  Int_t caller = GetCaller();
  if (caller == 6)  // valid only after usdraw
    return FINUC.np + FHEAVY.npheav;
  else
    return 0;
} // end of NSecondaries

//______________________________________________________________________________ 
void TFluka::GetSecondary(Int_t isec, Int_t& particleId,
		TLorentzVector& position, TLorentzVector& momentum)
{
  Int_t caller = GetCaller();
  if (caller == 6) {  // valid only after usdraw
    if (isec >= 0 && isec < FINUC.np) {
      particleId = PDGFromId(FINUC.kpart[isec]);
      position.SetX(fXsco);
      position.SetY(fYsco);
      position.SetZ(fZsco);
      position.SetT(TRACKR.atrack);
//    position.SetT(TRACKR.atrack+FINUC.agesec[isec]); //not yet implem.
      momentum.SetPx(FINUC.plr[isec]*FINUC.cxr[isec]);
      momentum.SetPy(FINUC.plr[isec]*FINUC.cyr[isec]);
      momentum.SetPz(FINUC.plr[isec]*FINUC.czr[isec]);
      momentum.SetE(FINUC.tki[isec] + PAPROP.am[FINUC.kpart[isec]+6]);
    }
    else if (isec >= FINUC.np && isec < FINUC.np + FHEAVY.npheav) {
      Int_t jsec = isec - FINUC.np;
      particleId = FHEAVY.kheavy[jsec]; // this is Fluka id !!!
      position.SetX(fXsco);
      position.SetY(fYsco);
      position.SetZ(fZsco);
      position.SetT(TRACKR.atrack);
//    position.SetT(TRACKR.atrack+FHEAVY.agheav[jsec]); //not yet implem.
      momentum.SetPx(FHEAVY.pheavy[jsec]*FHEAVY.cxheav[jsec]);
      momentum.SetPy(FHEAVY.pheavy[jsec]*FHEAVY.cyheav[jsec]);
      momentum.SetPz(FHEAVY.pheavy[jsec]*FHEAVY.czheav[jsec]);
      if (FHEAVY.tkheav[jsec] >= 3 && FHEAVY.tkheav[jsec] <= 6) 
        momentum.SetE(FHEAVY.tkheav[jsec] + PAPROP.am[jsec+6]);
      else if (FHEAVY.tkheav[jsec] > 6)
        momentum.SetE(FHEAVY.tkheav[jsec] + FHEAVY.amnhea[jsec]); // to be checked !!!
    }
    else
      Warning("GetSecondary","isec out of range");
  }
  else
    Warning("GetSecondary","no secondaries available");
} // end of GetSecondary

//______________________________________________________________________________ 
TMCProcess TFluka::ProdProcess(Int_t) const
// Name of the process that has produced the secondary particles
// in the current step
{
    const TMCProcess kIpNoProc = kPNoProcess;
    const TMCProcess kIpPDecay = kPDecay;
    const TMCProcess kIpPPair = kPPair;
//  const TMCProcess kIpPPairFromPhoton = kPPairFromPhoton;
//  const TMCProcess kIpPPairFromVirtualPhoton = kPPairFromVirtualPhoton;
    const TMCProcess kIpPCompton = kPCompton;
    const TMCProcess kIpPPhotoelectric = kPPhotoelectric;
    const TMCProcess kIpPBrem = kPBrem;
//  const TMCProcess kIpPBremFromHeavy = kPBremFromHeavy;
//  const TMCProcess kIpPBremFromElectronOrPositron = kPBremFromElectronOrPositron;
    const TMCProcess kIpPDeltaRay = kPDeltaRay;
//  const TMCProcess kIpPMoller = kPMoller;
//  const TMCProcess kIpPBhabha = kPBhabha;
    const TMCProcess kIpPAnnihilation = kPAnnihilation;
//  const TMCProcess kIpPAnnihilInFlight = kPAnnihilInFlight;
//  const TMCProcess kIpPAnnihilAtRest = kPAnnihilAtRest;
    const TMCProcess kIpPHadronic = kPHadronic;
    const TMCProcess kIpPMuonNuclear = kPMuonNuclear;
    const TMCProcess kIpPPhotoFission = kPPhotoFission;
    const TMCProcess kIpPRayleigh = kPRayleigh;
//  const TMCProcess kIpPCerenkov = kPCerenkov;
//  const TMCProcess kIpPSynchrotron = kPSynchrotron;

    Int_t mugamma = TRACKR.jtrack == 7 || TRACKR.jtrack == 10 || TRACKR.jtrack == 11;
    if (fIcode == 102) return kIpPDecay;
    else if (fIcode == 104 || fIcode == 217) return kIpPPair;
//  else if (fIcode == 104) return kIpPairFromPhoton;
//  else if (fIcode == 217) return kIpPPairFromVirtualPhoton;
    else if (fIcode == 219) return kIpPCompton;
    else if (fIcode == 221) return kIpPPhotoelectric;
    else if (fIcode == 105 || fIcode == 208) return kIpPBrem;
//  else if (fIcode == 105) return kIpPBremFromHeavy;
//  else if (fIcode == 208) return kPBremFromElectronOrPositron;
    else if (fIcode == 103 || fIcode == 400) return kIpPDeltaRay;
    else if (fIcode == 210 || fIcode == 212) return kIpPDeltaRay;
//  else if (fIcode == 210) return kIpPMoller;
//  else if (fIcode == 212) return kIpPBhabha;
    else if (fIcode == 214 || fIcode == 215) return kIpPAnnihilation;
//  else if (fIcode == 214) return kIpPAnnihilInFlight;
//  else if (fIcode == 215) return kIpPAnnihilAtRest;
    else if (fIcode == 101) return kIpPHadronic;
    else if (fIcode == 101) {
      if (!mugamma) return kIpPHadronic;
      else if (TRACKR.jtrack == 7) return kIpPPhotoFission;
      else return kIpPMuonNuclear;
    }
    else if (fIcode == 225) return kIpPRayleigh;
// Fluka codes 100, 300 and 400 still to be investigasted
    else return kIpNoProc;
}

//Int_t StepProcesses(TArrayI &proc) const
// Return processes active in the current step
//{
//ck = total energy of the particl ???????????????? 
//}


//______________________________________________________________________________ 
Int_t TFluka::VolId2Mate(Int_t id) const
{
//
// Returns the material number for a given volume ID
//
   return fGeom->VolId2Mate(id);
}

//______________________________________________________________________________ 
const char* TFluka::VolName(Int_t id) const
{
//
// Returns the volume name for a given volume ID
//
   return fGeom->VolName(id);
}

//______________________________________________________________________________ 
Int_t TFluka::VolId(const Text_t* volName) const
{
//
// Converts from volume name to volume ID.
// Time consuming. (Only used during set-up)
// Could be replaced by hash-table
//
   return fGeom->VolId(volName);
}

//______________________________________________________________________________ 
Int_t TFluka::CurrentVolID(Int_t& copyNo) const
{
//
// Return the logical id and copy number corresponding to the current fluka region
//
   return fGeom->CurrentVolID(copyNo);
} 

//______________________________________________________________________________ 
Int_t TFluka::CurrentVolOffID(Int_t off, Int_t& copyNo) const
{
//
// Return the logical id and copy number of off'th mother 
// corresponding to the current fluka region
//
   return fGeom->CurrentVolOffID(off, copyNo);
}

//______________________________________________________________________________ 
const char* TFluka::CurrentVolName() const
{
//
// Return the current volume name
//
   return fGeom->CurrentVolName();
}

//______________________________________________________________________________ 
const char* TFluka::CurrentVolOffName(Int_t off) const
{
//
// Return the volume name of the off'th mother of the current volume
//
   return fGeom->CurrentVolOffName(off);
}

//______________________________________________________________________________ 
Int_t TFluka::CurrentMaterial(Float_t & /*a*/, Float_t & /*z*/, 
		      Float_t & /*dens*/, Float_t & /*radl*/, Float_t & /*absl*/) const
{
//
//  Return the current medium number  ??? what about material properties
//
  Int_t copy;
  Int_t id  =  TFluka::CurrentVolID(copy);
  Int_t med =  TFluka::VolId2Mate(id);
  return med;
}

//______________________________________________________________________________ 
void TFluka::Gmtod(Float_t* xm, Float_t* xd, Int_t iflag)
{
// Transforms a position from the world reference frame
// to the current volume reference frame.
//
//  Geant3 desription:
//  ==================
//       Computes coordinates XD (in DRS) 
//       from known coordinates XM in MRS 
//       The local reference system can be initialized by
//         - the tracking routines and GMTOD used in GUSTEP
//         - a call to GMEDIA(XM,NUMED)
//         - a call to GLVOLU(NLEVEL,NAMES,NUMBER,IER) 
//             (inverse routine is GDTOM) 
//
//        If IFLAG=1  convert coordinates 
//           IFLAG=2  convert direction cosinus
//
// ---
   fGeom->Gmtod(xm,xd,iflag);   
}
  
//______________________________________________________________________________ 
void TFluka::Gmtod(Double_t* xm, Double_t* xd, Int_t iflag)
{
// Transforms a position from the world reference frame
// to the current volume reference frame.
//
//  Geant3 desription:
//  ==================
//       Computes coordinates XD (in DRS) 
//       from known coordinates XM in MRS 
//       The local reference system can be initialized by
//         - the tracking routines and GMTOD used in GUSTEP
//         - a call to GMEDIA(XM,NUMED)
//         - a call to GLVOLU(NLEVEL,NAMES,NUMBER,IER) 
//             (inverse routine is GDTOM) 
//
//        If IFLAG=1  convert coordinates 
//           IFLAG=2  convert direction cosinus
//
// ---
   fGeom->Gmtod(xm,xd,iflag);   
}

//______________________________________________________________________________ 
void TFluka::Gdtom(Float_t* xd, Float_t* xm, Int_t iflag)
{
// Transforms a position from the current volume reference frame
// to the world reference frame.
//
//  Geant3 desription:
//  ==================
//  Computes coordinates XM (Master Reference System
//  knowing the coordinates XD (Detector Ref System)
//  The local reference system can be initialized by
//    - the tracking routines and GDTOM used in GUSTEP
//    - a call to GSCMED(NLEVEL,NAMES,NUMBER)
//        (inverse routine is GMTOD)
// 
//   If IFLAG=1  convert coordinates
//      IFLAG=2  convert direction cosinus
//
// ---
   fGeom->Gdtom(xd,xm,iflag);   
}

//______________________________________________________________________________ 
void TFluka::Gdtom(Double_t* xd, Double_t* xm, Int_t iflag)
{
// Transforms a position from the current volume reference frame
// to the world reference frame.
//
//  Geant3 desription:
//  ==================
//  Computes coordinates XM (Master Reference System
//  knowing the coordinates XD (Detector Ref System)
//  The local reference system can be initialized by
//    - the tracking routines and GDTOM used in GUSTEP
//    - a call to GSCMED(NLEVEL,NAMES,NUMBER)
//        (inverse routine is GMTOD)
// 
//   If IFLAG=1  convert coordinates
//      IFLAG=2  convert direction cosinus
//
// ---
   fGeom->Gdtom(xd,xm,iflag);   
}

// ===============================================================
void TFluka::FutoTest() 
{
    Int_t icode, mreg, newreg, particleId;
    Double_t rull, xsco, ysco, zsco;
    TLorentzVector position, momentum;
    icode = GetIcode();
    if (icode == 0) {
	if (fVerbosityLevel >=3)
	    cout << " icode=" << icode << endl;
    } else if (icode > 0 && icode <= 5) {
// mgdraw
	mreg = GetMreg();
	if (fVerbosityLevel >=3)
	    cout << " icode=" << icode
		 << " mreg=" << mreg
		 << endl;
	TrackPosition(position);
	TrackMomentum(momentum);
	if (fVerbosityLevel >=3) {
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
	
	Float_t x = position.X();
	Float_t y = position.Y();
	Float_t z = position.Z();
	Float_t xm[3];
	Float_t xd[3];
	xm[0] = x; xm[1] = y; xm[2] = z;
	if (fVerbosityLevel >= 3)
	    printf("Global trackPosition: %f %f %f \n", x, y, z);
	Gmtod(xm, xd, 1);
	if (fVerbosityLevel >= 3)
	    printf("Local trackPosition: %f %f %f \n", xd[0], xd[1], xd[2]);
	Gdtom(xd, xm, 1);
	if (fVerbosityLevel >= 3)
	    printf("New trackPosition: %f %f %f \n", xm[0], xm[1], xm[2]);
    } else if((icode >= 10 && icode <= 15) ||
	      (icode >= 20 && icode <= 24) ||
	      (icode >= 30 && icode <= 33) ||
	      (icode >= 40 && icode <= 41) ||
	      (icode >= 50 && icode <= 52)) {
// endraw
	mreg = GetMreg();
	rull = GetRull();
	xsco = GetXsco();
	ysco = GetYsco();
	zsco = GetZsco();
	
	if (fVerbosityLevel >=3) {     
	    cout << " icode=" << icode
		 << " mreg=" << mreg
		 << " rull=" << rull
		 << " xsco=" << xsco
		 << " ysco=" << ysco
		 << " zsco=" << zsco << endl;
	}
	TrackPosition(position);
	TrackMomentum(momentum);
	if (fVerbosityLevel >=3) {
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
    } else if((icode >= 100 && icode <= 105) ||
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
	  xsco = GetXsco();
	  ysco = GetYsco();
	  zsco = GetZsco();
	  
	  if (fVerbosityLevel >=3) {
	      cout << " icode=" << icode
		   << " mreg=" << mreg
		   << " xsco=" << xsco
		   << " ysco=" << ysco
		   << " zsco=" << zsco << endl;
	      cout << "TrackPid=" << TrackPid() << endl;
	      cout << "NSecondaries=" << NSecondaries() << endl;
	  }
	  
	  for (Int_t isec=0; isec< NSecondaries(); isec++) {
	      TFluka::GetSecondary(isec, particleId, position, momentum);
	      if (fVerbosityLevel >=3) {
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
    } else if((icode == 19) ||
	      (icode == 29) ||
	      (icode == 39) ||
	      (icode == 49) ||
	      (icode == 59)) {
	mreg = GetMreg();
	newreg = GetNewreg();
	xsco = GetXsco();
	ysco = GetYsco();
	zsco = GetZsco();
	if (fVerbosityLevel >=3) {
	    cout << " icode=" << icode
		 << " mreg=" << mreg
		 << " newreg=" << newreg
		 << " xsco=" << xsco
		 << " ysco=" << ysco
		 << " zsco=" << zsco << endl;
	}
    }
} // end of FutoTest

