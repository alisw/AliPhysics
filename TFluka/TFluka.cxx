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

//
// Realisation of the TVirtualMC interface for the FLUKA code
// (See official web side http://www.fluka.org/).
//
// This implementation makes use of the TGeo geometry modeller.
// User configuration is via automatic generation of FLUKA input cards.
//
// Authors:
// A. Fasso
// E. Futo
// A. Gheata
// A. Morsch
//

#include <Riostream.h>

#include "TFluka.h"
#include "TCallf77.h"      //For the fortran calls
#include "Fdblprc.h"       //(DBLPRC) fluka common
#include "Fepisor.h"       //(EPISOR) fluka common
#include "Ffinuc.h"        //(FINUC)  fluka common
#include "Fiounit.h"       //(IOUNIT) fluka common
#include "Fpaprop.h"       //(PAPROP) fluka common
#include "Fpart.h"         //(PART)   fluka common
#include "Ftrackr.h"       //(TRACKR) fluka common
#include "Fpaprop.h"       //(PAPROP) fluka common
#include "Ffheavy.h"       //(FHEAVY) fluka common
#include "Fopphst.h"       //(OPPHST) fluka common
#include "Fstack.h"        //(STACK)  fluka common
#include "Fstepsz.h"       //(STEPSZ) fluka common
#include "Fopphst.h"       //(OPPHST) fluka common

#include "TVirtualMC.h"
#include "TMCProcess.h"
#include "TGeoManager.h"
#include "TGeoMaterial.h"
#include "TGeoMedium.h"
#include "TFlukaMCGeometry.h"
#include "TGeoMCGeometry.h"
#include "TFlukaCerenkov.h"
#include "TFlukaConfigOption.h"
#include "TFlukaScoringOption.h"
#include "TLorentzVector.h"
#include "TArrayI.h"

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
   fInputFileName(""),
   fProcesses(0), 
   fCuts(0),
   fUserScore(0)
{ 
  //
  // Default constructor
  //
   fGeneratePemf = kFALSE;
   fNVolumes = 0;
   fCurrentFlukaRegion = -1;
   fGeom = 0;
   fMCGeo = 0;
   fMaterials = 0;
   fDummyBoundary = 0;
   fFieldFlag = 1;
   fStopped   = 0;
   fStopEvent = 0;
   fStopRun   = 0;
   fNEvent    = 0;
} 
 
//______________________________________________________________________________ 
TFluka::TFluka(const char *title, Int_t verbosity, Bool_t isRootGeometrySupported)
  :TVirtualMC("TFluka",title, isRootGeometrySupported),
   fVerbosityLevel(verbosity),
   fInputFileName(""),
   fTrackIsEntering(0),
   fTrackIsExiting(0),
   fTrackIsNew(0),
   fProcesses(new TObjArray(100)),
   fCuts(new TObjArray(100)), 
   fUserScore(new TObjArray(100)) 
{
  // create geometry interface
   if (fVerbosityLevel >=3)
       cout << "<== TFluka::TFluka(" << title << ") constructor called." << endl;
   SetCoreInputFileName();
   SetInputFileName();
   SetGeneratePemf(kFALSE);
   fNVolumes      = 0;
   fCurrentFlukaRegion = -1;
   fDummyBoundary = 0;
   fFieldFlag = 1;
   fGeneratePemf = kFALSE;
   fMCGeo = new TGeoMCGeometry("MCGeo", "TGeo Implementation of VirtualMCGeometry", kTRUE);
   fGeom = new TFlukaMCGeometry("geom", "FLUKA VMC Geometry");
   if (verbosity > 2) fGeom->SetDebugMode(kTRUE);
   fMaterials = 0;
   fStopped   = 0;
   fStopEvent = 0;
   fStopRun   = 0;
   fNEvent    = 0;
}

//______________________________________________________________________________ 
TFluka::~TFluka() {
// Destructor
    if (fVerbosityLevel >=3)
	cout << "<== TFluka::~TFluka() destructor called." << endl;
    
    delete fGeom;
    delete fMCGeo;
    
    if (fCuts) {
	fCuts->Delete();
	delete fCuts;
    }

    if (fProcesses) {
	fProcesses->Delete();
	delete fProcesses;
    }


}

//
//______________________________________________________________________________
// TFluka control methods
//______________________________________________________________________________ 
void TFluka::Init() {
//
//  Geometry initialisation
//
    if (fVerbosityLevel >=3) cout << "==> TFluka::Init() called." << endl;
    
    if (!gGeoManager) new TGeoManager("geom", "FLUKA geometry");
    fApplication->ConstructGeometry();
    TGeoVolume *top = (TGeoVolume*)gGeoManager->GetListOfVolumes()->First();
    gGeoManager->SetTopVolume(top);
    gGeoManager->CloseGeometry("di");
    gGeoManager->DefaultColors();  // to be removed
    fNVolumes = fGeom->NofVolumes();
    fGeom->CreateFlukaMatFile("flukaMat.inp");   
    if (fVerbosityLevel >=3) {
       printf("== Number of volumes: %i\n ==", fNVolumes);
       cout << "\t* InitPhysics() - Prepare input file to be called" << endl; 
    }   
    // now we have TGeo geometry created and we have to patch FlukaVmc.inp
    // with the material mapping file FlukaMat.inp
}


//______________________________________________________________________________ 
void TFluka::FinishGeometry() {
//
// Build-up table with region to medium correspondance
//
  if (fVerbosityLevel >=3) {
    cout << "==> TFluka::FinishGeometry() called." << endl;
    printf("----FinishGeometry - nothing to do with TGeo\n");
    cout << "<== TFluka::FinishGeometry() called." << endl;
  }  
} 

//______________________________________________________________________________ 
void TFluka::BuildPhysics() {
//
//  Prepare FLUKA input files and call FLUKA physics initialisation
//
    
    if (fVerbosityLevel >=3)
	cout << "==> TFluka::BuildPhysics() called." << endl;
// Prepare input file with the current physics settings
    InitPhysics(); 
    cout << "\t* InitPhysics() - Prepare input file was called" << endl; 
    
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
    
    FinishGeometry();
    
    if (fVerbosityLevel >=3)
	cout << "<== TFluka::Init() called." << endl;
    
    
    if (fVerbosityLevel >=3)
	cout << "<== TFluka::BuildPhysics() called." << endl;
}  

//______________________________________________________________________________ 
void TFluka::ProcessEvent() {
//
// Process one event
//
    if (fStopRun) {
	printf("User Run Abortion: No more events handled !\n");
	fNEvent += 1;
	return;
    }

    if (fVerbosityLevel >=3)
	cout << "==> TFluka::ProcessEvent() called." << endl;
    fApplication->GeneratePrimaries();
    EPISOR.lsouit = true;
    flukam(1);
    if (fVerbosityLevel >=3)
	cout << "<== TFluka::ProcessEvent() called." << endl;
    //
    // Increase event number
    //
    fNEvent += 1;
}

//______________________________________________________________________________ 
Bool_t TFluka::ProcessRun(Int_t nevent) {
//
// Run steering
//

  if (fVerbosityLevel >=3)
    cout << "==> TFluka::ProcessRun(" << nevent << ") called." 
	 << endl;

  if (fVerbosityLevel >=2) {
    cout << "\t* GLOBAL.fdrtr = " << (GLOBAL.lfdrtr?'T':'F') << endl;
    cout << "\t* Calling flukam again..." << endl;
  }

  fApplication->InitGeometry();
  Int_t todo = TMath::Abs(nevent);
  for (Int_t ev = 0; ev < todo; ev++) {
      fApplication->BeginEvent();
      ProcessEvent();
      fApplication->FinishEvent();
  }

  if (fVerbosityLevel >=3)
    cout << "<== TFluka::ProcessRun(" << nevent << ") called." 
	 << endl;
  return kTRUE;
}

//_____________________________________________________________________________
// methods for building/management of geometry

// functions from GCONS 
//____________________________________________________________________________ 
void TFluka::Gfmate(Int_t imat, char *name, Float_t &a, Float_t &z,  
		    Float_t &dens, Float_t &radl, Float_t &absl,
		    Float_t* /*ubuf*/, Int_t& /*nbuf*/) {
//
   TGeoMaterial *mat;
   TIter next (gGeoManager->GetListOfMaterials());
   while ((mat = (TGeoMaterial*)next())) {
     if (mat->GetUniqueID() == (UInt_t)imat) break;
   }
   if (!mat) {
      Error("Gfmate", "no material with index %i found", imat);
      return;
   }
   sprintf(name, "%s", mat->GetName());
   a = mat->GetA();
   z = mat->GetZ();
   dens = mat->GetDensity();
   radl = mat->GetRadLen();
   absl = mat->GetIntLen();
} 

//______________________________________________________________________________ 
void TFluka::Gfmate(Int_t imat, char *name, Double_t &a, Double_t &z,  
		    Double_t &dens, Double_t &radl, Double_t &absl,
		    Double_t* /*ubuf*/, Int_t& /*nbuf*/) {
//
   TGeoMaterial *mat;
   TIter next (gGeoManager->GetListOfMaterials());
   while ((mat = (TGeoMaterial*)next())) {
     if (mat->GetUniqueID() == (UInt_t)imat) break;
   }
   if (!mat) {
      Error("Gfmate", "no material with index %i found", imat);
      return;
   }
   sprintf(name, "%s", mat->GetName());
   a = mat->GetA();
   z = mat->GetZ();
   dens = mat->GetDensity();
   radl = mat->GetRadLen();
   absl = mat->GetIntLen();
} 

// detector composition
//______________________________________________________________________________ 
void TFluka::Material(Int_t& kmat, const char* name, Double_t a, 
		      Double_t z, Double_t dens, Double_t radl, Double_t absl,
		      Float_t* buf, Int_t nwbuf) {
//
   Double_t* dbuf = fGeom->CreateDoubleArray(buf, nwbuf);  
   Material(kmat, name, a, z, dens, radl, absl, dbuf, nwbuf);
   delete [] dbuf;
} 

//______________________________________________________________________________ 
void TFluka::Material(Int_t& kmat, const char* name, Double_t a, 
		      Double_t z, Double_t dens, Double_t radl, Double_t absl,
		      Double_t* /*buf*/, Int_t /*nwbuf*/) {
//
  TGeoMaterial *mat;
  kmat = gGeoManager->GetListOfMaterials()->GetSize();
  if ((z-Int_t(z)) > 1E-3) {
     mat = fGeom->GetMakeWrongMaterial(z);
     if (mat) {
        mat->SetRadLen(radl,absl);
        mat->SetUniqueID(kmat);
        return;
     }
  }      
  gGeoManager->Material(name, a, z, dens, kmat, radl, absl);
} 

//______________________________________________________________________________ 
void TFluka::Mixture(Int_t& kmat, const char *name, Float_t *a, 
		     Float_t *z, Double_t dens, Int_t nlmat, Float_t *wmat) {
//
  Double_t* da = fGeom->CreateDoubleArray(a, TMath::Abs(nlmat));  
  Double_t* dz = fGeom->CreateDoubleArray(z, TMath::Abs(nlmat));  
  Double_t* dwmat = fGeom->CreateDoubleArray(wmat, TMath::Abs(nlmat));  

  Mixture(kmat, name, da, dz, dens, nlmat, dwmat);
  for (Int_t i=0; i<nlmat; i++) {
    a[i] = da[i]; z[i] = dz[i]; wmat[i] = dwmat[i];
  }  

  delete [] da;
  delete [] dz;
  delete [] dwmat;
} 

//______________________________________________________________________________ 
void TFluka::Mixture(Int_t& kmat, const char *name, Double_t *a, 
		     Double_t *z, Double_t dens, Int_t nlmat, Double_t *wmat) {
//
  // Defines mixture OR COMPOUND IMAT as composed by 
  // THE BASIC NLMAT materials defined by arrays A,Z and WMAT
  // 
  // If NLMAT > 0 then wmat contains the proportion by
  // weights of each basic material in the mixture. 
  // 
  // If nlmat < 0 then WMAT contains the number of atoms 
  // of a given kind into the molecule of the COMPOUND
  // In this case, WMAT in output is changed to relative
  // weigths.
  //
  Int_t i,j;
  if (nlmat < 0) {
     nlmat = - nlmat;
     Double_t amol = 0;
     for (i=0;i<nlmat;i++) {
        amol += a[i]*wmat[i];
     }
     for (i=0;i<nlmat;i++) {
        wmat[i] *= a[i]/amol;
     }
  }
  kmat = gGeoManager->GetListOfMaterials()->GetSize();
  // Check if we have elements with fractional Z
  TGeoMaterial *mat = 0;
  TGeoMixture *mix = 0;
  Bool_t mixnew = kFALSE;
  for (i=0; i<nlmat; i++) {
     if (z[i]-Int_t(z[i]) < 1E-3) continue;
     // We have found an element with fractional Z -> loop mixtures to look for it
     for (j=0; j<kmat; j++) {
        mat = (TGeoMaterial*)gGeoManager->GetListOfMaterials()->At(j);
        if (!mat) break;
        if (!mat->IsMixture()) continue;
        mix = (TGeoMixture*)mat;
        if (TMath::Abs(z[i]-mix->GetZ()) >1E-3) continue;
//        printf(" FOUND component %i as mixture %s\n", i, mat->GetName());
        mixnew = kTRUE;
        break;
     }
     if (!mixnew) Warning("Mixture","%s : cannot find component %i with fractional Z=%f\n", name, i, z[i]);
     break;
  }   
  if (mixnew) {
     Int_t nlmatnew = nlmat+mix->GetNelements()-1;
     Double_t *anew = new Double_t[nlmatnew];
     Double_t *znew = new Double_t[nlmatnew];
     Double_t *wmatnew = new Double_t[nlmatnew];
     Int_t ind=0;
     for (j=0; j<nlmat; j++) {
        if (j==i) continue;
        anew[ind] = a[j];
        znew[ind] = z[j];
        wmatnew[ind] = wmat[j];
        ind++;
     }
     for (j=0; j<mix->GetNelements(); j++) {
        anew[ind] = mix->GetAmixt()[j];
        znew[ind] = mix->GetZmixt()[j];
        wmatnew[ind] = wmat[i]*mix->GetWmixt()[j];
        ind++;
     }
     Mixture(kmat, name, anew, znew, dens, nlmatnew, wmatnew);
     delete [] anew;
     delete [] znew;
     delete [] wmatnew;
     return;
  }   
  // Now we need to compact identical elements within the mixture
  // First check if this happens   
  mixnew = kFALSE;  
  for (i=0; i<nlmat-1; i++) {
     for (j=i+1; j<nlmat; j++) {
        if (z[i] == z[j]) {
           mixnew = kTRUE;
           break;
        }
     }   
     if (mixnew) break;
  }   
  if (mixnew) {
     Int_t nlmatnew = 0;
     Double_t *anew = new Double_t[nlmat];
     Double_t *znew = new Double_t[nlmat];
     memset(znew, 0, nlmat*sizeof(Double_t));
     Double_t *wmatnew = new Double_t[nlmat];
     Bool_t skipi;
     for (i=0; i<nlmat; i++) {
        skipi = kFALSE;
        for (j=0; j<nlmatnew; j++) {
           if (z[i] == z[j]) {
              wmatnew[j] += wmat[i];
              skipi = kTRUE;
              break;
           }
        }   
        if (skipi) continue;    
        anew[nlmatnew] = a[i];
        znew[nlmatnew] = z[i];
        wmatnew[nlmatnew] = wmat[i];
        nlmatnew++;
     }
     Mixture(kmat, name, anew, znew, dens, nlmatnew, wmatnew);
     delete [] anew;
     delete [] znew;
     delete [] wmatnew;
     return;     
   }
   gGeoManager->Mixture(name, a, z, dens, nlmat, wmat, kmat);
} 

//______________________________________________________________________________ 
void TFluka::Medium(Int_t& kmed, const char *name, Int_t nmat, 
		    Int_t isvol, Int_t ifield, Double_t fieldm, Double_t tmaxfd, 
		    Double_t stemax, Double_t deemax, Double_t epsil, 
		    Double_t stmin, Float_t* ubuf, Int_t nbuf) { 
  //
  kmed = gGeoManager->GetListOfMedia()->GetSize()+1;
  fMCGeo->Medium(kmed, name, nmat, isvol, ifield, fieldm, tmaxfd, stemax, deemax, 
	     epsil, stmin, ubuf, nbuf);
} 

//______________________________________________________________________________ 
void TFluka::Medium(Int_t& kmed, const char *name, Int_t nmat, 
		    Int_t isvol, Int_t ifield, Double_t fieldm, Double_t tmaxfd, 
		    Double_t stemax, Double_t deemax, Double_t epsil, 
		    Double_t stmin, Double_t* ubuf, Int_t nbuf) { 
  //
  kmed = gGeoManager->GetListOfMedia()->GetSize()+1;
  fMCGeo->Medium(kmed, name, nmat, isvol, ifield, fieldm, tmaxfd, stemax, deemax, 
	     epsil, stmin, ubuf, nbuf);
} 

//______________________________________________________________________________ 
void TFluka::Matrix(Int_t& krot, Double_t thetaX, Double_t phiX, 
		    Double_t thetaY, Double_t phiY, Double_t thetaZ, 
		    Double_t phiZ) {
//		     
  krot = gGeoManager->GetListOfMatrices()->GetEntriesFast();
  fMCGeo->Matrix(krot, thetaX, phiX, thetaY, phiY, thetaZ, phiZ); 
} 

//______________________________________________________________________________ 
void TFluka::Gstpar(Int_t itmed, const char* param, Double_t parval) {
//
//
// Check if material is used    
   if (fVerbosityLevel >=3) 
       printf("Gstpar called with %6d %5s %12.4e %6d\n", itmed, param, parval, fGeom->GetFlukaMaterial(itmed));
   Int_t* reglist;
   Int_t nreg;
   reglist = fGeom->GetMaterialList(itmed, nreg);
   if (nreg == 0) return;
//
   Bool_t process = kFALSE;
   if (strncmp(param, "DCAY",  4) == 0 ||
       strncmp(param, "PAIR",  4) == 0 ||
       strncmp(param, "COMP",  4) == 0 ||
       strncmp(param, "PHOT",  4) == 0 ||
       strncmp(param, "PFIS",  4) == 0 ||
       strncmp(param, "DRAY",  4) == 0 ||
       strncmp(param, "ANNI",  4) == 0 ||
       strncmp(param, "BREM",  4) == 0 ||
       strncmp(param, "MUNU",  4) == 0 ||
       strncmp(param, "CKOV",  4) == 0 ||
       strncmp(param, "HADR",  4) == 0 ||
       strncmp(param, "LOSS",  4) == 0 ||
       strncmp(param, "MULS",  4) == 0 ||
       strncmp(param, "RAYL",  4) == 0) 
   {
       process = kTRUE;
   } 
   if (process) {
       SetProcess(param, Int_t (parval), fGeom->GetFlukaMaterial(itmed));
   } else {
       SetCut(param, parval, fGeom->GetFlukaMaterial(itmed));
   }
}    

// functions from GGEOM 
//_____________________________________________________________________________
void TFluka::Gsatt(const char *name, const char *att, Int_t val)
{ 
  // Set visualisation attributes for one volume
  char vname[5];
  fGeom->Vname(name,vname);
  char vatt[5];
  fGeom->Vname(att,vatt);
  gGeoManager->SetVolumeAttribute(vname, vatt, val);
}

//______________________________________________________________________________ 
Int_t TFluka::Gsvolu(const char *name, const char *shape, Int_t nmed,  
		     Float_t *upar, Int_t np)  {
//
    return fMCGeo->Gsvolu(name, shape, nmed, upar, np); 
}

//______________________________________________________________________________ 
Int_t TFluka::Gsvolu(const char *name, const char *shape, Int_t nmed,  
		     Double_t *upar, Int_t np)  {
//
    return fMCGeo->Gsvolu(name, shape, nmed, upar, np); 
}
 
//______________________________________________________________________________ 
void TFluka::Gsdvn(const char *name, const char *mother, Int_t ndiv, 
		   Int_t iaxis) {
//
    fMCGeo->Gsdvn(name, mother, ndiv, iaxis); 
} 

//______________________________________________________________________________ 
void TFluka::Gsdvn2(const char *name, const char *mother, Int_t ndiv, 
		    Int_t iaxis, Double_t c0i, Int_t numed) {
//
    fMCGeo->Gsdvn2(name, mother, ndiv, iaxis, c0i, numed); 
} 

//______________________________________________________________________________ 
void TFluka::Gsdvt(const char *name, const char *mother, Double_t step, 
		   Int_t iaxis, Int_t numed, Int_t ndvmx) {
//	
    fMCGeo->Gsdvt(name, mother, step, iaxis, numed, ndvmx); 
} 

//______________________________________________________________________________ 
void TFluka::Gsdvt2(const char *name, const char *mother, Double_t step, 
		    Int_t iaxis, Double_t c0, Int_t numed, Int_t ndvmx) { 
//
    fMCGeo->Gsdvt2(name, mother, step, iaxis, c0, numed, ndvmx); 
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
  fMCGeo->Gspos(name, nr, mother, x, y, z, irot, konly); 
} 

//______________________________________________________________________________ 
void TFluka::Gsposp(const char *name, Int_t nr, const char *mother,  
		    Double_t x, Double_t y, Double_t z, Int_t irot,
		    const char *konly, Float_t *upar, Int_t np)  {
  //
  fMCGeo->Gsposp(name, nr, mother, x, y, z, irot, konly, upar, np); 
} 

//______________________________________________________________________________ 
void TFluka::Gsposp(const char *name, Int_t nr, const char *mother,  
		    Double_t x, Double_t y, Double_t z, Int_t irot,
		    const char *konly, Double_t *upar, Int_t np)  {
  //
  fMCGeo->Gsposp(name, nr, mother, x, y, z, irot, konly, upar, np); 
} 

//______________________________________________________________________________ 
void TFluka::Gsbool(const char* /*onlyVolName*/, const char* /*manyVolName*/) {
//
// Nothing to do with TGeo
}

//______________________________________________________________________________ 
void TFluka::SetCerenkov(Int_t itmed, Int_t npckov, Float_t* ppckov,
			 Float_t* absco, Float_t* effic, Float_t* rindex) {
//
// Set Cerenkov properties for medium itmed
//
// npckov: number of sampling points
// ppckov: energy values
// absco:  absorption length
// effic:  quantum efficiency
// rindex: refraction index
//
//
//  
//  Create object holding Cerenkov properties
//  
    TFlukaCerenkov* cerenkovProperties = new TFlukaCerenkov(npckov, ppckov, absco, effic, rindex);
//
//  Pass object to medium
    TGeoMedium* medium = gGeoManager->GetMedium(itmed);
    medium->SetCerenkovProperties(cerenkovProperties);
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
  //                      Alpha     He3       Triton    Deuteron  gen. ion  opt. photon   
    Int_t idSpecial[6] = {10020040, 10020030, 10010030, 10010020, 10000000, 50000050};
  // IPTOKP array goes from official to internal

    if (id == -1) {
// Cerenkov photon
	if (fVerbosityLevel >= 3)
	    printf("\n PDGFromId: Cerenkov Photon \n");
	return  50000050;
    }
// Error id    
    if (id == 0 || id < -6 || id > 250) {
	if (fVerbosityLevel >= 3)
	    printf("PDGFromId: Error id = 0\n");
	return -1;
    }
// Good id    
    if (id > 0) {
	Int_t intfluka = GetFlukaIPTOKP(id);
	if (intfluka == 0) {
	    if (fVerbosityLevel >= 3)
		printf("PDGFromId: Error intfluka = 0: %d\n", id);
	    return -1;
	} else if (intfluka < 0) {
	    if (fVerbosityLevel >= 3)
		printf("PDGFromId: Error intfluka < 0: %d\n", id);
	    return -1;
	}
	if (fVerbosityLevel >= 3)
	    printf("mpdgha called with %d %d \n", id, intfluka);
	// MPDGHA() goes from fluka internal to pdg.
	return mpdgha(intfluka);
    } else {
	// ions and optical photons
	return idSpecial[id + 6];
    }
}

void TFluka::StopTrack()
{
    // Set stopping conditions
    // Works for photons and charged particles
    fStopped = kTRUE;
}
  
//_____________________________________________________________________________
// methods for physics management
//____________________________________________________________________________ 
//
// set methods
//

void TFluka::SetProcess(const char* flagName, Int_t flagValue, Int_t imed)
{
//  Set process user flag for material imat
//
    TFlukaConfigOption* proc = new TFlukaConfigOption(flagName, flagValue, imed);
    fProcesses->Add(proc);
}

//______________________________________________________________________________ 
Bool_t TFluka::SetProcess(const char* flagName, Int_t flagValue)
{
//  Set process user flag 
//
//    
//  Update if already in the list
//

    TIter next(fProcesses);
    TFlukaConfigOption* proc;
    while((proc = (TFlukaConfigOption*)next()))
    { 
	if (strcmp(proc->GetName(), flagName) == 0) {
	    proc->SetFlag(flagValue);
	    proc->SetMedium(-1);
	    return kTRUE;
         }
    }
//
// If not create a new process
//    

    proc = new TFlukaConfigOption(flagName, flagValue);
    fProcesses->Add(proc);
    
    return kTRUE;  
}

//______________________________________________________________________________ 
void TFluka::SetCut(const char* cutName, Double_t cutValue, Int_t imed)
{
// Set user cut value for material imed
//
    printf("TFluka::SetCut %s %e %d \n", cutName, cutValue, imed);
    
    TFlukaConfigOption* cut = new TFlukaConfigOption(cutName, cutValue, imed);
    fCuts->Add(cut);
}

//______________________________________________________________________________ 
Bool_t TFluka::SetCut(const char* cutName, Double_t cutValue)
{
// Set user cut value 
//
//    
//  Update if already in the list
//

    TIter next(fCuts);
    TFlukaConfigOption* cut;
    while((cut = (TFlukaConfigOption*)next()))
    { 
	if (strcmp(cut->GetName(), cutName) == 0) {
	    cut->SetCut(cutValue);
	    return kTRUE;
         }
    }
//
// If not create a new process
//    

    cut = new TFlukaConfigOption(cutName, cutValue);
    fCuts->Add(cut);
    
    return kTRUE;  
}

void TFluka::SetUserScoring(const char* option, Int_t npar, Float_t what[12])
{
//
// Ads a user scoring option to th list
//
    TFlukaScoringOption* opt = new TFlukaScoringOption(option, "User Scoring", npar, what);
    fUserScore->Add(opt);
}


//______________________________________________________________________________ 
Double_t TFluka::Xsec(char*, Double_t, Int_t, Int_t)
{
  printf("WARNING: Xsec not yet implemented !\n"); return -1.;
}


//______________________________________________________________________________ 
void TFluka::InitPhysics()
{
//
// Physics initialisation with preparation of FLUKA input cards
//
  printf("=>InitPhysics\n");
  Int_t j, k;
  Double_t fCut;

  FILE *pFlukaVmcCoreInp, *pFlukaVmcFlukaMat, *pFlukaVmcInp;

  Double_t zero  = 0.0;
  Double_t one   = 1.0;
  Double_t two   = 2.0;
  Double_t three = 3.0;

  Float_t fLastMaterial = fGeom->GetLastMaterialIndex();
  if (fVerbosityLevel >= 3) printf("   last FLUKA material is %g\n", fLastMaterial);

  // Prepare  Cerenkov
  TObjArray *matList = GetFlukaMaterials();
  Int_t nmaterial =  matList->GetEntriesFast();
  fMaterials = new Int_t[nmaterial+3];
	      
// construct file names

  TString sFlukaVmcCoreInp = getenv("ALICE_ROOT");
  sFlukaVmcCoreInp +="/TFluka/input/";
  TString sFlukaVmcTmp = "flukaMat.inp";
  TString sFlukaVmcInp = GetInputFileName();
  sFlukaVmcCoreInp += GetCoreInputFileName();

// open files 

  if ((pFlukaVmcCoreInp = fopen(sFlukaVmcCoreInp.Data(),"r")) == NULL) {
      printf("\nCannot open file %s\n",sFlukaVmcCoreInp.Data());
      exit(1);
  }
  if ((pFlukaVmcFlukaMat = fopen(sFlukaVmcTmp.Data(),"r")) == NULL) {
      printf("\nCannot open file %s\n",sFlukaVmcTmp.Data());
      exit(1);
  }
  if ((pFlukaVmcInp = fopen(sFlukaVmcInp.Data(),"w")) == NULL) {
      printf("\nCannot open file %s\n",sFlukaVmcInp.Data());
      exit(1);
  }

// copy core input file 
  Char_t sLine[255];
  Float_t fEventsPerRun;
  
  while ((fgets(sLine,255,pFlukaVmcCoreInp)) != NULL) {
      if (strncmp(sLine,"GEOEND",6) != 0)
	  fprintf(pFlukaVmcInp,"%s",sLine); // copy until GEOEND card
      else {
	  fprintf(pFlukaVmcInp,"GEOEND\n");   // add GEOEND card
	  goto flukamat;
      }
  } // end of while until GEOEND card
  

 flukamat:
  while ((fgets(sLine,255,pFlukaVmcFlukaMat)) != NULL) { // copy flukaMat.inp file
      fprintf(pFlukaVmcInp,"%s\n",sLine);
  }
  
  while ((fgets(sLine,255,pFlukaVmcCoreInp)) != NULL) { 
      if (strncmp(sLine,"START",5) != 0)
	  fprintf(pFlukaVmcInp,"%s\n",sLine);
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
  fprintf(pFlukaVmcInp,"*----------------------------------------------------------------------------- \n");
  fprintf(pFlukaVmcInp,"*----- The following data are generated from SetProcess and SetCut calls ----- \n");
  fprintf(pFlukaVmcInp,"*----------------------------------------------------------------------------- \n");

// Outer loop over processes
  TIter next(fProcesses);
  TFlukaConfigOption *proc;
// Inner loop over processes
  TIter nextp(fProcesses);
  TFlukaConfigOption *procp;
// Loop over cuts
  TIter nextc(fCuts);
  TFlukaConfigOption *cut = 0x0;

  while((proc = (TFlukaConfigOption*)next())) {
      Float_t matMin = three;
      Float_t matMax = fLastMaterial;
      Bool_t  global = kTRUE;
      if (proc->Medium() != -1) {
	  Int_t mat;
	  if ((mat = proc->Medium()) >= GetFlukaMaterials()->GetEntries()) continue;
	  matMin = Float_t(mat);
	  matMax = matMin;
	  global = kFALSE;

	  fprintf(pFlukaVmcInp,"*\n*Material specific process setting for #%8d \n", mat);
	  printf("Process for %d \n", mat);
	  TGeoMaterial*    material =  (TGeoMaterial*) (GetFlukaMaterials())->At(GetMaterialIndex(mat));
	  printf("Process for %d %s\n", mat, material->GetName());
	  
      }
      
    // annihilation
    // G3 default value: 1
    // G4 processes: G4eplusAnnihilation/G4IeplusAnnihilation
    // Particles: e+
    // Physics:   EM
    // flag = 0 no annihilation
    // flag = 1 annihilation, decays processed
    // flag = 2 annihilation, no decay product stored
    // gMC ->SetProcess("ANNI",1); // EMFCUT   -1.   0.  0. 3. lastmat 0. ANNH-THR
      if (strncmp(proc->GetName(),"ANNI",4) == 0) {
	  if (proc->Flag() == 1 || proc->Flag() == 2) {
	      fprintf(pFlukaVmcInp,"*\n*Kinetic energy threshold (GeV) for e+ annihilation - resets to default=0.\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('ANNI',1) or SetProcess('ANNI',2)\n");
	      // -one = kinetic energy threshold (GeV) for e+ annihilation (resets to default=0)
	      // zero = not used
	      // zero = not used
	      // matMin = lower bound of the material indices in which the respective thresholds apply
	      // matMax = upper bound of the material indices in which the respective thresholds apply
	      // one = step length in assigning indices
	      // "ANNH-THR"; 
	      fprintf(pFlukaVmcInp,"EMFCUT    %10.1f%10.1f%10.1f%10.1f%10.1f%10.1fANNH-THR\n",-one,zero,zero,matMin,matMax,one);
	  }
	  else if (proc->Flag() == 0) {
	      fprintf(pFlukaVmcInp,"*\n*No annihilation - no FLUKA card generated\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('ANNI',0)\n");
	  }
	  else  {
	      fprintf(pFlukaVmcInp,"*\n*Illegal flag value in SetProcess('ANNI',?) call.\n");
	      fprintf(pFlukaVmcInp,"*No FLUKA card generated\n");
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
    else if ((strncmp(proc->GetName(),"PAIR",4) == 0) && (proc->Flag() == 1 || proc->Flag() == 2)) {

	nextp.Reset();
	
	while ((procp = (TFlukaConfigOption*)nextp())) {
	    if ((strncmp(procp->GetName(),"BREM",4) == 0) && 
		(proc->Flag() == 1 || procp->Flag() == 2) &&
		(procp->Medium() == proc->Medium())) {
		fprintf(pFlukaVmcInp,"*\n*Bremsstrahlung and pair production by muons and charged hadrons both activated\n");
		fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('BREM',1) and SetProcess('PAIR',1)\n");
		fprintf(pFlukaVmcInp,"*Energy threshold set by call SetCut('BCUTM',cut) or set to 0.\n");
		fprintf(pFlukaVmcInp,"*Energy threshold set by call SetCut('PPCUTM',cut) or set to 0.\n");
		// three = bremsstrahlung and pair production by muons and charged hadrons both are activated
		fprintf(pFlukaVmcInp,"PAIRBREM  %10.1f",three);
		// direct pair production by muons
		// G4 particles: "e-", "e+"
		// G3 default value: 0.01 GeV
		//gMC ->SetCut("PPCUTM",cut); // total energy cut for direct pair prod. by muons
		fCut = 0.0;
		nextc.Reset();
		while ((cut = (TFlukaConfigOption*)nextc())) {
		    if (strncmp(cut->GetName(), "PPCUTM", 6) == 0 &&
			(cut->Medium() == proc->Medium())) fCut = cut->Cut();
		}
		fprintf(pFlukaVmcInp,"%10.4g",fCut);
		// fCut; = e+, e- kinetic energy threshold (in GeV) for explicit pair production.
		// muon and hadron bremsstrahlung
		// G4 particles: "gamma"
		// G3 default value: CUTGAM=0.001 GeV
		//gMC ->SetCut("BCUTM",cut);  // cut for muon and hadron bremsstrahlung
		fCut = 0.0;
		nextc.Reset();
		while ((cut = (TFlukaConfigOption*)nextc())) {
		    if (strncmp(cut->GetName(), "BCUTM", 5) == 0 &&
			(cut->Medium() == proc->Medium())) fCut = cut->Cut();
		}
		fprintf(pFlukaVmcInp,"%10.4g%10.1f%10.1f\n",fCut,matMin,matMax);
		// fCut = photon energy threshold (GeV) for explicit bremsstrahlung production
		// matMin = lower bound of the material indices in which the respective thresholds apply
		// matMax = upper bound of the material indices in which the respective thresholds apply
		
		// for e+ and e-
		fprintf(pFlukaVmcInp,"*\n*Kinetic energy threshold (GeV) for e+/e- bremsstrahlung - resets to default=0.\n");
		fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('BREM',1);\n");
		fCut = -1.0;
		nextc.Reset();
		while ((cut = (TFlukaConfigOption*)nextc())) {
		    if (strncmp(cut->GetName(), "BCUTE", 5) == 0 &&
			(cut->Medium() == proc->Medium())) fCut = cut->Cut();
		}
		//fCut = kinetic energy threshold (GeV) for e+/e- bremsstrahlung (resets to default=0)
		// zero = not used
		// zero = not used
		// matMin = lower bound of the material indices in which the respective thresholds apply
		// matMax = upper bound of the material indices in which the respective thresholds apply
		// one = step length in assigning indices
		// "ELPO-THR"; 
		fprintf(pFlukaVmcInp,"EMFCUT    %10.4g%10.1f%10.1f%10.1f%10.1f%10.1fELPO-THR\n",fCut,zero,zero,matMin,matMax,one);
		
          // for e+ and e-
		fprintf(pFlukaVmcInp,"*\n*Pair production by electrons is activated\n");
		fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('PAIR',1);\n");
		fCut = -1.0;
		nextc.Reset();
		while ((cut = (TFlukaConfigOption*)nextc())) {
		    if (strncmp(cut->GetName(), "CUTGAM", 6) == 0 &&
			(cut->Medium() == proc->Medium())) fCut = cut->Cut();
		}
		// fCut = energy threshold (GeV) for gamma pair production (< 0.0 : resets to default, = 0.0 : ignored)
		// matMin = lower bound of the material indices in which the respective thresholds apply
		// matMax =  upper bound of the material indices in which the respective thresholds apply
		// one = step length in assigning indices
		fprintf(pFlukaVmcInp,"EMFCUT    %10.1f%10.1f%10.4g%10.1f%10.1f%10.1fPHOT-THR\n",zero,zero,fCut,matMin,matMax,one);
		goto BOTH;
	    } // end of if for BREM
	} // end of loop for BREM
	
	// only pair production by muons and charged hadrons is activated
	fprintf(pFlukaVmcInp,"*\n*Pair production by muons and charged hadrons is activated\n");
	fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('PAIR',1) or SetProcess('PAIR',2)\n");
	fprintf(pFlukaVmcInp,"*Energy threshold set by call SetCut('PPCUTM',cut) or set to 0.\n");
	// direct pair production by muons
	// G4 particles: "e-", "e+"
	// G3 default value: 0.01 GeV
	//gMC ->SetCut("PPCUTM",cut); // total energy cut for direct pair prod. by muons
	// one = pair production by muons and charged hadrons is activated
	// zero = e+, e- kinetic energy threshold (in GeV) for explicit pair production.
	// zero = no explicit bremsstrahlung production is simulated
	// matMin = lower bound of the material indices in which the respective thresholds apply
	// matMax = upper bound of the material indices in which the respective thresholds apply
	fprintf(pFlukaVmcInp,"PAIRBREM  %10.1f%10.1f%10.1f%10.1f%10.1f\n",one,zero,zero,matMin,matMax);
	
	// for e+ and e-
	fprintf(pFlukaVmcInp,"*\n*Pair production by electrons is activated\n");
	fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('PAIR',1) or SetProcess('PAIR',2)\n");
	fCut = -1.0;
	nextc.Reset();
	while ((cut = (TFlukaConfigOption*)nextc())) {
	    if (strncmp(cut->GetName(), "CUTGAM", 6) == 0 &&
		(cut->Medium() == proc->Medium())) fCut = cut->Cut();
	}
	// zero = energy threshold (GeV) for Compton scattering (= 0.0 : ignored)
	// zero = energy threshold (GeV) for Photoelectric (= 0.0 : ignored)
	// fCut = energy threshold (GeV) for gamma pair production (< 0.0 : resets to default, = 0.0 : ignored)
	// matMin = lower bound of the material indices in which the respective thresholds apply
	// matMax = upper bound of the material indices in which the respective thresholds apply
	// one = step length in assigning indices
	fprintf(pFlukaVmcInp,"EMFCUT    %10.1f%10.1f%10.4g%10.1f%10.1f%10.1fPHOT-THR\n",zero,zero,fCut,matMin,matMax,one);
      
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
      else if (strncmp(proc->GetName(),"BREM",4) == 0) {
	  nextp.Reset();
	  while((procp = (TFlukaConfigOption*)nextp())) {
	      if ((strncmp(procp->GetName(),"PAIR",4) == 0) && 
		  procp->Flag() == 1 &&
		  (procp->Medium() == proc->Medium())) goto NOBREM;
	  }
          if (proc->Flag() == 1 || proc->Flag() == 2) { 
	      fprintf(pFlukaVmcInp,"*\n*Bremsstrahlung by muons and charged hadrons is activated\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('BREM',1) or SetProcess('BREM',2)\n");
	      fprintf(pFlukaVmcInp,"*Energy threshold set by call SetCut('BCUTM',cut) or set to 0.\n");
	      // two = bremsstrahlung by muons and charged hadrons is activated
	      // zero = no meaning
	      // muon and hadron bremsstrahlung
	      // G4 particles: "gamma"
	      // G3 default value: CUTGAM=0.001 GeV
	      //gMC ->SetCut("BCUTM",cut);  // cut for muon and hadron bremsstrahlung
	      fCut = 0.0;
	      nextc.Reset();
	      while ((cut = (TFlukaConfigOption*)nextc())) {
		  if (strncmp(cut->GetName(), "BCUTM", 5) == 0 &&
		      (cut->Medium() == proc->Medium())) fCut = cut->Cut();
	      }
	      // fCut = photon energy threshold (GeV) for explicit bremsstrahlung production
	      // matMin = lower bound of the material indices in which the respective thresholds apply
	      // matMax = upper bound of the material indices in which the respective thresholds apply
	      fprintf(pFlukaVmcInp,"PAIRBREM  %10.1f%10.1f%10.4g%10.1f%10.1f\n",two,zero,fCut,matMin,matMax);
	      
	      // for e+ and e-
	      fprintf(pFlukaVmcInp,"*\n*Kinetic energy threshold (GeV) for e+/e- bremsstrahlung - resets to default=0.\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('BREM',1);");
	      // - one = kinetic energy threshold (GeV) for e+/e- bremsstrahlung (resets to default=0)
	      // zero = not used
	      // zero = not used
	      // matMin = lower bound of the material indices in which the respective thresholds apply
	      // matMax = upper bound of the material indices in which the respective thresholds apply
	      // one = step length in assigning indices
	      //"ELPO-THR"; 
	      fprintf(pFlukaVmcInp,"EMFCUT    %10.1f%10.1f%10.1f%10.1f%10.1f%10.1fELPO-THR\n",-one,zero,zero,matMin,matMax,one);
	  }
	  else if (proc->Flag() == 0) {
	      fprintf(pFlukaVmcInp,"*\n*No bremsstrahlung - no FLUKA card generated\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('BREM',0)\n");
	  }
	  else  {
	      fprintf(pFlukaVmcInp,"*\n*Illegal flag value in SetProcess('BREM',?) call.\n");
	      fprintf(pFlukaVmcInp,"*No FLUKA card generated\n");
	  }
      NOBREM:
	  j = 0;
      } // end of else if (strncmp(proc->GetName(),"BREM",4) == 0)
      
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
      
      else if (strncmp(proc->GetName(),"CKOV",4) == 0) {
	  if ((proc->Flag() == 1 || proc->Flag() == 2) && global) {
	      // Write comments
	      fprintf(pFlukaVmcInp, "* \n"); 
	      fprintf(pFlukaVmcInp, "*Cerenkov photon generation\n"); 
	      fprintf(pFlukaVmcInp, "*Generated from call: SetProcess('CKOV',1) or SetProcess('CKOV',2)\n"); 
	      // Loop over media 
	      for (Int_t im = 0; im < nmaterial; im++)
	      {
		  TGeoMaterial* material = dynamic_cast<TGeoMaterial*> (matList->At(im));
		  Int_t idmat = material->GetIndex();

		  if (!global && idmat != proc->Medium()) continue;
		  
		  fMaterials[idmat] = im;
		  // Skip media with no Cerenkov properties
		  TFlukaCerenkov* cerenkovProp;
		  if (!(cerenkovProp = dynamic_cast<TFlukaCerenkov*>(material->GetCerenkovProperties()))) continue;
		  //
		  // This medium has Cerenkov properties 
		  //
		  //
		  // Write OPT-PROD card for each medium 
		  Float_t  emin  = cerenkovProp->GetMinimumEnergy();
		  Float_t  emax  = cerenkovProp->GetMaximumEnergy();	      
		  fprintf(pFlukaVmcInp, "OPT-PROD  %10.4g%10.4g%10.4g%10.4g%10.4g%10.4gCERENKOV\n", emin, emax, 0., 
			  Float_t(idmat), Float_t(idmat), 0.); 
		  //
		  // Write OPT-PROP card for each medium 
		  // Forcing FLUKA to call user routines (queffc.cxx, rflctv.cxx, rfrndx.cxx)
		  //
		  fprintf(pFlukaVmcInp, "OPT-PROP  %10.4g%10.4g%10.4g%10.1f%10.1f%10.1fWV-LIMIT\n",  
			  cerenkovProp->GetMinimumWavelength(),
			  cerenkovProp->GetMaximumWavelength(), 
			  cerenkovProp->GetMaximumWavelength(), 
			  Float_t(idmat), Float_t(idmat), 0.0);
		  
		  if (cerenkovProp->IsMetal()) {
		      fprintf(pFlukaVmcInp, "OPT-PROP  %10.1f%10.1f%10.1f%10.1f%10.1f%10.1fMETAL\n",  
			      -100., -100., -100., 
			      Float_t(idmat), Float_t(idmat), 0.0);
		  } else {
		      fprintf(pFlukaVmcInp, "OPT-PROP  %10.1f%10.1f%10.1f%10.1f%10.1f%10.1f\n",  
			      -100., -100., -100., 
			      Float_t(idmat), Float_t(idmat), 0.0);
		  }
		  
		  
		  for (Int_t j = 0; j < 3; j++) {
		      fprintf(pFlukaVmcInp, "OPT-PROP  %10.1f%10.1f%10.1f%10.1f%10.1f%10.1f&\n",  
			      -100., -100., -100., 
			      Float_t(idmat), Float_t(idmat), 0.0);
		  }
		  // Photon detection efficiency user defined
		  
		  if (cerenkovProp->IsSensitive())
		      fprintf(pFlukaVmcInp, "OPT-PROP  %10.1f%10.1f%10.1f%10.1f%10.1f%10.1fSENSITIV\n",  
			      -100., -100., -100., 
			      Float_t(idmat), Float_t(idmat), 0.0);
		  
	      } // materials
	  } else if (proc->Flag() == 0) {
	      fprintf(pFlukaVmcInp,"*\n*No Cerenkov photon generation\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('CKOV',0)\n");
	      // zero = not used
	      // zero = not used
	      // zero = not used
	      // matMin = lower bound of the material indices in which the respective thresholds apply
	      // matMax = upper bound of the material indices in which the respective thresholds apply
	      // one = step length in assigning indices
	      //"CERE-OFF"; 
	      fprintf(pFlukaVmcInp,"OPT-PROD  %10.1f%10.1f%10.1f%10.1f%10.1f%10.1fCERE-OFF\n",zero,zero,zero,matMin,matMax,one);
	  }
	  else  {
	      fprintf(pFlukaVmcInp,"*\n*Illegal flag value in SetProcess('CKOV',?) call.\n");
	      fprintf(pFlukaVmcInp,"*No FLUKA card generated\n");
	  }
      } // end of else if (strncmp(proc->GetName(),"CKOV",4) == 0)
      
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
      else if (strncmp(proc->GetName(),"COMP",4) == 0) {
	  if (proc->Flag() == 1 || proc->Flag() == 2) { 
	      fprintf(pFlukaVmcInp,"*\n*Energy threshold (GeV) for Compton scattering - resets to default=0.\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('COMP',1);\n");
	      // - one = energy threshold (GeV) for Compton scattering - resets to default=0.
	      // zero = not used
	      // zero = not used
	      // matMin = lower bound of the material indices in which the respective thresholds apply
	      // matMax = upper bound of the material indices in which the respective thresholds apply
	      // one = step length in assigning indices
	      //"PHOT-THR"; 
	      fprintf(pFlukaVmcInp,"EMFCUT    %10.1f%10.1f%10.1f%10.1f%10.1f%10.1fPHOT-THR\n",-one,zero,zero,matMin,matMax,one);
	  }
	  else if (proc->Flag() == 0) {
	      fprintf(pFlukaVmcInp,"*\n*No Compton scattering - no FLUKA card generated\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('COMP',0)\n");
	  }
	  else  {
	      fprintf(pFlukaVmcInp,"*\n*Illegal flag value in SetProcess('COMP',?) call.\n");
	      fprintf(pFlukaVmcInp,"*No FLUKA card generated\n");
	  }
      } // end of else if (strncmp(proc->GetName(),"COMP",4) == 0)
      
      // decay
      // G3 default value: 1
      // G4 process: G4Decay
      // 
      // Particles: all which decay is applicable for
      // Physics:   General
      // flag = 0 no decays
      // flag = 1 decays, secondaries processed
      // flag = 2 decays, no secondaries stored
      //gMC ->SetProcess("DCAY",0); // not available
      else if ((strncmp(proc->GetName(),"DCAY",4) == 0) && proc->Flag() == 0) 
	  cout << "SetProcess for flag =" << proc->GetName() << " value=" << proc->Flag() << " not avaliable!" << endl;
      else if ((strncmp(proc->GetName(),"DCAY",4) == 0) && proc->Flag() == 1) {
          // Nothing to do decays are switched on by default
      }
      
      
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
      else if (strncmp(proc->GetName(),"DRAY",4) == 0) {
	  if (proc->Flag() == 0 || proc->Flag() == 4) {
	      fprintf(pFlukaVmcInp,"*\n*Kinetic energy threshold (GeV) for delta ray production\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('DRAY',0) or SetProcess('DRAY',4)\n");
	      fprintf(pFlukaVmcInp,"*No delta ray production by muons - threshold set artificially high\n");
	      Double_t emin = 1.0e+6; // kinetic energy threshold (GeV) for delta ray production (discrete energy transfer)
	      // zero = ignored
	      // zero = ignored
	      // matMin = lower bound of the material indices in which the respective thresholds apply
	      // matMax = upper bound of the material indices in which the respective thresholds apply
	      // one = step length in assigning indices
	      fprintf(pFlukaVmcInp,"DELTARAY  %10.4g%10.1f%10.1f%10.1f%10.1f%10.1f\n",emin,zero,zero,matMin,matMax,one);
	  }
	  else if (proc->Flag() == 1 || proc->Flag() == 2 || proc->Flag() == 3) {
	      fprintf(pFlukaVmcInp,"*\n*Kinetic energy threshold (GeV) for delta ray production\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('DRAY',flag), flag=1,2,3\n");
	      fprintf(pFlukaVmcInp,"*Delta ray production by muons switched on\n");
	      fprintf(pFlukaVmcInp,"*Energy threshold set by call SetCut('DCUTM',cut) or set to 1.0e+6.\n");
	      fCut = 1.0e+6;
	      nextc.Reset();
	      while ((cut = (TFlukaConfigOption*)nextc())) {
		  if (strncmp(cut->GetName(), "DCUTM", 5) == 0 &&
		      cut->Medium() == proc->Medium()) fCut = cut->Cut();
	      }
	      // fCut = kinetic energy threshold (GeV) for delta ray production (discrete energy transfer)
	      // zero = ignored
	      // zero = ignored
	      // matMin = lower bound of the material indices in which the respective thresholds apply
	      // matMax =  upper bound of the material indices in which the respective thresholds apply
	      // one = step length in assigning indices
	      fprintf(pFlukaVmcInp,"DELTARAY  %10.4g%10.1f%10.1f%10.1f%10.1f%10.1f\n",fCut,zero,zero,matMin,matMax,one);
	  }
	  else  {
	      fprintf(pFlukaVmcInp,"*\n*Illegal flag value in SetProcess('DRAY',?) call.\n");
	      fprintf(pFlukaVmcInp,"*No FLUKA card generated\n");
	  }
      } // end of else if (strncmp(proc->GetName(),"DRAY",4) == 0)
      
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
      else if (strncmp(proc->GetName(),"HADR",4) == 0) {
	  if (proc->Flag() == 1 || proc->Flag() == 2) {
	      fprintf(pFlukaVmcInp,"*\n*Hadronic interaction is ON by default in FLUKA\n");
	      fprintf(pFlukaVmcInp,"*No FLUKA card generated\n");
	  }
	  else if (proc->Flag() == 0) {
	      fprintf(pFlukaVmcInp,"*\n*Hadronic interaction is set OFF\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('HADR',0);\n");
	      fprintf(pFlukaVmcInp,"*Switching off hadronic interactions not foreseen in FLUKA\n");
	      fprintf(pFlukaVmcInp,"THRESHOL  %10.1f%10.1f%10.1f%10.1e%10.1f\n",zero, zero, zero, 1.e10, zero);
	  }
	  else  {
	      fprintf(pFlukaVmcInp,"*\n*Illegal flag value in SetProcess('HADR',?) call.\n");
	      fprintf(pFlukaVmcInp,"*No FLUKA card generated\n");
	  }
      } // end of else if (strncmp(proc->GetName(),"HADR",4) == 0)
      
      
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
      else if (strncmp(proc->GetName(),"LOSS",4) == 0) {
	  if (proc->Flag() == 2) { // complete energy loss fluctuations
	      fprintf(pFlukaVmcInp,"*\n*Complete energy loss fluctuations do not exist in FLUKA\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('LOSS',2);\n");
	      fprintf(pFlukaVmcInp,"*flag=2=complete energy loss fluctuations\n");
	      fprintf(pFlukaVmcInp,"*No FLUKA card generated\n");
	  }
	  else if (proc->Flag() == 1 || proc->Flag() == 3) { // restricted energy loss fluctuations
	      fprintf(pFlukaVmcInp,"*\n*Restricted energy loss fluctuations\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('LOSS',1) or SetProcess('LOSS',3)\n");
	      // one = restricted energy loss fluctuations (for hadrons and muons) switched on
	      // one = restricted energy loss fluctuations (for e+ and e-) switched on
	      // one = minimal accuracy
	      // matMin = lower bound of the material indices in which the respective thresholds apply
	      // upper bound of the material indices in which the respective thresholds apply
	      fprintf(pFlukaVmcInp,"IONFLUCT  %10.1f%10.1f%10.1f%10.1f%10.1f\n",one,one,one,matMin,matMax);
	  }
	  else if (proc->Flag() == 4) { // no energy loss fluctuations
	      fprintf(pFlukaVmcInp,"*\n*No energy loss fluctuations\n");
	      fprintf(pFlukaVmcInp,"*\n*Generated from call: SetProcess('LOSS',4)\n");
	      // - one = restricted energy loss fluctuations (for hadrons and muons) switched off
	      // - one = restricted energy loss fluctuations (for e+ and e-) switched off
	      // one = minimal accuracy
	      // matMin = lower bound of the material indices in which the respective thresholds apply
	      // matMax = upper bound of the material indices in which the respective thresholds apply
	      fprintf(pFlukaVmcInp,"IONFLUCT  %10.1f%10.1f%10.1f%10.1f%10.1f\n",-one,-one,one,matMin,matMax);
	  }
	  else  {
	      fprintf(pFlukaVmcInp,"*\n*Illegal flag value in SetProcess('LOSS',?) call.\n");
	      fprintf(pFlukaVmcInp,"*No FLUKA card generated\n");
	  }
      } // end of else if (strncmp(proc->GetName(),"LOSS",4) == 0)
      
      
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
      else if (strncmp(proc->GetName(),"MULS",4) == 0) {
	  if (proc->Flag() == 1 || proc->Flag() == 2 || proc->Flag() == 3) {
	      fprintf(pFlukaVmcInp,"*\n*Multiple scattering is ON by default for e+e- and for hadrons/muons\n");
	      fprintf(pFlukaVmcInp,"*No FLUKA card generated\n");
	  }
	  else if (proc->Flag() == 0) {
	      fprintf(pFlukaVmcInp,"*\n*Multiple scattering is set OFF\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('MULS',0);\n");
	      // zero = ignored
	      // three = multiple scattering for hadrons and muons is completely suppressed
	      // three = multiple scattering for e+ and e- is completely suppressed
	      // matMin = lower bound of the material indices in which the respective thresholds apply
	      // matMax = upper bound of the material indices in which the respective thresholds apply
	      fprintf(pFlukaVmcInp,"MULSOPT   %10.1f%10.1f%10.1f%10.1f%10.1f\n",zero,three,three,matMin,matMax);
	  }
	  else  {
	      fprintf(pFlukaVmcInp,"*\n*Illegal flag value in SetProcess('MULS',?) call.\n");
	      fprintf(pFlukaVmcInp,"*No FLUKA card generated\n");
	  }
      } // end of else if (strncmp(proc->GetName(),"MULS",4) == 0)
      

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
      else if (strncmp(proc->GetName(),"MUNU",4) == 0) {
	  if (proc->Flag() == 1) {
	      fprintf(pFlukaVmcInp,"*\n*Muon nuclear interactions with production of secondary hadrons\n");
	      fprintf(pFlukaVmcInp,"*\n*Generated from call: SetProcess('MUNU',1);\n");
	      // one = full simulation of muon nuclear interactions and production of secondary hadrons
	      // zero = ratio of longitudinal to transverse virtual photon cross-section - Default = 0.25.
	      // zero = fraction of rho-like interactions ( must be < 1) - Default = 0.75.
	      // matMin = lower bound of the material indices in which the respective thresholds apply
	      // matMax = upper bound of the material indices in which the respective thresholds apply
	      fprintf(pFlukaVmcInp,"MUPHOTON  %10.1f%10.1f%10.1f%10.1f%10.1f\n",one,zero,zero,matMin,matMax);
	  }
	  else if (proc->Flag() == 2) {
	      fprintf(pFlukaVmcInp,"*\n*Muon nuclear interactions without production of secondary hadrons\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('MUNU',2);\n");
	      // two = full simulation of muon nuclear interactions and production of secondary hadrons
	      // zero = ratio of longitudinal to transverse virtual photon cross-section - Default = 0.25.
	      // zero = fraction of rho-like interactions ( must be < 1) - Default = 0.75.
	      // matMin = lower bound of the material indices in which the respective thresholds apply
	      // matMax = upper bound of the material indices in which the respective thresholds apply
	      fprintf(pFlukaVmcInp,"MUPHOTON  %10.1f%10.1f%10.1f%10.1f%10.1f\n",two,zero,zero,matMin,matMax);
	  }
	  else if (proc->Flag() == 0) {
	      fprintf(pFlukaVmcInp,"*\n*No muon nuclear interaction - no FLUKA card generated\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('MUNU',0)\n");
	  }
	  else  {
	      fprintf(pFlukaVmcInp,"*\n*Illegal flag value in SetProcess('MUNU',?) call.\n");
	      fprintf(pFlukaVmcInp,"*No FLUKA card generated\n");
	  }
      } // end of else if (strncmp(proc->GetName(),"MUNU",4) == 0)
      
      
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
      else if (strncmp(proc->GetName(),"PFIS",4) == 0) {
	  if (proc->Flag() == 0) {
	      fprintf(pFlukaVmcInp,"*\n*No photonuclear interactions\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('PFIS',0);\n");
	      // - one = no photonuclear interactions
	      // zero = not used
	      // zero = not used
	      // matMin = lower bound of the material indices in which the respective thresholds apply
	      // matMax = upper bound of the material indices in which the respective thresholds apply
	      fprintf(pFlukaVmcInp,"PHOTONUC  %10.1f%10.1f%10.1f%10.1f%10.1f\n",-one,zero,zero,matMin,matMax);
	  }
	  else if (proc->Flag() == 1) {
	      fprintf(pFlukaVmcInp,"*\n*Photon nuclear interactions are activated at all energies\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('PFIS',1);\n");
	      // one = photonuclear interactions are activated at all energies
	      // zero = not used
	      // zero = not used
	      // matMin = lower bound of the material indices in which the respective thresholds apply
	      // matMax = upper bound of the material indices in which the respective thresholds apply
	      fprintf(pFlukaVmcInp,"PHOTONUC  %10.1f%10.1f%10.1f%10.1f%10.1f\n",one,zero,zero,matMin,matMax);
	  }
	  else if (proc->Flag() == 0) {
	      fprintf(pFlukaVmcInp,"*\n*No photofission - no FLUKA card generated\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('PFIS',0)\n");
	  }
	  else {
	      fprintf(pFlukaVmcInp,"*\n*Illegal flag value in SetProcess('PFIS',?) call.\n");
	      fprintf(pFlukaVmcInp,"*No FLUKA card generated\n");
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
      else if (strncmp(proc->GetName(),"PHOT",4) == 0) {
	  if (proc->Flag() == 1 || proc->Flag() == 2) {
	      fprintf(pFlukaVmcInp,"*\n*Photo electric effect is activated\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('PHOT',1);\n");
	      // zero = ignored
	      // - one = resets to default=0.
	      // zero = ignored
	      // matMin = lower bound of the material indices in which the respective thresholds apply
	      // matMax = upper bound of the material indices in which the respective thresholds apply
	      // one = step length in assigning indices
	      //"PHOT-THR"; 
	      fprintf(pFlukaVmcInp,"EMFCUT    %10.1f%10.1f%10.1f%10.1f%10.1f%10.1fPHOT-THR\n",zero,-one,zero,matMin,matMax,one);
	  }
	  else if (proc->Flag() == 0) {
	      fprintf(pFlukaVmcInp,"*\n*No photo electric effect - no FLUKA card generated\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('PHOT',0)\n");
	  }
	  else {
	      fprintf(pFlukaVmcInp,"*\n*Illegal flag value in SetProcess('PHOT',?) call.\n");
	      fprintf(pFlukaVmcInp,"*No FLUKA card generated\n");
	  }
      } // else if (strncmp(proc->GetName(),"PHOT",4) == 0)
      
      
      // Rayleigh scattering
      // G3 default value: 0
      // G4 process: G4OpRayleigh
      // 
      // Particles: optical photon
      // Physics:   Optical
      // flag = 0 Rayleigh scattering off
      // flag = 1 Rayleigh scattering on
      //xx gMC ->SetProcess("RAYL",1);
      else if (strncmp(proc->GetName(),"RAYL",4) == 0) {
	  if (proc->Flag() == 1) {
	      fprintf(pFlukaVmcInp,"*\n*Rayleigh scattering is ON by default in FLUKA\n");
	      fprintf(pFlukaVmcInp,"*No FLUKA card generated\n");
	  }
	  else if (proc->Flag() == 0) {
	      fprintf(pFlukaVmcInp,"*\n*Rayleigh scattering is set OFF\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('RAYL',0);\n");
	      // - one = no Rayleigh scattering and no binding corrections for Compton
	      // matMin = lower bound of the material indices in which the respective thresholds apply
	      // matMax = upper bound of the material indices in which the respective thresholds apply
	      fprintf(pFlukaVmcInp,"EMFRAY    %10.1f%10.1f%10.1f%10.1f\n",-one,three,matMin,matMax);
	  }
	  else  {
	      fprintf(pFlukaVmcInp,"*\n*Illegal flag value in SetProcess('RAYL',?) call.\n");
	      fprintf(pFlukaVmcInp,"*No FLUKA card generated\n");
	  }
      } // end of else if (strncmp(proc->GetName(),"RAYL",4) == 0)
      
      
      // synchrotron radiation in magnetic field
      // G3 default value: 0
      // G4 process: G4SynchrotronRadiation
      // 
      // Particles: ??
      // Physics:   Not set
      // flag = 0 no synchrotron radiation
      // flag = 1 synchrotron radiation
      //xx gMC ->SetProcess("SYNC",1); // synchrotron radiation generation
      else if (strncmp(proc->GetName(),"SYNC",4) == 0) {
	  fprintf(pFlukaVmcInp,"*\n*Synchrotron radiation generation is NOT implemented in FLUKA\n");
	  fprintf(pFlukaVmcInp,"*No FLUKA card generated\n");
      }
      
      
      // Automatic calculation of tracking medium parameters
      // flag = 0 no automatic calculation
      // flag = 1 automatic calculation
      //xx gMC ->SetProcess("AUTO",1); // ??? automatic computation of the tracking medium parameters
      else if (strncmp(proc->GetName(),"AUTO",4) == 0) {
	  fprintf(pFlukaVmcInp,"*\n*Automatic calculation of tracking medium parameters is always ON in FLUKA\n");
	  fprintf(pFlukaVmcInp,"*No FLUKA card generated\n");
      }
      
      
      // To control energy loss fluctuation model
      // flag = 0 Urban model
      // flag = 1 PAI model
      // flag = 2 PAI+ASHO model (not active at the moment)
      //xx gMC ->SetProcess("STRA",1); // ??? energy fluctuation model
      else if (strncmp(proc->GetName(),"STRA",4) == 0) {
	  if (proc->Flag() == 0 || proc->Flag() == 2 || proc->Flag() == 3) {
	      fprintf(pFlukaVmcInp,"*\n*Ionization energy losses calculation is activated\n");
	      fprintf(pFlukaVmcInp,"*Generated from call: SetProcess('STRA',n);, n=0,1,2\n");
	      // one = restricted energy loss fluctuations (for hadrons and muons) switched on
	      // one = restricted energy loss fluctuations (for e+ and e-) switched on
	      // one = minimal accuracy
	      // matMin = lower bound of the material indices in which the respective thresholds apply
	      // matMax = upper bound of the material indices in which the respective thresholds apply
	      fprintf(pFlukaVmcInp,"IONFLUCT  %10.1f%10.1f%10.1f%10.1f%10.1f\n",one,one,one,matMin,matMax);
	  }
	  else {
	      fprintf(pFlukaVmcInp,"*\n*Illegal flag value in SetProcess('STRA',?) call.\n");
	      fprintf(pFlukaVmcInp,"*No FLUKA card generated\n");
	  }
      } // else if (strncmp(proc->GetName(),"STRA",4) == 0)
      



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
	  


	  cout << "SetProcess for flag=" << proc->GetName() << " value=" << proc->Flag() << " not yet implemented!" << endl;
      }
  } //end of loop number of SetProcess calls

 
// Loop over number of SetCut calls  
	    
  nextc.Reset();
  while ((cut = (TFlukaConfigOption*)nextc())) {
      Float_t matMin = three;
      Float_t matMax = fLastMaterial;
      Bool_t global  = kTRUE;
      if (cut->Medium() != -1) {
	  Int_t mat;
	  if ((mat = cut->Medium()) >= GetFlukaMaterials()->GetEntries()) continue;
	  matMin = Float_t(mat);
	  matMax = matMin;
	  global = kFALSE;
	  TGeoMaterial*    material =  (TGeoMaterial*) (GetFlukaMaterials())->At(GetMaterialIndex(mat));
	  fprintf(pFlukaVmcInp,"*\n*Material specific cut setting for #%8d %s %s %13.3e\n", mat, material->GetName(), cut->GetName(), cut->Cut());  

      } 

      // cuts handled in SetProcess calls
      if (strncmp(cut->GetName(),"BCUTM",5) == 0) continue;
      else if (strncmp(cut->GetName(),"BCUTE",5) == 0) continue;
      else if (strncmp(cut->GetName(),"DCUTM",5) == 0) continue;
      else if (strncmp(cut->GetName(),"PPCUTM",6) == 0) continue;
      
      // delta-rays by electrons
      // G4 particles: "e-"
      // G3 default value: 10**4 GeV
      // gMC ->SetCut("DCUTE",cut);  // cut for deltarays by electrons 
      else if (strncmp(cut->GetName(),"DCUTE",5) == 0) {
	fprintf(pFlukaVmcInp,"*Cut for delta rays by electrons\n");
	// -cut->Cut();
	// zero = ignored
	// zero = ignored
	// matMin = lower bound of the material indices in which the respective thresholds apply
	// matMax = upper bound of the material indices in which the respective thresholds apply
        // loop over materials for EMFCUT FLUKA cards
        for (j=0; j < matMax-matMin+1; j++) {
          Int_t nreg, imat, *reglist;
          Float_t ireg;
          imat = (Int_t) matMin + j;
          reglist = fGeom->GetMaterialList(imat, nreg);
          // loop over regions of a given material
          for (k=0; k<nreg; k++) {
            ireg = reglist[k];
	    fprintf(pFlukaVmcInp,"EMFCUT    %10.4g%10.1f%10.1f%10.1f%10.1f\n",-cut->Cut(),zero,zero,ireg,ireg);
          }
        }
	fprintf(pFlukaVmcInp,"DELTARAY  %10.4g%10.3f%10.3f%10.1f%10.1f%10.1f\n",cut->Cut(), 100., 1.03, matMin, matMax, 1.0);
      } // end of if for delta-rays by electrons
    

      // gammas
      // G4 particles: "gamma"
      // G3 default value: 0.001 GeV
      // gMC ->SetCut("CUTGAM",cut); // cut for gammas
      
      else if (strncmp(cut->GetName(),"CUTGAM",6) == 0 && global) {
	fprintf(pFlukaVmcInp,"*\n*Cut for gamma\n");
	fprintf(pFlukaVmcInp,"*Generated from call: SetCut('CUTGAM',cut);\n");
	// -cut->Cut();
	// 7.0 = lower bound of the particle id-numbers to which the cut-off
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f\n",-cut->Cut(),7.0);
      }
      else if (strncmp(cut->GetName(),"CUTGAM",6) == 0 && !global) {
        // loop over materials for EMFCUT FLUKA cards
        for (j=0; j < matMax-matMin+1; j++) {
          Int_t nreg, imat, *reglist;
          Float_t ireg;
          imat = (Int_t) matMin + j;
          reglist = fGeom->GetMaterialList(imat, nreg);
          // loop over regions of a given material
          for (Int_t k = 0; k < nreg; k++) {
            ireg = reglist[k];
	    fprintf(pFlukaVmcInp,"EMFCUT    %10.4g%10.4g%10.1f%10.1f%10.1f%10.1f\n", zero, cut->Cut(), zero, ireg, ireg, one);
          }
        }
      } // end of else if for gamma


      // electrons
      // G4 particles: "e-"
      // ?? positrons
      // G3 default value: 0.001 GeV
      //gMC ->SetCut("CUTELE",cut); // cut for e+,e-
      else if (strncmp(cut->GetName(),"CUTELE",6) == 0 && global) {
	fprintf(pFlukaVmcInp,"*\n*Cut for electrons\n");
	fprintf(pFlukaVmcInp,"*Generated from call: SetCut('CUTELE',cut);\n");
	// -cut->Cut();
	// three = lower bound of the particle id-numbers to which the cut-off
	// 4.0 = upper bound of the particle id-numbers to which the cut-off
	// one = step length in assigning numbers
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f%10.1f\n",-cut->Cut(),three,4.0,one);
      }
      else if (strncmp(cut->GetName(),"CUTELE",6) == 0 && !global) {
        // loop over materials for EMFCUT FLUKA cards
        for (j=0; j < matMax-matMin+1; j++) {
          Int_t nreg, imat, *reglist;
          Float_t ireg;
          imat = (Int_t) matMin + j;
          reglist = fGeom->GetMaterialList(imat, nreg);
          // loop over regions of a given material
          for (k=0; k<nreg; k++) {
            ireg = reglist[k];
	    fprintf(pFlukaVmcInp,"EMFCUT    %10.4g%10.4g%10.1f%10.1f%10.1f%10.1f\n", -cut->Cut(), zero, zero, ireg, ireg, one);
          }
        }
      } // end of else if for electrons

    
      // neutral hadrons
      // G4 particles: of type "baryon", "meson", "nucleus" with zero charge
      // G3 default value: 0.01 GeV
      //gMC ->SetCut("CUTNEU",cut); // cut for neutral hadrons
      else if (strncmp(cut->GetName(),"CUTNEU",6) == 0 && global) {
	fprintf(pFlukaVmcInp,"*\n*Cut for neutral hadrons\n");
	fprintf(pFlukaVmcInp,"*Generated from call: SetCut('CUTNEU',cut);\n");
	  
	// 8.0 = Neutron
	// 9.0 = Antineutron
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f\n",-cut->Cut(),8.0,9.0);
	  
	// 12.0 = Kaon zero long
	// 12.0 = Kaon zero long
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f\n",-cut->Cut(),12.0,12.0);
	  
	// 17.0 = Lambda, 18.0 = Antilambda
	// 19.0 = Kaon zero short
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f\n",-cut->Cut(),17.0,19.0);
	  
	// 22.0 = Sigma zero, Pion zero, Kaon zero
	// 25.0 = Antikaon zero
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f\n",-cut->Cut(),22.0,25.0);
	  
	// 32.0 = Antisigma zero
	// 32.0 = Antisigma zero
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f\n",-cut->Cut(),32.0,32.0);
	  
	// 34.0 = Xi zero
	// 35.0 = AntiXi zero
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f\n",-cut->Cut(),34.0,35.0);
	  
	// 47.0 = D zero
	// 48.0 = AntiD zero
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f\n",-cut->Cut(),47.0,48.0);
	  
	// 53.0 = Xi_c zero
	// 53.0 = Xi_c zero
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f\n",-cut->Cut(),53.0,53.0);
	  
	// 55.0 = Xi'_c zero
	// 56.0 = Omega_c zero
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f\n",-cut->Cut(),55.0,56.0);
	  
	// 59.0 = AntiXi_c zero
	// 59.0 = AntiXi_c zero
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f\n",-cut->Cut(),59.0,59.0);
	  
	// 61.0 = AntiXi'_c zero
	// 62.0 = AntiOmega_c zero
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f\n",-cut->Cut(),61.0,62.0);
      }
      
      // charged hadrons
      // G4 particles: of type "baryon", "meson", "nucleus" with non-zero charge
      // G3 default value: 0.01 GeV
      //gMC ->SetCut("CUTHAD",cut); // cut for charged hadrons
      else if (strncmp(cut->GetName(),"CUTHAD",6) == 0 && global) {
	fprintf(pFlukaVmcInp,"*\n*Cut for charged hadrons\n");
	fprintf(pFlukaVmcInp,"*Generated from call: SetCut('CUTHAD',cut);\n");
	  
	// 1.0 = Proton
	// 2.0 = Antiproton
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f\n",-cut->Cut(),1.0,2.0);
	  
	// 13.0 = Positive Pion, Negative Pion, Positive Kaon
	// 16.0 = Negative Kaon
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f\n",-cut->Cut(),13.0,16.0);
	  
	// 20.0 = Negative Sigma
	// 21.0 = Positive Sigma
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f\n",-cut->Cut(),20.0,21.0);
	  
	// 31.0 = Antisigma minus
	// 33.0 = Antisigma plus
	// 2.0 = step length
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f%10.1f\n",-cut->Cut(),31.0,33.0,2.0);
	  
	// 36.0 = Negative Xi, Positive Xi, Omega minus
	// 39.0 = Antiomega
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f\n",-cut->Cut(),36.0,39.0);
	  
	// 45.0 = D plus
	// 46.0 = D minus
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f\n",-cut->Cut(),45.0,46.0);
	  
	// 49.0 = D_s plus, D_s minus, Lambda_c plus
	// 52.0 = Xi_c plus
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f\n",-cut->Cut(),49.0,52.0);
	  
	// 54.0 = Xi'_c plus
	// 60.0 = AntiXi'_c minus
	// 6.0 = step length
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f%10.1f\n",-cut->Cut(),54.0,60.0,6.0);
	  
	// 57.0 = Antilambda_c minus
	// 58.0 = AntiXi_c minus
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f\n",-cut->Cut(),57.0,58.0);
      }

      // muons
      // G4 particles: "mu+", "mu-"
      // G3 default value: 0.01 GeV
      //gMC ->SetCut("CUTMUO",cut); // cut for mu+, mu-
      else if (strncmp(cut->GetName(),"CUTMUO",6)== 0 && global) {
	fprintf(pFlukaVmcInp,"*\n*Cut for muons\n");
	fprintf(pFlukaVmcInp,"*Generated from call: SetCut('CUTMUO',cut);\n");
	// 10.0 = Muon+
	// 11.0 = Muon-
	fprintf(pFlukaVmcInp,"PART-THR  %10.4g%10.1f%10.1f\n",-cut->Cut(),10.0,11.0);
      }
      
      //
      // time of flight cut in seconds
      // G4 particles: all
      // G3 default value: 0.01 GeV
      //gMC ->SetCut("TOFMAX",tofmax); // time of flight cuts in seconds
      else if (strncmp(cut->GetName(),"TOFMAX",6) == 0) {
	fprintf(pFlukaVmcInp,"*\n*Time of flight cuts in seconds\n");
	fprintf(pFlukaVmcInp,"*Generated from call: SetCut('TOFMAX',tofmax);\n");
	// zero = ignored
	// zero = ignored
	// -6.0 = lower bound of the particle numbers for which the transport time cut-off and/or the start signal is to be applied
	// 64.0 = upper bound of the particle numbers for which the transport time cut-off and/or the start signal is to be applied
	fprintf(pFlukaVmcInp,"TIME-CUT  %10.4g%10.1f%10.1f%10.1f%10.1f\n",cut->Cut()*1.e9,zero,zero,-6.0,64.0);
      }
      
      else if (global){
	cout << "SetCut for flag=" << cut->GetName() << " value=" << cut->Cut() << " not yet implemented!" << endl;
      }
      else {
	cout << "SetCut for flag=" << cut->GetName() << " value=" << cut->Cut() << " (material specific) not yet implemented!" << endl;
      }
      
  } //end of loop over SetCut calls
  
// Add START and STOP card
  fprintf(pFlukaVmcInp,"START     %10.1f\n",fEventsPerRun);
  fprintf(pFlukaVmcInp,"STOP      \n");
   
  
// Close files
  
   fclose(pFlukaVmcCoreInp);
   fclose(pFlukaVmcFlukaMat);
   fclose(pFlukaVmcInp);
   
} // end of InitPhysics


//______________________________________________________________________________ 
void TFluka::SetMaxStep(Double_t step)
{
// Set the maximum step size
    if (step > 1.e4) return;
    
    Int_t mreg, latt;
    fGeom->GetCurrentRegion(mreg, latt);
    STEPSZ.stepmx[mreg - 1] = step;
}


Double_t TFluka::MaxStep() const
{
// Return the maximum for current medium
    Int_t mreg, latt;
    fGeom->GetCurrentRegion(mreg, latt);
    return (STEPSZ.stepmx[mreg - 1]);
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
  if (caller == 3 || caller == 6 || caller == 11 || caller == 12) { //bxdraw,endraw,usdraw
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
  if (caller == 3 || caller == 6 || caller == 11 || caller == 12) { //bxdraw,endraw,usdraw
    x = GetXsco();
    y = GetYsco();
    z = GetZsco();
  }
  else if (caller == 4 || caller == 5) { // mgdraw, sodraw
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
  if (caller == 11 || caller==12 || caller == 3 || caller == 6) //bxdraw,endraw,usdraw
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
  if (caller == 11 || caller==12 || caller == 3 || caller == 4 || caller == 6) //bxdraw,endraw,mgdraw,usdraw
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
  if (caller == 11 || caller==12 || caller == 3 || caller == 4 || caller == 6) //bxdraw,endraw,mgdraw,usdraw
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
// TRACKR.dtrack = energy deposition of the jth deposition event

  // If coming from bxdraw we have 2 steps of 0 length and 0 edep
  Int_t caller = GetCaller();
  if (caller == 11 || caller==12) return 0.0;
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
  if (caller != 2) { // not eedraw 
      return PDGFromId(TRACKR.jtrack);
  }
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
// Return true for the first call of Stepping()
   return fTrackIsNew;
}

void     TFluka::SetTrackIsNew(Bool_t flag)
{
// Return true for the first call of Stepping()
   fTrackIsNew = flag;

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
  if (caller == 11 || caller==12)  // bxdraw
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
// True if track is exiting volume
//
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
      fIcode == 103 || // delta ray generation by hadron
      fIcode == 104 || // direct pair production
      fIcode == 105 || // bremsstrahlung (muon)
      fIcode == 208 || // bremsstrahlung (electron)
      fIcode == 214 || // in-flight annihilation
      fIcode == 215 || // annihilation at rest
      fIcode == 217 || // pair production
      fIcode == 219 || // Compton scattering
      fIcode == 221 || // Photoelectric effect
      fIcode == 300 || // hadronic interaction
      fIcode == 400    // delta-ray
      ) return 1;
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

{
// Number of secondary particles generated in the current step
// FINUC.np = number of secondaries except light and heavy ions
// FHEAVY.npheav = number of secondaries for light and heavy secondary ions
    Int_t caller = GetCaller();
    if (caller == 6)  // valid only after usdraw
	return FINUC.np + FHEAVY.npheav;
    else if (caller == 50) {
	// Cerenkov Photon production
	return fNCerenkov;
    }
    return 0;
} // end of NSecondaries

//______________________________________________________________________________ 
void TFluka::GetSecondary(Int_t isec, Int_t& particleId,
		TLorentzVector& position, TLorentzVector& momentum)
{
// Copy particles from secondary stack to vmc stack
//

    Int_t caller = GetCaller();
    if (caller == 6) {  // valid only after usdraw
	if (FINUC.np > 0) {
	    // Hadronic interaction
	    if (isec >= 0 && isec < FINUC.np) {
		particleId = PDGFromId(FINUC.kpart[isec]);
		position.SetX(fXsco);
		position.SetY(fYsco);
		position.SetZ(fZsco);
		position.SetT(TRACKR.atrack);
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
    } else if (caller == 50) {
	Int_t index = OPPHST.lstopp - isec;
	position.SetX(OPPHST.xoptph[index]);
	position.SetY(OPPHST.yoptph[index]);
	position.SetZ(OPPHST.zoptph[index]);
	position.SetT(OPPHST.agopph[index]);
	Double_t p = OPPHST.poptph[index];
	
	momentum.SetPx(p * OPPHST.txopph[index]);
	momentum.SetPy(p * OPPHST.tyopph[index]);
	momentum.SetPz(p * OPPHST.tzopph[index]);
	momentum.SetE(p);
    }
    else
	Warning("GetSecondary","no secondaries available");
    
} // end of GetSecondary


//______________________________________________________________________________ 
TMCProcess TFluka::ProdProcess(Int_t) const

{
// Name of the process that has produced the secondary particles
// in the current step

    Int_t mugamma = (TRACKR.jtrack == 7 || TRACKR.jtrack == 10 || TRACKR.jtrack == 11);

    if      (fIcode == 102)                  return kPDecay;
    else if (fIcode == 104 || fIcode == 217) return kPPair;
    else if (fIcode == 219)                  return kPCompton;
    else if (fIcode == 221)                  return kPPhotoelectric;
    else if (fIcode == 105 || fIcode == 208) return kPBrem;
    else if (fIcode == 103 || fIcode == 400) return kPDeltaRay;
    else if (fIcode == 210 || fIcode == 212) return kPDeltaRay;
    else if (fIcode == 214 || fIcode == 215) return kPAnnihilation;
    else if (fIcode == 101)                  return kPHadronic;
    else if (fIcode == 101) {
	if (!mugamma)                        return kPHadronic;
	else if (TRACKR.jtrack == 7)         return kPPhotoFission;
	else return kPMuonNuclear;
    }
    else if (fIcode == 225)                  return kPRayleigh;
// Fluka codes 100, 300 and 400 still to be investigasted
    else                                     return kPNoProcess;
}


Int_t TFluka::StepProcesses(TArrayI &proc) const
{
  //
  // Return processes active in the current step
  //
    proc.Set(1);
    TMCProcess iproc;
    switch (fIcode) {
    case 15:
    case 24:
    case 33:
    case 41:
    case 52:
	iproc =  kPTOFlimit;
	break;
    case 12:
    case 14:
    case 21:
    case 22:
    case 23:
    case 31:
    case 32:
    case 40:
    case 51:
	iproc =  kPStop;
	break;
    case 50:
	iproc = kPLightAbsorption;
	break;
    case 20: 
	iproc = kPPhotoelectric;
	break;
    default:
	iproc = ProdProcess(0);
    }
    proc[0] = iproc;
    return 1;
}
//______________________________________________________________________________ 
Int_t TFluka::VolId2Mate(Int_t id) const
{
//
// Returns the material number for a given volume ID
//
   return fMCGeo->VolId2Mate(id);
}

//______________________________________________________________________________ 
const char* TFluka::VolName(Int_t id) const
{
//
// Returns the volume name for a given volume ID
//
   return fMCGeo->VolName(id);
}

//______________________________________________________________________________ 
Int_t TFluka::VolId(const Text_t* volName) const
{
//
// Converts from volume name to volume ID.
// Time consuming. (Only used during set-up)
// Could be replaced by hash-table
//
   return fMCGeo->VolId(volName);
}

//______________________________________________________________________________ 
Int_t TFluka::CurrentVolID(Int_t& copyNo) const
{
//
// Return the logical id and copy number corresponding to the current fluka region
//
  if (gGeoManager->IsOutside()) return 0;
  TGeoNode *node = gGeoManager->GetCurrentNode();
  copyNo = node->GetNumber();
  Int_t id = node->GetVolume()->GetNumber();
  return id;
} 

//______________________________________________________________________________ 
Int_t TFluka::CurrentVolOffID(Int_t off, Int_t& copyNo) const
{
//
// Return the logical id and copy number of off'th mother 
// corresponding to the current fluka region
//
  if (off<0 || off>gGeoManager->GetLevel()) return 0;
  if (off==0) return CurrentVolID(copyNo);
  TGeoNode *node = gGeoManager->GetMother(off);
  if (!node) return 0;
  copyNo = node->GetNumber();
  return node->GetVolume()->GetNumber();
}

//______________________________________________________________________________ 
const char* TFluka::CurrentVolName() const
{
//
// Return the current volume name
//
  if (gGeoManager->IsOutside()) return 0;
  return gGeoManager->GetCurrentVolume()->GetName();
}

//______________________________________________________________________________ 
const char* TFluka::CurrentVolOffName(Int_t off) const
{
//
// Return the volume name of the off'th mother of the current volume
//
  if (off<0 || off>gGeoManager->GetLevel()) return 0;
  if (off==0) return CurrentVolName();
  TGeoNode *node = gGeoManager->GetMother(off);
  if (!node) return 0;
  return node->GetVolume()->GetName();
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
   Double_t xmL[3], xdL[3];
   Int_t i;
   for (i=0;i<3;i++) xmL[i]=xm[i];
   if (iflag == 1) gGeoManager->MasterToLocal(xmL,xdL);
   else            gGeoManager->MasterToLocalVect(xmL,xdL);
   for (i=0;i<3;i++) xd[i] = xdL[i];
}
  
//______________________________________________________________________________ 
void TFluka::Gmtod(Double_t* xm, Double_t* xd, Int_t iflag)
{
   if (iflag == 1) gGeoManager->MasterToLocal(xm,xd);
   else            gGeoManager->MasterToLocalVect(xm,xd);
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
   Double_t xmL[3], xdL[3];
   Int_t i;
   for (i=0;i<3;i++) xdL[i] = xd[i];
   if (iflag == 1) gGeoManager->LocalToMaster(xdL,xmL);
   else            gGeoManager->LocalToMasterVect(xdL,xmL);
   for (i=0;i<3;i++) xm[i]=xmL[i];
}

//______________________________________________________________________________ 
void TFluka::Gdtom(Double_t* xd, Double_t* xm, Int_t iflag)
{
   if (iflag == 1) gGeoManager->LocalToMaster(xd,xm);
   else            gGeoManager->LocalToMasterVect(xd,xm);
}

//______________________________________________________________________________
TObjArray *TFluka::GetFlukaMaterials()
{
   return fGeom->GetMatList();
}   

//______________________________________________________________________________
void TFluka::SetMreg(Int_t l) 
{
// Set current fluka region
   fCurrentFlukaRegion = l;
   fGeom->SetMreg(l);
}




TString TFluka::ParticleName(Int_t pdg) const
{
    // Return particle name for particle with pdg code pdg.
    Int_t ifluka = IdFromPDG(pdg);
    return TString((CHPPRP.btype[ifluka+6]), 8);
}
 

Double_t TFluka::ParticleMass(Int_t pdg) const
{
    // Return particle mass for particle with pdg code pdg.
    Int_t ifluka = IdFromPDG(pdg);
    return (PAPROP.am[ifluka+6]);
}

Double_t TFluka::ParticleCharge(Int_t pdg) const
{
    // Return particle charge for particle with pdg code pdg.
    Int_t ifluka = IdFromPDG(pdg);
    return Double_t(PAPROP.ichrge[ifluka+6]);
}

Double_t TFluka::ParticleLifeTime(Int_t pdg) const
{
    // Return particle lifetime for particle with pdg code pdg.
    Int_t ifluka = IdFromPDG(pdg);
    return (PAPROP.thalf[ifluka+6]);
}

void TFluka::Gfpart(Int_t pdg, char* name, Int_t& type, Float_t& mass, Float_t& charge, Float_t& tlife)
{
    // Retrieve particle properties for particle with pdg code pdg.
    
    strcpy(name, ParticleName(pdg).Data());
    type   = ParticleMCType(pdg);
    mass   = ParticleMass(pdg);
    charge = ParticleCharge(pdg);
    tlife  = ParticleLifeTime(pdg);
}



#define pushcerenkovphoton pushcerenkovphoton_
#define usersteppingckv    usersteppingckv_


extern "C" {
    void pushcerenkovphoton(Double_t & px, Double_t & py, Double_t & pz, Double_t & e,
			    Double_t & vx, Double_t & vy, Double_t & vz, Double_t & tof,
			    Double_t & polx, Double_t & poly, Double_t & polz, Double_t & wgt, Int_t& ntr)
    {
	//
	// Pushes one cerenkov photon to the stack
	//
	
	TFluka* fluka =  (TFluka*) gMC;
	TVirtualMCStack* cppstack = fluka->GetStack();
	Int_t parent =  TRACKR.ispusr[mkbmx2-1];
	cppstack->PushTrack(0, parent, 50000050,
			    px, py, pz, e,
                            vx, vy, vz, tof,
			    polx, poly, polz,
			    kPCerenkov, ntr, wgt, 0); 
    }

    void usersteppingckv(Int_t & nphot, Int_t & mreg, Double_t & x, Double_t & y, Double_t & z)
    {
	//
	// Calls stepping in order to signal cerenkov production
	//
	TFluka *fluka = (TFluka*)gMC;
	fluka->SetMreg(mreg);
	fluka->SetXsco(x);
	fluka->SetYsco(y);
	fluka->SetZsco(z);
	fluka->SetNCerenkov(nphot);
	fluka->SetCaller(50);
	printf("userstepping ckv: %10d %10d %13.3f %13.3f %13.2f\n", nphot, mreg, x, y, z);
	(TVirtualMCApplication::Instance())->Stepping();
    }
}

