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
# define newplo newplo_
#else 
# define flukam  FLUKAM
# define fluka_openinp FLUKA_OPENINP
# define fluka_closeinp FLUKA_CLOSEINP
# define mcihad MCIHAD
# define mpdgha MPDGHA
# define newplo NEWPLO
#endif

extern "C" 
{
  //
  // Prototypes for FLUKA functions
  //
  void type_of_call flukam(const int&);
  void type_of_call newplo();
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
   fUserConfig(0), 
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
   fUserConfig(new TObjArray(100)),
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
   fGeom  = new TFlukaMCGeometry("geom", "FLUKA VMC Geometry");
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
    
    if (fUserConfig) {
	fUserConfig->Delete();
	delete fUserConfig;
    }
    
    if (fUserScore) {
	fUserScore->Delete();
	delete fUserScore;
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
    
    // Now we have TGeo geometry created and we have to patch FlukaVmc.inp
    // with the material mapping file FlukaMat.inp
 
    fNVolumes = fGeom->NofVolumes();
    fGeom->CreateFlukaMatFile("flukaMat.inp");   
    if (fVerbosityLevel >=3) {
       printf("== Number of volumes: %i\n ==", fNVolumes);
       cout << "\t* InitPhysics() - Prepare input file to be called" << endl; 
    }
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

    
    if (fVerbosityLevel >=3) {
	TList *medlist = gGeoManager->GetListOfMedia();
	TIter next(medlist);
	TGeoMedium*   med = 0x0;
	TGeoMaterial* mat = 0x0;
	Int_t ic = 0;
	
	while((med = (TGeoMedium*)next()))
	{
	    mat = med->GetMaterial();
	    printf("Medium %5d %12s %5d %5d\n", ic, (med->GetName()), med->GetId(), mat->GetIndex());
	    ic++;
	}
    }
    
    //
    // At this stage we have the information on materials and cuts available.
    // Now create the pemf file
    
    if (fGeneratePemf) fGeom->CreatePemfFile();
    
    //
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
  // Write fluka specific scoring output
  newplo();
  
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
// Define a material
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
// Define a material mixture
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
  // Define a medium
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
  // Define a medium
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
   if (fVerbosityLevel >= 3) 
       printf("Gstpar called with %6d %5s %12.4e %6d\n", itmed, param, parval, fGeom->GetFlukaMaterial(itmed));
   Int_t* reglist;
   Int_t nreg;
   reglist = fGeom->GetMaterialList(fGeom->GetFlukaMaterial(itmed), nreg);
   if (nreg == 0) {
       return;
   }
   
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
//    
//  Update if already in the list
//
    TIter next(fUserConfig);
    TFlukaConfigOption* proc;
    while((proc = (TFlukaConfigOption*)next()))
    { 
	if (proc->Medium() == imed) {
	    proc->SetProcess(flagName, flagValue);
	    return;
	}
    }
    proc = new TFlukaConfigOption(imed);
    proc->SetProcess(flagName, flagValue);
    fUserConfig->Add(proc);
}

//______________________________________________________________________________ 
Bool_t TFluka::SetProcess(const char* flagName, Int_t flagValue)
{
//  Set process user flag 
//
//    
    SetProcess(flagName, flagValue, -1);
    return kTRUE;  
}

//______________________________________________________________________________ 
void TFluka::SetCut(const char* cutName, Double_t cutValue, Int_t imed)
{
// Set user cut value for material imed
//
    TIter next(fUserConfig);
    TFlukaConfigOption* proc;
    while((proc = (TFlukaConfigOption*)next()))
    { 
	if (proc->Medium() == imed) {
	    proc->SetCut(cutName, cutValue);
	    return;
	}
    }

    proc = new TFlukaConfigOption(imed);
    proc->SetCut(cutName, cutValue);
    fUserConfig->Add(proc);
}

//______________________________________________________________________________ 
Bool_t TFluka::SetCut(const char* cutName, Double_t cutValue)
{
// Set user cut value 
//
//    
    SetCut(cutName, cutValue, -1);
    return kTRUE;
}


void TFluka::SetUserScoring(const char* option, Int_t npr, char* outfile, Float_t* what)
{
//
// Adds a user scoring option to the list
//
    TFlukaScoringOption* opt = new TFlukaScoringOption(option, "User Scoring", npr,outfile,what);
    fUserScore->Add(opt);
}
//______________________________________________________________________________
void TFluka::SetUserScoring(const char* option, Int_t npr, char* outfile, Float_t* what, const char* det1, const char* det2, const char* det3)
{
//
// Adds a user scoring option to the list
//
    TFlukaScoringOption* opt = new TFlukaScoringOption(option, "User Scoring", npr, outfile, what, det1, det2, det3);
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

// Construct file names
    FILE *pFlukaVmcCoreInp, *pFlukaVmcFlukaMat, *pFlukaVmcInp;
    TString sFlukaVmcCoreInp = getenv("ALICE_ROOT");
    sFlukaVmcCoreInp +="/TFluka/input/";
    TString sFlukaVmcTmp = "flukaMat.inp";
    TString sFlukaVmcInp = GetInputFileName();
    sFlukaVmcCoreInp += GetCoreInputFileName();
    
// Open files 
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

// Copy core input file 
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

    
// Pass information to configuration objects
    
    Float_t fLastMaterial = fGeom->GetLastMaterialIndex();
    TFlukaConfigOption::SetStaticInfo(pFlukaVmcInp, 3, fLastMaterial, fGeom);
    
    TIter next(fUserConfig);
    TFlukaConfigOption* proc;
    while((proc = dynamic_cast<TFlukaConfigOption*> (next()))) proc->WriteFlukaInputCards();
//
// Process Fluka specific scoring options
//
    TFlukaScoringOption::SetStaticInfo(pFlukaVmcInp, fGeom);
    Float_t loginp        = 49.0;
    Int_t inp             = 0;
    Int_t nscore          = fUserScore->GetEntries();
    
    TFlukaScoringOption *mopo = 0x0;
    TFlukaScoringOption *mopi = 0x0;

    for (Int_t isc = 0; isc < nscore; isc++) 
    {
	mopo = dynamic_cast<TFlukaScoringOption*> (fUserScore->At(isc));
	char*    fileName = mopo->GetFileName();
	Int_t    size     = strlen(fileName);
	Float_t  lun      = -1.;
//
// Check if new output file has to be opened
	for (Int_t isci = 0; isci < isc; isci++) {
	    mopi = dynamic_cast<TFlukaScoringOption*> (fUserScore->At(isc));
	    if(strncmp(mopi->GetFileName(), fileName, size)==0) {
		// 
		// No, the file already exists
		lun = mopi->GetLun();
		mopo->SetLun(lun);
		break;
	    }
	} // inner loop

	if (lun == -1.) {
	    // Open new output file
	    inp++;
	    mopo->SetLun(loginp + inp);
	    mopo->WriteOpenFlukaFile();
	}
	mopo->WriteFlukaInputCards();
    }
    
// Add START and STOP card
    fprintf(pFlukaVmcInp,"START     %10.1f\n",fEventsPerRun);
    fprintf(pFlukaVmcInp,"STOP      \n");
   
  
// Close files
   fclose(pFlukaVmcCoreInp);
   fclose(pFlukaVmcFlukaMat);
   fclose(pFlukaVmcInp);


//
// Initialisation needed for Cerenkov photon production and transport
    TObjArray *matList = GetFlukaMaterials();
    Int_t nmaterial =  matList->GetEntriesFast();
    fMaterials = new Int_t[nmaterial+3];
    
    for (Int_t im = 0; im < nmaterial; im++)
    {
	TGeoMaterial* material = dynamic_cast<TGeoMaterial*> (matList->At(im));
	Int_t idmat = material->GetIndex();
	fMaterials[idmat] = im;
    }
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
	iproc = kPStop;
	break;
    case 50:
	iproc = kPLightAbsorption;
	break;
    case 59:
	iproc = kPLightRefraction;
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
    char sname[20];
    Int_t len;
    strncpy(sname, volName, len = strlen(volName));
    sname[len] = 0;
    while (sname[len - 1] == ' ') sname[--len] = 0;
    return fMCGeo->VolId(sname);
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
Int_t TFluka::CurrentMaterial(Float_t & a, Float_t & z, 
		      Float_t & dens, Float_t & radl, Float_t & absl) const
{
//
//  Return the current medium number and material properties
//
  Int_t copy;
  Int_t id  =  TFluka::CurrentVolID(copy);
  Int_t med =  TFluka::VolId2Mate(id);
  TGeoVolume*     vol = gGeoManager->GetCurrentVolume();
  TGeoMaterial*   mat = vol->GetMaterial();
  a    = mat->GetA();
  z    = mat->GetZ();
  dens = mat->GetDensity();
  radl = mat->GetRadLen();
  absl = mat->GetIntLen();
  
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

