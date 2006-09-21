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
#include "TFlukaCodes.h"
#include "TCallf77.h"      //For the fortran calls
#include "Fdblprc.h"       //(DBLPRC) fluka common
#include "Fsourcm.h"       //(SOURCM) fluka common
#include "Fgenstk.h"       //(GENSTK)  fluka common
#include "Fiounit.h"       //(IOUNIT) fluka common
#include "Fpaprop.h"       //(PAPROP) fluka common
#include "Fpart.h"         //(PART)   fluka common
#include "Ftrackr.h"       //(TRACKR) fluka common
#include "Fpaprop.h"       //(PAPROP) fluka common
#include "Ffheavy.h"       //(FHEAVY) fluka common
#include "Fopphst.h"       //(OPPHST) fluka common
#include "Fflkstk.h"       //(FLKSTK) fluka common
#include "Fstepsz.h"       //(STEPSZ) fluka common
#include "Fopphst.h"       //(OPPHST) fluka common
#include "Fltclcm.h"       //(LTCLCM) fluka common
#include "Falldlt.h"       //(ALLDLT) fluka common

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
#include "TArrayD.h"
#include "TDatabasePDG.h"
#include "TStopwatch.h"


// Fluka methods that may be needed.
#ifndef WIN32 
# define flukam  flukam_
# define fluka_openinp fluka_openinp_
# define fluka_openout fluka_openout_
# define fluka_closeinp fluka_closeinp_
# define mcihad mcihad_
# define mpdgha mpdgha_
# define newplo newplo_
# define genout genout_
# define flkend flkend_
#else 
# define flukam  FLUKAM
# define fluka_openinp FLUKA_OPENINP
# define fluka_openout FLUKA_OPENOUT
# define fluka_closeinp FLUKA_CLOSEINP
# define mcihad MCIHAD
# define mpdgha MPDGHA
# define newplo NEWPLO
# define genout GENOUT
# define flkend FLKEND
#endif

extern "C" 
{
  //
  // Prototypes for FLUKA functions
  //
  void type_of_call flukam(const int&);
  void type_of_call newplo();
  void type_of_call genout();
  void type_of_call flkend();
  void type_of_call fluka_openinp(const int&, DEFCHARA);
  void type_of_call fluka_openout(const int&, DEFCHARA);
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
   fNEvent(0),
   fInputFileName(""),
   fCoreInputFileName(""),
   fCaller(kNoCaller),
   fIcode(kNoProcess),
   fNewReg(-1),
   fRull(0),
   fXsco(0),
   fYsco(0),
   fZsco(0),
   fTrackIsEntering(kFALSE),
   fTrackIsExiting(kFALSE),
   fTrackIsNew(kFALSE),
   fFieldFlag(kTRUE),
   fGeneratePemf(kFALSE),
   fDummyBoundary(kFALSE),
   fStopped(kFALSE),
   fStopEvent(kFALSE),
   fStopRun(kFALSE),
   fMaterials(0),
   fNVolumes(0),
   fCurrentFlukaRegion(-1),
   fNCerenkov(0),
   fGeom(0),
   fMCGeo(0),
   fUserConfig(0), 
   fUserScore(0)
{ 
  //
  // Default constructor
  //
} 
 
//______________________________________________________________________________ 
TFluka::TFluka(const char *title, Int_t verbosity, Bool_t isRootGeometrySupported)
  :TVirtualMC("TFluka",title, isRootGeometrySupported),
   fVerbosityLevel(verbosity),
   fNEvent(0),
   fInputFileName(""),
   fCoreInputFileName(""),
   fCaller(kNoCaller),
   fIcode(kNoProcess),
   fNewReg(-1),
   fRull(0),
   fXsco(0),
   fYsco(0),
   fZsco(0),
   fTrackIsEntering(kFALSE),
   fTrackIsExiting(kFALSE),
   fTrackIsNew(kFALSE),
   fFieldFlag(kTRUE),
   fGeneratePemf(kFALSE),
   fDummyBoundary(kFALSE),
   fStopped(kFALSE),
   fStopEvent(kFALSE),
   fStopRun(kFALSE),
   fMaterials(0),
   fNVolumes(0),
   fCurrentFlukaRegion(-1),
   fNCerenkov(0),
   fGeom(0),
   fMCGeo(0),
   fUserConfig(new TObjArray(100)),
   fUserScore(new TObjArray(100)) 
{
  // create geometry interface
   if (fVerbosityLevel >=3)
       cout << "<== TFluka::TFluka(" << title << ") constructor called." << endl;
   SetCoreInputFileName();
   SetInputFileName();
   fMCGeo = new TGeoMCGeometry("MCGeo", "TGeo Implementation of VirtualMCGeometry", kFALSE);
   fGeom  = new TFlukaMCGeometry("geom", "FLUKA VMC Geometry");
   if (verbosity > 2) fGeom->SetDebugMode(kTRUE);
   PrintHeader();
}

//______________________________________________________________________________ 
TFluka::~TFluka()
{
    // Destructor
    if (fVerbosityLevel >=3)
        cout << "<== TFluka::~TFluka() destructor called." << endl;
    if (fMaterials) delete [] fMaterials;
    
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
    if (!gGeoManager->IsClosed()) {
       TGeoVolume *top = (TGeoVolume*)gGeoManager->GetListOfVolumes()->First();
       gGeoManager->SetTopVolume(top);
       gGeoManager->CloseGeometry("di");
    } else {
       TGeoNodeCache *cache = gGeoManager->GetCache();
       if (!cache->HasIdArray()) {
          Warning("Init", "Node ID tracking must be enabled with TFluka: enabling...\n");
          cache->BuildIdArray();
       }   
    }           
    fNVolumes = fGeom->NofVolumes();
    fGeom->CreateFlukaMatFile("flukaMat.inp");   
    if (fVerbosityLevel >=3) {
       printf("== Number of volumes: %i\n ==", fNVolumes);
       cout << "\t* InitPhysics() - Prepare input file to be called" << endl; 
    }

    fApplication->InitGeometry();

    //
    // Add ions to PDG Data base
    //
     AddParticlesToPdgDataBase();
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
//  Open fortran files    
    const char* fname = fInputFileName;
    fluka_openinp(lunin, PASSCHARA(fname));
    fluka_openout(11, PASSCHARA("fluka.out"));
//  Read input cards    
    cout << "==> TFluka::BuildPhysics() Read input cards." << endl;
    TStopwatch timer;
    timer.Start();
    GLOBAL.lfdrtr = true;
    flukam(1);
    cout << "<== TFluka::BuildPhysics() Read input cards End"
         << Form(" R:%.2fs C:%.2fs", timer.RealTime(),timer.CpuTime()) << endl;
//  Close input file
    fluka_closeinp(lunin);
//  Finish geometry    
    FinishGeometry();
}  

//______________________________________________________________________________ 
void TFluka::ProcessEvent() {
//
// Process one event
//
    if (fStopRun) {
        Warning("ProcessEvent", "User Run Abortion: No more events handled !\n");
        fNEvent += 1;
        return;
    }

    if (fVerbosityLevel >=3)
        cout << "==> TFluka::ProcessEvent() called." << endl;
    fApplication->GeneratePrimaries();
    SOURCM.lsouit = true;
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

  Int_t todo = TMath::Abs(nevent);
  for (Int_t ev = 0; ev < todo; ev++) {
      TStopwatch timer;
      timer.Start();
      fApplication->BeginEvent();
      ProcessEvent();
      fApplication->FinishEvent();
      cout << "Event: "<< ev
           << Form(" R:%.2fs C:%.2fs", timer.RealTime(),timer.CpuTime()) << endl;
  }

  if (fVerbosityLevel >=3)
    cout << "<== TFluka::ProcessRun(" << nevent << ") called." 
         << endl;
  
  // Write fluka specific scoring output
  genout();
  newplo();
  flkend();
  
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
//
   Bool_t process = kFALSE;
   Bool_t modelp  = kFALSE;
   
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
   
   if (strncmp(param, "PRIMIO_N",  8) == 0 ||
       strncmp(param, "PRIMIO_E",  8) == 0)
   {
       modelp = kTRUE;
   }
   
   if (process) {
       // Process switch
       SetProcess(param, Int_t (parval), itmed);
   } else if (modelp) {
       // Model parameters
       SetModelParameter(param, parval, itmed);
   } else {
       // Cuts
       SetCut(param, parval, itmed);
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

//______________________________________________________________________
Bool_t TFluka::GetTransformation(const TString &volumePath,TGeoHMatrix &mat)
{
    // Returns the Transformation matrix between the volume specified
    // by the path volumePath and the Top or mater volume. The format
    // of the path volumePath is as follows (assuming ALIC is the Top volume)
    // "/ALIC_1/DDIP_1/S05I_2/S05H_1/S05G_3". Here ALIC is the top most
    // or master volume which has only 1 instance of. Of all of the daughter
    // volumes of ALICE, DDIP volume copy #1 is indicated. Similarly for
    // the daughter volume of DDIP is S05I copy #2 and so on.
    // Inputs:
    //   TString& volumePath  The volume path to the specific volume
    //                        for which you want the matrix. Volume name
    //                        hierarchy is separated by "/" while the
    //                        copy number is appended using a "_".
    // Outputs:
    //  TGeoHMatrix &mat      A matrix with its values set to those
    //                        appropriate to the Local to Master transformation
    // Return:
    //   A logical value if kFALSE then an error occurred and no change to
    //   mat was made.

   // We have to preserve the modeler state
   return fMCGeo->GetTransformation(volumePath, mat);
}   
   
//______________________________________________________________________
Bool_t TFluka::GetShape(const TString &volumePath,TString &shapeType,
                        TArrayD &par)
{
    // Returns the shape and its parameters for the volume specified
    // by volumeName.
    // Inputs:
    //   TString& volumeName  The volume name
    // Outputs:
    //   TString &shapeType   Shape type
    //   TArrayD &par         A TArrayD of parameters with all of the
    //                        parameters of the specified shape.
    // Return:
    //   A logical indicating whether there was an error in getting this
    //   information
   return fMCGeo->GetShape(volumePath, shapeType, par);
}
   
//______________________________________________________________________
Bool_t TFluka::GetMaterial(const TString &volumeName,
                            TString &name,Int_t &imat,
                            Double_t &a,Double_t &z,Double_t &dens,
                            Double_t &radl,Double_t &inter,TArrayD &par)
{
    // Returns the Material and its parameters for the volume specified
    // by volumeName.
    // Note, Geant3 stores and uses mixtures as an element with an effective
    // Z and A. Consequently, if the parameter Z is not integer, then
    // this material represents some sort of mixture.
    // Inputs:
    //   TString& volumeName  The volume name
    // Outputs:
    //   TSrting   &name       Material name
    //   Int_t     &imat       Material index number
    //   Double_t  &a          Average Atomic mass of material
    //   Double_t  &z          Average Atomic number of material
    //   Double_t  &dens       Density of material [g/cm^3]
    //   Double_t  &radl       Average radiation length of material [cm]
    //   Double_t  &inter      Average interaction length of material [cm]
    //   TArrayD   &par        A TArrayD of user defined parameters.
    // Return:
    //   kTRUE if no errors
   return fMCGeo->GetMaterial(volumeName,name,imat,a,z,dens,radl,inter,par);
}

//______________________________________________________________________
Bool_t TFluka::GetMedium(const TString &volumeName,TString &name,
                         Int_t &imed,Int_t &nmat,Int_t &isvol,Int_t &ifield,
                         Double_t &fieldm,Double_t &tmaxfd,Double_t &stemax,
                         Double_t &deemax,Double_t &epsil, Double_t &stmin,
                         TArrayD &par)
{
    // Returns the Medium and its parameters for the volume specified
    // by volumeName.
    // Inputs:
    //   TString& volumeName  The volume name.
    // Outputs:
    //   TString  &name       Medium name
    //   Int_t    &nmat       Material number defined for this medium
    //   Int_t    &imed       The medium index number
    //   Int_t    &isvol      volume number defined for this medium
    //   Int_t    &iflield    Magnetic field flag
    //   Double_t &fieldm     Magnetic field strength
    //   Double_t &tmaxfd     Maximum angle of deflection per step
    //   Double_t &stemax     Maximum step size
    //   Double_t &deemax     Maximum fraction of energy allowed to be lost
    //                        to continuous process.
    //   Double_t &epsil      Boundary crossing precision
    //   Double_t &stmin      Minimum step size allowed
    //   TArrayD  &par        A TArrayD of user parameters with all of the
    //                        parameters of the specified medium.
    // Return:
    //   kTRUE if there where no errors
   return fMCGeo->GetMedium(volumeName,name,imed,nmat,isvol,ifield,fieldm,tmaxfd,stemax,deemax,epsil,stmin,par);
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

void TFluka::SetCerenkov(Int_t itmed, Int_t npckov, Float_t* ppckov,
                         Float_t* absco, Float_t* effic, Float_t* rindex, Float_t* rfl) {
//
// Set Cerenkov properties for medium itmed
//
// npckov: number of sampling points
// ppckov: energy values
// absco:  absorption length
// effic:  quantum efficiency
// rindex: refraction index
// rfl:    reflectivity for boundary to medium itmed
//
//  
//  Create object holding Cerenkov properties
//  
    TFlukaCerenkov* cerenkovProperties = new TFlukaCerenkov(npckov, ppckov, absco, effic, rindex, rfl);
//
//  Pass object to medium
    TGeoMedium* medium = gGeoManager->GetMedium(itmed);
    medium->SetCerenkovProperties(cerenkovProperties);
}  


//______________________________________________________________________________ 
void TFluka::SetCerenkov(Int_t /*itmed*/, Int_t /*npckov*/, Double_t * /*ppckov*/,
                         Double_t * /*absco*/, Double_t * /*effic*/, Double_t * /*rindex*/) {
//
//  Double_t version not implemented
}  

void TFluka::SetCerenkov(Int_t /*itmed*/, Int_t /*npckov*/, Double_t* /*ppckov*/,
                         Double_t* /*absco*/, Double_t* /*effic*/, Double_t* /*rindex*/, Double_t* /*rfl*/) {
//
// //  Double_t version not implemented
}

// Euclid
//______________________________________________________________________________ 
void TFluka::WriteEuclid(const char* /*fileName*/, const char* /*topVol*/, 
                          Int_t /*number*/, Int_t /*nlevel*/) {
//
// Not with TGeo
   Warning("WriteEuclid", "Not implemented !");
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
Int_t TFluka::GetDummyRegion() const
{
// Returns index of the dummy region.
   return fGeom->GetDummyRegion();
}   

//____________________________________________________________________________ 
Int_t TFluka::GetDummyLattice() const
{
// Returns index of the dummy lattice.
   return fGeom->GetDummyLattice();
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
    if (pdg == 50000051) return (kFLUKAoptical);
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

    if (id == kFLUKAoptical) {
// Cerenkov photon
//        if (fVerbosityLevel >= 3)
//            printf("\n PDGFromId: Cerenkov Photon \n");
        return  50000050;
    }
// Error id    
    if (id == 0 || id < kFLUKAcodemin || id > kFLUKAcodemax) {
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
//        if (fVerbosityLevel >= 3)
//            printf("mpdgha called with %d %d \n", id, intfluka);
        return mpdgha(intfluka);
    } else {
        // ions and optical photons
        return idSpecial[id - kFLUKAcodemin];
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
void TFluka::SetModelParameter(const char* parName, Double_t parValue, Int_t imed)
{
// Set model parameter for material imed
//
    TIter next(fUserConfig);
    TFlukaConfigOption* proc;
    while((proc = (TFlukaConfigOption*)next()))
    { 
        if (proc->Medium() == imed) {
            proc->SetModelParameter(parName, parValue);
            return;
        }
    }

    proc = new TFlukaConfigOption(imed);
    proc->SetModelParameter(parName, parValue);
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
  Warning("Xsec", "Not yet implemented.!\n"); return -1.;
}


//______________________________________________________________________________ 
void TFluka::InitPhysics()
{
//
// Physics initialisation with preparation of FLUKA input cards
//
// Construct file names
    FILE *pFlukaVmcCoreInp, *pFlukaVmcFlukaMat, *pFlukaVmcInp;
    TString sFlukaVmcCoreInp = getenv("ALICE_ROOT");
    sFlukaVmcCoreInp +="/TFluka/input/";
    TString sFlukaVmcTmp = "flukaMat.inp";
    TString sFlukaVmcInp = GetInputFileName();
    sFlukaVmcCoreInp += GetCoreInputFileName();
    
// Open files 
    if ((pFlukaVmcCoreInp = fopen(sFlukaVmcCoreInp.Data(),"r")) == NULL) {
        Warning("InitPhysics", "\nCannot open file %s\n",sFlukaVmcCoreInp.Data());
        exit(1);
    }
    if ((pFlukaVmcFlukaMat = fopen(sFlukaVmcTmp.Data(),"r")) == NULL) {
        Warning("InitPhysics", "\nCannot open file %s\n",sFlukaVmcTmp.Data());
        exit(1);
    }
    if ((pFlukaVmcInp = fopen(sFlukaVmcInp.Data(),"w")) == NULL) {
        Warning("InitPhysics", "\nCannot open file %s\n",sFlukaVmcInp.Data());
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
    
    TFlukaScoringOption *mopo = 0;
    TFlukaScoringOption *mopi = 0;

    for (Int_t isc = 0; isc < nscore; isc++) 
    {
        mopo = dynamic_cast<TFlukaScoringOption*> (fUserScore->At(isc));
        char*    fileName = mopo->GetFileName();
        Int_t    size     = strlen(fileName);
        Float_t  lun      = -1.;
//
// Check if new output file has to be opened
        for (Int_t isci = 0; isci < isc; isci++) {

        
            mopi = dynamic_cast<TFlukaScoringOption*> (fUserScore->At(isci));
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

// Add RANDOMIZ card
    fprintf(pFlukaVmcInp,"RANDOMIZ  %10.1f%10.0f\n", 1., Float_t(gRandom->GetSeed()));
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
//    if (step > 1.e4) return;
    
//    Int_t mreg=0, latt=0;
//    fGeom->GetCurrentRegion(mreg, latt);
    Int_t mreg = fGeom->GetCurrentRegion();
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
  FlukaCallerCode_t caller = GetCaller();
  if (caller == kENDRAW    || caller == kUSDRAW || 
      caller == kBXExiting || caller == kBXEntering || 
      caller == kUSTCKV) { 
    position.SetX(GetXsco());
    position.SetY(GetYsco());
    position.SetZ(GetZsco());
    position.SetT(TRACKR.atrack);
  }
  else if (caller == kMGDRAW) { 
    position.SetX(TRACKR.xtrack[TRACKR.ntrack]);
    position.SetY(TRACKR.ytrack[TRACKR.ntrack]);
    position.SetZ(TRACKR.ztrack[TRACKR.ntrack]);
    position.SetT(TRACKR.atrack);
  }
  else if (caller == kSODRAW) { 
    position.SetX(TRACKR.xtrack[TRACKR.ntrack]);
    position.SetY(TRACKR.ytrack[TRACKR.ntrack]);
    position.SetZ(TRACKR.ztrack[TRACKR.ntrack]);
    position.SetT(0);
  } else if (caller == kMGResumedTrack) { 
    position.SetX(TRACKR.spausr[0]);
    position.SetY(TRACKR.spausr[1]);
    position.SetZ(TRACKR.spausr[2]);
    position.SetT(TRACKR.spausr[3]);
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
  FlukaCallerCode_t caller = GetCaller();
  if (caller == kENDRAW    || caller == kUSDRAW || 
      caller == kBXExiting || caller == kBXEntering || 
      caller == kUSTCKV) { 
    x = GetXsco();
    y = GetYsco();
    z = GetZsco();
  }
  else if (caller == kMGDRAW || caller == kSODRAW) { 
    x = TRACKR.xtrack[TRACKR.ntrack];
    y = TRACKR.ytrack[TRACKR.ntrack];
    z = TRACKR.ztrack[TRACKR.ntrack];
  }
  else if (caller == kMGResumedTrack) {
    x = TRACKR.spausr[0];
    y = TRACKR.spausr[1];
    z = TRACKR.spausr[2];
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
  FlukaCallerCode_t  caller = GetCaller();
  FlukaProcessCode_t icode  = GetIcode();
  
  if (caller != kEEDRAW && caller != kMGResumedTrack && 
      (caller != kENDRAW || (icode != kEMFSCOstopping1 && icode != kEMFSCOstopping2))) {
    if (TRACKR.ptrack >= 0) {
      momentum.SetPx(TRACKR.ptrack*TRACKR.cxtrck);
      momentum.SetPy(TRACKR.ptrack*TRACKR.cytrck);
      momentum.SetPz(TRACKR.ptrack*TRACKR.cztrck);
      momentum.SetE(TRACKR.etrack);
      return;
    }
    else {
      Double_t p = sqrt(TRACKR.etrack * TRACKR.etrack - ParticleMassFPC(TRACKR.jtrack) * ParticleMassFPC(TRACKR.jtrack));
      momentum.SetPx(p*TRACKR.cxtrck);
      momentum.SetPy(p*TRACKR.cytrck);
      momentum.SetPz(p*TRACKR.cztrck);
      momentum.SetE(TRACKR.etrack);
      return;
    }
  } else if  (caller == kMGResumedTrack) {
    momentum.SetPx(TRACKR.spausr[4]);
    momentum.SetPy(TRACKR.spausr[5]);
    momentum.SetPz(TRACKR.spausr[6]);
    momentum.SetE (TRACKR.spausr[7]);
    return;
  } else if (caller == kENDRAW && (icode == kEMFSCOstopping1 || icode == kEMFSCOstopping2)) {
      momentum.SetPx(0.);
      momentum.SetPy(0.);
      momentum.SetPz(0.);
      momentum.SetE(TrackMass());
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
  FlukaCallerCode_t   caller = GetCaller();
  FlukaProcessCode_t  icode  = GetIcode();
  if (caller != kEEDRAW && caller != kMGResumedTrack && 
      (caller != kENDRAW || (icode != kEMFSCOstopping1 && icode != kEMFSCOstopping2))) {
    if (TRACKR.ptrack >= 0) {
      px = TRACKR.ptrack*TRACKR.cxtrck;
      py = TRACKR.ptrack*TRACKR.cytrck;
      pz = TRACKR.ptrack*TRACKR.cztrck;
      e  = TRACKR.etrack;
      return;
    }
    else {
      Double_t p = sqrt(TRACKR.etrack * TRACKR.etrack - ParticleMassFPC(TRACKR.jtrack) *  ParticleMassFPC(TRACKR.jtrack));
      px = p*TRACKR.cxtrck;
      py = p*TRACKR.cytrck;
      pz = p*TRACKR.cztrck;
      e  = TRACKR.etrack;
      return;
    }
  } else if (caller == kMGResumedTrack) {
      px = TRACKR.spausr[4];
      py = TRACKR.spausr[5];
      pz = TRACKR.spausr[6];
      e  = TRACKR.spausr[7];
      return;
  } else if (caller == kENDRAW && (icode == kEMFSCOstopping1 || icode == kEMFSCOstopping2)) {
      px = 0.;
      py = 0.;
      pz = 0.;
      e  = TrackMass();
  }
  else
    Warning("TrackMomentum","momentum not available");
}

//______________________________________________________________________________ 
Double_t TFluka::TrackStep() const
{
// Return the length in centimeters of the current step
// TRACKR.ctrack = total curved path
  FlukaCallerCode_t caller = GetCaller();
  if (caller == kBXEntering || caller == kBXExiting || 
      caller == kENDRAW     || caller == kUSDRAW || 
      caller == kUSTCKV     || caller == kMGResumedTrack)
    return 0.0;
  else if (caller == kMGDRAW)
    return TRACKR.ctrack;
  else {
    Warning("TrackStep", "track step not available");
    return 0.0;
  }  
}

//______________________________________________________________________________ 
Double_t TFluka::TrackLength() const
{
// TRACKR.cmtrck = cumulative curved path since particle birth
  FlukaCallerCode_t caller = GetCaller();
  if (caller == kBXEntering || caller == kBXExiting || 
      caller == kENDRAW || caller == kUSDRAW || caller == kMGDRAW || 
      caller == kUSTCKV) 
    return TRACKR.cmtrck;
  else if (caller == kMGResumedTrack) 
    return TRACKR.spausr[8];
  else {
    Warning("TrackLength", "track length not available");
    return 0.0;
  } 
}

//______________________________________________________________________________ 
Double_t TFluka::TrackTime() const
{
// Return the current time of flight of the track being transported
// TRACKR.atrack = age of the particle
  FlukaCallerCode_t caller = GetCaller();
  if (caller == kBXEntering || caller == kBXExiting || 
      caller == kENDRAW     || caller == kUSDRAW    || caller == kMGDRAW || 
      caller == kUSTCKV)
    return TRACKR.atrack;
  else if (caller == kMGResumedTrack)
    return TRACKR.spausr[3];
  else {
    Warning("TrackTime", "track time not available");
    return 0.0;
  }   
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
  // If coming from usdraw we just signal particle production - no edep
  // If just first time after resuming, no edep for the primary
  FlukaCallerCode_t caller = GetCaller();
  if (caller == kBXExiting || caller == kBXEntering || 
      caller == kUSDRAW    || caller == kMGResumedTrack) return 0.0;
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
Int_t TFluka::CorrectFlukaId() const
{
   // since we don't put photons and e- created bellow transport cut on the vmc stack
   // and there is a call to endraw for energy deposition for each of them
   // and they have the track number of their parent, but different identity (pdg)
   // so we want to assign also their parent identity.
   if( (IsTrackStop() )
        && TRACKR.ispusr[mkbmx2 - 4] == TRACKR.ispusr[mkbmx2 - 1]
        && TRACKR.jtrack != TRACKR.ispusr[mkbmx2 - 3] ) {
      if (fVerbosityLevel >=3)
         cout << "CorrectFlukaId() for icode=" << GetIcode()
               << " track=" << TRACKR.ispusr[mkbmx2 - 1]
               << " current PDG=" << PDGFromId(TRACKR.jtrack)
               << " assign parent PDG=" << PDGFromId(TRACKR.ispusr[mkbmx2 - 3]) << endl;
      return TRACKR.ispusr[mkbmx2 - 3]; // assign parent identity
   }
   return TRACKR.jtrack;
}


//______________________________________________________________________________ 
Int_t TFluka::TrackPid() const
{
// Return the id of the particle transported
// TRACKR.jtrack = identity number of the particle
  FlukaCallerCode_t caller = GetCaller();
  if (caller != kEEDRAW) {
     return PDGFromId( CorrectFlukaId() );
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
  FlukaCallerCode_t caller = GetCaller();
  if (caller != kEEDRAW) 
     return PAPROP.ichrge[CorrectFlukaId()+6];
  else
    return -1000.0;
}

//______________________________________________________________________________ 
Double_t TFluka::TrackMass() const
{
// PAPROP.am = particle mass in GeV
// TRACKR.jtrack = identity number of the particle
  FlukaCallerCode_t caller = GetCaller();
  if (caller != kEEDRAW)
     return PAPROP.am[CorrectFlukaId()+6];
  else
    return -1000.0;
}

//______________________________________________________________________________ 
Double_t TFluka::Etot() const
{
// TRACKR.etrack = total energy of the particle
  FlukaCallerCode_t caller = GetCaller();
  if (caller != kEEDRAW)
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
  FlukaCallerCode_t caller = GetCaller();
  if (caller == kBXEntering || caller == kBXExiting)
    return 0;
  else
    return 1;
}

//______________________________________________________________________________ 
Bool_t   TFluka::IsTrackEntering() const
{
// True if this is the first step of the track in the current volume

  FlukaCallerCode_t caller = GetCaller();
  if (caller == kBXEntering)
    return 1;
  else return 0;
}

//______________________________________________________________________________ 
Bool_t   TFluka::IsTrackExiting() const
{
// True if track is exiting volume
//
  FlukaCallerCode_t caller = GetCaller();
  if (caller == kBXExiting)
    return 1;
  else return 0;
}

//______________________________________________________________________________ 
Bool_t   TFluka::IsTrackOut() const
{
// True if the track is out of the setup
// means escape
  FlukaProcessCode_t icode = GetIcode();
    
  if (icode == kKASKADescape ||
      icode == kEMFSCOescape ||
      icode == kKASNEUescape ||
      icode == kKASHEAescape ||
      icode == kKASOPHescape) 
       return 1;
  else return 0;
}

//______________________________________________________________________________ 
Bool_t   TFluka::IsTrackDisappeared() const
{
// All inelastic interactions and decays
// fIcode from usdraw
  FlukaProcessCode_t icode = GetIcode();
  if (icode == kKASKADinelint    || // inelastic interaction
      icode == kKASKADdecay      || // particle decay
      icode == kKASKADdray       || // delta ray generation by hadron
      icode == kKASKADpair       || // direct pair production
      icode == kKASKADbrems      || // bremsstrahlung (muon)
      icode == kEMFSCObrems      || // bremsstrahlung (electron)
      icode == kEMFSCOmoller     || // Moller scattering
      icode == kEMFSCObhabha     || // Bhaba scattering
      icode == kEMFSCOanniflight || // in-flight annihilation
      icode == kEMFSCOannirest   || // annihilation at rest
      icode == kEMFSCOpair       || // pair production
      icode == kEMFSCOcompton    || // Compton scattering
      icode == kEMFSCOphotoel    || // Photoelectric effect
      icode == kKASNEUhadronic   || // hadronic interaction
      icode == kKASHEAdray          // delta-ray
      ) return 1;
  else return 0;
}

//______________________________________________________________________________ 
Bool_t   TFluka::IsTrackStop() const
{
// True if the track energy has fallen below the threshold
// means stopped by signal or below energy threshold
  FlukaProcessCode_t icode = GetIcode();
  if (icode == kKASKADstopping  || // stopping particle
      icode == kKASKADtimekill  || // time kill 
      icode == kEMFSCOstopping1 || // below user-defined cut-off
      icode == kEMFSCOstopping2 || // below user cut-off
      icode == kEMFSCOtimekill  || // time kill
      icode == kKASNEUstopping  || // neutron below threshold
      icode == kKASNEUtimekill  || // time kill
      icode == kKASHEAtimekill  || // time kill
      icode == kKASOPHtimekill) return 1; // time kill
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
// GENSTK.np = number of secondaries except light and heavy ions
// FHEAVY.npheav = number of secondaries for light and heavy secondary ions
    FlukaCallerCode_t caller = GetCaller();
    if (caller == kUSDRAW)  // valid only after usdraw
        return GENSTK.np + FHEAVY.npheav;
    else if (caller == kUSTCKV) {
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

    FlukaCallerCode_t caller = GetCaller();
    if (caller == kUSDRAW) {  // valid only after usdraw
        if (GENSTK.np > 0) {
            // Hadronic interaction
            if (isec >= 0 && isec < GENSTK.np) {
                particleId = PDGFromId(GENSTK.kpart[isec]);
                position.SetX(fXsco);
                position.SetY(fYsco);
                position.SetZ(fZsco);
                position.SetT(TRACKR.atrack);
                momentum.SetPx(GENSTK.plr[isec]*GENSTK.cxr[isec]);
                momentum.SetPy(GENSTK.plr[isec]*GENSTK.cyr[isec]);
                momentum.SetPz(GENSTK.plr[isec]*GENSTK.czr[isec]);
                momentum.SetE(GENSTK.tki[isec] + PAPROP.am[GENSTK.kpart[isec]+6]);
            }
            else if (isec >= GENSTK.np && isec < GENSTK.np + FHEAVY.npheav) {
                Int_t jsec = isec - GENSTK.np;
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
    } else if (caller == kUSTCKV) {
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

    Int_t mugamma = (TRACKR.jtrack == kFLUKAphoton || 
                     TRACKR.jtrack == kFLUKAmuplus ||
                     TRACKR.jtrack == kFLUKAmuminus);
    FlukaProcessCode_t icode = GetIcode();

    if      (icode == kKASKADdecay)                                   return kPDecay;
    else if (icode == kKASKADpair || icode == kEMFSCOpair)            return kPPair;
    else if (icode == kEMFSCOcompton)                                 return kPCompton;
    else if (icode == kEMFSCOphotoel)                                 return kPPhotoelectric;
    else if (icode == kKASKADbrems      || icode == kEMFSCObrems)     return kPBrem;
    else if (icode == kKASKADdray       || icode == kKASHEAdray)      return kPDeltaRay;
    else if (icode == kEMFSCOmoller     || icode == kEMFSCObhabha)    return kPDeltaRay;
    else if (icode == kEMFSCOanniflight || icode == kEMFSCOannirest)  return kPAnnihilation;
    else if (icode == kKASKADinelint) {
        if (!mugamma)                                                 return kPHadronic;
        else if (TRACKR.jtrack == kFLUKAphoton)                       return kPPhotoFission;
        else                                                          return kPMuonNuclear;
    }
    else if (icode == kEMFSCOrayleigh)                                return kPRayleigh;
// Fluka codes 100, 300 and 400 still to be investigasted
    else                                                              return kPNoProcess;
}


Int_t TFluka::StepProcesses(TArrayI &proc) const
{
  //
  // Return processes active in the current step
  //
    FlukaProcessCode_t icode = GetIcode();
    proc.Set(1);
    TMCProcess iproc;
    switch (icode) {
    case kKASKADtimekill:
    case kEMFSCOtimekill:
    case kKASNEUtimekill:
    case kKASHEAtimekill:
    case kKASOPHtimekill:
        iproc =  kPTOFlimit;
        break;
    case kKASKADstopping:
    case kKASKADescape:
    case kEMFSCOstopping1:
    case kEMFSCOstopping2:
    case kEMFSCOescape:
    case kKASNEUstopping:
    case kKASNEUescape:
    case kKASHEAescape:
    case kKASOPHescape:
        iproc = kPStop;
        break;
    case kKASOPHabsorption:
        iproc = kPLightAbsorption;
        break;
    case kKASOPHrefraction:
        iproc = kPLightRefraction;
    case kEMSCOlocaledep : 
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

const char* TFluka::CurrentVolPath() {
  // Return the current volume path
  return gGeoManager->GetPath(); 
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
//
// See Gmtod(Float_t*, Float_t*, Int_t)
//
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
//
// See Gdtom(Float_t*, Float_t*, Int_t)
//
   if (iflag == 1) gGeoManager->LocalToMaster(xd,xm);
   else            gGeoManager->LocalToMasterVect(xd,xm);
}

//______________________________________________________________________________
TObjArray *TFluka::GetFlukaMaterials()
{
//
// Get array of Fluka materials
   return fGeom->GetMatList();
}   

//______________________________________________________________________________
void TFluka::SetMreg(Int_t l, Int_t lttc) 
{
// Set current fluka region
   fCurrentFlukaRegion = l;
   fGeom->SetMreg(l,lttc);
}




//______________________________________________________________________________
TString TFluka::ParticleName(Int_t pdg) const
{
    // Return particle name for particle with pdg code pdg.
    Int_t ifluka = IdFromPDG(pdg);
    return TString((CHPPRP.btype[ifluka - kFLUKAcodemin]), 8);
}
 

//______________________________________________________________________________
Double_t TFluka::ParticleMass(Int_t pdg) const
{
    // Return particle mass for particle with pdg code pdg.
    Int_t ifluka = IdFromPDG(pdg);
    return (PAPROP.am[ifluka - kFLUKAcodemin]);
}

//______________________________________________________________________________
Double_t TFluka::ParticleMassFPC(Int_t fpc) const
{
    // Return particle mass for particle with Fluka particle code fpc
    return (PAPROP.am[fpc - kFLUKAcodemin]);
}

//______________________________________________________________________________
Double_t TFluka::ParticleCharge(Int_t pdg) const
{
    // Return particle charge for particle with pdg code pdg.
    Int_t ifluka = IdFromPDG(pdg);
    return Double_t(PAPROP.ichrge[ifluka - kFLUKAcodemin]);
}

//______________________________________________________________________________
Double_t TFluka::ParticleLifeTime(Int_t pdg) const
{
    // Return particle lifetime for particle with pdg code pdg.
    Int_t ifluka = IdFromPDG(pdg);
    return (PAPROP.tmnlf[ifluka - kFLUKAcodemin]);
}

//______________________________________________________________________________
void TFluka::Gfpart(Int_t pdg, char* name, Int_t& type, Float_t& mass, Float_t& charge, Float_t& tlife)
{
    // Retrieve particle properties for particle with pdg code pdg.
    
    strcpy(name, ParticleName(pdg).Data());
    type   = ParticleMCType(pdg);
    mass   = ParticleMass(pdg);
    charge = ParticleCharge(pdg);
    tlife  = ParticleLifeTime(pdg);
}

//______________________________________________________________________________
void TFluka::PrintHeader()
{
    //
    // Print a header
    printf("\n");
    printf("\n");    
    printf("------------------------------------------------------------------------------\n");
    printf("- You are using the TFluka Virtual Monte Carlo Interface to FLUKA.           -\n");    
    printf("- Please see the file fluka.out for FLUKA output and licensing information.  -\n");    
    printf("------------------------------------------------------------------------------\n");
    printf("\n");
    printf("\n");    
}


#define pshckp pshckp_
#define ustckv ustckv_


extern "C" {
  void pshckp(Double_t & px, Double_t & py, Double_t & pz, Double_t & e,
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
    if (fluka->GetVerbosityLevel() >= 3)
            printf("pshckp: track=%d parent=%d lattc=%d %s\n", ntr, parent, TRACKR.lt1trk, fluka->CurrentVolName());
  }
    
    void ustckv(Int_t & nphot, Int_t & mreg, Double_t & x, Double_t & y, Double_t & z)
    {
        //
        // Calls stepping in order to signal cerenkov production
        //
        TFluka *fluka = (TFluka*)gMC;
        fluka->SetMreg(mreg, TRACKR.lt1trk); //LTCLCM.mlatm1);
        fluka->SetXsco(x);
        fluka->SetYsco(y);
        fluka->SetZsco(z);
        fluka->SetNCerenkov(nphot);
        fluka->SetCaller(kUSTCKV);
        if (fluka->GetVerbosityLevel() >= 3)
            printf("ustckv: %10d mreg=%d lattc=%d  newlat=%d (%f, %f, %f) edep=%f vol=%s\n",
                    nphot, mreg, TRACKR.lt1trk, LTCLCM.newlat, x, y, z, fluka->Edep(), fluka->CurrentVolName());
   
    // check region lattice consistency (debug Ernesto)
    // *****************************************************
   Int_t nodeId;
   Int_t volId = fluka->CurrentVolID(nodeId);
   Int_t crtlttc = gGeoManager->GetCurrentNodeId()+1;

   if( mreg != volId  && !gGeoManager->IsOutside() ) {
       cout << "  ustckv:   track=" << TRACKR.ispusr[mkbmx2-1] << " pdg=" << fluka->PDGFromId(TRACKR.jtrack)
            << " icode=" << fluka->GetIcode() << " gNstep=" << fluka->GetNstep() << endl
            << "               fluka   mreg=" << mreg << " mlttc=" << TRACKR.lt1trk << endl
            << "               TGeo   volId=" << volId << " crtlttc=" << crtlttc << endl
            << "     common TRACKR   lt1trk=" << TRACKR.lt1trk << " lt2trk=" << TRACKR.lt2trk << endl
            << "     common LTCLCM   newlat=" << LTCLCM.newlat << " mlatld=" <<  LTCLCM.mlatld << endl
            << "                     mlatm1=" << LTCLCM.mlatm1 << " mltsen=" <<  LTCLCM.mltsen << endl
            << "                     mltsm1=" << LTCLCM.mltsm1 << " mlattc=" << LTCLCM.mlattc << endl;
        if( TRACKR.lt1trk == crtlttc ) cout << "   *************************************************************" << endl;
    }
    // *****************************************************



        (TVirtualMCApplication::Instance())->Stepping();
    }
}

//______________________________________________________________________________
void TFluka::AddParticlesToPdgDataBase() const
{

//
// Add particles to the PDG data base

    TDatabasePDG *pdgDB = TDatabasePDG::Instance();

    const Int_t kion=10000000;

    const Double_t kAu2Gev   = 0.9314943228;
    const Double_t khSlash   = 1.0545726663e-27;
    const Double_t kErg2Gev  = 1/1.6021773349e-3;
    const Double_t khShGev   = khSlash*kErg2Gev;
    const Double_t kYear2Sec = 3600*24*365.25;
//
// Ions
//

  pdgDB->AddParticle("Deuteron","Deuteron",2*kAu2Gev+8.071e-3,kTRUE,
                     0,3,"Ion",kion+10020);
  pdgDB->AddParticle("Triton","Triton",3*kAu2Gev+14.931e-3,kFALSE,
                     khShGev/(12.33*kYear2Sec),3,"Ion",kion+10030);
  pdgDB->AddParticle("Alpha","Alpha",4*kAu2Gev+2.424e-3,kTRUE,
                     khShGev/(12.33*kYear2Sec),6,"Ion",kion+20040);
  pdgDB->AddParticle("HE3","HE3",3*kAu2Gev+14.931e-3,kFALSE,
                     0,6,"Ion",kion+20030);
}

//
// Info about primary ionization electrons
//

//______________________________________________________________________________
Int_t TFluka::GetNPrimaryElectrons()
{
    // Get number of primary electrons
    return ALLDLT.nalldl;
}

//______________________________________________________________________________
Double_t GetPrimaryElectronKineticEnergy(Int_t i)
{
    Double_t ekin = -1.;
    // Returns kinetic energy of primary electron i
    if (i >= 0 && i < ALLDLT.nalldl) {
        ekin =  ALLDLT.talldl[i];
    } else {
        Warning("GetPrimaryElectronKineticEnergy",
                "Primary electron index out of range %d %d \n",
                i, ALLDLT.nalldl);
    }
    return ekin;
}
