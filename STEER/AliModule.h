#ifndef ALIMODULE_H
#define ALIMODULE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//
// This is the basic class for any
// ALICE detector module, whether it is 
// sensitive or not. Detector classes depend
// on this.
//

#include <Riostream.h>
#include <TNamed.h>

#include "AliTriggerDetector.h"

class TClonesArray;
class TBrowser;
class TArrayI;
class TFile;
class TTree;
class AliLoader;
class AliTrackReference;
class AliDigitizer;
class AliRunDigitizer;
class AliVertexer;
class AliTracker;
class AliESD;
class AliRunLoader;
class AliRawReader;


class AliModule : public TNamed {
public:

  // Creators - distructors
  AliModule(const char* name, const char *title);
  AliModule();
  virtual ~AliModule();

  // Inline functions
  virtual  int           GetNdigits() const {return 0;}
  virtual  int           GetNhits()  const {return 0;}
  virtual  TArrayI      *GetIdtmed()   const {return fIdtmed;}
  virtual  TList        *Histograms() const {return fHistograms;}
  virtual  TList        *Nodes()  const {return fNodes;}
  virtual  TClonesArray *Digits() const {return 0;}
  virtual  TClonesArray *Hits()   const {return 0;}
  virtual  TObjArray    *Points() const {return 0;}
  virtual  Int_t         GetIshunt() const {return 0;}
  virtual  void          SetIshunt(Int_t) {}
  virtual  Bool_t        IsActive() const {return fActive;}
  virtual  Bool_t        IsFolder() const {return kTRUE;}
  virtual  Int_t&        LoMedium() {return fLoMedium;}
  virtual  Int_t&        HiMedium() {return fHiMedium;}

  // Module composition
  virtual void AliMaterial(Int_t imat, const char* name, Float_t a, 
			   Float_t z, Float_t dens, Float_t radl,
			   Float_t absl, Float_t *buf=0, Int_t nwbuf=0) const;
  virtual void AliGetMaterial(Int_t imat, char* name, Float_t &a, 
			      Float_t &z, Float_t &dens, Float_t &radl,
			      Float_t &absl) const;
  virtual void AliMixture(Int_t imat, const char *name, Float_t *a,
                          Float_t *z, Float_t dens, Int_t nlmat,
                          Float_t *wmat) const;
  virtual void AliMedium(Int_t numed, const char *name, Int_t nmat,
                          Int_t isvol, Int_t ifield, Float_t fieldm,
                          Float_t tmaxfd, Float_t stemax, Float_t deemax,
                          Float_t epsil, Float_t stmin, Float_t *ubuf=0,
                          Int_t nbuf=0) const;
  virtual void AliMatrix(Int_t &nmat, Float_t theta1, Float_t phi1,
                          Float_t theta2, Float_t phi2, Float_t theta3,
                          Float_t phi3) const;
  
  // Virtual methods
  virtual void  BuildGeometry() {};
  virtual Int_t IsVersion() const =0;
  
  
  // Other methods
  virtual void        AddDigit(Int_t*, Int_t*){
  Error("AddDigit","Digits cannot be added to module %s\n",fName.Data());}
  virtual void        AddHit(Int_t, Int_t*, Float_t *) {
  Error("AddDigit","Hits cannot be added to module %s\n",fName.Data());}
  virtual void        Hits2SDigits() {}
  virtual AliDigitizer* CreateDigitizer(AliRunDigitizer* /*manager*/) const 
    {return NULL;}
  virtual AliTriggerDetector* CreateTriggerDetector() const
    { AliTriggerDetector* det = new AliTriggerDetector(); det->SetName(GetName()); return det;}
  virtual void        SDigits2Digits() {}
  virtual void        Hits2Digits() {}
  virtual void        Digits2Reco() {}
  virtual void        Digits2Raw();
  virtual void        Raw2Digits()  {}
  virtual Bool_t      Raw2SDigits(AliRawReader*) {return kFALSE;}
  virtual void        Browse(TBrowser *) {} //PH Do we need it?
  virtual void        CreateGeometry() {}
  virtual void        CreateMaterials() {}
  virtual void        AddAlignableVolumes() const;
  virtual void        Disable();
  virtual void        Enable();
  virtual void        PreTrack(){}
  virtual void        PostTrack(){}
  virtual void        FinishEvent() {}
  virtual void        FinishRun() {}
  virtual void        FinishPrimary() {}
  virtual void        Init() {}
  virtual void        LoadPoints(Int_t ) {}
  virtual void        UpdateInternalGeometry() {}



  virtual void        MakeBranch(Option_t * /*opt =" "*/) {}
  virtual void        MakeTree(Option_t *) {}//skowron 

  virtual AliLoader*  MakeLoader(const char* topfoldername);  
  virtual AliLoader*  GetLoader() const {return 0x0;} //skowron
  

  virtual void        Paint(Option_t *) {} //PH Do we need it?
  virtual void        ResetDigits() {}
  virtual void        ResetSDigits() {}
  virtual void        ResetHits() {}
  virtual void        ResetPoints() {}
  virtual void        SetTreeAddress();
  virtual void        SetTimeGate(Float_t) {}
  virtual Float_t     GetTimeGate() const {return 1.e10;}
  virtual void        StepManager() {}
  virtual void        DisableStepManager() {fEnable = kFALSE;}
  virtual Bool_t      StepManagerIsEnabled() const {return fEnable;}
  virtual void        SetBufferSize(Int_t) {}  
  virtual Float_t     ZMin() const;
  virtual Float_t     ZMax() const;
  virtual void        SetEuclidFile(char *material,char *geometry=0);
  virtual void        ReadEuclid(const char *filnam, char *topvol);
  virtual void        ReadEuclidMedia(const char *filnam);
// Track reference related
  TClonesArray *TrackReferences()   const {return fTrackReferences;}
  virtual void        RemapTrackHitIDs(Int_t *) {}
  virtual void        RemapTrackReferencesIDs(Int_t *map); //remaping track references MI
  virtual void        ResetTrackReferences();
  virtual  AliTrackReference * AddTrackReference(Int_t label);
  virtual  AliTrackReference * FirstTrackReference(Int_t track);
  virtual  AliTrackReference * NextTrackReference();
  virtual void        MakeBranchTR(Option_t *opt=" ");
  TTree* TreeTR();  //shorcut method for accessing treeTR from folder

  void                SetRunLoader(AliRunLoader* runLoader) 
    {fRunLoader = runLoader;}
  
//
 
protected:      

  // Data members
  
  TString       fEuclidMaterial;  //!Name of the Euclid file for materials (if any)
  TString       fEuclidGeometry;  //!Name of the Euclid file for geometry (if any)
  
  TArrayI      *fIdtmed;      //!List of tracking medium numbers
  TArrayI      *fIdmate;      //!List of material numbers
  Int_t         fLoMedium;   //!Minimum tracking medium ID for this Module
  Int_t         fHiMedium;   //!Maximum tracking medium ID for this Module

  Bool_t        fActive;      //Detector activity flag
  TList        *fHistograms;  //List of histograms
  TList        *fNodes;       //List of geometry nodes
  Bool_t        fEnable;      //StepManager enabling flag
  TClonesArray *fTrackReferences;     //!list of track references - for one primary track only -MI
  Int_t         fMaxIterTrackRef;     //!for track refernce iterator routines
  Int_t         fCurrentIterTrackRef; //!for track refernce iterator routines

  AliRunLoader* fRunLoader;   //!local pointer to run loader

 private:
  AliModule(const AliModule &mod);
  AliModule& operator=(const AliModule &mod);

  ClassDef(AliModule,7)  //Base class for ALICE Modules
};
#endif
