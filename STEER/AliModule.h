#ifndef ALIMODULE_H
#define ALIMODULE_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "TNamed.h"
#include "TAttLine.h"
#include "TAttMarker.h"
#include "AliRndm.h"

//#include <TSystem.h>
class TClonesArray;
class TBrowser;
class TArrayI;

class AliModule : public TNamed , public TAttLine, public TAttMarker,
                  public AliRndm {
public:

  // Creators - distructors
  AliModule(const char* name, const char *title);
  AliModule();
  AliModule(const AliModule &mod);
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
  virtual  Int_t         GetDebug() const {return fDebug;}
  virtual  void          SetDebug(Int_t deb=0) {fDebug=deb;}

  // Module composition
  virtual void AliMaterial(Int_t imat, const char* name, Float_t a, 
			   Float_t z, Float_t dens, Float_t radl,
			   Float_t absl, Float_t *buf=0, Int_t nwbuf=0) const;
  virtual void AliGetMaterial(Int_t imat, char* name, Float_t &a, 
				Float_t &z, Float_t &dens, Float_t &radl,
				Float_t &absl);
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
  virtual void        Browse(TBrowser *) {}
  virtual void        CreateGeometry() {}
  virtual void        CreateMaterials() {}
  virtual void        Disable();
  virtual Int_t       DistancetoPrimitive(Int_t px, Int_t py);
  virtual void        Enable();
  virtual void        PreTrack(){}
  virtual void        PostTrack(){}
  virtual void        FinishEvent() {}
  virtual void        FinishRun() {}
  virtual void        FinishPrimary() {}
  virtual void        RemapTrackHitIDs(Int_t *map) {}

  //virtual void        Hits2Digits() {}
  virtual void        Init() {}
  virtual void        LoadPoints(Int_t ) {}
  virtual void        MakeBranch(Option_t *) {}
  virtual void        Paint(Option_t *) {}
  virtual void        ResetDigits() {}
  virtual void        ResetHits() {}
  virtual void        ResetPoints() {}
  virtual void        SetTreeAddress() {}
  virtual void        SetTimeGate(Float_t) {}
  virtual Float_t     GetTimeGate() const {return 1.e10;}
  virtual void        StepManager() {}
  virtual void        SetBufferSize(Int_t) {}  
  virtual Float_t     ZMin() const;
  virtual Float_t     ZMax() const;
  virtual void        SetEuclidFile(char *material,char *geometry=0);
  virtual void ReadEuclid(const char *filnam, char *topvol);
  virtual void ReadEuclidMedia(const char *filnam);
  AliModule& operator=(const AliModule &mod);
  void Copy(AliModule &mod) const;
 
protected:      
  // Data members
  
  TString       fEuclidMaterial;  //Name of the Euclid file for materials (if any)
  TString       fEuclidGeometry;  //Name of the Euclid file for geometry (if any)
  
  TArrayI      *fIdtmed;      //List of tracking medium numbers
  TArrayI      *fIdmate;      //List of material numbers
  Int_t         fLoMedium;   //Minimum tracking medium ID for this Module
  Int_t         fHiMedium;   //Maximum tracking medium ID for this Module

  Bool_t        fActive;      //Detector activity flag
  TList        *fHistograms;  //List of histograms
  TList        *fNodes;       //List of geometry nodes
  Int_t         fDebug;       //Debug flag

  ClassDef(AliModule,1)  //Base class for ALICE Modules
};
#endif
