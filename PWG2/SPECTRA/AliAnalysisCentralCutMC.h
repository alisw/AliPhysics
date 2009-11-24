#ifndef ALIANALYSISCENTRALCUTMC_H
#define ALIANALYSISCENTRALCUTMC_H
/*
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved.
 * See cxx source for full Copyright notice
 * $Id$
 */

//  ***************************************************
//  * MC particle level cuts for azimuthal isotropic  *
//  * expansion in highly central collisions analysis *
//  * author: Cristian Andrei                         *
//  *         acristian@niham.nipne.ro                *
//  * *************************************************

#include "AliAnalysisCuts.h"

class TObject;
class TList;

class AliStack;
class AliMCEvent;
class AliMCParticle;

class AliAnalysisCentralCutMC: public AliAnalysisCuts{
 public:
  AliAnalysisCentralCutMC(const char *name="AliAnalysisCentralCutMC", const char *title="MC_cuts");
  virtual ~AliAnalysisCentralCutMC();

Bool_t  IsSelected(TObject* const obj);
Bool_t  IsSelected(TList* const /*list*/) {return kTRUE;}

  void SetOnlyPrimaries(Bool_t val = kFALSE) {fOnlyPrim = val;}

  void SetPDGCode(Int_t pdg = 0) {fPDGCode = pdg;}

  void SetPtRange(Float_t r1=0, Float_t r2=1e10)      { fPtMin=r1;  fPtMax=r2;}
  void SetEtaRange(Float_t r1=-1e10, Float_t r2=1e10) { fEtaMin=r1; fEtaMax=r2;}
  
  virtual void ReceiveEvt(TObject* mcEvent);

private:
  AliAnalysisCentralCutMC(const AliAnalysisCentralCutMC &ref);
  AliAnalysisCentralCutMC& operator=(const AliAnalysisCentralCutMC &ref);

  Bool_t fOnlyPrim; //kTRUE -> select only primary particles
  
  AliMCEvent *fMCEvt; //MC event is needed in order too get the Stack (for IsPrimary)

  Int_t fPDGCode; //the PDG code of the wanted particle

  Double_t   fPtMin,  fPtMax;       // definition of the range of the Pt
  Double_t   fEtaMin, fEtaMax;      // definition of the range of the eta

  
  Bool_t IsPrimary(AliMCParticle* const part, AliStack* const Stack);


  Bool_t CheckPDG(AliMCParticle* const mcPart, Int_t const pdg = 0);


  ClassDef(AliAnalysisCentralCutMC, 1);
};

#endif
