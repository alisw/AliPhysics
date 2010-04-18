#ifndef ALIRDHFCUTS_H
#define ALIRDHFCUTS_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
// Class AliRDHFCuts
// base class for cuts on AOD reconstructed heavy-flavour decays
// Author: A.Dainese, andrea.dainese@pd.infn.it
//***********************************************************

#include "AliAnalysisCuts.h"
#include "AliESDtrackCuts.h"
#include "AliVEvent.h"

class AliAODRecoDecayHF;

class AliRDHFCuts : public AliAnalysisCuts 
{
 public:

  AliRDHFCuts(const Char_t* name="RDHFCuts", const Char_t* title="");
  
  virtual ~AliRDHFCuts();
  
  AliRDHFCuts(const AliRDHFCuts& source);
  AliRDHFCuts& operator=(const AliRDHFCuts& source); 
 
  void SetMinVtxType(Int_t type=3) {fMinVtxType=type;}  
  void SetMinVtxContr(Int_t contr=1) {fMinVtxContr=contr;}  
  void SetMaxVtxRdChi2(Float_t chi2=1e6) {fMaxVtxRedChi2=chi2;}  
  void SetMinSPDMultiplicity(Int_t mult=0) {fMinSPDMultiplicity=mult;}  
  void SetTriggerMask(ULong64_t mask=0) {fTriggerMask=mask;} 
  void SetVarsForOpt(Int_t nVars,Bool_t *forOpt);
  void SetGlobalIndex(){fGlobalIndex=fnVars*fnPtBins;}
  void SetGlobalIndex(Int_t nVars,Int_t nptBins){fnVars=nVars; fnPtBins=nptBins; SetGlobalIndex();}
  void SetVarNames(Int_t nVars,TString *varNames,Bool_t *isUpperCut);  
  void SetPtBins(Int_t nPtBinLimits,Float_t *ptBinLimits);
  void SetCuts(Int_t nVars,Int_t nPtBins,Float_t** cutsRD);
  void SetCuts(Int_t glIndex, Float_t* cutsRDGlob);
  void AddTrackCuts(const AliESDtrackCuts *cuts) 
         {fTrackCuts=new AliESDtrackCuts(*cuts); return;}

  Float_t *GetPtBinLimits() const {return fPtBinLimits;}
  Int_t   GetNPtBins() const {return fnPtBins;}
  Int_t   GetNVars() const {return fnVars;} 
  TString *GetVarNames() const {return fVarNames;} 
  Bool_t  *GetVarsForOpt() const {return fVarsForOpt;} 
  Int_t   GetNVarsForOpt() const {return fnVarsForOpt;}
  const Float_t *GetCuts() const {return fCutsRD;} 
  void    GetCuts(Float_t**& cutsRD) const;
  Float_t GetCutValue(Int_t iVar,Int_t iPtBin) const;
  Bool_t  *GetIsUpperCut() const {return fIsUpperCut;}
  AliESDtrackCuts *GetTrackCuts() const {return fTrackCuts;}
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters) = 0;
  Int_t   GetGlobalIndex(Int_t iVar,Int_t iPtBin) const;
  void    GetVarPtIndex(Int_t iGlob, Int_t& iVar, Int_t& iPtBin) const;

  Bool_t IsSelected(TObject *obj) {return IsSelected(obj,AliRDHFCuts::kAll);}
  Bool_t IsSelected(TList *list) {if(!list) return kTRUE; return kFALSE;}
  Bool_t IsEventSelected(AliVEvent *event) const;
  Bool_t AreDaughtersSelected(AliAODRecoDecayHF *rd) const;

  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel) = 0;

  Int_t PtBin(Double_t pt) const;
  void PrintAll()const;

  enum{kAll,kTracks,kCandidate};

 protected:

  void SetNPtBins(Int_t nptBins){fnPtBins=nptBins;}
  void SetNVars(Int_t nVars){fnVars=nVars;}

  // cuts on the event
  Int_t fMinVtxType; // 0: not cut; 1: SPDZ; 2: SPD3D; 3: Tracks
  Int_t fMinVtxContr;   // minimum vertex contributors
  Float_t fMaxVtxRedChi2; // maximum chi2/ndf
  Int_t fMinSPDMultiplicity; // SPD multiplicity
  ULong64_t fTriggerMask; // trigger mask
  // quality cuts on the daughter tracks
  AliESDtrackCuts *fTrackCuts; // tracks for daughter tracks (AOD converted to ESD on the flight!)
  // cuts on the candidate
  Int_t fnPtBins;  // number of pt bins for cuts
  Int_t fnPtBinLimits; // "number of limits", that is fnPtBins+1
  Float_t* fPtBinLimits; //[fnPtBinLimits]  pt bins
  Int_t fnVars;    // number of cut vars for candidates
  TString *fVarNames; //[fnVars] names of the variables
  Int_t fnVarsForOpt;    // number of cut vars to be optimized for candidates
  Bool_t *fVarsForOpt; //[fnVars] kTRUE for vars to be used in optimization
  Int_t fGlobalIndex; // fnVars*fnPtBins
  Float_t *fCutsRD; //[fGlobalIndex] the cuts values
  Bool_t  *fIsUpperCut; //[fnVars] use > or < to select

  ClassDef(AliRDHFCuts,1);  // base class for cuts on AOD reconstructed 
                            // heavy-flavour decays
};

#endif

