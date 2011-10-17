#ifndef ALIRDHFCUTS_H
#define ALIRDHFCUTS_H
/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//***********************************************************
// Class AliRDHFCuts
// base class for cuts on AOD reconstructed heavy-flavour decays
// Author: A.Dainese, andrea.dainese@pd.infn.it
//***********************************************************

#include <TString.h>

#include "AliAnalysisCuts.h"
#include "AliESDtrackCuts.h"
#include "AliAODPidHF.h"
#include "AliAODEvent.h"
#include "AliVEvent.h"

class AliAODTrack;
class AliAODRecoDecayHF;
class AliESDVertex;

class AliRDHFCuts : public AliAnalysisCuts 
{
 public:

  enum ECentrality {kCentOff,kCentV0M,kCentTRK,kCentTKL,kCentCL1,kCentInvalid};
  enum ESelLevel {kAll,kTracks,kPID,kCandidate};
  enum EPileup {kNoPileupSelection,kRejectPileupEvent,kRejectTracksFromPileupVertex};
  enum ESele {kD0toKpiCuts,kD0toKpiPID,kD0fromDstarCuts,kD0fromDstarPID,kDplusCuts,kDplusPID,kDsCuts,kDsPID,kLcCuts,kLcPID,kDstarCuts,kDstarPID};
  enum ERejBits {kNotSelTrigger,kNoVertex,kTooFewVtxContrib,kZVtxOutFid,kPileupSPD,kOutsideCentrality,kPhysicsSelection};
  AliRDHFCuts(const Char_t* name="RDHFCuts", const Char_t* title="");
  
  virtual ~AliRDHFCuts();
  
  AliRDHFCuts(const AliRDHFCuts& source);
  AliRDHFCuts& operator=(const AliRDHFCuts& source); 

  virtual void SetStandardCutsPP2010() {return;}  
  virtual void SetStandardCutsPbPb2010() {return;}  


  void SetMinCentrality(Float_t minCentrality=0.) {fMinCentrality=minCentrality;} 
  void SetMaxCentrality(Float_t maxCentrality=100.) {fMaxCentrality=maxCentrality;} 
  void SetMinVtxType(Int_t type=3) {fMinVtxType=type;}  
  void SetMinVtxContr(Int_t contr=1) {fMinVtxContr=contr;}  
  void SetMaxVtxRdChi2(Float_t chi2=1e6) {fMaxVtxRedChi2=chi2;}  
  void SetMaxVtxZ(Float_t z=1e6) {fMaxVtxZ=z;}  
  void SetMinSPDMultiplicity(Int_t mult=0) {fMinSPDMultiplicity=mult;}  
  void SetTriggerMask(ULong64_t mask=0) {fTriggerMask=mask;} 
  void SetTriggerClass(TString trclass) {fTriggerClass=trclass;} 
  void SetVarsForOpt(Int_t nVars,Bool_t *forOpt);
  void SetGlobalIndex(){fGlobalIndex=fnVars*fnPtBins;}
  void SetGlobalIndex(Int_t nVars,Int_t nptBins){fnVars=nVars; fnPtBins=nptBins; SetGlobalIndex();}
  void SetVarNames(Int_t nVars,TString *varNames,Bool_t *isUpperCut);  
  void SetPtBins(Int_t nPtBinLimits,Float_t *ptBinLimits);
  void SetCuts(Int_t nVars,Int_t nPtBins,Float_t** cutsRD);
  void SetCuts(Int_t glIndex, Float_t* cutsRDGlob);
  void AddTrackCuts(const AliESDtrackCuts *cuts) 
          {delete fTrackCuts; fTrackCuts=new AliESDtrackCuts(*cuts); return;}
  void SetUsePID(Bool_t flag=kTRUE) {fUsePID=flag; return;}
  void SetUseAOD049(Bool_t flag=kTRUE) {fUseAOD049=flag; return;}
  void SetUseCentrality(Int_t flag=1);    // see enum below
  void SetPidHF(AliAODPidHF* pidObj) {
    if(fPidHF) delete fPidHF;
    fPidHF=new AliAODPidHF(*pidObj);
  }
  void SetRemoveDaughtersFromPrim(Bool_t removeDaughtersPrim) {fRemoveDaughtersFromPrimary=removeDaughtersPrim;}
  void SetMinPtCandidate(Double_t ptCand=-1.) {fMinPtCand=ptCand; return;}
  void SetMaxPtCandidate(Double_t ptCand=1000.) {fMaxPtCand=ptCand; return;}
  void SetOptPileup(Int_t opt=0){
    // see enum below
    fOptPileup=opt;
  }
  void ConfigurePileupCuts(Int_t minContrib=3, Float_t minDz=0.6){
    fMinContrPileup=minContrib;
    fMinDzPileup=minDz;
  }


  AliAODPidHF* GetPidHF() const {return fPidHF;}
  Float_t *GetPtBinLimits() const {return fPtBinLimits;}
  Int_t   GetNPtBins() const {return fnPtBins;}
  Int_t   GetNVars() const {return fnVars;} 
  TString *GetVarNames() const {return fVarNames;} 
  Bool_t  *GetVarsForOpt() const {return fVarsForOpt;} 
  Int_t   GetNVarsForOpt() const {return fnVarsForOpt;}
  const Float_t *GetCuts() const {return fCutsRD;} 
  void    GetCuts(Float_t**& cutsRD) const;
  Float_t GetCutValue(Int_t iVar,Int_t iPtBin) const;
  Double_t GetMaxVtxZ() const {return fMaxVtxZ;}  
  Float_t GetCentrality(AliAODEvent* aodEvent){return GetCentrality(aodEvent,(AliRDHFCuts::ECentrality)fUseCentrality);}
  Float_t GetCentrality(AliAODEvent* aodEvent, AliRDHFCuts::ECentrality estimator);
  Bool_t  *GetIsUpperCut() const {return fIsUpperCut;}
  AliESDtrackCuts *GetTrackCuts() const {return fTrackCuts;}
  virtual AliESDtrackCuts *GetTrackCutsSoftPi() const {return 0;}
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters) = 0;
  virtual void GetCutVarsForOpt(AliAODRecoDecayHF *d,Float_t *vars,Int_t nvars,Int_t *pdgdaughters,AliAODEvent * /*aod*/)
            {return GetCutVarsForOpt(d,vars,nvars,pdgdaughters);}
  Int_t   GetGlobalIndex(Int_t iVar,Int_t iPtBin) const;
  void    GetVarPtIndex(Int_t iGlob, Int_t& iVar, Int_t& iPtBin) const;
  Bool_t  GetIsUsePID() const {return fUsePID;}
  Bool_t  GetUseAOD049() const {return fUseAOD049;}
  Bool_t  GetIsPrimaryWithoutDaughters() const {return fRemoveDaughtersFromPrimary;}
  Bool_t GetOptPileUp() const {return fOptPileup;}
  Int_t GetUseCentrality() const {return fUseCentrality;}
  Float_t GetMinCentrality() const {return fMinCentrality;}
  Float_t GetMaxCentrality() const {return fMaxCentrality;}
  Double_t GetMinPtCandidate() const {return fMinPtCand;}
  Double_t GetMaxPtCandidate() const {return fMaxPtCand;}
  Bool_t IsSelected(TObject *obj) {return IsSelected(obj,AliRDHFCuts::kAll);}
  Bool_t IsSelected(TList *list) {if(!list) return kTRUE; return kFALSE;}
  Int_t  IsEventSelectedInCentrality(AliVEvent *event);
  Bool_t IsEventSelected(AliVEvent *event);
  Bool_t AreDaughtersSelected(AliAODRecoDecayHF *rd) const;
  Bool_t IsDaughterSelected(AliAODTrack *track,const AliESDVertex *primary,AliESDtrackCuts *cuts) const;
  virtual Int_t IsSelectedPID(AliAODRecoDecayHF * /*rd*/) {return 1;}

  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel) = 0;
  virtual Int_t IsSelected(TObject* obj,Int_t selectionLevel,AliAODEvent* /*aod*/)
                {return IsSelected(obj,selectionLevel);}
  Int_t PtBin(Double_t pt) const;
  void PrintAll()const;

  virtual Bool_t IsInFiducialAcceptance(Double_t /*pt*/,Double_t /*y*/) const {return kTRUE;}

  void SetWhyRejection(Int_t why) {fWhyRejection=why; return;}
  Int_t GetWhyRejection() const {return fWhyRejection;}
  UInt_t GetEventRejectionBitMap() const {return fEvRejectionBits;}
  Bool_t IsEventRejectedDueToTrigger() const {
    return fEvRejectionBits&(1<<kNotSelTrigger);
  }
  Bool_t IsEventRejectedDueToNotRecoVertex() const {
    return fEvRejectionBits&(1<<kNoVertex);
  }
  Bool_t IsEventRejectedDueToVertexContributors() const {
    return fEvRejectionBits&(1<<kTooFewVtxContrib);
  }
  Bool_t IsEventRejectedDueToZVertexOutsideFiducialRegion() const {
    return fEvRejectionBits&(1<<kZVtxOutFid);
  }
  Bool_t IsEventRejectedDueToPileupSPD() const {
    return fEvRejectionBits&(1<<kPileupSPD);
  }
  Bool_t IsEventRejectedDueToCentrality() const {
    return fEvRejectionBits&(1<<kOutsideCentrality);
  }
  Bool_t IsEventRejectedDuePhysicsSelection() const {
    return fEvRejectionBits&(1<<kPhysicsSelection);
  }

  void SetFixRefs(Bool_t fix=kTRUE) {fFixRefs=fix; return;}

  Bool_t CompareCuts(const AliRDHFCuts *obj) const;
  void MakeTable()const;

  Int_t GetIsSelectedCuts() const {return fIsSelectedCuts;}
  Int_t GetIsSelectedPID() const  {return fIsSelectedPID;}

  void SetUseMCVertex() { fUseMCVertex=kTRUE; }
  Bool_t GetUseMCVertex() const { return fUseMCVertex; }

  Bool_t RecalcOwnPrimaryVtx(AliAODRecoDecayHF *d,AliAODEvent *aod) const;
  Bool_t SetMCPrimaryVtx(AliAODRecoDecayHF *d,AliAODEvent *aod) const;
  void   CleanOwnPrimaryVtx(AliAODRecoDecayHF *d,AliAODEvent *aod,AliAODVertex *origownvtx) const;

  Bool_t CountEventForNormalization() const 
  { if(fWhyRejection==0) {return kTRUE;} else {return kFALSE;} }

  void SetKeepSignalMC() {fKeepSignalMC=kTRUE; return;}


 protected:

  void SetNPtBins(Int_t nptBins){fnPtBins=nptBins;}
  void SetNVars(Int_t nVars){fnVars=nVars;}

  Bool_t IsSignalMC(AliAODRecoDecay *d,AliAODEvent *aod,Int_t pdg) const;

  // cuts on the event
  Int_t fMinVtxType; // 0: not cut; 1: SPDZ; 2: SPD3D; 3: Tracks
  Int_t fMinVtxContr;   // minimum vertex contributors
  Float_t fMaxVtxRedChi2; // maximum chi2/ndf
  Float_t fMaxVtxZ; // maximum |z| of primary vertex
  Int_t fMinSPDMultiplicity; // SPD multiplicity
  ULong64_t fTriggerMask; // trigger mask
  TString  fTriggerClass; // trigger class
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
  Bool_t fUsePID; // enable PID usage (off by default)
  Bool_t fUseAOD049; // enable AOD049 centrality cleanup
  AliAODPidHF *fPidHF; // PID for heavy flavours manager
  Int_t fWhyRejection; // used to code the step at which candidate was rejected
  UInt_t fEvRejectionBits; //bit map storing the full info about event rejection
  Bool_t fRemoveDaughtersFromPrimary; // flag to switch on the removal of duaghters from the primary vertex computation
  Bool_t fUseMCVertex; // use MC primary vertex 
  Bool_t fUsePhysicsSelection; // use Physics selection criteria
  Int_t  fOptPileup;      // option for pielup selection
  Int_t  fMinContrPileup; // min. n. of tracklets in pileup vertex
  Float_t fMinDzPileup;   // min deltaz between main and pileup vertices
  Int_t   fUseCentrality; // off =0 (default)
                          // 1 = V0 
                          // 2 = Tracks
                          // 3 = Tracklets
                          // 4 = SPD clusters outer 
  Float_t fMinCentrality; // minimum centrality for selected events
  Float_t fMaxCentrality; // maximum centrality for selected events
  Bool_t  fFixRefs;       // fix the daughter track references 
  Int_t  fIsSelectedCuts; // outcome of cuts selection
  Int_t  fIsSelectedPID;  // outcome of PID selection
  Double_t fMinPtCand; // minimum pt of the candidate
  Double_t fMaxPtCand; // minimum pt of the candidate
  Bool_t  fKeepSignalMC; // IsSelected returns always kTRUE for MC signal

  ClassDef(AliRDHFCuts,19);  // base class for cuts on AOD reconstructed heavy-flavour decays
};

#endif
