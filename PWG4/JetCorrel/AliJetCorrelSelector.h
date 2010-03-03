#ifndef __ALIJETCORRELSELECTOR_H__
#define __ALIJETCORRELSELECTOR_H__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//__________________________________________
// Class for user selections.
// Object is created by ConfigJetCorrel.C macro and passed to 
// AddTaskJetCorrel.C running macro via the main class AliAnalysisTaskJetCorrel
//-- Author: Paul Constantin

#include "CorrelTrack.h"
#include "CorrelKFTrack.h"
#include "CorrelRecoParent.h"
#include "CorrelList.h"

class AliJetCorrelSelector : public TObject {
  
 public:
  AliJetCorrelSelector();
  ~AliJetCorrelSelector();
  
  // Selection getters:
  UInt_t PoolDepth() {return fPoolDepth;}
  UInt_t NoOfCorrel() {return fNumCorrel;}
  UInt_t* CorrelTypes() {return fCorrelType;}
  UInt_t NoOfBins(BinType_t cType) {return fNumBins[cType]-1;}  
  Float_t BinBorder(BinType_t cType, UInt_t ii);
  Int_t GetBin(BinType_t cType, Float_t val);
  Float_t MinAssocPt() {return minAssocPt;}
  Float_t MaxAssocPt() {return maxAssocPt;}
  UInt_t  NumAssocPt();
  Float_t MinTriggPt() {return minTriggPt;}
  Float_t MaxTriggPt() {return maxTriggPt;}
  UInt_t  NumTriggPt();
  Bool_t GenQA()       {return fGenQA;}
  Bool_t UseAliKF()    {return fUseAliKF;}
  UInt_t DPhiNumBins() {return fDPhiNumBins;}
  UInt_t DEtaNumBins() {return fDEtaNumBins;}
  void Show();
  // Selection Setters:
  void SetPoolDepth(UInt_t v)      {fPoolDepth=v;}
  void SetCorrelTypes(UInt_t s,  UInt_t * const v);
  void SetBinningCentr(UInt_t s, Float_t * const v);
  void SetBinningZvert(UInt_t s, Float_t * const v);
  void SetBinningTrigg(Float_t min, Float_t max, Float_t bw);
  void SetBinningAssoc(Float_t min, Float_t max, Float_t bw);
  void SetTriggers(UInt_t s, TString * const v);
  void SetITSRefit(Bool_t v)    {fITSRefit=v;}
  void SetTPCRefit(Bool_t v)    {fTPCRefit=v;}
  void SetTRDRefit(Bool_t v)    {fTRDRefit=v;}
  void SetMaxEta(Float_t v)     {fMaxEta=v;}
  void SetMaxITSChi2(Float_t v) {fMaxITSChi2=v;}
  void SetMaxTPCChi2(Float_t v) {fMaxTPCChi2=v;}
  void SetMinNClusITS(UInt_t v) {fMinNClusITS=v;}
  void SetMinNClusTPC(UInt_t v) {fMinNClusTPC=v;}
  void SetMaxNsigmaVtx(Float_t v) {fMaxNsigmaVtx=v;}
  void SetMaxTrkVtx(Float_t v)  {fMaxTrkVtx=v;}
  void SetRejectKinkChild(Bool_t v) {fRejectKinkChild=v;}
  void SetQA(Bool_t v) {fGenQA=v;}
  void SetDPhiNumBins(UInt_t v) {fDPhiNumBins=v;}
  void SetDEtaNumBins(UInt_t v) {fDEtaNumBins=v;}
  void SetTrkProximityCut(Float_t v) {trkMinProx=v;}
  void SetUseAliKF(Bool_t v)    {fUseAliKF=v;}
  // Cutting methods:
  Int_t AssocBin(Float_t pT);
  Int_t TriggBin(Float_t pT);
  Bool_t SelectedEvtTrigger(AliESDEvent * const jcESD);
  Bool_t CloseTrackPair(Float_t dist);
  Bool_t LowQualityTrack(AliESDtrack* t);  
  Bool_t PassPID(AliESDtrack* t, PartType_t pType);
  Float_t GetSigmaToVertex(AliESDtrack* trk);
  void GetPID(AliESDtrack* trk, Stat_t& fpid, Stat_t& fweight);
  
 private: 
  // Generic Selections:
  Bool_t fGenQA;                 // generate QA histos
  UInt_t fDPhiNumBins, fDEtaNumBins; // number of bins in DeltaPhi, DeltaEta histos
  UInt_t fNumCorrel, nEvtTriggs, fPoolDepth; // number of correlations, event triggers, pool depth
  UInt_t *fCorrelType;           // array of correlation types
  TString *fEvtTriggs;           // array of event triggers
  UInt_t fNumBins[2];            // number of bins: centr, zvert
  Float_t* fBinning[2];          // bin margins: centr, zvert
  Float_t minTriggPt, maxTriggPt, bwTriggPt; // trigg Pt binning
  Float_t minAssocPt, maxAssocPt, bwAssocPt; // assoc Pt binning
  // Track Selections:
  Bool_t fITSRefit, fTPCRefit, fTRDRefit, fRejectKinkChild; // on/off cuts
  Float_t fMaxEta;                   // single-particle eta cut
  Float_t fMaxNsigmaVtx;             // track-primary vertex cut (sigma)
  Float_t fMaxTrkVtx;                // track-primary vertex cut (value)
  Float_t fMaxITSChi2, fMaxTPCChi2;  // ITS/TPC Chi2/cluster cut
  UInt_t fMinNClusITS, fMinNClusTPC; // ITS/TPC number of clusters cut
  Float_t trkMinProx;                // two-track proximity cut (dist at TPC entrance)
  Bool_t fUseAliKF;                  // use AliKF or TLorentzVector for parent reconstruction
  
  // disable (make private) copy constructor, and assignment operator:
  AliJetCorrelSelector(const AliJetCorrelSelector&);
  AliJetCorrelSelector& operator=(const AliJetCorrelSelector&);
  
  ClassDef(AliJetCorrelSelector, 1);
};

#endif 
