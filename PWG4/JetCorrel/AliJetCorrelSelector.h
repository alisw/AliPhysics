#ifndef ALIJETCORRELSELECTOR_H
#define ALIJETCORRELSELECTOR_H
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
  UInt_t PoolDepth() const {return fPoolDepth;}
  UInt_t NoOfCorrel() const {return fNumCorrel;}
  UInt_t* CorrelTypes() const {return fCorrelType;}
  UInt_t NoOfBins(cBinType_t cType) const {return fNumBins[cType]-1;}  
  Float_t BinBorder(cBinType_t cType, UInt_t ii) const;
  Float_t MinLowBin(cBinType_t cType) const {return fBinning[cType][0];}
  Float_t MaxHighBin(cBinType_t cType) const {return fBinning[cType][fNumBins[cType]-1];}
  Int_t GetBin(cBinType_t cType, Float_t val) const;
  Bool_t GenQA() const       {return fGenQA;}
  Bool_t UseAliKF() const    {return fUseAliKF;}
  UInt_t DPhiNumBins() const {return fDPhiNumBins;}
  UInt_t DEtaNumBins() const {return fDEtaNumBins;}
  Float_t PoutBW() const     {return fPoutBW;}
  void Show() const;
  // Selection Setters:
  void SetPoolDepth(UInt_t v)      {fPoolDepth=v;}
  void SetCorrelTypes(UInt_t s,  UInt_t * const v);
  void SetBinningCentr(UInt_t s, Float_t * const v);
  void SetBinningZvert(UInt_t s, Float_t * const v);
  void SetBinningTrigg(UInt_t s, Float_t * const v);
  void SetBinningAssoc(UInt_t s, Float_t * const v);
  void SetTriggers(UInt_t s, TString * const v);
  void SetITSRefit(Bool_t v)    {fITSRefit=v;}
  void SetTPCRefit(Bool_t v)    {fTPCRefit=v;}
  void SetTRDRefit(Bool_t v)    {fTRDRefit=v;}
  void SetMaxEta(Float_t v)     {fMaxEta=v;}
  void SetPoutBinWidth(Float_t v) {fPoutBW=v;}
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
  void SetTrkProximityCut(Float_t v) {fTrkMinProx=v;}
  void SetUseAliKF(Bool_t v)    {fUseAliKF=v;}
  // Cutting methods:
  Bool_t SelectedEvtTrigger(AliESDEvent * const jcESD) const;
  Bool_t CloseTrackPair(Float_t dist) const;
  Bool_t LowQualityTrack(AliESDtrack* t) const;
  Bool_t PassPID(AliESDtrack* t, cPartType_t pType) const;
  Float_t GetSigmaToVertex(AliESDtrack* trk) const;
  void GetPID(AliESDtrack* trk, Stat_t& fpid, Stat_t& fweight) const;

  enum {kMaxCorrel = 1};   // Maximum no of correlations
  enum {kMaxVert = 3};    // Maximum no of vertex bins
  enum {kMaxCent = 3};     // Maximum no of centrality bins
  enum {kMaxTrig = 13};    // Maximum no of trigger bins
  enum {kMaxAsso = 7};    // Maximum no of associated bins
  
 private: 
  // Generic Selections:
  Bool_t fGenQA;                     //! generate QA histos
  UInt_t fDPhiNumBins, fDEtaNumBins; //! number of bins in DeltaPhi, DeltaEta histos
  UInt_t fNumCorrel, fNumEvtTriggs, fPoolDepth; //! number of correlations, event triggers, pool depth
  UInt_t* fCorrelType;           //! array of correlation types
  TString* fEvtTriggs;           //! array of event triggers
  UInt_t fNumBins[4];            //! number of bins: centr, zvert, trigg, assoc
  Float_t* fBinning[4];          //! bin margins: centr, zvert, trigg, assoc
  // Track Selections:
  Bool_t fITSRefit, fTPCRefit, fTRDRefit, fRejectKinkChild; //! on/off cuts
  Float_t fMaxEta;                   //! single-particle eta cut
  Float_t fPoutBW;                   //! Pout bin width
  Float_t fMaxNsigmaVtx;             //! track-primary vertex cut (sigma)
  Float_t fMaxTrkVtx;                //! track-primary vertex cut (value)
  Float_t fMaxITSChi2, fMaxTPCChi2;  //! ITS/TPC Chi2/cluster cut
  UInt_t fMinNClusITS, fMinNClusTPC; //! ITS/TPC number of clusters cut
  Float_t fTrkMinProx;               //! two-track proximity cut (dist at TPC entrance)
  Bool_t fUseAliKF;                  //! use AliKF or TLorentzVector for parent reconstruction
  
  // disable (make private) copy constructor, and assignment operator:
  AliJetCorrelSelector(const AliJetCorrelSelector&);
  AliJetCorrelSelector& operator=(const AliJetCorrelSelector&);
  
  ClassDef(AliJetCorrelSelector, 1); //AliJetCorrelSelector
};

#endif 
