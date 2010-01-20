#ifndef __ALIJETCORRELSELECTOR_H__
#define __ALIJETCORRELSELECTOR_H__
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//__________________________________________
// Class for user selections. Object is created by ConfigJetCorrel.C macro
// and passed to AddTaskJetCorrel.C running macro.
//-- Author: Paul Constantin

#include "CorrelTrack.h"
#include "CorrelKFTrack.h"
#include "CorrelRecoParent.h"
#include "CorrelList.h"

namespace JetCorrelHD {
  
  class AliJetCorrelSelector : public TObject {
    
  public:
    AliJetCorrelSelector();
    ~AliJetCorrelSelector();

    // Selection getters:
    UInt_t PoolDepth() {return fPoolDepth;}
    UInt_t NoOfCorrel() {return fNumCorrel;}
    UInt_t* CorrelTypes() {return fCorrelType;}
    const UInt_t NoOfBins(BinType_t cType) {return fNumBins[cType]-1;}  
    Float_t BinBorder(BinType_t cType, UInt_t ii);
    Int_t GetBin(BinType_t cType, Float_t val);
    Float_t MinAssocPt() {return minAssocPt;}
    Float_t MaxAssocPt() {return maxAssocPt;}
    UInt_t  NumAssocPt() {return UInt_t((maxAssocPt-minAssocPt)/bwAssocPt);}
    Float_t MinTriggPt() {return minTriggPt;}
    Float_t MaxTriggPt() {return maxTriggPt;}
    UInt_t  NumTriggPt() {return UInt_t((maxTriggPt-minTriggPt)/bwTriggPt);}
    Bool_t GenQA() {return fGenQA;};
    void Show();
    // Selection Setters:
    void SetPoolDepth(UInt_t v)      {fPoolDepth=v;}
    void SetCorrelTypes(UInt_t s,  UInt_t * const v);
    void SetBinningCentr(UInt_t s, Float_t * const v);
    void SetBinningZvert(UInt_t s, Float_t * const v);
    void SetBinningTrigg(Float_t min, Float_t max, Float_t bw);
    void SetBinningAssoc(Float_t min, Float_t max, Float_t bw);
    void SetITSRefit(Bool_t v)    {fITSRefit=v;}
    void SetTPCRefit(Bool_t v)    {fTPCRefit=v;}
    void SetTRDRefit(Bool_t v)    {fTRDRefit=v;}
    void SetMaxITSChi2(Float_t v) {fMaxITSChi2=v;}
    void SetMaxTPCChi2(Float_t v) {fMaxTPCChi2=v;}
    void SetMinNClusITS(UInt_t v) {fMinNClusITS=v;}
    void SetMinNClusTPC(UInt_t v) {fMinNClusTPC=v;}
    void SetMaxNsigmaVtx(Float_t v) {fMaxNsigmaVtx=v;}
    void SetRejectKinkChild(Bool_t v) {fRejectKinkChild=v;}
    void SetQA(Bool_t v) {fGenQA=v;}
    // Cutting methods:
    Bool_t IsAssoc(Float_t pT) {return (pT>=minAssocPt && pT<=maxAssocPt);}
    Bool_t IsTrigg(Float_t pT) {return (pT>=minTriggPt && pT<=maxTriggPt);}
    Bool_t GoodTrackPair(CorrelTrack_t* t1, CorrelTrack_t* t2);
    Bool_t LowQualityTrack(AliESDtrack* t);  
    Bool_t PassPID(AliESDtrack* t, PartType_t pType);
    Float_t GetSigmaToVertex(AliESDtrack* trk);
    void GetPID(AliESDtrack* trk, Stat_t& fpid, Stat_t& fweight);

  private: 
    // Generic Selections:
    UInt_t fNumCorrel, fPoolDepth; // number of correlations, pool depth
    UInt_t *fCorrelType;           // array of correlation types
    UInt_t fNumBins[2];            // number of bins: centr, zvert
    Float_t* fBinning[2];          // bin margins: centr, zvert
    Bool_t fGenQA;                 // generate QA histos
    Float_t minTriggPt, maxTriggPt, bwTriggPt; // trigg Pt binning
    Float_t minAssocPt, maxAssocPt, bwAssocPt; // assoc Pt binning
    // Track Selections:
    Bool_t fITSRefit, fTPCRefit, fTRDRefit, fRejectKinkChild; // on/off cuts
    Float_t fMaxNsigmaVtx;
    Float_t fMaxITSChi2, fMaxTPCChi2;
    UInt_t fMinNClusITS, fMinNClusTPC;
    
    // disable (make private) copy constructor, and assignment operator:
    AliJetCorrelSelector(const AliJetCorrelSelector&);
    AliJetCorrelSelector& operator=(const AliJetCorrelSelector&);

    ClassDef(AliJetCorrelSelector, 1);
  };

} // namespace

#endif 
