#ifndef  ALIJETCORRELWRITER_H
#define  ALIJETCORRELWRITER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
/* $Id:  $ */

//______________________________________________________
// Class for output (histograms) definition and filling.
// Contains also the methods for calculation of correlation parameters.
//-- Author: Paul Constantin

#include "AliJetCorrelMaker.h"

class AliJetCorrelWriter : public TObject {
  
 public:
  AliJetCorrelWriter();
  ~AliJetCorrelWriter();    
  void Init(AliJetCorrelSelector * const s, AliJetCorrelMaker * const m);
  
  void CreateGeneric(TList* histosContainer);
  void CreateQA(TList *histosContainer);
  void CreateCorrelations(TList* histosContainer);
  
  void FillGlobal(Float_t cent, Float_t vert);
  void FillSingleHistos(UInt_t cBin, CorrelList_t * const TriggList, UInt_t tIdx,
			CorrelList_t * const AssocList, UInt_t aIdx);
  void FillTrackQA(AliESDtrack * const track, UInt_t idx);
  void FillParentNtuple(CorrelList_t * const ParentList);
  void FillCorrelations(UInt_t fTyp, UInt_t iCorr, UInt_t cBin, UInt_t vBin,
			CorrelParticle_t * const Trigger, CorrelParticle_t * const Associated);
  
  Float_t DeltaPhi(Float_t phi1, Float_t phi2);
  
 private:
  AliJetCorrelSelector *fSelector;  //! user selection object
  AliJetCorrelMaker *fMaker;        //! correlation maker object
  TString fHname, fHtit;            //! histos name&title
  Bool_t fRecoTrigg;                //! is trigger reconstructed
  TRandom2 fRndm;                   //! random number generator
  
  // Output Histograms
  TNtuple *fNtuParent;                                       // reconstructed parent ntuple
  TH1F *fHBinsCentr, *fHBinsZVert, *fHBinsTrigg, *fHBinsAssoc, *fHCentr, *fHZVert;  // binning histos
  TH2F *fHTrkITSQA[2], *fHTrkTPCQA[2], *fHTrkVTXQA[2];       // track QA histos
  TH3F *fHTrkProx[2][AliJetCorrelSelector::kMaxCent];        // distance at TPC entrance between tracks
  TH3F *fHTriggAcc[AliJetCorrelSelector::kMaxCorrel];        // trigg PhixEta acceptance
  TH3F *fHAssocAcc[AliJetCorrelSelector::kMaxCorrel];        // assoc PhixEta acceptance
  TH1F *fHTriggPt[AliJetCorrelSelector::kMaxCorrel][AliJetCorrelSelector::kMaxCent]; // trigg Pt
  TH1F *fHAssocPt[AliJetCorrelSelector::kMaxCorrel][AliJetCorrelSelector::kMaxCent]; // assoc Pt
  TH3F *fHReal[AliJetCorrelSelector::kMaxCorrel][AliJetCorrelSelector::kMaxCent][AliJetCorrelSelector::kMaxVert][AliJetCorrelSelector::kMaxTrig][AliJetCorrelSelector::kMaxAsso];  // real DPhixDEtaxPout
  TH3F *fHMix[AliJetCorrelSelector::kMaxCorrel][AliJetCorrelSelector::kMaxCent][AliJetCorrelSelector::kMaxVert][AliJetCorrelSelector::kMaxTrig][AliJetCorrelSelector::kMaxAsso];   // mix  DPhixDEtaxPout
  
  // disable (make private) copy constructor and assignment operator:
  AliJetCorrelWriter(const AliJetCorrelWriter&);
  AliJetCorrelWriter& operator=(const AliJetCorrelWriter&);
  
  ClassDef(AliJetCorrelWriter, 1);
};

#endif
