#ifndef AliTransverseEventShape_H
#define AliTransverseEventShape_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//   Class AliTransverseEventShape
//   author: Antonio Ortiz Velasquez ( antonio.ortiz@nucleares,unam.mx )
//           Eleazar Cuautle Flores  ( ecuautle@nucleares.unam.mx  )
//*****************************************************

#include "TObject.h"

#include <AliAnalysisFilter.h>
#include <vector>

class AliVEvent;
class AliVVertex;
class AliESDEvent;
class AliAODEvent;
class AliStack;

class AliTransverseEventShape : public TObject {
  
 public:  
  AliTransverseEventShape();
  virtual ~AliTransverseEventShape() {};
  
  //Extra const
  AliTransverseEventShape(const AliTransverseEventShape& pd);
  AliTransverseEventShape &operator=(const AliTransverseEventShape &c);
  virtual void Init();

  void  SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}  
  void  SetUseHybridESA(Bool_t usehyb) {fUseHybrid = usehyb;}
//Filters for ESDs
  void  SetTrackFilterESAHyb1(AliAnalysisFilter* trackH1F) {fTrackFilterHybrid1 = trackH1F;}
  void  SetTrackFilterESAHyb2(AliAnalysisFilter* trackH2F) {fTrackFilterHybrid2 = trackH2F;}
  void  SetTrackFilterESA(AliAnalysisFilter* trackF) {fTrackFilterGlobal = trackF;}
//Filters for AODs
  void  SetAODTrackFilterESAHyb1(Int_t aodtrackH1F) {fAODFilterHybrid1 = aodtrackH1F;}
  void  SetAODTrackFilterESAHyb2(Int_t aodtrackH2F) {fAODFilterHybrid2 = aodtrackH2F;}
  void  SetAODTrackFilterESA(Int_t aodtrackF) {fAODFilterGlobal = aodtrackF;}

  void  SetMinMultForESA(Int_t minnch)     {fMinMultESA = minnch;}
  void  SetStepSizeESA(Float_t sizestep)   {fSizeStepESA = sizestep;}
  void  SetIsEtaAbsESA(Bool_t isabseta)    {fIsAbsEtaESA = isabseta;}
  void  SetTrackEtaMinESA(Float_t etaminF) {fEtaMinCutESA = etaminF;}
  void  SetTrackEtaMaxESA(Float_t etamaxF) {fEtaMaxCutESA = etamaxF;}
  void  SetTrackPtMinESA(Float_t ptminF)   {fPtMinCutESA = ptminF;}
  void  SetTrackPtMaxESA(Float_t ptmaxF)   {fPtMaxCutESA = ptmaxF;}
  
  void SaveHistosSo(const char* folder = 0 );
  void SaveHistosSt(const char* folder = 0 );  

  TH1D * GetHistData( Int_t bin_histo, TString lMethod  );
  
  //Utility functions
  //for the base virtual event class: all methods are common
  Float_t GetEventShape(AliVEvent *event, TString lMethod = "V0M", Bool_t fillHist = kTRUE);
  Float_t GetEventShapeTrue(AliStack *event, TString lMethod = "V0M", Bool_t fillHist = kTRUE);
  Float_t GetSpherocity(Bool_t fillHist);
  Float_t GetSphericity(Bool_t fillHist);
  Float_t GetSphericityMC(Bool_t fillHist);
  Float_t GetSpherocityMC(Bool_t fillHist);

  Int_t   ReadAODEvent(std::vector<Float_t> &pt, std::vector<Float_t> &eta, std::vector<Float_t> &phi);
  Int_t   ReadESDEvent(std::vector<Float_t> &pt, std::vector<Float_t> &eta, std::vector<Float_t> &phi);
  Int_t   ReadMC(std::vector<Float_t> &pt, std::vector<Float_t> &eta, std::vector<Float_t> &phi);

  Float_t AnalyseGetSpherocity(Bool_t fillHist, const std::vector<Float_t> &pt,
		  const std::vector<Float_t> &eta,
		  const std::vector<Float_t> &phi);

  Float_t AnalyseGetSphericity(Bool_t fillHist, const std::vector<Float_t> &pt,
		  const std::vector<Float_t> &eta,
		  const std::vector<Float_t> &phi);


  //EvSel Snippets
  Float_t MinVal( Float_t A, Float_t B ); 

 private:
  Bool_t fAnalysisMC;          //  Real(kFALSE) or MC(kTRUE) flag
  AliESDEvent * fESDEvent;
  AliAODEvent * fAODEvent;
  AliStack * fMCStack;
  Int_t fNrec;
  Bool_t fUseHybrid;
  AliAnalysisFilter *fTrackFilterHybrid1;
  AliAnalysisFilter *fTrackFilterHybrid2;
  AliAnalysisFilter *fTrackFilterGlobal; // Track filter for Event Shapes

  Int_t fAODFilterHybrid1;
  Int_t fAODFilterHybrid2;
  Int_t fAODFilterGlobal; // Track filter for Event Shapes

  Int_t   fMinMultESA;
  Float_t fSizeStepESA;
  Bool_t  fIsAbsEtaESA;
  Float_t fEtaMaxCutESA;
  Float_t fEtaMinCutESA;
  Float_t fPtMaxCutESA;
  Float_t fPtMinCutESA;
  Int_t   fRunNumber; // for control of run changes
  TH1D    *fhetaSo;
  TH1D    *fhphiSo;
  TH1D    *fhptSo;
  TH1D    *fhetaSt;
  TH1D    *fhphiSt;
  TH1D    *fhptSt;
  TH1D    *fhetaSoMC;
  TH1D    *fhphiSoMC;
  TH1D    *fhptSoMC;
  TH1D    *fhetaStMC;
  TH1D    *fhphiStMC;
  TH1D    *fhptStMC;


  ClassDef(AliTransverseEventShape,2) // base helper class
};
#endif

