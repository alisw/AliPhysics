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
  
  void  SetUseHybridESA(Bool_t usehyb) {fUseHybrid = usehyb;}
  void  SetTrackFilterESAHyb1(AliAnalysisFilter* trackH1F) {fTrackFilterHybrid1 = trackH1F;}
  void  SetTrackFilterESAHyb2(AliAnalysisFilter* trackH2F) {fTrackFilterHybrid2 = trackH2F;}
  void  SetTrackFilterESA(AliAnalysisFilter* trackF) {fTrackFilterESA = trackF;}
  void  SetMinMultForESA(Int_t minnch)     {fMinMultESA = minnch;}
  void  SetStepSizeESA(Float_t sizestep)   {fSizeStepESA = sizestep;}
  void  SetIsEtaAbsESA(Bool_t isabseta)    {fIsAbsEtaESA = isabseta;}
  void  SetTrackEtaMinESA(Float_t etaminF) {fEtaMinCutESA = etaminF;}
  void  SetTrackEtaMaxESA(Float_t etamaxF) {fEtaMaxCutESA = etamaxF;}
  void  SetTrackPtMinESA(Float_t ptminF)   {fPtMinCutESA = ptminF;}
  void  SetTrackPtMaxESA(Float_t ptmaxF)   {fPtMaxCutESA = ptmaxF;}
  
  void SaveHistos(const char* folder = 0);
  
  TH1D * GetHistData( Int_t bin_histo );
  
  //Utility functions
  //for the base virtual event class: all methods are common
  Float_t GetEventShape(AliVEvent *event, TString lMethod = "V0M", Bool_t fillHist = kTRUE);
  Float_t GetEventShapeTrue(AliStack *event, TString lMethod = "V0M", Bool_t fillHist = kTRUE);
  Float_t GetSpherocity(AliVEvent *event, Bool_t fillHist);
  Float_t GetSphericity(AliVEvent *event, Bool_t fillHist);
  Float_t GetSphericityMC(AliStack *eventMC, Bool_t fillHist);
  Float_t GetSpherocityMC(AliStack *eventMC, Bool_t fillHist);
  
  
  //EvSel Snippets
  Float_t MinVal( Float_t A, Float_t B ); 
  
 private:
  Bool_t fUseHybrid;
  AliAnalysisFilter *fTrackFilterHybrid1;
  AliAnalysisFilter *fTrackFilterHybrid2;
  AliAnalysisFilter *fTrackFilterESA; // Track filter for Event Shapes
  Int_t   fMinMultESA;
  Float_t fSizeStepESA;
  Bool_t  fIsAbsEtaESA;
  Float_t fEtaMaxCutESA;
  Float_t fEtaMinCutESA;
  Float_t fPtMaxCutESA;
  Float_t fPtMinCutESA;
  Int_t   fRunNumber; // for control of run changes
  TH1D    *fheta;
  TH1D    *fhphi;
  TH1D    *fhpt;
  TH1D    *fhetaMC;
  TH1D    *fhphiMC;
  TH1D    *fhptMC;
  TH1D    *fAverageAmplitudes; 
  
  
  ClassDef(AliTransverseEventShape,2) // base helper class
    };
#endif

