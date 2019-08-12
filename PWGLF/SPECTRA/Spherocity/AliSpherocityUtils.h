#ifndef AliSpherocityUtils_H
#define AliSpherocityUtils_H
/* This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

//*****************************************************
//   Class AliSpherocityUtils
//   author: Antonio Ortiz Velasquez ( antonio.ortiz@nucleares,unam.mx )
//*****************************************************

#include "TObject.h"

#include <AliAnalysisFilter.h>
#include <vector>

class AliVEvent;
class AliVVertex;
class AliESDEvent;
class AliAODEvent;
class AliStack;

class AliSpherocityUtils : public TObject {
  
 public:  
  AliSpherocityUtils();
  virtual ~AliSpherocityUtils() {};
  
  //Extra const
  AliSpherocityUtils(const AliSpherocityUtils& pd);
  AliSpherocityUtils &operator=(const AliSpherocityUtils &c);
  virtual void Init();

  void  SetAnalysisMC(Bool_t isMC) {fAnalysisMC = isMC;}  
//Filters for ESDs
  void  SetTrackFilter(AliAnalysisFilter* trackF) {fTrackFilterGlobal = trackF;}
//Filters for AODs
  void  SetAODTrackFilter(Int_t aodtrackF) {fAODFilterGlobal = aodtrackF;}

  void  SetMinMult(Int_t minnch)        {fMinMult    = minnch;}
  void  SetStepSize(Float_t sizestep)   {fSizeStep   = sizestep;}
  void  SetIsEtaAbs(Bool_t isabseta)    {fIsAbsEta   = isabseta;}
  void  SetTrackEtaMin(Float_t etaminF) {fEtaMinCut  = etaminF;}
  void  SetTrackEtaMax(Float_t etamaxF) {fEtaMaxCut  = etamaxF;}
  void  SetTrackPtMin(Float_t ptminF)   {fPtMinCut   = ptminF;}
  void  SetTrackPtMax(Float_t ptmaxF)   {fPtMaxCut   = ptmaxF;}
  
  void SaveHistosSo(const char* folder = 0 );

  TH1D * GetHistData( Int_t bin_histo );
  
  //Utility functions
  //for the base virtual event class: all methods are common
  Float_t GetEventShape(AliVEvent *event, TH1D * hphi =0, TH1D * heta =0);
  Float_t GetEventShapeTrue(AliStack *event, TH1D * hphi =0, TH1D * heta =0);
  Float_t GetSpherocity(TH1D * hphi =0, TH1D * heta =0);
  Float_t GetSpherocityMC(TH1D * hphi =0, TH1D * heta =0);

  Int_t   ReadAODEvent(std::vector<Float_t> &pt, std::vector<Float_t> &eta, std::vector<Float_t> &phi, TH1D * hphi =0, TH1D * heta =0);
  Int_t   ReadESDEvent(std::vector<Float_t> &pt, std::vector<Float_t> &eta, std::vector<Float_t> &phi, TH1D * hphi =0, TH1D * heta =0);
  Int_t   ReadMC(std::vector<Float_t> &pt, std::vector<Float_t> &eta, std::vector<Float_t> &phi, TH1D * hphi =0, TH1D * heta =0);

  Float_t AnalyseGetSpherocity(const std::vector<Float_t> &pt,
		  const std::vector<Float_t> &eta,
		  const std::vector<Float_t> &phi);


  //EvSel Snippets
  Float_t MinVal( Float_t A, Float_t B ); 

 private:
  AliESDEvent * fESDEvent;
  AliAODEvent * fAODEvent;
  AliStack * fMCStack;
  Int_t fNrec;
  Bool_t fAnalysisMC;          //  Real(kFALSE) or MC(kTRUE) flag
  AliAnalysisFilter *fTrackFilterGlobal; // Track filter for Event Shapes
  Int_t fAODFilterGlobal; // Track filter for Event Shapes

  Int_t   fMinMult;
  Float_t fSizeStep;
  Bool_t  fIsAbsEta;
  Float_t fEtaMaxCut;
  Float_t fEtaMinCut;
  Float_t fPtMaxCut;
  Float_t fPtMinCut;
  Int_t   fRunNumber; // for control of run changes

  ClassDef(AliSpherocityUtils,2) // base helper class
};
#endif

