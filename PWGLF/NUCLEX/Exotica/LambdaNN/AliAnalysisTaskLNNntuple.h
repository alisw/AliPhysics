#ifndef ALIANALYSISTASKNTUPLE_H
#define ALIANALYSISTASKNTUPLE_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskLNNntuple class
//-----------------------------------------------------------------

class TList;
class TH1F;
class TH2F;
class TH3F;
class TTree;
class TNtupleD;
class TVector3;
class AliESDtrack;
class AliESDv0;

#include <AliPIDResponse.h>

#include "TString.h"

#include "AliAnalysisTaskSE.h"
#include "AliEventCuts.h"

class AliAnalysisTaskLNNntuple : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskLNNntuple();
  AliAnalysisTaskLNNntuple(const char *name, Bool_t mc=kFALSE);
  virtual ~AliAnalysisTaskLNNntuple();
  
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);
  
  Bool_t PassTrackCuts(AliESDtrack *tr);
  Bool_t Passv0Cuts(AliESDv0 *v0, Double_t decayLength);
  Double_t GetTOFmass(AliESDtrack *tr);

  Int_t SelectCentrality(Float_t perc, Bool_t isPrimVtx=kFALSE);
  Bool_t IsPionCandidate(AliESDtrack *tr);
  Bool_t IsTritonCandidate(AliESDtrack *tr);

  void SetPtLimits(Double_t pionMax, Double_t tritonMin) {fMaxPtPion=pionMax; fMinPtTriton=tritonMin;}
  void SetYear(UInt_t year){fYear=year;}
  Double_t fMaxPtPion;
  Double_t fMinPtTriton; 

 private:
 
  Double_t Chi2perNDF(AliESDtrack *tr); 


  UInt_t fYear;
  Bool_t fMC;
  AliEventCuts fEventCuts;

  TList	*fListHist;	             //! List of  histograms
  TH1F *fHistEventMultiplicity;
  TH2F *fHistTrackMultiplicity;
  TH2F *fhBB;
  TH2F *fhTOF;
  TH2F *fhMassTOF;
  TH2F *fhBBPions;
  TH2F *fhBBH3;
  TH2F *fhBBH3TofSel;
  TH2F *fhTestNsigma;
  TH2F *fTPCclusPID;
  TH2F *fhTestQ;
 
  
  TNtupleD *fNt; //! 
  AliPIDResponse *fPIDResponse;     //! pointer to PID response

  AliAnalysisTaskLNNntuple(const AliAnalysisTaskLNNntuple&);            // not implemented
  AliAnalysisTaskLNNntuple& operator=(const AliAnalysisTaskLNNntuple&); // not implemented
  
  ClassDef(AliAnalysisTaskLNNntuple, 1);
};

#endif
