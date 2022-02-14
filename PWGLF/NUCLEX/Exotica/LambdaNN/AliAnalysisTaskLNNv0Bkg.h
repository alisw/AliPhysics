#ifndef ALIANALYSISTASKLNNV0BKG_H
#define ALIANALYSISTASKLNNV0BKG_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskLNNv0Bkg class
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

class AliAnalysisTaskLNNv0Bkg : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskLNNv0Bkg();
  AliAnalysisTaskLNNv0Bkg(const char *name);
  virtual ~AliAnalysisTaskLNNv0Bkg();
  
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
  void SetBkgType(UInt_t type=1) {fBkgType=type;} // default is track rotation
  Double_t fMaxPtPion;
  Double_t fMinPtTriton; 
  void SetMC(){fMC=kTRUE;}
  Bool_t fMC;
 private:
 
  Double_t Chi2perNDF(AliESDtrack *tr); 


  UInt_t fYear;
  UInt_t fBkgType; // 0 is for opposite sign, 1 is for like sign, 2 is for track rotation
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

  AliAnalysisTaskLNNv0Bkg(const AliAnalysisTaskLNNv0Bkg&);            // not implemented
  AliAnalysisTaskLNNv0Bkg& operator=(const AliAnalysisTaskLNNv0Bkg&); // not implemented
  
  ClassDef(AliAnalysisTaskLNNv0Bkg, 1);
};

#endif
