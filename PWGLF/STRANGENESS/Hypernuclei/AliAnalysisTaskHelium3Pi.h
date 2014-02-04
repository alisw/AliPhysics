#ifndef ALIANALYSISTASKHELIUM3PI_H
#define ALIANALYSISTASKHELIUM3PI_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskHelium3Pion class
//-----------------------------------------------------------------

class TList;
class TH1F;
class TH2F;
class TH3F;
class TNtuple;
class AliESDcascade;
//class AliCascadeVertexer; 

#include "TString.h"

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskHelium3Pi : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskHelium3Pi();
  AliAnalysisTaskHelium3Pi(const char *name);
  virtual ~AliAnalysisTaskHelium3Pi();
  
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);
  
  void SetCollidingSystems(Short_t collidingSystems = 0)     {fCollidingSystems = collidingSystems;}
  void SetAnalysisType    (const char* analysisType = "ESD") {fAnalysisType = analysisType;}
  void SetDataType    (const char* dataType = "REAL") {fDataType = dataType;}
  
  Double_t BetheBloch(Double_t bg,Double_t Charge,Bool_t isPbPb);

  
 private:
  
  TString fAnalysisType;	     //! "ESD" or "AOD" analysis type	
  
  Short_t fCollidingSystems;	     //! 0 = pp collisions or 1 = AA collisions
  TString fDataType;		     //! "REAL" or "SIM" data type	
  TList	*fListHistCascade;	     //! List of Cascade histograms
  TH1F *fHistEventMultiplicity;
  TH2F *fHistTrackMultiplicity;
  TH2F *fHistTrackMultiplicityCent;
  TH2F *fHistTrackMultiplicitySemiCent;
  TH2F *fHistTrackMultiplicityMB;
  TH2F *fHistTrackMultiplicityPVCent;
  TH2F *fHistTrackMultiplicityPVSemiCent;
  TH2F *fHistTrackMultiplicityPVMB;
  TH1F *fHistMult;
  TH2F *fhBB;
  TH2F *fhTOF;
  TH1F *fhMassTOF;
  TH2F *fhBBPions;
  TH2F *fhBBHe;
  TH2F *fhNaPos;
  TH2F *fhNaNeg;
  TH2F *fBetavsTPCsignalPos;
  TH2F *fBetavsTPCsignalNeg;
  TH2F *fHelium3TOF;
   
  TNtuple *fNtuple1;                  //! Ntupla Pairs Pi/Proton "standard"
  TNtuple *fNtuple4;                  //! He carateristiche

  static const Int_t fgNrot;
 

  AliAnalysisTaskHelium3Pi(const AliAnalysisTaskHelium3Pi&);            // not implemented
  AliAnalysisTaskHelium3Pi& operator=(const AliAnalysisTaskHelium3Pi&); // not implemented
  
  ClassDef(AliAnalysisTaskHelium3Pi, 0);
};

#endif
