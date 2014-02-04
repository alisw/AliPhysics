#ifndef ALIANALYSISTASKHELIUM3PIMC_H
#define ALIANALYSISTASKHELIUM3PIMC_H

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

class AliAnalysisTaskHelium3PiMC : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskHelium3PiMC();
  AliAnalysisTaskHelium3PiMC(const char *name);
  virtual ~AliAnalysisTaskHelium3PiMC();
  
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);
  
  void SetCollidingSystems(Short_t collidingSystems = 0)     {fCollidingSystems = collidingSystems;}
  void SetAnalysisType    (const char* analysisType = "ESD") {fAnalysisType = analysisType;}
  void SetDataType    (const char* dataType = "SIM") {fDataType = dataType;}

  Double_t BetheBloch(Double_t bg,Double_t Charge,Bool_t isPbPb);

  
 private:

  TString fAnalysisType;	     //! "ESD" or "AOD" analysis type	

  Short_t fCollidingSystems;	     //! 0 = pp collisions or 1 = AA collisions
  TString fDataType;		     //! "REAL" or "SIM" data type	
  TList	*fListHistCascade;	     //! List of Cascade histograms

  TH1F *fHistEventMultiplicity;      //! event multiplicity
  TH1F *fHistTrackMultiplicity;      //! track multiplicity  
  TH1F *fHistMCMultiplicityTracks;
  TH1F *fHistMCEta; 
  TH1F *fHistMCPt; 
  TH1F *fHistMCTheta; 
  TH1F *fHistMCDecayPosition;
  TH1F *fHistMCDecayRho; 
  TH2F *fhRigidityHevsMomPiMC;
  TH2F *fhRigidityHevsMomPiRec;
  TH1F *fhInvMassMC;
  TH1F *fhInvMassMum;
  TH1F *fhInvMassRec;

  TH1F *fhInvMassRec1;
  TH1F *fhInvMassRec2;
  TH1F *fhInvMassRec3;
  TH1F *fhInvMassRec4;
  TH1F *fhInvMassRec5;
  TH1F *fhInvMassRec6;
  TH1F *fhInvMassRec7;

  TH1F *fhHeMCRigidity;
  TH1F *fhPioneMC;
  TH2F *hBBTPCnoCuts;
  TH2F *fhBBTPC;
  TH2F *fhBBTPCNegativePions;
  TH2F *fhBBTPCPositivePions;
  TH2F *fhBBTPCHe3;
  TH2F *fHistProvaDCA;
  TH2F *fHistPercentileVsTrackNumber;
  TH1F *fhHeDCAXY; 
  TH1F *fhHeDCAZ; 
  TH1F *fhPiDCAXY; 
  TH1F *fhPiDCAZ; 
  TH1F *hITSClusterMap;
  
 
  TNtuple *fNtuple1;                  //! Ntupla Pairs Pi/Proton "standard"
  TNtuple *fNtuple2;                  //! Ntupla Pairs PiPos/Proton "background"


 static const Int_t fgNrot;
 

  AliAnalysisTaskHelium3PiMC(const AliAnalysisTaskHelium3PiMC&);            // not implemented
  AliAnalysisTaskHelium3PiMC& operator=(const AliAnalysisTaskHelium3PiMC&); // not implemented
  
  ClassDef(AliAnalysisTaskHelium3PiMC, 0);
};

#endif
