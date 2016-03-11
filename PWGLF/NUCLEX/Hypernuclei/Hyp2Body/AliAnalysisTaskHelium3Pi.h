#ifndef ALIANALYSISTASKHELIUM3PI_H
#define ALIANALYSISTASKHELIUM3PI_H

/*  See cxx source for full Copyright notice */

//-----------------------------------------------------------------
//                 AliAnalysisTaskHelium3Pion class
//-----------------------------------------------------------------

#include <AliPIDResponse.h>
#include "TString.h"
#include "AliAnalysisTaskSE.h"

class TList;
class TH1F;
class TH2F;
class TH3F;
class TTree;
class AliESDtrackCuts;

class AliAnalysisTaskHelium3Pi : public AliAnalysisTaskSE {
 public:

  AliAnalysisTaskHelium3Pi(); 
  AliAnalysisTaskHelium3Pi(TString name);
  virtual ~AliAnalysisTaskHelium3Pi();
   
  virtual void  UserCreateOutputObjects();
  virtual void  UserExec(Option_t *option);
  virtual void  Terminate(Option_t *);
 
  Double_t BetheBloch(Double_t betaGamma,Double_t charge,Bool_t isPbPb);
  Bool_t  Flatten(Float_t cent);
  
  void SetCollidingSystems(Short_t collidingSystems = 1){fCollidingSystems = collidingSystems;}; 
  void SetAnalysisType(TString analysisType = "ESD"){fAnalysisType= analysisType;};
  void SetDataType(TString dataType = "PbPb"){fDataType = dataType;};
  void SetYear(Int_t year  = 2011){fYear= year;};
  void SetVzMax(Float_t Vzmax = 10){fVzmax= Vzmax;};
  void SetApplyFlatten(Bool_t  applyFlatten = kFALSE){fApplyFlatten= applyFlatten;};
  void SetFill3Htree(Bool_t  fill3hetree = kFALSE){fFill3Hetree= fill3hetree;};
  void ComputeFlow(Bool_t  doFlow = kFALSE){fDoFlow= doFlow;};
 
 private:

  TString fAnalysisType;	     //  "ESD" analysis type	
  Short_t fCollidingSystems;	     //  0 = pp collisions or 1 = AA collisions
  TString fDataType;	             //  pp, pPb or PbPb
  Int_t   fYear;                     //  2010, 2011, 2015
  Float_t fVzmax;                    //  Vz max
  Bool_t  fApplyFlatten;             //  Apply flatter
  Bool_t  fFill3Hetree ;             //  Store the 3He tree
  Bool_t  fDoFlow;                   //  Compute flow-related (SP) variables
 				        
  AliESDEvent *fESDevent;            //  
  AliVEvent   *fevent;               //  
  				        
  				        
  TList	*fListHist;	             //  List of  histograms
  
  TH1F *fHistEventMultiplicity;
  
  TH2F *fHistTrackMultiplicity;
  TH2F *fHistTrackMultiplicityCent;
  TH2F *fHistTrackMultiplicitySemiCent;
  TH2F *fHistTrackMultiplicityMB;
  TH2F *fHistTrackMultiplicityINT7;
  TH2F *fHistTrackMultiplicityPVCent;
  TH2F *fHistTrackMultiplicityPVSemiCent;
  TH2F *fHistTrackMultiplicityPVMB;
  TH2F *fHistTrackMultiplicityPVINT7;

  TH2F *fhBB;
  TH2F *fhTOF;
  TH1F *fhMassTOF;
  TH2F *fhBBPions;
  TH2F *fhBBHe;

  // For SP resolution
  TH2F *hQVzAQVzCvsCentrality;

  // Controll Histograms

  TH2F *hqEPCvsCentrality; 
  TH2F *hqEPAvsCentrality;
  TH2F *hqEPvsCentrality;

  TTree *fNtuple1;                  // Tree Pairs Pi/Proton "standard"
  
  Float_t teventtype           ;
  Float_t tTrackNumber         ;
  Float_t tpercentile          ;
  Float_t txPrimaryVertex      ;
  Float_t tyPrimaryVertex      ;
  Float_t tzPrimaryVertex      ;
  Float_t txSecondaryVertex    ;
  Float_t tySecondaryVertex    ;
  Float_t tzSecondaryVertex    ;
  Float_t tdcaTracks           ;
  Float_t tCosPointingAngle    ;
  Float_t tDCAV0toPrimaryVertex;
  Float_t tHeSign              ;
  Float_t tHepInTPC            ;
  Float_t tHeTPCsignal         ;
  Float_t tDcaHeToPrimVertex   ;
  Float_t tHeEta               ;
  Float_t tmomHex              ;
  Float_t tmomHey              ;
  Float_t tmomHez              ;
  Float_t tmomHeAtSVx          ;
  Float_t tmomHeAtSVy          ;
  Float_t tmomHeAtSVz          ;
  Float_t tHeTPCNcls           ;
  Float_t tHeimpactXY          ;
  Float_t tHeimpactZ           ;
  Float_t tHeITSClusterMap     ;
  Float_t tIsHeITSRefit        ;
  Float_t tPionSign            ;
  Float_t tPionpInTPC          ;
  Float_t tPionTPCsignal       ;
  Float_t tDcaPionToPrimVertex ;
  Float_t tPionEta             ;
  Float_t tmomPionx            ;
  Float_t tmomPiony            ;
  Float_t tmomPionz            ;
  Float_t tmomNegPionAtSVx     ;
  Float_t tmomNegPionAtSVy     ;
  Float_t tmomNegPionAtSVz     ;
  Float_t tPionTPCNcls         ;
  Float_t tPionimpactXY        ;
  Float_t tPionimpactZ         ;
  Float_t tPionITSClusterMap   ;
  Float_t tIsPiITSRefit        ;
  Float_t txn                  ;
  Float_t txp                  ;
  Float_t tuqV0A               ;
  Float_t tuqV0C               ;
  
  TTree *fNtuple4;                  // Tree He 

  Float_t  tHeleventtype   ;
  Float_t  tHelpercentile  ;
  Float_t  tHelSign	   ;
  Float_t  tHelpinTPC	   ;
  Float_t  tHelGetTPCsignal;
  Float_t  tHelPx	   ;
  Float_t  tHelPy	   ;
  Float_t  tHelPz	   ;
  Float_t  tHelEta	   ;
  Float_t  tHelisTOF	   ;
  Float_t  tHelTOFpull	   ;
  Float_t  tHeMass         ;
  Float_t  tHelimpactXY	   ;
  Float_t  tHelimpactZ	   ;
  Float_t  tHelmapITS      ;
  Float_t  tHelBetaTOF     ;
  Float_t  tHelIsITSrefit  ;
 
  //---------------------------------------------------------------------------
  AliESDtrackCuts *fESDtrackCuts; 
  AliPIDResponse  *fPIDResponse;      // pointer to PID response
  //_______________________________________________________________________


  AliAnalysisTaskHelium3Pi(const AliAnalysisTaskHelium3Pi&); // not implemented
  AliAnalysisTaskHelium3Pi& operator=(const AliAnalysisTaskHelium3Pi&); // not implemented

  ClassDef(AliAnalysisTaskHelium3Pi, 2);
};

#endif
