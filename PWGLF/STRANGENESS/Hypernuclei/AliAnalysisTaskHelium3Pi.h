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
class TTree;
class AliESDtrackCuts;

#include <AliPIDResponse.h>

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

  //  Bool_t IsTrackAccepted(AliVTrack *track);
  
  private:
  
  TString fAnalysisType;	     //! "ESD" or "AOD" analysis type	
  

  Short_t fCollidingSystems;	     //! 0 = pp collisions or 1 = AA collisions
  
  AliESDtrackCuts *fESDtrackCuts; 

  TString fDataType;		     //! "REAL" or "SIM" data type	

  TList	*fListHist;	             //! List of  histograms

  TH1F *fHistEventMultiplicity;
  TH2F *fHistTrackMultiplicity;
  TH2F *fHistTrackMultiplicityCent;
  TH2F *fHistTrackMultiplicitySemiCent;
  TH2F *fHistTrackMultiplicityMB;
  TH2F *fHistTrackMultiplicityPVCent;
  TH2F *fHistTrackMultiplicityPVSemiCent;
  TH2F *fHistTrackMultiplicityPVMB;
  TH2F *fhBB;
  TH2F *fhTOF;
  TH1F *fhMassTOF;
  TH2F *fhBBPions;
  TH2F *fhBBHe;
  TH2F *fhNaPos;
  TH2F *fhNaNeg;
  TH2F *fBetavsTPCsignalPos;
  TH2F *fBetavsTPCsignalNeg;
   
  TTree *fNtuple1;                  //! Tree Pairs Pi/Proton "standard"
  
  Float_t trunNumber;
  Float_t tbunchcross;
  Float_t torbit;
  Float_t tperiod;
  Float_t teventtype;
  Float_t tTrackNumber;
  Float_t tpercentile;
  Float_t txPrimaryVertex;
  Float_t tyPrimaryVertex;
  Float_t tzPrimaryVertex;
  Float_t txSecondaryVertex;
  Float_t tySecondaryVertex;
  Float_t tzSecondaryVertex;
  Float_t tdcaTracks;
  Float_t tCosPointingAngle;
  Float_t tDCAV0toPrimaryVertex;
  Float_t tHeSign;
  Float_t tHepInTPC;
  Float_t tHeTPCsignal;
  Float_t tDcaHeToPrimVertex;
  Float_t tHeEta;
  Float_t tmomHex;
  Float_t tmomHey;
  Float_t tmomHez;
  Float_t tmomHeAtSVx;
  Float_t tmomHeAtSVy;
  Float_t tmomHeAtSVz;
  Float_t tHeTPCNcls;
  Float_t tHeimpactXY;
  Float_t tHeimpactZ;
  Float_t tHeITSClusterMap;
  Float_t tIsHeITSRefit;
  Float_t tPionSign;
  Float_t tPionpInTPC;
  Float_t tPionTPCsignal;
  Float_t tDcaPionToPrimVertex;
  Float_t tPionEta;
  Float_t tmomPionx;
  Float_t tmomPiony;
  Float_t tmomPionz;
  Float_t tmomNegPionAtSVx;
  Float_t tmomNegPionAtSVy;
  Float_t tmomNegPionAtSVz;
  Float_t tPionTPCNcls;
  Float_t tPionimpactXY;
  Float_t tPionimpactZ;
  Float_t tPionITSClusterMap;
  Float_t tIsPiITSRefit;
  Float_t txn;
  Float_t txp;
  Float_t tchi2He;
  Float_t tchi2Pi;
  
  TTree *fNtuple4;                  //! Tree He 

  Float_t tHelrunNumber;
  Float_t tHelBCNumber;
  Float_t tHelOrbitNumber;
  Float_t tHelPeriodNumber;
  Float_t tHeleventtype;
  Float_t tHelisHeITSrefit;
  Float_t tHelpercentile;
  Float_t tHelSign;
  Float_t tHelpinTPC;
  Float_t tHelGetTPCsignal;
  Float_t tHelPx;
  Float_t tHelPy;
  Float_t tHelPz;
  Float_t tHelEta;
  Float_t tHelisTOF;
  Float_t tHelpoutTPC;
  Float_t tHeltimeTOF;
  Float_t tHeltrackLenghtTOF;
  Float_t tHelimpactXY;
  Float_t tHelimpactZ;
  Float_t tHelmapITS;
  Float_t tHelTPCNcls;
  Float_t tHelTRDsignal;
  Float_t tHelxPrimaryVertex;
  Float_t tHelyPrimaryVertex;
  Float_t tHelzPrimaryVertex;
  Float_t tHelchi2PerClusterTPC;
  
  static const Int_t fgNrot;
 
  AliPIDResponse *fPIDResponse;     //! pointer to PID response

  AliAnalysisTaskHelium3Pi(const AliAnalysisTaskHelium3Pi&);            // not implemented
  AliAnalysisTaskHelium3Pi& operator=(const AliAnalysisTaskHelium3Pi&); // not implemented
  
  ClassDef(AliAnalysisTaskHelium3Pi, 0);
};

#endif
