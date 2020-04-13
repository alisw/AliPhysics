#ifndef ALIANALYSISTASKHYPV2PBPB18_H
#define ALIANALYSISTASKHYPV2PBPB18_H

// ROOT includes
#include <TList.h>
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TF1.h>
#include <TString.h>
#include <THnSparse.h>

// AliRoot includes
#include "AliAnalysisTaskSE.h"
#include "AliAODEvent.h"
#include "AliVHeader.h"
#include "AliVVertex.h"
#include "AliVEvent.h"
#include "AliVTrack.h"
#include "AliPIDResponse.h"
#include "AliEventCuts.h"


class AliAnalysisTaskHypv2PbPb18 : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskHypv2PbPb18();
  AliAnalysisTaskHypv2PbPb18(const char *name);

  virtual ~AliAnalysisTaskHypv2PbPb18();

  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *); 

  virtual void Initialize();

  Bool_t   PassedBasicTrackQualityCuts_Pos (AliESDtrack *track);
  Bool_t   PassedBasicTrackQualityCuts_Neg (AliESDtrack *track);
  Double_t GetTransverseDCA                (AliESDtrack *track);
  Double_t GetDCA                          (AliESDtrack *track);
  Bool_t   PassedMinimalQualityCutsV0      (AliESDv0 *V0);
  Bool_t   IsHyperTritonCandidate          (AliESDv0 *V0);
  Double_t GetDecayLengthV0                (AliESDv0 *V0);
  Bool_t   Is3HeCandidate                  (AliESDtrack *track);
  Bool_t   IsPionCandidate                 (AliESDtrack *track);
  Double_t InvariantMassHypertriton        (TVector3 P1, TVector3 P2);


  Float_t GetPhi0Pi(Float_t phi);
 
  void SetVzMax(Float_t Vzmax)                   {fVzmax     = Vzmax;    };
  void SetCentralityEstimator(Short_t centEst)   {fCenCalV0  = centEst;  };
  void SetHarmonic(Int_t  harmonic)              {fNHarm     = harmonic; };
  void SetPeriod(Bool_t period)                  {fPeriod    = period;};
  
  AliEventCuts fEventCuts;

 private:
  
  void Analyze(AliVEvent* esd, Double_t vz, Int_t evttype);
  //  void Analyze(AliVEvent* aod);
  void OpenInfoCalbration(Int_t run);
  
  AliESDEvent *fESDevent;                         // 
  AliVEvent   *fevent;                            // 
  
  Int_t        fRun;                // run number - for calibration
  TH1D*        fMultV0;             // profile from V0 multiplicity
  
  TH1D*        fQxnmV0A;            // <Qx2> V0A
  TH1D*        fQynmV0A;            // <Qy2> V0A
  TH1D*        fQxnsV0A;            // sigma Qx2 V0A
  TH1D*        fQynsV0A;            // sigma Qy2 V0A
  TH1D*        fQxnmV0C;            // <Qx2> V0C
  TH1D*        fQynmV0C;            // <Qy2> V0C
  TH1D*        fQxnsV0C;            // sigma Qx2 V0C
  TH1D*        fQynsV0C;            // sigma Qy2 V0C
    
  Double_t     fNHarm;              // harmonic number: 2, 3
  Short_t      fRecPass;            // flag for reconstruction pass: 0->Fast, 1->NoSDD, 2->SDD
  Short_t      fCenCalV0;           // flag for centrality estimators used for V0 recentering: 0->V0A, 1->V0, 2->V0AEq, 3-> CL1
  Short_t      fFilterBit;          // flag for AOD filterbit
  
  Int_t fptc ;
  Float_t fVzmax;
  Bool_t fPeriod;
  
  //output hist
  TList	*fListHist;	           // List of  histograms
  
  TH1F  *fHistEventMultiplicity;           // event multiplicity
  TH2F  *fHistTrackMultiplicity;           // track multiplicity
  
  TH2F  *fhBB;                             // ScatterPlot Total
  TH2F  *fhBBHyp;                          // ScatterPlot Total
  TH2F  *fhBBAHyp;                          // ScatterPlot Total
  
  // Event Plane vs Centrality
  
  TH2D *EPVzAvsCentrality  ; 
  TH2D *EPVzCvsCentrality  ; 
  
  // For SP resolution
  
  TH2F *hQVzAQVzCvsCentrality;
  TH2F *hQVzAQTPCvsCentrality;
  TH2F *hQVzCQTPCvsCentrality;
  // For NUA correction
  
  TH2F *hQxVzAvsCentrality;
  TH2F *hQyVzAvsCentrality;
  TH2F *hQxVzCvsCentrality;
  TH2F *hQyVzCvsCentrality;
  
  // for EP
  TH2F *hCos2DeltaTPCVzAvsCentrality;
  TH2F *hCos2DeltaTPCVzCvsCentrality;
  TH2F *hCos2DeltaVzAVzCvsCentrality;
  TH2F *hCos2DeltaVzATPCvsCentrality;
  TH2F *hCos2DeltaVzCTPCvsCentrality;
  TH2F *hCos2DeltaVzCVzAvsCentrality;

  Int_t eventtype;
  
  // TTree
  TTree *ftree;//!
  
  //Global Variables
  Int_t    iEvent;//
  Double_t zVertex;//
  Double_t centrality;//
  
  //Variables for HyperTriton - First Daughter
  Double_t px_Daughter1;//
  Double_t py_Daughter1;//
  Double_t pz_Daughter1;//
  Int_t    q_Daughter1;//
  Double_t dcaxy_Daughter1;//
  Int_t    nTPC_Clusters_Daughter1;//
  Int_t    nTPC_Clusters_dEdx_Daughter1;//
  Double_t chi2_TPC_Daughter1;//
  Double_t nSigmaTPC_He3_Daughter1;//
  Double_t nSigmaTPC_Pion_Daughter1;//
  
  //Variables for HyperTriton - Second Daughter
  Double_t px_Daughter2;//
  Double_t py_Daughter2;//
  Double_t pz_Daughter2;//
  Int_t    q_Daughter2;//
  Double_t dcaxy_Daughter2;//
  Int_t    nTPC_Clusters_Daughter2;//
  Int_t    nTPC_Clusters_dEdx_Daughter2;//
  Double_t chi2_TPC_Daughter2;//
  Double_t nSigmaTPC_He3_Daughter2;//
  Double_t nSigmaTPC_Pion_Daughter2;//
    
  //Pair Variables
  Int_t    isOnTheFlyV0;//
  Double_t cosPointingAngle;//
  Double_t dcaV0Daughters;//
  Double_t radius;//
  Double_t chi2V0;//
  Double_t decayLength;//

  //flow variables
  Double_t deltaphiV0A;//
  Double_t deltaphiV0C;//
  Double_t uqV0A;//
  Double_t uqV0C;//
    
  //-------------------------
  AliPIDResponse  *fPIDResponse;   //! pointer to PID response
  AliESDtrackCuts *fESDtrackCuts_Pos;
  AliESDtrackCuts *fESDtrackCuts_Neg;
  AliESDtrackCuts *fESDtrackCutsEP;
  //--------------------------

  AliAnalysisTaskHypv2PbPb18(const AliAnalysisTaskHypv2PbPb18&);
  AliAnalysisTaskHypv2PbPb18& operator=(const AliAnalysisTaskHypv2PbPb18&);
  
  ClassDef(AliAnalysisTaskHypv2PbPb18, 1);    //Analysis task for high pt analysis
  
};

#endif
