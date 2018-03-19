#ifndef AliAnalysisTaskTOFppSpectra_cxx
#define AliAnalysisTaskTOFppSpectra_cxx


class AliVTrack;
class TH1F;
class AliESDEvent;
class AliVEvent;
class TH3F;
class TH2F;
class AliPIDResponse;
class AliESDtrackCuts;

#include <AliPID.h>
#include <AliPIDResponse.h>
#include "AliAnalysisTaskSE.h"
#include "AliTrackerBase.h"

#include "AliAnalysisUtils.h"

#include "AliMultiplicity.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliMultSelectionTask.h"
//#include "AliMultSelectionBase.h"


class AliAnalysisTaskTOFppSpectra : public AliAnalysisTaskSE {
 public:


  AliAnalysisTaskTOFppSpectra();
  AliAnalysisTaskTOFppSpectra(const char *Periodname, Int_t nTPC_CR, Int_t Chi2_TPCcluser, Int_t DCAz, Int_t DCAxy, Int_t value_Sigma, Float_t value_Slope);
  virtual ~AliAnalysisTaskTOFppSpectra() {}

	virtual void   UserCreateOutputObjects();
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t *);
	
//	virtual Float_t GetVertex(AliESDEvent* fESD) const;
	Bool_t TPCPID(AliESDtrack *track);
	Bool_t TOFPID(AliESDtrack *track);
	Double_t Rapidity(AliESDtrack *track , Double_t mass);

//	Double_t SecondaryTail(AliESDtrack *track);


//  Bool_t AcceptEvent(AliESDEvent *fESD) const;
	void	Exit(const char *msg);

	Bool_t	selectVertex2015pp(AliESDEvent *esd, Bool_t checkSPDres=kTRUE, Bool_t requireSPDandTrk=kFALSE, Bool_t checkProximity=kTRUE);
	Bool_t  IsGoodSPDvertexRes(const AliESDVertex *spdVertex);



	void SetMinNCrossedRowsTPC(Int_t MinTPCcr) {fMinTPCcr=MinTPCcr;}
        void SetMaxChi2PerClusterTPC(Int_t MaxChi2PerTPC) {fMaxChi2PerTPC=MaxChi2PerTPC;}
        void SetMaxDCAToVertexZ(Int_t MaxDCAz) {fMaxDCAz=MaxDCAz;}
        void SetDCAtoVertexXYPtDep(const char *MaxDCAxy) {fMaxDCAxy=MaxDCAxy;}
        
//	void SetParameter(2, Double_t Sigma) {fSigma=Sigma;}
//	void SetParameter(4, Double_t Slope) {fSlope=Slope;}

 private:
	



        AliPIDResponse *fPIDResponse;     //! PID response object
	AliMultSelection *fMultSelection;

	Int_t           fMinTPCcr;
        Int_t           fMaxChi2PerTPC;
        Int_t           fMaxDCAz;
        TString           fMaxDCAxy;
        Int_t           fSigma;
        Float_t           fSlope;


        AliESDtrackCuts *fesdTrackCuts;
        AliESDtrackCuts *fesdTrackCuts_no_dca;
	AliESDEvent	*fESD;    //! ESD object
	TList       	*fOutputList; //! Output list
	TH1F          	*fEventCounter;
	TH1F          	*fEventPS;
	TH1F          	*fEventVtx;
	TH1F          	*fEventVtx10;
	TH1F          	*fZVertex;
	
	TH3F          	*fTOFTimeV0MPtPi;
	TH3F          	*fTOFTimeV0MPtPosPi;
	TH3F          	*fTOFTimeV0MPtNegPi;

	TH3F          	*fTOFSigmaV0MPtPi;
	TH3F          	*fTOFSigmaV0MPtPosPi;
	TH3F          	*fTOFSigmaV0MPtNegPi;

	TH3F          	*fTOFResolutionV0MPtPi;
	TH3F          	*fTOFResolutionV0MPtPosPi;
	TH3F          	*fTOFResolutionV0MPtNegPi;

	TH3F          	*fTOFExpTimeV0MPtPi_El;
	TH3F          	*fTOFExpTimeV0MPtPosPi_El;
	TH3F          	*fTOFExpTimeV0MPtNegPi_El;
	TH3F          	*fTOFExpSigmaV0MPtPi_El;
	TH3F          	*fTOFExpSigmaV0MPtPosPi_El;
	TH3F          	*fTOFExpSigmaV0MPtNegPi_El;

	TH3F          	*fTOFExpTimeV0MPtPi_Mu;
	TH3F          	*fTOFExpTimeV0MPtPosPi_Mu;
	TH3F          	*fTOFExpTimeV0MPtNegPi_Mu;
	TH3F          	*fTOFExpSigmaV0MPtPi_Mu;
	TH3F          	*fTOFExpSigmaV0MPtPosPi_Mu;
	TH3F          	*fTOFExpSigmaV0MPtNegPi_Mu;

	TH3F          	*fTOFExpTimeV0MPtPi_Pi;
	TH3F          	*fTOFExpTimeV0MPtPosPi_Pi;
	TH3F          	*fTOFExpTimeV0MPtNegPi_Pi;
	TH3F          	*fTOFExpSigmaV0MPtPi_Pi;
	TH3F          	*fTOFExpSigmaV0MPtPosPi_Pi;
	TH3F          	*fTOFExpSigmaV0MPtNegPi_Pi;

	TH3F          	*fTOFExpTimeV0MPtPi_K;
	TH3F          	*fTOFExpTimeV0MPtPosPi_K;
	TH3F          	*fTOFExpTimeV0MPtNegPi_K;
	TH3F          	*fTOFExpSigmaV0MPtPi_K;
	TH3F          	*fTOFExpSigmaV0MPtPosPi_K;
	TH3F          	*fTOFExpSigmaV0MPtNegPi_K;

	TH3F          	*fTOFExpTimeV0MPtPi_P;
	TH3F          	*fTOFExpTimeV0MPtPosPi_P;
	TH3F          	*fTOFExpTimeV0MPtNegPi_P;
	TH3F          	*fTOFExpSigmaV0MPtPi_P;
	TH3F          	*fTOFExpSigmaV0MPtPosPi_P;
	TH3F          	*fTOFExpSigmaV0MPtNegPi_P;

	TH3F          	*fTOFExpTimeV0MPtPi_D;
	TH3F          	*fTOFExpTimeV0MPtPosPi_D;
	TH3F          	*fTOFExpTimeV0MPtNegPi_D;
	TH3F          	*fTOFExpSigmaV0MPtPi_D;
	TH3F          	*fTOFExpSigmaV0MPtPosPi_D;
	TH3F          	*fTOFExpSigmaV0MPtNegPi_D;

	TH3F          	*fTOFMismatchTimeV0MPtPi;
	TH3F          	*fTOFMismatchTimeV0MPtPosPi;
	TH3F          	*fTOFMismatchTimeV0MPtNegPi;

	TH3F          	*fTOFMismatchSigmaV0MPtPi;
	TH3F          	*fTOFMismatchSigmaV0MPtPosPi;
	TH3F          	*fTOFMismatchSigmaV0MPtNegPi;

	TH3F            *fTOFDCAxyV0MPtPi;
	TH3F            *fTOFDCAxyV0MPtPosPi;
	TH3F            *fTOFDCAxyV0MPtNegPi;

	TH2F            *fEventV0MPS;
        TH2F            *fEventV0MVtx;
	TH2F            *fEventV0M;

	TH1F            *fT0Resolution;
	TH1F            *fTimeOfFlightRes;
	TH1F            *fTimeOfFlightTOFRes;
	TH1F            *fTimeOfFlightGoodRes;

	TH2F            *fBetaP;
	TH2F            *fBetaPNoMismatch;
	TH2F            *fBetaPNoMismatchEtaCut;
	TH2F            *fBetaPt;
	TH2F            *fBetaPtNoMismatch;
	TH2F            *fBetaPtNoMismatchEtaCut;

	TH2F            *fTPCdEdxP;
	TH2F            *fTPCdEdxPt;
	
	TH2F            *fTPCTOFnSigmaPi;
	TH2F            *fTOFChannelVsTime;

	TH1F            *fGausTime;
	TH1F            *fTOFGausTime;



	TH3F          	*fTOFTimeV0MPtK;
	TH3F          	*fTOFTimeV0MPtPosK;
	TH3F          	*fTOFTimeV0MPtNegK;

	TH3F          	*fTOFSigmaV0MPtK;
	TH3F          	*fTOFSigmaV0MPtPosK;
	TH3F          	*fTOFSigmaV0MPtNegK;

	TH3F          	*fTOFResolutionV0MPtK;
	TH3F          	*fTOFResolutionV0MPtPosK;
	TH3F          	*fTOFResolutionV0MPtNegK;

	TH3F          	*fTOFExpTimeV0MPtK_El;
	TH3F          	*fTOFExpTimeV0MPtPosK_El;
	TH3F          	*fTOFExpTimeV0MPtNegK_El;
	TH3F          	*fTOFExpSigmaV0MPtK_El;
	TH3F          	*fTOFExpSigmaV0MPtPosK_El;
	TH3F          	*fTOFExpSigmaV0MPtNegK_El;

	TH3F          	*fTOFExpTimeV0MPtK_Mu;
	TH3F          	*fTOFExpTimeV0MPtPosK_Mu;
	TH3F          	*fTOFExpTimeV0MPtNegK_Mu;
	TH3F          	*fTOFExpSigmaV0MPtK_Mu;
	TH3F          	*fTOFExpSigmaV0MPtPosK_Mu;
	TH3F          	*fTOFExpSigmaV0MPtNegK_Mu;

	TH3F          	*fTOFExpTimeV0MPtK_Pi;
	TH3F          	*fTOFExpTimeV0MPtPosK_Pi;
	TH3F          	*fTOFExpTimeV0MPtNegK_Pi;
	TH3F          	*fTOFExpSigmaV0MPtK_Pi;
	TH3F          	*fTOFExpSigmaV0MPtPosK_Pi;
	TH3F          	*fTOFExpSigmaV0MPtNegK_Pi;

	TH3F          	*fTOFExpTimeV0MPtK_K;
	TH3F          	*fTOFExpTimeV0MPtPosK_K;
	TH3F          	*fTOFExpTimeV0MPtNegK_K;
	TH3F          	*fTOFExpSigmaV0MPtK_K;
	TH3F          	*fTOFExpSigmaV0MPtPosK_K;
	TH3F          	*fTOFExpSigmaV0MPtNegK_K;

	TH3F          	*fTOFExpTimeV0MPtK_P;
	TH3F          	*fTOFExpTimeV0MPtPosK_P;
	TH3F          	*fTOFExpTimeV0MPtNegK_P;
	TH3F          	*fTOFExpSigmaV0MPtK_P;
	TH3F          	*fTOFExpSigmaV0MPtPosK_P;
	TH3F          	*fTOFExpSigmaV0MPtNegK_P;

	TH3F          	*fTOFExpTimeV0MPtK_D;
	TH3F          	*fTOFExpTimeV0MPtPosK_D;
	TH3F          	*fTOFExpTimeV0MPtNegK_D;
	TH3F          	*fTOFExpSigmaV0MPtK_D;
	TH3F          	*fTOFExpSigmaV0MPtPosK_D;
	TH3F          	*fTOFExpSigmaV0MPtNegK_D;

	TH3F          	*fTOFMismatchTimeV0MPtK;
	TH3F          	*fTOFMismatchTimeV0MPtPosK;
	TH3F          	*fTOFMismatchTimeV0MPtNegK;

	TH3F          	*fTOFMismatchSigmaV0MPtK;
	TH3F          	*fTOFMismatchSigmaV0MPtPosK;
	TH3F          	*fTOFMismatchSigmaV0MPtNegK;

	TH3F            *fTOFDCAxyV0MPtK;
	TH3F            *fTOFDCAxyV0MPtPosK;
	TH3F            *fTOFDCAxyV0MPtNegK;

	TH2F            *fTPCTOFnSigmaK;

	

	TH3F          	*fTOFTimeV0MPtP;
	TH3F          	*fTOFTimeV0MPtPosP;
	TH3F          	*fTOFTimeV0MPtNegP;

	TH3F          	*fTOFSigmaV0MPtP;
	TH3F          	*fTOFSigmaV0MPtPosP;
	TH3F          	*fTOFSigmaV0MPtNegP;

	TH3F          	*fTOFResolutionV0MPtP;
	TH3F          	*fTOFResolutionV0MPtPosP;
	TH3F          	*fTOFResolutionV0MPtNegP;

	TH3F          	*fTOFExpTimeV0MPtP_El;
	TH3F          	*fTOFExpTimeV0MPtPosP_El;
	TH3F          	*fTOFExpTimeV0MPtNegP_El;
	TH3F          	*fTOFExpSigmaV0MPtP_El;
	TH3F          	*fTOFExpSigmaV0MPtPosP_El;
	TH3F          	*fTOFExpSigmaV0MPtNegP_El;

	TH3F          	*fTOFExpTimeV0MPtP_Mu;
	TH3F          	*fTOFExpTimeV0MPtPosP_Mu;
	TH3F          	*fTOFExpTimeV0MPtNegP_Mu;
	TH3F          	*fTOFExpSigmaV0MPtP_Mu;
	TH3F          	*fTOFExpSigmaV0MPtPosP_Mu;
	TH3F          	*fTOFExpSigmaV0MPtNegP_Mu;

	TH3F          	*fTOFExpTimeV0MPtP_Pi;
	TH3F          	*fTOFExpTimeV0MPtPosP_Pi;
	TH3F          	*fTOFExpTimeV0MPtNegP_Pi;
	TH3F          	*fTOFExpSigmaV0MPtP_Pi;
	TH3F          	*fTOFExpSigmaV0MPtPosP_Pi;
	TH3F          	*fTOFExpSigmaV0MPtNegP_Pi;

	TH3F          	*fTOFExpTimeV0MPtP_K;
	TH3F          	*fTOFExpTimeV0MPtPosP_K;
	TH3F          	*fTOFExpTimeV0MPtNegP_K;
	TH3F          	*fTOFExpSigmaV0MPtP_K;
	TH3F          	*fTOFExpSigmaV0MPtPosP_K;
	TH3F          	*fTOFExpSigmaV0MPtNegP_K;

	TH3F          	*fTOFExpTimeV0MPtP_P;
	TH3F          	*fTOFExpTimeV0MPtPosP_P;
	TH3F          	*fTOFExpTimeV0MPtNegP_P;
	TH3F          	*fTOFExpSigmaV0MPtP_P;
	TH3F          	*fTOFExpSigmaV0MPtPosP_P;
	TH3F          	*fTOFExpSigmaV0MPtNegP_P;

	TH3F          	*fTOFExpTimeV0MPtP_D;
	TH3F          	*fTOFExpTimeV0MPtPosP_D;
	TH3F          	*fTOFExpTimeV0MPtNegP_D;
	TH3F          	*fTOFExpSigmaV0MPtP_D;
	TH3F          	*fTOFExpSigmaV0MPtPosP_D;
	TH3F          	*fTOFExpSigmaV0MPtNegP_D;

	TH3F          	*fTOFMismatchTimeV0MPtP;
	TH3F          	*fTOFMismatchTimeV0MPtPosP;
	TH3F          	*fTOFMismatchTimeV0MPtNegP;

	TH3F          	*fTOFMismatchSigmaV0MPtP;
	TH3F          	*fTOFMismatchSigmaV0MPtPosP;
	TH3F          	*fTOFMismatchSigmaV0MPtNegP;

	TH3F            *fTOFDCAxyV0MPtP;
	TH3F            *fTOFDCAxyV0MPtPosP;
	TH3F            *fTOFDCAxyV0MPtNegP;


	TH2F            *fTPCTOFnSigmaP;
	TH1F            *fGausTime_K;
	TH1F            *fTOFGausTime_K;
	TH1F            *fGausTime_P;
	TH1F            *fTOFGausTime_P;

	
	TH3F          	*fTOFNoMismatchTimeV0MPtPi;
	TH3F          	*fTOFNoMismatchTimeV0MPtPosPi;
	TH3F          	*fTOFNoMismatchTimeV0MPtNegPi;
	TH3F          	*fTOFNoMismatchSigmaV0MPtPi;
	TH3F          	*fTOFNoMismatchSigmaV0MPtPosPi;
	TH3F          	*fTOFNoMismatchSigmaV0MPtNegPi;

	TH3F          	*fTOFNoMismatchTimeV0MPtK;
	TH3F          	*fTOFNoMismatchTimeV0MPtPosK;
	TH3F          	*fTOFNoMismatchTimeV0MPtNegK;
	TH3F          	*fTOFNoMismatchSigmaV0MPtK;
	TH3F          	*fTOFNoMismatchSigmaV0MPtPosK;
	TH3F          	*fTOFNoMismatchSigmaV0MPtNegK;


	TH3F          	*fTOFNoMismatchTimeV0MPtP;
	TH3F          	*fTOFNoMismatchTimeV0MPtPosP;
	TH3F          	*fTOFNoMismatchTimeV0MPtNegP;
	TH3F          	*fTOFNoMismatchSigmaV0MPtP;
	TH3F          	*fTOFNoMismatchSigmaV0MPtPosP;
	TH3F          	*fTOFNoMismatchSigmaV0MPtNegP;

	TH1F *fSec;
	TH1F *fSec_k;
	TH1F *fSec_p;
	TFile *fSecondary;
	TH1F *fV0MPC;
	TH1F *fV0MPC_vertexcut;

	TH1F *ftail;
	TH1F *ftail_Random;
	TH1F *fPtTPC_AllP;
	TH1F *fPtTPC_AllN;
	TH1F *fPtTOF_AllP;
	TH1F *fPtTOF_AllN;

	TH1F *fTPC_CR;
        TH1F *fChi2TPCcluster;
        TH1F *fDCAZ;
        TH1F *fDCAxy;

UInt_t   fTrigSel;


  AliAnalysisTaskTOFppSpectra(const AliAnalysisTaskTOFppSpectra&); // not implemented
  AliAnalysisTaskTOFppSpectra& operator=(const AliAnalysisTaskTOFppSpectra&); // not implemented

  ClassDef(AliAnalysisTaskTOFppSpectra, 1); // example of analysis
};

#endif
