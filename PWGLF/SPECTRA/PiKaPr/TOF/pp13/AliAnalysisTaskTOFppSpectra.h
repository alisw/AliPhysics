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
  AliAnalysisTaskTOFppSpectra(const char *Periodname, Int_t nTPC_CR, Int_t Chi2_TPCcluser, Int_t DCAz, Int_t DCAxy, Int_t value_Sigma, Int_t value_Slope);
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
        Int_t           fSlope;


        AliESDtrackCuts *fesdTrackCuts;
        AliESDtrackCuts *fesdTrackCuts_no_dca;
	AliESDEvent	*fESD;    //! ESD object
	TList       	*fOutputList; //! Output list
	TH1F          	*fEventCounter;
	TH1F          	*fEventPS;
	TH1F          	*fEventVtx;
	TH1F          	*fEventVtx10;
	TH1F          	*fZVertex;
	
	TH3F          	*fTOFTimeV0MPtPosPi;
	TH3F          	*fTOFTimeV0MPtNegPi;

	TH3F          	*fTOFExpTimeV0MPtPosPi_El;
	TH3F          	*fTOFExpTimeV0MPtNegPi_El;

	TH3F          	*fTOFExpTimeV0MPtPosPi_Mu;
	TH3F          	*fTOFExpTimeV0MPtNegPi_Mu;

	TH3F          	*fTOFExpTimeV0MPtPosPi_Pi;
	TH3F          	*fTOFExpTimeV0MPtNegPi_Pi;

	TH3F          	*fTOFExpTimeV0MPtPosPi_K;
	TH3F          	*fTOFExpTimeV0MPtNegPi_K;

	TH3F          	*fTOFExpTimeV0MPtPosPi_P;
	TH3F          	*fTOFExpTimeV0MPtNegPi_P;

	TH3F          	*fTOFExpTimeV0MPtPosPi_D;
	TH3F          	*fTOFExpTimeV0MPtNegPi_D;

	TH3F          	*fTOFMismatchTimeV0MPtPosPi;
	TH3F          	*fTOFMismatchTimeV0MPtNegPi;

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

	TH3F          	*fTOFTimeV0MPtPosK;
	TH3F          	*fTOFTimeV0MPtNegK;

	TH3F          	*fTOFExpTimeV0MPtPosK_El;
	TH3F          	*fTOFExpTimeV0MPtNegK_El;

	TH3F          	*fTOFExpTimeV0MPtPosK_Mu;
	TH3F          	*fTOFExpTimeV0MPtNegK_Mu;

	TH3F          	*fTOFExpTimeV0MPtPosK_Pi;
	TH3F          	*fTOFExpTimeV0MPtNegK_Pi;

	TH3F          	*fTOFExpTimeV0MPtPosK_K;
	TH3F          	*fTOFExpTimeV0MPtNegK_K;

	TH3F          	*fTOFExpTimeV0MPtPosK_P;
	TH3F          	*fTOFExpTimeV0MPtNegK_P;

	TH3F          	*fTOFExpTimeV0MPtPosK_D;
	TH3F          	*fTOFExpTimeV0MPtNegK_D;

	TH3F          	*fTOFMismatchTimeV0MPtPosK;
	TH3F          	*fTOFMismatchTimeV0MPtNegK;


	TH3F            *fTOFDCAxyV0MPtPosK;
	TH3F            *fTOFDCAxyV0MPtNegK;

	TH2F            *fTPCTOFnSigmaK;
	

	TH3F          	*fTOFTimeV0MPtPosP;
	TH3F          	*fTOFTimeV0MPtNegP;

	TH3F          	*fTOFExpTimeV0MPtPosP_El;
	TH3F          	*fTOFExpTimeV0MPtNegP_El;

	TH3F          	*fTOFExpTimeV0MPtPosP_Mu;
	TH3F          	*fTOFExpTimeV0MPtNegP_Mu;

	TH3F          	*fTOFExpTimeV0MPtPosP_Pi;
	TH3F          	*fTOFExpTimeV0MPtNegP_Pi;

	TH3F          	*fTOFExpTimeV0MPtPosP_K;
	TH3F          	*fTOFExpTimeV0MPtNegP_K;

	TH3F          	*fTOFExpTimeV0MPtPosP_P;
	TH3F          	*fTOFExpTimeV0MPtNegP_P;

	TH3F          	*fTOFExpTimeV0MPtPosP_D;
	TH3F          	*fTOFExpTimeV0MPtNegP_D;

	TH3F          	*fTOFMismatchTimeV0MPtPosP;
	TH3F          	*fTOFMismatchTimeV0MPtNegP;


	TH3F            *fTOFDCAxyV0MPtPosP;
	TH3F            *fTOFDCAxyV0MPtNegP;


	TH2F            *fTPCTOFnSigmaP;
	
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
