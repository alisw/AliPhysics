#ifndef AliAnalysisTaskTOFMC_cxx
#define AliAnalysisTaskTOFMC_cxx

// example of an analysis task creating a p_t spectrum
class AliVTrack;
class TH1F;
class AliESDEvent;
class AliVEvent;
class TH2F;
class TH3F;
class AliPIDResponse;
class AliPIDCombined;
class AliESDtrackCuts;


class AliPPVsMultUtils;

#include <AliPID.h>
#include <AliPIDResponse.h>
#include <AliPIDCombined.h>
#include "AliAnalysisTaskSE.h"
#include "AliTrackerBase.h"

#include "AliAnalysisUtils.h"
#include "AliPPVsMultUtils.h"

#include "AliMultiplicity.h"
#include "AliOADBContainer.h"
#include "AliOADBMultSelection.h"
#include "AliMultSelectionTask.h"
#include "AliMultEstimator.h"
#include "AliMultVariable.h"
#include "AliMultInput.h"
#include "AliMultSelection.h"
#include "AliCentrality.h"
#include "AliMCParticle.h"


class AliAnalysisTaskTOFMC : public AliAnalysisTaskSE {
 public:


/*  AliAnalysisTaskTOFMC() : AliAnalysisTaskSE(),fESD(0), fOutputList(0),fPIDResponse(0),fesdTrackCuts(0x0),fesdTrackCuts_no_dca(0x0),
fTrigSel(AliVEvent::kINT7),fMultSelection(0x0),
fdEdxP(0),fdEdxPt(0),fdEdxPq(0),fdEdxPtq(0),fbetaAllPt(0),fbetaAllP(0),fbetaAllPtq(0),fbetaAllPq(0),
fPtVsTPion(0),fPtVsTKaon(0),fPtVsTProton(0),
fEventCounter(0),fEventPS(0),fEventVtx(0),fEventVtx10(0),fZVertex(0),fZVertexRec(0),fZVertexEff(0),fZVertex10(0),
fCorrRefMultVsNch(0),fCorrRefMultVsV0M(0),fV0MPC(0),

fPtV0MGenPion(0),fPtV0MTPCRecPion(0),fPtV0MTOFRecPion(0),
fPtV0MGenPionP(0),fPtV0MTPCRecPionP(0),fPtV0MTOFRecPionP(0),
fPtV0MGenPionM(0),fPtV0MTPCRecPionM(0),fPtV0MTOFRecPionM(0),

fPtV0MGenKaon(0),fPtV0MTPCRecKaon(0),fPtV0MTOFRecKaon(0),
fPtV0MGenKaonP(0),fPtV0MTPCRecKaonP(0),fPtV0MTOFRecKaonP(0),
fPtV0MGenKaonM(0),fPtV0MTPCRecKaonM(0),fPtV0MTOFRecKaonM(0),

fPtV0MGenProton(0),fPtV0MTPCRecProton(0),fPtV0MTOFRecProton(0),
fPtV0MGenProtonP(0),fPtV0MTPCRecProtonP(0),fPtV0MTOFRecProtonP(0),
fPtV0MGenProtonM(0),fPtV0MTPCRecProtonM(0),fPtV0MTOFRecProtonM(0),

fPtV0MDCAxyTOFPriPion(0),fPtV0MDCAxyTOFPriPionP(0),fPtV0MDCAxyTOFPriPionM(0),
fPtV0MDCAxyTOFWeakPion(0),fPtV0MDCAxyTOFWeakPionP(0),fPtV0MDCAxyTOFWeakPionM(0),
fPtV0MDCAxyTOFMatPion(0),fPtV0MDCAxyTOFMatPionP(0),fPtV0MDCAxyTOFMatPionM(0),

fPtV0MDCAxyTOFPriProton(0),fPtV0MDCAxyTOFPriProtonP(0),fPtV0MDCAxyTOFPriProtonM(0),
fPtV0MDCAxyTOFWeakProton(0),fPtV0MDCAxyTOFWeakProtonP(0),fPtV0MDCAxyTOFWeakProtonM(0),
fPtV0MDCAxyTOFMatProton(0),fPtV0MDCAxyTOFMatProtonP(0),fPtV0MDCAxyTOFMatProtonM(0),


fPPVsMultUtils(new AliPPVsMultUtils()),
fPtGenPion_kINT7(0),fPtGenPion_inel(0),fPtGenPion_signal_loss(0),
fPtGenKaon_kINT7(0),fPtGenKaon_inel(0),fPtGenKaon_signal_loss(0),
fPtGenProton_kINT7(0),fPtGenProton_inel(0),fPtGenProton_signal_loss(0),


fPtV0MTOFRecPion_nSigma(0),fPtV0MTOFRecPionP_nSigma(0),fPtV0MTOFRecPionM_nSigma(0),
fPtV0MTOFRecKaon_nSigma(0),fPtV0MTOFRecKaonP_nSigma(0),fPtV0MTOFRecKaonM_nSigma(0),
fPtV0MTOFRecProton_nSigma(0),fPtV0MTOFRecProtonP_nSigma(0),fPtV0MTOFRecProtonM_nSigma(0),

fPtV0MTOFRecPion_nSigma_excl(0),fPtV0MTOFRecPionP_nSigma_excl(0),fPtV0MTOFRecPionM_nSigma_excl(0),
fPtV0MTOFRecKaon_nSigma_excl(0),fPtV0MTOFRecKaonP_nSigma_excl(0),fPtV0MTOFRecKaonM_nSigma_excl(0),
fPtV0MTOFRecProton_nSigma_excl(0),fPtV0MTOFRecProtonP_nSigma_excl(0),fPtV0MTOFRecProtonM_nSigma_excl(0),

fTOFTimeV0MPtPi(0),fTOFTimeV0MPtK(0),fTOFTimeV0MPtP(0),
fTOFTimeV0MPtMismatchDecayPi(0),fTOFTimeV0MPtMismatchDecayK(0),fTOFTimeV0MPtMismatchDecayP(0), 

fEventV0MPS(0),fEventV0MVtx(0),fEventV0M(0),
fTOFLabel(),

fV0MPC_vertexcut(0),
fPtTPC_AllP(0),fPtTOF_AllP(0), fPtTPC_AllN(0),fPtTOF_AllN(0)




{}

*/
 AliAnalysisTaskTOFMC();
 AliAnalysisTaskTOFMC(const char *PeriodName, Int_t nTPC_CR,Int_t Chi2_TPCcluser, Int_t DCAz, Int_t DCAxy);
 //AliAnalysisTaskTOFMC(const char *name);
  virtual ~AliAnalysisTaskTOFMC() {}

	virtual void   UserCreateOutputObjects();
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t *);
	virtual Float_t GetVertex(AliESDEvent* fESD) const;
	Bool_t TPCPID(AliESDtrack *track);
	Bool_t TOFPID(AliESDtrack *track);
	Double_t Rapidity(AliESDtrack *track , Double_t mass);

	Bool_t  selectVertex2015pp(AliESDEvent *esd, Bool_t checkSPDres=kTRUE, Bool_t requireSPDandTrk=kFALSE, Bool_t checkProximity=kTRUE);
        Bool_t  IsGoodSPDvertexRes(const AliESDVertex *spdVertex);


//	AliESDtrackCuts *GetTrackCuts();// const {return *fesdTrackCuts;}; // getter
//	AliESDtrackCuts *GetTrackCuts() const {return *fesdTrackCuts_no_dca;}; // getter
//	void SetTrackCuts(AliESDtrackCuts *value) {fesdTrackCuts = value;}; // setter
//	void SetTrackCuts2(AliESDtrackCuts *value) {fesdTrackCuts_no_dca = value;}; // setter

	void SetMinNCrossedRowsTPC(Int_t MinTPCcr) {fMinTPCcr=MinTPCcr;}
	void SetMaxChi2PerClusterTPC(Int_t MaxChi2PerTPC) {fMaxChi2PerTPC=MaxChi2PerTPC;}
	void SetMaxDCAToVertexZ(Int_t MaxDCAz) {fMaxDCAz=MaxDCAz;}
	void SetDCAtoVertexXYPtDep(const char *MaxDCAxy) {fMaxDCAxy=MaxDCAxy;}
        

 private:


	AliPIDResponse *fPIDResponse;     //! PID response object
        AliMultSelection *fMultSelection;
        AliMultSelection *fMultSelection_INEL;
	AliESDtrackCuts *fesdTrackCuts;
        AliESDtrackCuts *fesdTrackCuts_no_dca;

	AliPPVsMultUtils *fPPVsMultUtils;


	Int_t           fMinTPCcr;
	Int_t           fMaxChi2PerTPC;
	Int_t           fMaxDCAz;
	TString           fMaxDCAxy;

	AliESDEvent *fESD;    //! ESD object
	TList       *fOutputList; //! Output list
  	TH1F        *fEventCounter;
	TH1F            *fEventPS;
        TH1F            *fEventVtx;
        TH1F            *fEventVtx10;
        TH1F            *fZVertex;
        TH1F            *fZVertexRec;
        TH1F            *fZVertexEff;
        TH1F            *fZVertex10;
        TH1F            *fV0MPC;
        TH2F            *fCorrRefMultVsNch;
        TH2F            *fCorrRefMultVsV0M;
	TH2F        *fdEdxP;
        TH2F        *fdEdxPt;
        TH2F        *fdEdxPq;
        TH2F        *fdEdxPtq;
	TH2F        *fbetaAllPt;
        TH2F        *fbetaAllP;
	TH2F        *fbetaAllPtq;
        TH2F        *fbetaAllPq;
	TH2F            *fPtVsTPion;
        TH2F            *fPtVsTKaon;
        TH2F            *fPtVsTProton;


        TH2F            *fPtV0MGenPion;
        TH2F            *fPtV0MTPCRecPion;
        TH2F            *fPtV0MTOFRecPion;
        TH2F            *fPtV0MGenPionP;
        TH2F            *fPtV0MTPCRecPionP;
        TH2F            *fPtV0MTOFRecPionP;
        TH2F            *fPtV0MGenPionM;
        TH2F            *fPtV0MTPCRecPionM;
        TH2F            *fPtV0MTOFRecPionM;

        TH2F            *fPtV0MGenKaon;
        TH2F            *fPtV0MTPCRecKaon;
        TH2F            *fPtV0MTOFRecKaon;
        TH2F            *fPtV0MGenKaonP;
        TH2F            *fPtV0MTPCRecKaonP;
        TH2F            *fPtV0MTOFRecKaonP;
        TH2F            *fPtV0MGenKaonM;
        TH2F            *fPtV0MTPCRecKaonM;
        TH2F            *fPtV0MTOFRecKaonM;

        TH2F            *fPtV0MGenProton;
        TH2F            *fPtV0MTPCRecProton;
        TH2F            *fPtV0MTOFRecProton;
        TH2F            *fPtV0MGenProtonP;
        TH2F            *fPtV0MTPCRecProtonP;
        TH2F            *fPtV0MTOFRecProtonP;
        TH2F            *fPtV0MGenProtonM;
        TH2F            *fPtV0MTPCRecProtonM;
        TH2F            *fPtV0MTOFRecProtonM;

        TH3F            *fPtV0MDCAxyTOFPriPion;
        TH3F            *fPtV0MDCAxyTOFWeakPion;
        TH3F            *fPtV0MDCAxyTOFMatPion;
        TH3F            *fPtV0MDCAxyTOFPriPionP;
        TH3F            *fPtV0MDCAxyTOFWeakPionP;
        TH3F            *fPtV0MDCAxyTOFMatPionP;
        TH3F            *fPtV0MDCAxyTOFPriPionM;
        TH3F            *fPtV0MDCAxyTOFWeakPionM;
        TH3F            *fPtV0MDCAxyTOFMatPionM;

        TH3F            *fPtV0MDCAxyTOFPriProton;
        TH3F            *fPtV0MDCAxyTOFWeakProton;
        TH3F            *fPtV0MDCAxyTOFMatProton;
        TH3F            *fPtV0MDCAxyTOFPriProtonP;
        TH3F            *fPtV0MDCAxyTOFWeakProtonP;
        TH3F            *fPtV0MDCAxyTOFMatProtonP;
        TH3F            *fPtV0MDCAxyTOFPriProtonM;
        TH3F            *fPtV0MDCAxyTOFWeakProtonM;
        TH3F            *fPtV0MDCAxyTOFMatProtonM;


        TH2F            *fPtV0MGenPion_kINT7;
        TH2F            *fPtV0MGenPion_inel;
        TH2F            *fPtV0MGenPion_signal_loss;
        TH2F            *fPtV0MGenKaon_kINT7;
        TH2F            *fPtV0MGenKaon_inel;
        TH2F            *fPtV0MGenKaon_signal_loss;
        TH2F            *fPtV0MGenProton_kINT7;
        TH2F            *fPtV0MGenProton_inel;
        TH2F            *fPtV0MGenProton_signal_loss;


        TH2F            *fPtV0MTOFRecPion_nSigma;
        TH2F            *fPtV0MTOFRecPionP_nSigma;
        TH2F            *fPtV0MTOFRecPionM_nSigma;
        TH2F            *fPtV0MTOFRecKaon_nSigma;
        TH2F            *fPtV0MTOFRecKaonP_nSigma;
        TH2F            *fPtV0MTOFRecKaonM_nSigma;
        TH2F            *fPtV0MTOFRecProton_nSigma;
        TH2F            *fPtV0MTOFRecProtonP_nSigma;
        TH2F            *fPtV0MTOFRecProtonM_nSigma;

        TH2F            *fPtV0MTOFRecPion_nSigma_excl;
        TH2F            *fPtV0MTOFRecPionP_nSigma_excl;
        TH2F            *fPtV0MTOFRecPionM_nSigma_excl;
        TH2F            *fPtV0MTOFRecKaon_nSigma_excl;
        TH2F            *fPtV0MTOFRecKaonP_nSigma_excl;
        TH2F            *fPtV0MTOFRecKaonM_nSigma_excl;
        TH2F            *fPtV0MTOFRecProton_nSigma_excl;
        TH2F            *fPtV0MTOFRecProtonP_nSigma_excl;
        TH2F            *fPtV0MTOFRecProtonM_nSigma_excl;


	TH3F            *fTOFTimeV0MPtPi;
	TH3F            *fTOFTimeV0MPtK;
	TH3F            *fTOFTimeV0MPtP;
	TH3F            *fTOFTimeV0MPtMismatchDecayPi;
	TH3F            *fTOFTimeV0MPtMismatchDecayK;
	TH3F            *fTOFTimeV0MPtMismatchDecayP;
	

	TH2F            *fEventV0MPS;
        TH2F            *fEventV0MVtx;
        TH2F            *fEventV0M;
	
        TH1F *fV0MPC_vertexcut;
        TH1F *fPtTPC_AllP;
        TH1F *fPtTPC_AllN;
        TH1F *fPtTOF_AllP;
        TH1F *fPtTOF_AllN;


	TH1F *fTPC_CR;
        TH1F *fChi2TPCcluster;
        TH1F *fDCAZ;
        TH1F *fDCAxy;

	Int_t fTOFLabel[3]; // TOF label


UInt_t   fTrigSel;


  AliAnalysisTaskTOFMC(const AliAnalysisTaskTOFMC&); // not implemented
  AliAnalysisTaskTOFMC& operator=(const AliAnalysisTaskTOFMC&); // not implemented

  ClassDef(AliAnalysisTaskTOFMC, 1); // example of analysis
};

#endif

