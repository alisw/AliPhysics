#ifndef ALIANALYSISTASKFLOWTPCTOFQCSP_H
#define ALIANALYSISTASKFLOWTPCTOFQCSP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  Task for Heavy Flavour Electron Flow  TPC plus TOF                //
//                                                                    //
//  Author: Andrea Dubla (Utrecht University)                         //
//                                                                    //
////////////////////////////////////////////////////////////////////////

class THnSparse;
class TH2F;
class TLorentzVector;
class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliAnalysisFilter;
class AliESDtrackCuts;
class AliESDtrack;
class AliAODtrack;
class AliHFEcontainer;
class AliHFEcuts;
//class AliHFEpid;
//class AliHFEpidQAmanager;
class AliCFManager;
class AliFlowTrackCuts;
class AliFlowTrack;
class AliFlowEvent;
class AliFlowCandidateTrack;
class AliFlowEventSimple;
class AliCentrality;
class AliPIDResponse;
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskFlowTPCTOFQCSP : public AliAnalysisTaskSE {

  public:
    AliAnalysisTaskFlowTPCTOFQCSP();
    AliAnalysisTaskFlowTPCTOFQCSP(const char *name);
    virtual ~AliAnalysisTaskFlowTPCTOFQCSP();

    void                                 SetEnableDebugMode() {fDebug = kTRUE; };
    void                                 SetCentralityParameters(Double_t CentralityMin, Double_t CentralityMax, const char* CentralityMethod); //select centrality
    void                                 CheckCentrality(AliAODEvent *event,Bool_t &centralitypass); //to use only events with the correct centrality....
    void                                 SelectPhotonicElectron(Int_t itrack,const AliAODTrack *track, Bool_t &fFlagPhotonicElec);
    void                                 SetInvariantMassCut(Double_t invmass) {fInvmassCut = invmass;};
    void                                 SetTrigger(Int_t trig) {fTrigger = trig;};
    void                                 SetTPCS(Int_t sig) {fTPCS = sig;};
    void                                 SetVz(Int_t ver) {fVz = ver;};
    void                                 SetQAPIDSparse(Bool_t qapidsparse) {fQAPIDSparse = qapidsparse;};
    void                                 SetPhiCut(Bool_t phicut){fPhiCut = phicut;};
    template <typename T> void           PlotVZeroMultiplcities(const T* event) const;
    template <typename T> void           SetNullCuts(T* aod);
    void                                 PrepareFlowEvent(Int_t iMulti, AliFlowEvent *FlowEv) const;
    virtual void                         UserCreateOutputObjects();
    virtual void                         UserExec(Option_t *option);
    virtual void                         Terminate(Option_t *);
    void                                 SetRPCuts(AliFlowTrackCuts *cutsRP) { fCutsRP = cutsRP; }
    void                                 SetHFECuts(AliHFEcuts * const cuts) { fCuts = cuts; };
    void                                 SetIDCuts(Double_t minTPCnsigma, Double_t maxTPCnsigma, Double_t minTOFnSigma, Double_t maxTOFnSigma);



  //  AliHFEpid *GetPID() const { return fPID; };

  private:

    Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);

    Bool_t               fDebug; //! enable debug mode
    AliAODEvent 	     *fAOD;			//AOD object
    TList       	   	 *fOutputList;		//output list
    AliHFEcuts   		 *fCuts;                 //Cut Collection
    Bool_t 			     fIdentifiedAsOutInz;    //Out Of Range in z
    Bool_t 	   	   	     fPassTheEventCut;       //Pass The Event Cut
    AliCFManager 	   	 *fCFM;                  //!Correction Framework Manager
 //   AliHFEpid 		     *fPID;                  //PID
 //   AliHFEpidQAmanager   *fPIDqa;		//! PID QA manager
    AliFlowTrackCuts     *fCutsRP; // track cuts for reference particles
    AliFlowTrackCuts     *fNullCuts; // dummy cuts for flow event tracks
    AliFlowEvent         *fFlowEvent; //! flow events Inclusive e 
    const char           *fkCentralityMethod; // method used to determine centrality (V0 by default)
    Double_t             fCentrality; // event centrality for QA
    Double_t             fCentralityMin; // lower bound of cenrality bin
    Double_t             fCentralityMax; // upper bound of centrality bin
    Double_t	         fInvmassCut;		//invariant mass cut value
    Int_t	             fTrigger;		//invariant mass cut value
    TH1F                 *fPhi; //! QA plot of azimuthal distribution of tracks used for event plane estimation
    TH1F                 *fEta; //! QA plot of eta distribution of tracks used for event plane estimation
    TH1F                 *fVZEROA; //! QA plot vzeroa multiplicity (all tracks in event)
    TH1F                 *fVZEROC; //! QA plot vzeroc multiplicity (all tracks in event)
    TH1F                 *fTPCM; //! QA plot TPC multiplicity (tracks used for event plane estimation)
    TH1F		         *fNoEvents;		//no of events
    TH1F                 *fInclusiveElecPt; // Inclusive elec pt
//    TH2F			     *fTPCnsigma;		//TPC n sigma vs p
    TH2F		         *fTPCnsigmaAft;		//TPC n sigma vs p after HFE pid
//    TH2F			     *fTOFns;		//track TOFnSigma
    TH2F			     *fTOFnsAft;		//track TOFnSigma after PID
    TH2F			     *fTOFBetaAft;
    TH1F                 *fCentralityPass; // ! QA histogram of events that pass centrality cut
    TH1F                 *fCentralityNoPass; //! QA histogram of events that do not pass centrality cut
    TH1F                 *fInvmassLS1; //LS Invmass for all rec par
    TH1F                 *fInvmassULS1;//ULS Invmass for all rec par
    TH1F			     *fPhotoElecPt;		//photonic elec pt
    TH1F			     *fSemiInclElecPt;	//Semi inclusive ele pt
    TH1F                 *fULSElecPt; //ULS elec Pt
    TH1F                 *fLSElecPt;// LS elec pt
    TH2F                 *fMultCorAfterCuts; //! QA profile global and tpc multiplicity after outlier cut
    TH2F                 *fMultvsCentr; //! QA profile of centralty vs multiplicity
	TProfile			 *fSubEventDPhiv2;
	TH1D 		         *EPVzA;//ep v0a
	TH1D 		         *EPVzC;//ep v0c
	TH1D 		         *EPTPC;//ep tpc
    THnSparseF           *fV2Phi;//! v2 analysis of EP-V0
    TH1D                 *fvertex;//vertex
    TH2F                 *fMultCorBeforeCuts; //! QA profile global and tpc multiplicity after outlier cut
    AliPIDResponse       *fPIDResponse;//!  pidresponse
    THnSparseF           *fQAPid;             //! v2 analysis of EP-V0
    Double_t              fminTPCnsigma;  //ID cuts tpc
    Double_t              fmaxTPCnsigma;  //ID cuts tpc
    Bool_t 	   	   	      fQAPIDSparse;       //QAPIDSPARSE
    Double_t              fminTOFnSigma;  //ID cuts tof
    Double_t              fmaxTOFnSigma;//ID cuts tof
    THnSparseF           *fQAPidSparse;             //! v2 analysis of EP-V0
    Int_t                 fTPCS; //tpc signal cluster
    Int_t                 fVz; //vertex range
    Bool_t 	   	   	      fPhiCut;       //Phi cut to simulate emcal acc

    AliAnalysisTaskFlowTPCTOFQCSP(const AliAnalysisTaskFlowTPCTOFQCSP&); // not implemented
    AliAnalysisTaskFlowTPCTOFQCSP& operator=(const AliAnalysisTaskFlowTPCTOFQCSP&); // not implemented

    ClassDef(AliAnalysisTaskFlowTPCTOFQCSP, 2); //!example of analysis
};

#endif


