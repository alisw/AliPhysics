#ifndef ALIANALYSISTASKFLOWTPCEMCALQCSP_H
#define ALIANALYSISTASKFLOWTPCEMCALQCSP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  Task for Heavy Flavour Electron Flow                              //
//                                                                    //
//  Author: Andrea Dubla (Utrecht University)                         //
//                                                                    //
////////////////////////////////////////////////////////////////////////

class THnSparse;
class TH2F;
class TLorentzVector;

class AliEMCALTrack;
class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliEMCALGeometry;
class AliEMCALRecoUtils;
class AliAnalysisFilter;
class AliESDtrackCuts;
class AliESDtrack;
class AliAODtrack;
class AliHFEcontainer;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliCFManager;
class AliFlowTrackCuts;
class AliFlowTrack;
class AliFlowEvent;
class AliFlowCandidateTrack;
class AliFlowEventSimple;
class AliCentrality;
#include "AliAnalysisTaskSE.h"
//#include "AliCentrality.h"

class AliAnalysisTaskFlowTPCEMCalQCSP : public AliAnalysisTaskSE {

  public:
    AliAnalysisTaskFlowTPCEMCalQCSP();
    AliAnalysisTaskFlowTPCEMCalQCSP(const char *name);
    virtual ~AliAnalysisTaskFlowTPCEMCalQCSP();

    void                                 SetEnableDebugMode() {fDebug = kTRUE; };
    void                                 SetCentralityParameters(Double_t CentralityMin, Double_t CentralityMax, const char* CentralityMethod); //select centrality
    void                                 CheckCentrality(AliAODEvent *event,Bool_t &centralitypass); //to use only events with the correct centrality....
    void                                 SelectPhotonicElectron(Int_t itrack,const AliAODTrack *track, Bool_t &fFlagPhotonicElec);
    void                                 SetInvariantMassCut(Double_t invmass) {fInvmassCut = invmass;};
    void                                 SetTrigger(Int_t trig) {fTrigger = trig;};
    template <typename T> void           PlotVZeroMultiplcities(const T* event) const;
    template <typename T> void           SetNullCuts(T* aod);
    void                                 PrepareFlowEvent(Int_t iMulti, AliFlowEvent *FlowEv) const;
    virtual void                         UserCreateOutputObjects();
    virtual void                         UserExec(Option_t *option);
    virtual void                         Terminate(Option_t *);
    void                                 SetRPCuts(AliFlowTrackCuts *cutsRP) { fCutsRP = cutsRP; }
    void                                 SetHFECuts(AliHFEcuts * const cuts) { fCuts = cuts; };
    void                                 SetIDCuts(Double_t minTPC, Double_t maxTPC, Double_t minEovP, Double_t maxEovP, Double_t minM20, Double_t maxM20, Double_t minM02, Double_t maxM02, Double_t Dispersion);



    AliHFEpid *GetPID() const { return fPID; };

  private:

    Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);

    Bool_t               fDebug; //! enable debug mode
    AliAODEvent 	     *fAOD;			//AOD object
    AliEMCALGeometry     *fGeom; 		// emcal geometry
    TList       	   	 *fOutputList;		//output list
    AliHFEcuts   		 *fCuts;                 //Cut Collection
    Bool_t 			     fIdentifiedAsOutInz;    //Out Of Range in z
    Bool_t 	   	   	     fPassTheEventCut;       //Pass The Event Cut
    Double_t 	         fVz;                    //z position of the primary vertex
    AliCFManager 	   	 *fCFM;                  //!Correction Framework Manager
    AliHFEpid 		     *fPID;                  //PID
    AliHFEpidQAmanager   *fPIDqa;		//! PID QA manager
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
    TH2F		         *fTrkEovPBef;		//track E/p before HFE pid
    TH2F			     *fdEdxBef;		//track dEdx vs p before HFE pid
    TH1F                 *fInclusiveElecPt; // Inclusive elec pt
    TH2F			     *fTPCnsigma;		//TPC n sigma vs p
    TH2F		         *fTPCnsigmaAft;		//TPC n sigma vs p after HFE pid
    TH1F                 *fCentralityPass; // ! QA histogram of events that pass centrality cut
    TH1F                 *fCentralityNoPass; //! QA histogram of events that do not pass centrality cut
    THnSparse            *fSparseElectron;//! Trk matching sparse for v1 clusterizer
    Double_t             *fvalueElectron;//! 6dim sparse fill for trk matching quantities
    TH1F                 *fInvmassLS1; //LS Invmass for all rec par
    TH1F                 *fInvmassULS1;//ULS Invmass for all rec par
    TH1F			     *fPhotoElecPt;		//photonic elec pt
    TH1F			     *fSemiInclElecPt;	//Semi inclusive ele pt
    TH1F                 *fULSElecPt; //ULS elec Pt
    TH1F                 *fLSElecPt;// LS elec pt
	
    Double_t fminTPC;  //ID cuts tpc
    Double_t fmaxTPC;  //ID cuts tpc
    Double_t fminEovP;  //ID cuts eovp
    Double_t fmaxEovP;//ID cuts eovp
    Double_t fminM20;//ID cuts SS
    Double_t fmaxM20;//ID cuts SS
    Double_t fminM02;//ID cuts SS
    Double_t fmaxM02;//ID cuts SS
    Double_t fDispersion;//ID cuts SS
    
	
    AliAnalysisTaskFlowTPCEMCalQCSP(const AliAnalysisTaskFlowTPCEMCalQCSP&); // not implemented
    AliAnalysisTaskFlowTPCEMCalQCSP& operator=(const AliAnalysisTaskFlowTPCEMCalQCSP&); // not implemented

    ClassDef(AliAnalysisTaskFlowTPCEMCalQCSP, 2); //!example of analysis
};

#endif


