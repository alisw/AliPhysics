#ifndef ALIANALYSISTASKFLOWTPCEMCALQCSP_H
#define ALIANALYSISTASKFLOWTPCEMCALQCSP_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  Task for Heavy Flavour Electron Flow TPC plus EMCal               //
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
class AliAnalysisFilter;
class AliESDtrackCuts;
class AliESDtrack;
class AliAODtrack;
class AliHFEcontainer;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidTPC;
class AliHFEpidQAmanager;
class AliCFManager;
class AliFlowTrackCuts;
class AliFlowTrack;
class AliFlowEvent;
class AliFlowCandidateTrack;
class AliFlowEventSimple;
class AliCentrality;
class AliSelectNonHFE;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskFlowTPCEMCalQCSP : public AliAnalysisTaskSE {
    
public:
    AliAnalysisTaskFlowTPCEMCalQCSP();
    AliAnalysisTaskFlowTPCEMCalQCSP(const char *name);
    virtual ~AliAnalysisTaskFlowTPCEMCalQCSP();
    
    void                                 SetEnableDebugMode() {fDebug = kTRUE; };
    void                                 SetCentralityParameters(Double_t CentralityMin, Double_t CentralityMax, const char* CentralityMethod); //select centrality
    void                                 CheckCentrality(AliAODEvent *event,Bool_t &centralitypass); //to use only events with the correct centrality....
    void                                 SelectPhotonicElectron(Int_t itrack,const AliAODTrack *track, Double_t eopinc, Double_t evPlAnglV0, Bool_t &fFlagPhotonicElec, Double_t weightEPflat, Double_t multev);
    void                                 SetInvariantMassCut(Double_t invmass) {fInvmassCut = invmass;};
    void                                 SetpTCuttrack(Double_t ptcut) {fpTCut = ptcut;};
    void                                 SetTrigger(Int_t trig) {fTrigger = trig;};
    void                                 SetAssoTPCCluster(Int_t tpc_clust) {fAssoTPCCluster = tpc_clust;};
    void                                 SetMultCorrelationCut(Bool_t multcut) {fMultCut = multcut;};
    void                                 SetAssoITSRefit(Bool_t itsref) {fAssoITSRefit = itsref;};
    void                                 SetFlowSideBands(Bool_t sidebandsflow){fSideBandsFlow = sidebandsflow;}
    void                                 SelectPhotonicElectronMethod(Bool_t dca){fDCA = dca;}
    void                                 Setphiminuspsi(Bool_t phipsi){fPhiminusPsi = phipsi;}
    void                                 SetPurity(Bool_t Purityel){fpurity = Purityel;}
    template <typename T> void           PlotVZeroMultiplcities(const T* event) const;
    template <typename T> void           SetNullCuts(T* aod);
    void                                 PrepareFlowEvent(Int_t iMulti, AliFlowEvent *FlowEv) const;
    virtual void                         UserCreateOutputObjects();
    virtual void                         UserExec(Option_t *option);
    virtual void                         Terminate(Option_t *);
    void                                 SetRPCuts(AliFlowTrackCuts *cutsRP) { fCutsRP = cutsRP; }
    void                                 SetHFECuts(AliHFEcuts * const cuts) { fCuts = cuts; };
    void                                 SetIDCuts(Double_t minTPC, Double_t maxTPC, Double_t minEovP, Double_t maxEovP, Double_t minM20, Double_t maxM20, Double_t minM02, Double_t maxM02, Double_t Dispersion);
    void                                 SetOpeningAngleflag(Bool_t opang){fOP_angle = opang;};
    void                                 SetOpeningAngleCut(Double_t opanglecut) {fOpeningAngleCut = opanglecut;};
    void                                 SetHistoForCentralityFlattening(TH1F *h,Double_t minCentr,Double_t maxCentr,Double_t centrRef=0.,Int_t switchTRand=0);
    Bool_t                               IsEventSelectedForCentrFlattening(Float_t centvalue);
    
    void                                 SetHistoForEPFlattWeights(TH1D *h);
    Double_t                             GiveMeWeight(Double_t EP);
    void                                 SetEPWeight(Bool_t epw){EPweights = epw;};
    
    void                                 SetTPCPID(AliHFEpidTPC *pidcorr){ftpcpid = pidcorr;};
    void                                 SetMultCorrectionTheo(Bool_t mulcorr){multCorrection = mulcorr;}
    
    void                                 SetPtMinAssoCut(Double_t ptminimumasso) {fptminAsso = ptminimumasso;};
    
    
    void                                 SetEtaMinPos(Double_t minetapos){fEtaMinimumPositive = minetapos;}
    void                                 SetEtaMinNeg(Double_t minetaneg){fEtaMinimumNegative = minetaneg;}
    
    
    void                                 SetHistoForCentralityFlattening_Bis(TH1F *h,Double_t minCentr,Double_t maxCentr,Double_t centrRef);
    Bool_t                               IsEventSelectedForCentrFlattening_Bis(Double_t centvalue);
    
    void                                 SetCentralityMine(Bool_t CentFlatMine){fCentFlatMine = CentFlatMine;}

    
    AliHFEpid *GetPID() const { return fPID; };
    
private:
    
    Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);
    
    Bool_t               fDebug; //! enable debug mode
    AliAODEvent          *fAOD;                     //AOD object
    AliVEvent            *fVevent;			//ESD object
    AliEMCALGeometry     *fGeom;                // emcal geometry
    TList                *fOutputList;          //output list
    AliHFEcuts           *fCuts;                 //Cut Collection
    Bool_t               fIdentifiedAsOutInz;    //Out Of Range in z
    Bool_t               fPassTheEventCut;       //Pass The Event Cut
    AliCFManager         *fCFM;                  //!Correction Framework Manager
    AliHFEpid            *fPID;                  //PID
    AliHFEpidTPC         *ftpcpid;              // for TPC mult/eta correction, are done in the HFE class
    AliHFEpidQAmanager   *fPIDqa;               //! PID QA manager
    AliFlowTrackCuts     *fCutsRP; // track cuts for reference particles
    AliFlowTrackCuts     *fNullCuts; // dummy cuts for flow event tracks
    AliFlowEvent         *fFlowEvent; //! flow events Inclusive e
    const char           *fkCentralityMethod; // method used to determine centrality (V0 by default)
    Double_t             fCentrality; // event centrality for QA
    Double_t             fCentralityMin; // lower bound of cenrality bin
    Double_t             fCentralityMax; // upper bound of centrality bin
    Double_t             fInvmassCut;           //invariant mass cut value
    Double_t             fpTCut;    //pt cut value
    Int_t                fTrigger;          //invariant mass cut value
    TH1F                 *fPhi; //! QA plot of azimuthal distribution of tracks used for event plane estimation
    TH1F                 *fEta; //! QA plot of eta distribution of tracks used for event plane estimation
    TH1F                 *fVZEROA; //! QA plot vzeroa multiplicity (all tracks in event)
    TH1F                 *fVZEROC; //! QA plot vzeroc multiplicity (all tracks in event)
    TH1F                 *fTPCM; //! QA plot TPC multiplicity (tracks used for event plane estimation)
    TH1F                 *fNoEvents;            //!no of events
    TH2F                 *fTrkEovPBef;          //!track E/p before HFE pid
    //   TH2F                           *fdEdxBef;         //track dEdx vs p before HFE pid
    TH1F                 *fInclusiveElecPt; //! Inclusive elec pt
    TH2F                 *fTPCnsigma;               //!TPC n sigma vs p
    TH2F                 *fTPCnsigmaAft;                //!TPC n sigma vs p after HFE pid
    TH1F                 *fCentralityPass; // ! QA histogram of events that pass centrality cut
    TH1F                 *fCentralityNoPass; //! QA histogram of events that do not pass centrality cut
    TH1F                 *fInvmassLS1; //!LS Invmass for all rec par
    TH1F                 *fInvmassULS1;//!ULS Invmass for all rec par
    TH1F                 *fPhotoElecPt;             //!photonic elec pt
    TH1F                 *fSemiInclElecPt;  //!Semi inclusive ele pt
    TH1F                 *fULSElecPt; //!ULS elec Pt
    TH1F                 *fLSElecPt;//! LS elec pt
    Double_t             fminTPC;  //ID cuts tpc
    Double_t             fmaxTPC;  //ID cuts tpc
    Double_t             fminEovP;  //ID cuts eovp
    Double_t             fmaxEovP;//ID cuts eovp
    Double_t             fminM20;//ID cuts SS
    Double_t             fmaxM20;//ID cuts SS
    Double_t             fminM02;//ID cuts SS
    Double_t             fmaxM02;//ID cuts SS
    Double_t             fDispersion;//ID cuts SS
    TH2F                 *fMultCorAfterCuts; //! QA profile global and tpc multiplicity after outlier cut
    TH2F                 *fMultvsCentr; //! QA profile of centralty vs multiplicity
    TProfile             *fSubEventDPhiv2;//!lll
    TH1D                 *EPVzA;//!v0aep
    TH1D                 *EPVzC;//!v0cep
    TH1D                 *EPTPC;//!tpcep
    THnSparseF           *fV2Phi;//! v2 analysis of EP-V0
    THnSparseD           *fSparseElectronHadron;//! Trk matching sparse for v1 clusterizer
    TH1D                 *fvertex;//!huge ThNsparse
    TH2F                 *fMultCorBeforeCuts; //! QA profile global and tpc multiplicity after outlier cut
    Bool_t               fSideBandsFlow;//flow from side bands for contamination
    Bool_t               fPhiminusPsi;//flow from phi minus psi method
    AliFlowEvent         *fFlowEventCont; // flow events for elect Contamination
    Bool_t               fpurity; //for purity evaluation
    THnSparseD           *fSparseElectronpurity;//! Trk matching sparse for v1 clusterizer
    TH1F                 *fOpeningAngleLS;  //!opening angle for LS pairs
    TH1F                 *fOpeningAngleULS; //!opening angle for ULS pairs
    AliSelectNonHFE      *fNonHFE;//new elienos stuff
    Bool_t               fDCA;//selection PHelectron
    Double_t             fOpeningAngleCut; //openingAngle cut value
    Bool_t               fOP_angle; //to shitch on and off the op_angle cut
    Int_t                fAssoTPCCluster;//asso tpc cluster
    Bool_t               fAssoITSRefit;//asso its refit
    Bool_t               fMultCut;//for mult correlationcut
    TH2F                 *fMultCorAfterCentrBeforeCuts; //! QA profile global and tpc multiplicity after outlier cut
    TH2F                 *fMultCorAfterVZTRKComp;//!after cent comp
    TH1F                 *fCentralityBeforePileup;//!cent chneck
    TH1F                 *fCentralityAfterVZTRK;//!cent chneck2
    TH1F                 *fCentralityAfterCorrCut;//!cent chneck2
    TH2F                 *fMultCorAfterCorrCut;//!after cent comp
    TH1D                 *EPVz;//!v0cep
    TH1D                 *EPTPCp;//!tpcep
    TH1D                 *EPTPCn;//!tpcep
    TProfile             *fSubEventDPhiv2new;
    THnSparseF           *fV2Phivzerotot;//! v2 analysis of EP-V0
    TH1F                 *fHistCentrDistr;//-> isto for Centr Flat
    TH1F                 *fCentralityNoPassForFlattening; //! QA histogram of events that do not pass centrality cut for flattening
    TH1F                 *fInvmassLS1highpt; //!LS Invmass for all rec par high pt
    TH1F                 *fInvmassULS1highpt;//!ULS Invmass for all rec par high pt
    THnSparseF           *fSparsephipsiULS;//!ssss
    THnSparseF           *fSparsephipsiLS;//!ssss
    THnSparseF           *fSparseMassULS;//!ssss
    THnSparseF           *fSparseMassLS;//!ssssss
    
    TH1D                 *fHistEPDistrWeight;// isto for Centr Flat
    Bool_t               EPweights;//for mult correlationcut
    TH1D                 *EPVzAftW;//!v0cep
    
    Bool_t                multCorrection;//Flag to activate mult/etacorrection
    Double_t              fptminAsso;//minassopt
    
    Double_t              fEtaMinimumPositive;//for reso
    Double_t              fEtaMinimumNegative;//for reso
    
    TH1F                 *fCentralityAll;//!centall
    Bool_t               fCentFlatMine; //for purity evaluation

    
    AliAnalysisTaskFlowTPCEMCalQCSP(const AliAnalysisTaskFlowTPCEMCalQCSP&); // not implemented
    AliAnalysisTaskFlowTPCEMCalQCSP& operator=(const AliAnalysisTaskFlowTPCEMCalQCSP&); // not implemented
    
    ClassDef(AliAnalysisTaskFlowTPCEMCalQCSP, 2); //!example of analysis
};

#endif

