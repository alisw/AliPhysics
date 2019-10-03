/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIANALYSISTASKFlowTPCEMCalEP_H
#define ALIANALYSISTASKFlowTPCEMCalEP_H AliPIDResponse.h

class THnSparse;
class TH2F;
class TLorentzVector;

class AliEMCALTrack;
class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliAODMCParticle;
class AliAODMCHeader;
class AliGenEventHeader;
class AliVEvent;
class AliEMCALGeometry;
class AliAnalysisFilter;
class AliESDtrackCuts;
class AliESDtrack;
class AliHFEcontainer;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliCFManager;
class AliPIDResponse;
class AliMultSelection;

#include "AliAnalysisTaskFlowVectorCorrections.h"
#include "AliQnCorrectionsQnVector.h"
#include "AliMCEventHandler.h"
#include "AliMCEvent.h"
#include "AliMCParticle.h"
#include "AliStack.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskFlowTPCEMCalEP : public AliAnalysisTaskSE {
public:
    AliAnalysisTaskFlowTPCEMCalEP();
    AliAnalysisTaskFlowTPCEMCalEP(const char *name);
    virtual ~AliAnalysisTaskFlowTPCEMCalEP();
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    
    void SetHFECuts(AliHFEcuts * const cuts) { fCuts = cuts; };
    AliHFEpid *GetPID() const { return fPID; }
    
    void SetCorrectionPassV0(const char *passV0) {fpassV0 = passV0;};
    void SetCorrectionPassTPC(const char *passTPC){fpassTPC = passTPC;};
    void SetEP(Bool_t UseNewEP) {fUseNewEP = UseNewEP;};
    void SetTender(Bool_t UseTender) {fUseTender = UseTender;};
    void SetPeriod (Double_t period) {fWhichPeriod = period;};
    void SetAssPtCut (Double_t AssPtCut) {fAssPtCut = AssPtCut;};
    void SetITSnCut (Int_t ITSncut) {fITSncut = ITSncut;};
    void SetAssTPCnCut (Int_t AssTPCnCut) {fAssTPCnCut = AssTPCnCut;};
    void SetTPCnCut(Int_t TPCnCut) {fTPCnCut = TPCnCut;};
    void SetSigmaITScut(Double_t SigmaITScut) {fSigmaITScut = SigmaITScut;};
    void SetSigmaTOFcut(Double_t SigmaTOFcut) {fSigmaTOFcut = SigmaTOFcut;};
    void SetSigmaTPCcut(Double_t SigmaTPCcut) {fSigmaTPCcut = SigmaTPCcut;};
    void SetTimeCut(Bool_t TimeCut) {fTimeCut = TimeCut;};
    void SetWeightSyst(Bool_t WeightSyst) {fWeightSyst = WeightSyst;};
    void SetSystTOFcut(Bool_t SystTOFcut) {fSystTOFcut = SystTOFcut;};
    void SetCutM02(Double_t CutM02) {fCutM02 = CutM02;};
    void SetCutM20(Double_t CutM20) {fCutM20 = CutM20;};
    void SetSScut(Bool_t SScut) {fSScut = SScut;};
    void SetAssITSrefitCut(Bool_t AssITSrefitCut) {fAssITSrefitCut = AssITSrefitCut;};
    void SetRejectKinkMother(Bool_t rejectKinkMother = kFALSE) { fRejectKinkMother = rejectKinkMother; };
    void SetPileUpCut(Bool_t EnablePileupRejVZEROTPCout){fEnablePileupRejVZEROTPCout = EnablePileupRejVZEROTPCout;};
    void SelectPhotonicElectron(Int_t iTracks,AliAODTrack *track,Bool_t &fFlagPhotonicElec, Bool_t &fFlagPhotonicElecBCG,Double_t weight, Int_t iCent, Int_t iHijing, Int_t iDecay, Double_t EovP, Double_t fTPCnSigma, Double_t evPlaneV0);
    void GetWeightAndDecay(AliAODMCParticle *particle, Int_t iCent, Int_t &decay, Double_t &weight);
    Bool_t InclElecTrackCuts(AliAODTrack *ietrack);
    Bool_t AssElecTrackCuts(AliAODTrack *aetrack);
    const AliQnCorrectionsQnVector *GetQnVectorFromList( const TList *list,const char *subdetector,const char *expectedstep,const char *altstep);
    void InitParameters();
    
    
    Double_t GetCos2DeltaPhi(Double_t phiA,Double_t phiB)		const;
    Double_t GetDeltaPhi(Double_t phiA,Double_t phiB)	const;
    Double_t GetPi0weight(Double_t mcPi0pT,Int_t iCent) const;
    Double_t GetEtaweight(Double_t mcEtapT,Int_t iCent) const;
    Double_t GetSigmaEMCal(Double_t EoverP, Double_t pt, Int_t iCent) const;
    Double_t GetSigmaEMCalMC(Double_t EoverP, Double_t pt, Int_t iCent) const;
    
    Double_t GetCentWeight(Int_t centbin);
    Double_t GetEPweight(Int_t bin, Int_t iCent);
    Bool_t   RejectEvent(Double_t cent, Int_t centbin);
    Bool_t   RejectEventPlane(Double_t EP, Int_t EPbin);
    Bool_t   IsFromHFdecay(AliAODMCParticle *particle);
    Bool_t   IsFromLMdecay(AliAODMCParticle *particle);
    Bool_t   IsPrimary(AliAODMCParticle *particle);
    
private:
    
    Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);
    
    Int_t                 fWhichPeriod;       // period
    Double_t              fAssPtCut;          // pt cut for associated electron
    Int_t                 fITSncut;             // ITC number of clusters for tagged electrons
    Int_t                 fAssTPCnCut;		// TPC number of clusters for associated electron
    Int_t                 fTPCnCut;         // TPC number of clusters for tagged electron
    Bool_t                fAssITSrefitCut;	// ITS refit for associated electron
    Bool_t                fUseNewEP;          // Use new EP framework
    Bool_t                fUseTender;          // Use tender
    Double_t              fSigmaITScut;
    Double_t              fSigmaTOFcut;
    Double_t              fSigmaTPCcut;
    Bool_t                fTimeCut;
    Bool_t                fWeightSyst;
    Bool_t                fSystTOFcut;
    Double_t              fCutM02;
    Double_t              fCutM20;
    Bool_t                fSScut;
    Bool_t                fEnablePileupRejVZEROTPCout;
    
    // AliESDEvent        	*fESD;	            	 //! ESD object
    AliAODEvent           *fAOD;                  //! AOD object
    AliVEvent             *fVevent;               //! VEvent
    AliPIDResponse        *fpidResponse;          //! PID response
    AliMultSelection      *fMultSelection;
    AliCentrality         *fCentrality;
    
    AliMCEvent            *fMC;                   //! MC object
    AliStack              *fStack;              //! stack
    AliAODMCParticle      *fMCparticle;         //! MC particle
    TClonesArray          *fMCarray;            //! MC array
    AliAODMCHeader        *fMCheader;

    
    TClonesArray  *fTracks_tender;              //Tender tracks
    TClonesArray  *fCaloClusters_tender;        //Tender clusters
    
    
    TList              	*fOutputList;		 //! output list
    
    //  AliESDtrackCuts     	*fTrackCuts;      	 //! ESD track cuts
    //  AliESDtrackCuts     	*fAssTrackCuts;      	 //! ESD track cuts
    AliHFEcuts            *fCuts;                  //! Cut Collection
    
    Bool_t                fIdentifiedAsOutInz;    //! Out Of Range in z
    Bool_t                fPassTheEventCut;       //! Pass The Event Cut
    Bool_t                fRejectKinkMother;      //! Reject Kink Mother
    Bool_t                fIsMC;                  //! flag for MC analysis
    Bool_t                fIsAOD;                //! flag for AOD analysis
    Bool_t                fSetMassConstraint;	//! set mass constraint
    
    Double_t              fVz;                    //! z position of the primary vertex
    AliCFManager          *fCFM;                  //! Correction Framework Manager
    AliHFEpid             *fPID;                  //! PID
    AliHFEpidQAmanager 	*fPIDqa;		//! PID QA manager
    AliAnalysisTaskFlowVectorCorrections *flowQnVectorTask; //! new Qn vector framework
    AliQnCorrectionsManager *fFlowQnVectorMgr; //! new ep
    Double_t              fOpeningAngleCut;	//! openingAngle cut for non-HFE selection
    Double_t              fInvmassCut;		//! invariant mass cut  for non-HFE selection
    Double_t              fChi2Cut;              //! Chi2 cut  for non-HFE selection
    Double_t              fDCAcut;               //! DCA cut  for non-HFE selection
    
    Int_t                 fWhichDecay;		//! which decay
    Double_t              fPi0EtaWeight;		//! weight for pi0/eta in MC with enhanced signal
    
    TH1F                  *fCentAftThr;		//! centrality for events with electron pt  > 8 GeV/c
    TH1F                  *fevPlaneV0AftThr[3];	//! V0 event plane distribution for events with electron pt  > 8 GeV/c
    
    TH2F                  *fTrigger;		//! check trigger
    TH1F                  *fNoEvents;		//! no of events
    TH1F                  *fTrkpt;		//! track pt
    
    TH2F                  *fElecPtULSInvmassCut[3];//! electron pt, ULS pair, invariant mass cut
    TH2F                  *fElecPtLSInvmassCut[3]; //! electron pt, LS pair, invariant mass cut
    TH2F                  *fElecPtInvmassCut[3];   //! electron pt, invariant mass cut
    TH2F                  *fInclElec[3];   	 //! inclusive electron pt
    
    TH2F                  *fInvmassLS[3];	//! Inv mass of LS (e,e)
    TH2F                  *fInvmassULS[3];	//! Inv mass of ULS (e,e)
    TH2F                  *fOpeningAngleLS[3];	//! opening angle for LS pairs
    TH2F                  *fOpeningAngleULS[3];	//! opening angle for ULS pairs
    
    TH1F                  *fTrackPtBefTrkCuts;	//! Track pt before track cuts
    TH1F                  *fTrackPtAftTrkCuts;	//! Track pt after track cuts
    TH1F                  *fChargedParticlePhi;         //! Track phi
    TH1F                  *fElectronPhi;                //! electron phi
    
    TH1F                  *fCent;			//! centrality distribution
    TH1F                  *fCentAftFlt;		//! centrality distribution after centrality flattening
    
    TH1F                  *fevPlaneV0[3];		//! V0 event plane distribution
    TH1F                  *fTPCsubEPres;		//! TPC event plane resolution
    TH2F                  *fEPres[3];		//! histograms for event plane resolution calculation
    THnSparse             *fCorr;			//! correlations
    THnSparse             *fElecMC;		//! electron background from MC
    TH2F                  *feTPCV2[3];		//! CosDeltaPhi vs pt of inclusive eletron (only TPC PID)
    TH2F                  *feV2[3];		//! CosDeltaPhi vs pt of inclusive eletron (TPC + EMCAL PID)
    TH2F                  *fChargPartV2[3];	//! CosDeltaPhi vs pt of charged particle for trigger correction
    TH2F                  *fMtcPartV2[3];		//! CosDeltaPhi vs pt of matched particle for trigger correction
    
    TH2F                  *fPi0Pt[3];		//! primary pi0 pt to compute the weight
    TH2F                  *fEtaPt[3];		//! primary eta pt to compute the weight
    
    TH2F                  *fHistITSnSig[3];      //! ITS sigma vs p
    TH2F                  *fHistTOFnSig[3];      //! TOF sigma vs p
    TH2F                  *fHistTPCnSig[3];      //! TPC sigma vs p
    TH2F                  *fHistTPCnSigITScut[3];      //! TPC sigma vs p (ITS cut)
    TH2F                  *fHistTPCnSigTOFcut[3];      //! TPC sigma vs p (TOF cut)
    TH2F                  *fHistTPCnSigITSTOFcut[3];      //! TPC sigma vs p (ITS+TOF cuts)
    TH2F                  *fHistTPCnSigEop[3];      //! TPC sigma vs E/p
    TH2F                  *fHistTPCnSigEMCalnSig[3];      //! TPC sigma vs EMCal Sig
    TH2F                  *fHistM02sig[3]; //! m02 vs pt for electrons
    TH2F                  *fHistM20sig[3]; //! m20 vs pt for electrons
    TH2F                  *fHistM02backg[3]; //! m02 vs pt for hadrons
    TH2F                  *fHistM20backg[3]; //! m02 vs pt for hadrons
    TH2F                  *fHistM02EoverP[3]; //! m02 vs E/p
    TH2F                  *fHistM20EoverP[3]; //! m20 vs E/p
    
    TH2F                  *fEoverPsignalTPC[3];    //! E/p for electrons (TPC cut)
    TH2F                  *fEoverPsignalTPCM02[3];    //! E/p for electrons (TPCcut + M02 cuts)
    TH2F                  *fEoverPbackg[3];    //! E/p for background
    TH2F                  *fcpV2_3040;          //! charged particle v2 in 30-40%
    TH2F                  *fcpV2_4050;//! charged particle v2 in 30-40%
    
    TString fpassV0;
    TString fpassTPC;
    
    AliAnalysisTaskFlowTPCEMCalEP(const AliAnalysisTaskFlowTPCEMCalEP&); // not implemented
    AliAnalysisTaskFlowTPCEMCalEP& operator=(const AliAnalysisTaskFlowTPCEMCalEP&); // not implemented
    
    ClassDef(AliAnalysisTaskFlowTPCEMCalEP, 1); //!example of analysis
};

#endif
