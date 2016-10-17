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
class AliVEvent;
class AliEMCALGeometry;
class AliEMCALRecoUtils;
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
    
  void SetEP(Bool_t UseNewEP) {fUseNewEP = UseNewEP;};
  void SetPeriod (Double_t period) {fWhichPeriod = period;};
  void SetAssPtCut (Double_t AssPtCut) {fAssPtCut = AssPtCut;};
  void SetAssTPCnCut (Int_t AssTPCnCut) {fAssTPCnCut = AssTPCnCut;};
  void SetAssITSrefitCut(Bool_t AssITSrefitCut) {fAssITSrefitCut = AssITSrefitCut;};
  void SetRejectKinkMother(Bool_t rejectKinkMother = kFALSE) { fRejectKinkMother = rejectKinkMother; };
  void SelectPhotonicElectron(Int_t iTracks,AliAODTrack *track,Bool_t &fFlagPhotonicElec, Bool_t &fFlagPhotonicElecBCG,Double_t weight, Int_t iCent, Int_t iHijing, Int_t iDecay, Double_t fEMCalnSigma, Double_t fTPCnSigma);
  void GetWeightAndDecay(TParticle *particle, Int_t iCent, Int_t &decay, Double_t &weight);
  Bool_t InclElecTrackCuts(AliAODTrack *ietrack);
  Bool_t AssElecTrackCuts(AliAODTrack *aetrack);
  void InitParameters();

  
  Double_t GetCos2DeltaPhi(Double_t phiA,Double_t phiB)		const;
  Double_t GetDeltaPhi(Double_t phiA,Double_t phiB)	const;
  Double_t GetPi0weight(Double_t mcPi0pT,Int_t iCent) const;
  Double_t GetEtaweight(Double_t mcEtapT,Int_t iCent) const;
  Double_t GetSigmaEMCal(Double_t EoverP, Double_t pt, Int_t iCent) const;
  Double_t GetSigmaEMCalMC(Double_t EoverP, Double_t pt, Int_t iCent) const;
  
  Double_t GetCentWeight(Int_t centbin);
  Double_t GetEPweight(Int_t bin);
  Bool_t   RejectEvent(Double_t cent, Int_t centbin);
  Bool_t   RejectEventPlane(Double_t EP, Int_t EPbin);
  Bool_t   IsFromHFdecay(TParticle *particle);
  Bool_t   IsFromLMdecay(TParticle *particle);
  Bool_t   IsPrimary(TParticle *particle);

 private:
  
  Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);
  
  Int_t                 fWhichPeriod;       // period
  Double_t              fAssPtCut;          // pt cut for associated electron
  Int_t                 fAssTPCnCut;		// TPC number of clusters for associated electron
  Bool_t                fAssITSrefitCut;	// ITS refit for associated electron
  Bool_t                fUseNewEP;          // Use new EP framework
  
 // AliESDEvent        	*fESD;	            	 //! ESD object
  AliAODEvent           *fAOD;                  //! AOD object
  AliVEvent             *fVevent;               //! VEvent
  AliPIDResponse        *fpidResponse;          //! PID response
  AliMultSelection      *fMultSelection;
  AliCentrality         *fCentrality;
    
  AliMCEvent            *fMC;                   //! MC object
  AliStack              *fStack;		//! stack
    
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
  TH2F                  *fTrkEovPBef;		//! track E/p before HFE pid
  TH2F                  *fTrkEovPAft;		//! track E/p after HFE pid
  TH2F                  *fdEdxBef;		//! track dEdx vs p before HFE pid
  TH2F                  *fdEdxAft;		//! track dEdx vs p after HFE pid
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
  TH2F                  *fTPCnsigma;		//! TPC n sigma vs p
  
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
    
  TH1F                  *fPi0Pt[3];		//! primary pi0 pt to compute the weight
  TH1F                  *fEtaPt[3];		//! primary eta pt to compute the weight
  
  TH1F                  *fEoverPsig[3][8][4];	//! E/p distribution for electrons
  TH1F                  *fEoverPuls[3][8][4];	//! E/p distribution for electrons from unlike-sign pairs
  TH1F                  *fEoverPls[3][8][4];	//! E/p distribution for electrons from like-sign pairs
  TH1F                  *fEoverPbcg[3][8][4];	//! E/p distribution for hadrons

  TH1F                  *fDe[6];
  TH1F                  *fD0e[6];
  TH1F                  *fDpluse[6];
  TH1F                  *fDminuse[6];
  
  TH2F                  *fD0_e;
  
  TH1F                  *fTot_pi0e;		//! inclusive electron
  TH1F                  *fPhot_pi0e;		//! ULS pair
  TH1F                  *fPhotBCG_pi0e;		//! LS pair
  TH1F                  *fTot_etae;		//! inclusive electron
  TH1F                  *fPhot_etae;		//! ULS pair
  TH1F                  *fPhotBCG_etae;		//! LS pair

  TH1F                  *fInvMass;		//! Invariant mass of ULS pairs
  TH1F                  *fInvMassBack;		//! Invariant mass if LS pairs
  TH1F                  *fDCA;		        //! DCA of ULS pairs
  TH1F                  *fDCABack;		//! DCA of LS pairs
  TH1F                  *fOpAngle;		//! Opening angle of ULS pairs
  TH1F                  *fOpAngleBack;		//! Opening angle of LS pairs


  AliAnalysisTaskFlowTPCEMCalEP(const AliAnalysisTaskFlowTPCEMCalEP&); // not implemented
  AliAnalysisTaskFlowTPCEMCalEP& operator=(const AliAnalysisTaskFlowTPCEMCalEP&); // not implemented
  
  ClassDef(AliAnalysisTaskFlowTPCEMCalEP, 1); //!example of analysis
};

#endif
