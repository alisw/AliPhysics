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
class AliSelectNonHFE;
class AliPIDResponse;

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
  void SetOpeningAngleCut(Double_t openingAngle) {fOpeningAngleCut = openingAngle;};
  void SetInvariantMassCut(Double_t invMass) {fInvmassCut = invMass;};
  void SetNonHFEalgorithm(TString nonHFEalgorithm)  {fnonHFEalgorithm = nonHFEalgorithm;};

  AliHFEpid *GetPID() const { return fPID; }
  void SetRejectKinkMother(Bool_t rejectKinkMother = kFALSE) { fRejectKinkMother = rejectKinkMother; };
  void SelectPhotonicElectron(Int_t iTracks,AliESDtrack *track,Bool_t &fFlagPhotonicElec, Bool_t &fFlagPhotonicElecBCG,Double_t weight, Int_t iCent, Int_t iHijing, Int_t iDecay);
  void InitParameters();

  
  Double_t GetCos2DeltaPhi(Double_t phiA,Double_t phiB)		const;
  Double_t GetDeltaPhi(Double_t phiA,Double_t phiB)	const;
  Double_t GetPi0weight(Double_t mcPi0pT,Int_t iCent) const;
  Double_t GetEtaweight(Double_t mcEtapT,Int_t iCent) const;
  Double_t GetSigmaEMCal(Double_t EoverP, Double_t pt, Int_t iCent) const;
  Double_t GetWeight(TParticle *particle, Int_t iCent); 
  Bool_t   IsElectronFromPi0(TParticle *particle,  Double_t &weight, Int_t iCent);
  Bool_t   IsElectronFromEta(TParticle *particle,  Double_t &weight, Int_t iCent);
  Bool_t   IsElectronFromGamma(TParticle *particle, Double_t &weight, Int_t iCent);
  Bool_t   IsPi0EtaFromHFdecay(TParticle *particle);
  Bool_t   IsPi0EtaFromLMdecay(TParticle *particle);
  Bool_t   IsPrimary(TParticle *particle);

 private:
  
  Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);
  

  AliESDEvent        	*fESD;	            	 //! ESD object
  AliAODEvent           *fAOD;                  //! AOD object
  AliVEvent             *fVevent;               //! VEvent
  AliPIDResponse        *fpidResponse;          //! PID response

  AliMCEvent            *fMC;                   //! MC object
  AliStack		*fStack;		//! stack
    
  TList              	*fOutputList;		 //! output list
  
  AliESDtrackCuts     	*fTrackCuts;      	 //! ESD track cuts
  AliHFEcuts 		*fCuts;                  //! Cut Collection
  AliSelectNonHFE       *fNonHFE;               //! Select non heavy flavour electrons

  Bool_t 		fIdentifiedAsOutInz;    //! Out Of Range in z
  Bool_t 		fPassTheEventCut;       //! Pass The Event Cut
  Bool_t 		fRejectKinkMother;      //! Reject Kink Mother
  Bool_t 		fIsMC;                  //! flag for MC analysis   
  Bool_t                fIsAOD;                //! flag for AOD analysis
  Bool_t		fSetMassConstraint;	//! set mass constraint

  Double_t 		fVz;                    //! z position of the primary vertex
  AliCFManager 		*fCFM;                  //! Correction Framework Manager
  AliHFEpid 		*fPID;                  //! PID
  AliHFEpidQAmanager 	*fPIDqa;		//! PID QA manager
  Double_t 		fOpeningAngleCut;	//! openingAngle cut for non-HFE selection
  Double_t	        fInvmassCut;		//! invariant mass cut  for non-HFE selection
  Double_t	        fChi2Cut;              //! Chi2 cut  for non-HFE selection
  Double_t	        fDCAcut;               //! DCA cut  for non-HFE selection
  TString               fnonHFEalgorithm;      //! algorithm to select non-HFE pairs (KF or DCA) 

  TH1F		        *fNoEvents;		//! no of events
  TH1F		        *fTrkpt;		//! track pt
  TH2F		        *fTrkEovPBef;		//! track E/p before HFE pid
  TH2F		        *fTrkEovPAft;		//! track E/p after HFE pid
  TH2F		        *fdEdxBef;		//! track dEdx vs p before HFE pid
  TH2F		        *fdEdxAft;		//! track dEdx vs p after HFE pid
  TH1F		        *fElecPtULSInvmassCut[3];//! electron pt, ULS pair, invariant mass cut
  TH1F		        *fElecPtLSInvmassCut[3]; //! electron pt, LS pair, invariant mass cut
  TH1F		        *fElecPtInvmassCut[3];   //! electron pt, invariant mass cut
  TH1F		        *fInclElec[3];   	 //! inclusive electron pt
  
  TH2F                  *fInvmassLS[3];	//! Inv mass of LS (e,e)
  TH2F                  *fInvmassULS[3];	//! Inv mass of ULS (e,e)
  TH2F                  *fOpeningAngleLS[3];	//! opening angle for LS pairs
  TH2F                  *fOpeningAngleULS[3];	//! opening angle for ULS pairs

  TH1F		        *fTrackPtBefTrkCuts;	//! Track pt before track cuts	
  TH1F		        *fTrackPtAftTrkCuts;	//! Track pt after track cuts
  TH2F		        *fTPCnsigma;		//! TPC n sigma vs p	
  
  TH1F		        *fCent;			//! centrality
  TH1F		        *fevPlaneV0[3];		//! V0 event plane distribution
  TH1F		        *fTPCsubEPres;		//! TPC event plane resolution
  THnSparse	        *fEPres;		//! event plane resolution
  THnSparse	        *fCorr;			//! correlations
  THnSparse	        *fElecMC;		//! electron background from MC 
  TH2F		        *feTPCV2[3];		//! CosDeltaPhi vs pt of inclusive eletron (only TPC PID)
  TH2F		        *feV2[3];		//! CosDeltaPhi vs pt of inclusive eletron (TPC + EMCAL PID)
  TH2F		        *fChargPartV2[3];	//! CosDeltaPhi vs pt of charged particle for trigger correction
  TH2F		        *fMtcPartV2[3];		//! CosDeltaPhi vs pt of matched particle for trigger correction
    
  TH1F		        *fPi0Pt[3];		//! primary pi0 pt to compute the weight
  TH1F		        *fEtaPt[3];		//! primary eta pt to compute the weight
  
  TH1F			*fEoverPsig[3][8][4];	//! E/p distribution for electrons
  TH1F			*fEoverPuls[3][8][4];	//! E/p distribution for electrons from unlike-sign pairs
  TH1F			*fEoverPls[3][8][4];	//! E/p distribution for electrons from like-sign pairs
  TH1F			*fEoverPbcg[3][8][4];	//! E/p distribution for hadrons

  TH1F		        *fDe[6];
  TH1F		        *fD0e[6];
  TH1F		        *fDpluse[6];
  TH1F		        *fDminuse[6];
  
  TH2F		        *fD0_e;
  
  TH1F		        *fTot_pi0e;		//! inclusive electron
  TH1F		        *fPhot_pi0e;		//! ULS pair 
  TH1F		        *fPhotBCG_pi0e;		//! LS pair
  TH1F		        *fTot_etae;		//! inclusive electron
  TH1F		        *fPhot_etae;		//! ULS pair 
  TH1F		        *fPhotBCG_etae;		//! LS pair

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
