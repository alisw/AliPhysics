#ifndef ALIANALYSISTASKELECHADRONCORREL_H
#define ALIANALYSISTASKELECHADRONCORREL_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//  Task for Heavy Flavour Electron-Hadron DeltaPhi Correlation       //
//                                                                    //
//  Author: Deepa Thomas (Utrecht University)                         //
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
class AliHFEcontainer;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliCFManager;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskElecHadronCorrel : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskElecHadronCorrel();
  AliAnalysisTaskElecHadronCorrel(const char *name);
  virtual ~AliAnalysisTaskElecHadronCorrel();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void SetHFECuts(AliHFEcuts * const cuts) { fCuts = cuts; };
  void SetOpeningAngleCut (Double_t openingAngle) {fOpeningAngleCut = openingAngle;};
  void SetInvariantMassCut (Double_t invmass) {fInvmassCut = invmass;};
  AliHFEpid *GetPID() const { return fPID; }
  void SetRejectKinkMother(Bool_t rejectKinkMother = kFALSE) { fRejectKinkMother = rejectKinkMother; };
  void SelectPhotonicElectron(Int_t itrack, AliESDtrack *track, Bool_t &fFlagPhotonicElec);
  void ElectronHadCorrel(Int_t itrack, AliESDtrack *track, TH2F *DphiPt);	
  Bool_t IsLEDEvent() const;	
 private:
  
  Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);
  
  AliESDEvent 		*fESD;			//!ESD object
  AliEMCALGeometry	*fGeom; 		//! emcal geometry 
  
  TList       		*fOutputList;		//! output list
  
  AliESDtrackCuts	*fTrackCuts1;		//! ESD track cuts
  AliESDtrackCuts	*fTrackCuts2;		//! ESD track cuts
  AliHFEcuts 		*fCuts;                 //! Cut Collection
  Bool_t 		fIdentifiedAsOutInz;    //Out Of Range in z
  Bool_t 		fPassTheEventCut;       //Pass The Event Cut
  Bool_t 		fRejectKinkMother;      //Reject Kink Mother
  Double_t 		fVz;                    //z position of the primary vertex
  AliCFManager 		*fCFM;                  //! Correction Framework Manager
  AliHFEpid 		*fPID;                  //! PID
  AliHFEpidQAmanager 	*fPIDqa;		//! PID QA manager
  Double_t 		fOpeningAngleCut;	//openingAngle cut value
  Double_t		fInvmassCut;		//invariant mass cut value
  
  TH1F			*fNoEvents;		//! no of events
  TH1F			*fTrkpt;		//! track pt
  TH2F			*fTrkEovPBef;		//! track E/p before HFE pid
  TH2F			*fTrkEovPAft;		//! track E/p after HFE pid
  TH2F			*fdEdxBef;		//! track dEdx vs p before HFE pid
  TH2F			*fdEdxAft;		//! track dEdx vs p before HFE pid
  TH1F			*fInvmassLS;		//! Inv mass of LS (e,e)
  TH1F			*fInvmassULS;		//! Inv mass of ULS (e,e)
  TH1F			*fOpeningAngleLS;	//! opening angle for LS pairs
  TH1F			*fOpeningAngleULS;	//! opening angle for ULS pairs
  TH2F			*fSemiIncElecDphi;  	//! Semi Inclusive elec - had DPhi
  TH2F			*fPhotElecDphi;  	//! Photon elec - had DPhi
  TH2F			*fInclusiveElecDphi;  	//! Inclusive elec - had DPhi
  TH2F			*fDphiMassHigh;		//! Dphi - LS+ULS, mass>0.01
  TH2F			*fDphiULSMassLow;	//! Dphi - ULS, mass<0.01
  TH2F			*fDphiLSMassLow;	//! Dphi - LS, mass<0.01
  TH1F			*fPhotoElecPt;		//! photonic elec pt 
  TH1F			*fSemiInclElecPt;	//! Semi inclusive ele pt
  
  TH1F			*fTrackPtBefTrkCuts;	//! Track pt before track cuts	
  TH1F			*fTrackPtAftTrkCuts;	//! Track pt after track cuts
  TH2F			*fTPCnsigma;		//! TPC n sigma vs p	
  
  AliAnalysisTaskElecHadronCorrel(const AliAnalysisTaskElecHadronCorrel&); // not implemented
  AliAnalysisTaskElecHadronCorrel& operator=(const AliAnalysisTaskElecHadronCorrel&); // not implemented
  
  ClassDef(AliAnalysisTaskElecHadronCorrel, 1); //!example of analysis
};

#endif


