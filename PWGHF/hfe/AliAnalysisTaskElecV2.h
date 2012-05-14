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

#ifndef ALIANALYSISTASKELECV2_H
#define ALIANALYSISTASKELECV2_H

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

class AliAnalysisTaskElecV2 : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskElecV2();
  AliAnalysisTaskElecV2(const char *name);
  virtual ~AliAnalysisTaskElecV2();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void SetHFECuts(AliHFEcuts * const cuts) { fCuts = cuts; };
  void SetOpeningAngleCut (Double_t openingAngle) {fOpeningAngleCut = openingAngle;};
  void SetInvariantMassCut (Double_t invmass) {fInvmassCut = invmass;};
  AliHFEpid *GetPID() const { return fPID; }
  void SetRejectKinkMother(Bool_t rejectKinkMother = kFALSE) { fRejectKinkMother = rejectKinkMother; };
  void SelectPhotonicElectron(Int_t itrack, AliESDtrack *track, Bool_t &fFlagPhotonicElec);
  Double_t GetCos2DeltaPhi(Double_t phiA,Double_t phiB)		const;
  Double_t GetDeltaPhi(Double_t phiA,Double_t phiB)	const;
 private:
  
  Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);
  
  AliESDEvent 		*fESD;			//!ESD object
    
  TList       		*fOutputList;		//! output list
  
  AliESDtrackCuts	*fTrackCuts;		//! ESD track cuts
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
  TH2F			*fdEdxAft;		//! track dEdx vs p after HFE pid
  TH1F			*fInvmassLS;		//! Inv mass of LS (e,e)
  TH1F			*fInvmassULS;		//! Inv mass of ULS (e,e)
  TH1F			*fOpeningAngleLS;	//! opening angle for LS pairs
  TH1F			*fOpeningAngleULS;	//! opening angle for ULS pairs
  TH1F			*fPhotoElecPt;		//! photonic elec pt 
  TH1F			*fSemiInclElecPt;	//! Semi inclusive ele pt
  
  TH1F			*fTrackPtBefTrkCuts;	//! Track pt before track cuts	
  TH1F			*fTrackPtAftTrkCuts;	//! Track pt after track cuts
  TH2F			*fTPCnsigma;		//! TPC n sigma vs p	
  
  TH1F			*fCent;			//! centrality
  TH1F			*fevPlaneV0A;		//! V0A event plane distribution
  TH1F			*fevPlaneV0C;		//! V0C event plane distribution
  TH1F			*fevPlaneTPC;		//! TPC event plane distribution
  TH2F			*fTPCsubEPres;		//! TPC event plane resolution
  THnSparse		*fEPres;		//! event plane resolution
  THnSparse		*fCorr;			//! correlations
  THnSparse		*feTPCV2;		//! inclusive eletron v2 (only TPC PID)
  THnSparse		*feV2;			//! inclusive eletron v2 (TPC + EMCAL PID)
  THnSparse		*fphoteV2;		//! photonic electron v2 (TPC + EMCAL PID)
  THnSparse		*fChargPartV2;		//! charged particle v2
  
  AliAnalysisTaskElecV2(const AliAnalysisTaskElecV2&); // not implemented
  AliAnalysisTaskElecV2& operator=(const AliAnalysisTaskElecV2&); // not implemented
  
  ClassDef(AliAnalysisTaskElecV2, 1); //!example of analysis
};

#endif


