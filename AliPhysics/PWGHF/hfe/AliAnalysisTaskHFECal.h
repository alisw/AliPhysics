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

#ifndef ALIANALYSISTASKHFECAL_H
#define ALIANALYSISTASKHFECAL_H

class THnSparse;
class TH2F;
class TLorentzVector;
class TGraphErrors;

class AliEMCALTrack;
class AliMagF;
class AliESDEvent;
class AliAODEvent;
class AliEMCALGeometry;
class AliAnalysisFilter;
class AliESDtrackCuts;
class AliESDtrack;
class AliHFEcontainer;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliCFManager;

#include "AliAnalysisTaskSE.h"
#include "AliStack.h"

class AliAnalysisTaskHFECal : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskHFECal();
  AliAnalysisTaskHFECal(const char *name);
  virtual ~AliAnalysisTaskHFECal();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  void SetHFECuts(AliHFEcuts * const cuts) { fCuts = cuts; };
  void SetOpeningAngleCut (Double_t openingAngle) {fOpeningAngleCut = openingAngle;};
  void SetMimpTassCut (Double_t MimpTassCut) {fMimpTassCut = MimpTassCut;};
  void SetMimNsigassCut (Double_t MimNsigassCut) {fMimNsigassCut = MimNsigassCut;};
  void SetInvariantMassCut (Double_t invmass) {fInvmassCut = invmass;};
  void SetMassConstraint	(Bool_t MassConstraint)		{ fSetMassConstraint	= MassConstraint; };
  void SetMassWidthCut  	(Bool_t MassWidthCut)		{ fSetMassWidthCut	= MassWidthCut; };
  void SetMassNonlinear  	(Bool_t MassNonlinear)		{ fSetMassNonlinear	= MassNonlinear; };
  void SetMassCalMethod  	(Bool_t KFpart)		{ fSetKFpart = KFpart; };
  void SetQAHist (int qahist) {fqahist = qahist;};
  AliHFEpid *GetPID() const { return fPID; }
  void SetRejectKinkMother(Bool_t rejectKinkMother = kFALSE) { fRejectKinkMother = rejectKinkMother; };
  void SelectPhotonicElectron(Int_t itrack, Double_t cent, AliESDtrack *track, Bool_t &fFlagPhotonic, Bool_t &fFlagConvinat, Double_t nSig, Double_t shower, Double_t ep, Double_t mce, Double_t w, Int_t ibgevent, Bool_t tagpi0, Bool_t tageta, Int_t iCal);
  void FindMother(TParticle* part, int &label, int &pid);
  double GetMCweight(double mcPi0pT);
  double GetMCweightEta(double mcEtapT);
  void FindTriggerClusters();
  double MCEopMeanCorrection(double pTmc, float central);
  double NsigmaCorrection(double tmpeta, float central);
 private:
  
  Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);
  
  AliESDEvent 		*fESD;			//!ESD object
  AliMCEvent 		*fMC;			//!MC object
  AliStack 		*stack;			//!MC object
  AliEMCALGeometry  	*fGeom; 		// emcal geometry 
    
  TList       		*fOutputList;		//! output list
  Int_t                 fqahist;  

  AliESDtrackCuts	*fTrackCuts;		//! ESD track cuts
  AliHFEcuts 		*fCuts;                 //! Cut Collection
  Bool_t 		fIdentifiedAsOutInz;    //Out Of Range in z
  Bool_t 		fPassTheEventCut;       //Pass The Event Cut
  Bool_t 		fRejectKinkMother;      //Reject Kink Mother
  Bool_t                fmcData;
  Double_t 		fVz;                    //z position of the primary vertex
  AliCFManager 		*fCFM;                  //! Correction Framework Manager
  AliHFEpid 		*fPID;                  //! PID
  AliHFEpidQAmanager 	*fPIDqa;		//! PID QA manager
  Double_t 		fOpeningAngleCut;	//openingAngle cut value
  Double_t 		fMimpTassCut;	//openingAngle cut value
  Double_t 		fMimNsigassCut;	//openingAngle cut value
  Double_t		fInvmassCut;		//invariant mass cut value
  Bool_t		 fSetMassConstraint;		// Set mass constraint
  Bool_t		 fSetMassWidthCut;		// Set mass constraint
  Bool_t		 fSetMassNonlinear;		// Set mass constraint
  Bool_t		 fSetKFpart;		// Set mass constraint
 
  int ftriggers[48][60];//!
  int ftriggersCut[48][60];//!
  int ftriggersTime[48][60];//!
 

  TH1F			*fNoEvents;		//! no of events
  THnSparseD		*fEMCAccE;		//! EMC acc
  TH2F  		*hEMCAccE;		//! EMC acc
  TH1F			*fTrkpt;		//! track pt
  TH2F			*fTrkEovPBef;		//! track E/p before HFE pid
  TH2F			*fTrkEovPAft;		//! track E/p after HFE pid
  TH2F			*fdEdxBef;		//! track dEdx vs p before HFE pid
  TH2F			*fdEdxAft;		//! track dEdx vs p after HFE pid
  TH2F			*fIncpT;		//! HFE pid electron vs centrality
  TH2F			*fIncpTM20;		//! HFE pid electron vs centrality
  THnSparseD		*fInvmassLS;		//! Inv mass of LS (e,e)
  THnSparseD		*fInvmassULS;		//! Inv mass of ULS (e,e)
  THnSparseD		*fInvmassLSmc;		//! Inv mass of LS (e,e)
  THnSparseD		*fInvmassULSmc;		//! Inv mass of ULS (e,e)
  TH2D		*fInvmassLSreco;		//! Inv mass of LS (e,e)
  TH2D		*fInvmassULSreco;		//! Inv mass of ULS (e,e)
  TH2D		*fInvmassLSmc0;		//! Inv mass of ULS (e,e)
  TH2D		*fInvmassLSmc1;		//! Inv mass of ULS (e,e)
  TH2D		*fInvmassLSmc2;		//! Inv mass of ULS (e,e)
  TH2D		*fInvmassLSmc3;		//! Inv mass of ULS (e,e)
  TH2D		*fInvmassULSmc0;		//! Inv mass of ULS (e,e)
  TH2D		*fInvmassULSmc1;		//! Inv mass of ULS (e,e)
  TH2D		*fInvmassULSmc2;		//! Inv mass of ULS (e,e)
  TH2D		*fInvmassULSmc3;		//! Inv mass of ULS (e,e)
  TH1F			*fOpeningAngleLS;	//! opening angle for LS pairs
  TH1F			*fOpeningAngleULS;	//! opening angle for ULS pairs
  TH1F			*fPhotoElecPt;		//! photonic elec pt 
  TH2F			*fPhoElecPt;	        //! Pho inclusive ele pt
  TH2F			*fPhoElecPtM20;	        //! Pho inclusive ele pt
  TH2F			*fPhoElecPtM20Mass;	        //! Pho inclusive ele pt
  TH2F			*fSameElecPt;	        //! Same inclusive ele pt
  TH2F			*fSameElecPtM20;	        //! Same inclusive ele pt
  TH2F			*fSameElecPtM20Mass;	        //! Same inclusive ele pt
  TH2F			*fSemiElecPtM20;	        //! Same inclusive ele pt

  TH1F			*fTrackPtBefTrkCuts;	//! Track pt before track cuts	
  TH1F			*fTrackPtAftTrkCuts;	//! Track pt after track cuts
  TH2F			*fTPCnsigma;		//! TPC n sigma vs p	
  
  TH1F			*fCent;			//! centrality
  THnSparseD		*fEleInfo;		//! EMC acc
  THnSparseD		*fElenSigma;		//! EMC acc
  /*
  //<---- trigger info
  TH1F	      *fClsEBftTrigCut;	//Cluster E before trigger selection
  TH1F        *fClsEAftTrigCut;	//Cluster E after trigger selection
  TH1F	      *fClsEAftTrigCut1;	//Cluster E after trigger selection
  TH1F	      *fClsEAftTrigCut2;	//Cluster E after trigger selection
  TH1F	      *fClsEAftTrigCut3;	//Cluster E after trigger selection
  TH1F	      *fClsEAftTrigCut4;	//Cluster E after trigger selection
  TH2F        *fClsETime; //ClsE vs time distribution
  TH2F        *fClsETime1; //ClsE vs time distribution
  TH1F        *fTrigTimes;// trigger time
  TH2F        *fCellCheck;// trigger time
  */
  //<------ MC
  TH2F                  *fInputHFEMC;
  TH2F                  *fInputAlle;
  TH2F			*fIncpTMChfe;		//! MC HFE pid electron vs centrality
  TH2F			*fIncpTMChfeAll;		//! MC HFE pid electron vs centrality
  TH2F			*fIncpTMCM20hfe;	//! MC HFE pid electron vs centrality
  TH2F			*fIncpTMCM20hfeAll;	//! MC HFE pid electron vs centrality
  TH2F			*fIncpTMCM20hfeCheck;	//! MC HFE pid electron vs centrality
 THnSparseD		*fInputHFEMC_weight;		//! MC HFE pid electron vs centrality
 THnSparseD		*fIncpTMCM20hfeCheck_weight;		//! MC HFE pid electron vs centrality
 THnSparseD		*fIncpTMCpho;		//! MC HFE pid electron vs centrality
 THnSparseD 		*fIncpTMCM20pho;	//! MC HFE pid electron vs centrality
 THnSparseD 		*fPhoElecPtMC;	        //! Pho inclusive ele pt
 THnSparseD 		*fPhoElecPtMCM20;	        //! Pho inclusive ele pt
 TH2D 		        *fPhoElecPtMCM20Mass;	        //! Pho inclusive ele pt
 THnSparseD 		*fSameElecPtMC;	        //! Same inclusive ele pt
 THnSparseD 		*fSameElecPtMCM20;	        //! Same inclusive ele pt
 TH2D 		        *fSameElecPtMCM20Mass;	        //! Same inclusive ele pt
 THnSparseD 		*fIncpTMCM20pho_pi0e;	//! MC HFE pid electron vs centrality
 THnSparseD 		*fPhoElecPtMCM20_pi0e;	        //! Pho inclusive ele pt
 THnSparseD 		*fSameElecPtMCM20_pi0e;	        //! Same inclusive ele pt
 THnSparseD 		*fIncpTMCM20pho_pi0Dalitz;	//! MC HFE pid electron vs centrality
 THnSparseD 		*fPhoElecPtMCM20_pi0Dalitz;	        //! Pho inclusive ele pt
 THnSparseD 		*fSameElecPtMCM20_pi0Dalitz;	        //! Same inclusive ele pt
 THnSparseD 		*fIncpTMCM20pho_eta;	//! MC HFE pid electron vs centrality
 THnSparseD 		*fPhoElecPtMCM20_eta;	        //! Pho inclusive ele pt
 THnSparseD 		*fSameElecPtMCM20_eta;	        //! Same inclusive ele pt
 THnSparseD 		*fIncpTMCpho_pi0e_TPC;	//! MC HFE pid electron vs centrality
 THnSparseD 		*fPhoElecPtMC_pi0e_TPC;	        //! Pho inclusive ele pt
 THnSparseD 		*fSameElecPtMC_pi0e_TPC;	        //! Same inclusive ele pt
 THnSparseD 		*fIncpTMCpho_eta_TPC;	//! MC HFE pid electron vs centrality
 THnSparseD 		*fPhoElecPtMC_eta_TPC;	        //! Pho inclusive ele pt
 THnSparseD 		*fSameElecPtMC_eta_TPC;	        //! Same inclusive ele pt
 TH1D                   *CheckNclust;  
 TH1D                   *CheckNits;  
 TH2D                   *CheckDCA;  
 THnSparseD             *Hpi0pTcheck; 
 THnSparseD             *HETApTcheck; 
 TH2D                   *HphopTcheck; 
 TH2D                   *HDpTcheck; 
 TH2D                   *HBpTcheck; 
 TH1D                   *fpTCheck; 
 TH2D                   *fMomDtoE; 
 TH2D                   *fLabelCheck;
 TH2D                   *fgeoFake;
 TH2D                   *fFakeTrk0;
 TH2D                   *fFakeTrk1;
 TH2D                   *ftimingEle;
 TH2D                   *fIncMaxE;
 TH2D                   *fIncReco;
 TH2D                   *fPhoReco;
 TH2D                   *fSamReco; 
 TH2D                   *fIncRecoMaxE;
 TH2D                   *fPhoRecoMaxE;
 TH2D                   *fSamRecoMaxE; 
 TH2D                   *fPhoVertexReco_TPC;
 TH2D                   *fPhoVertexReco_TPC_Invmass;
 TH2D                   *fPhoVertexReco_EMCal;
 TH2D                   *fPhoVertexReco_Invmass;
 TH2D                   *fPhoVertexReco_step0;
 TH2D                   *fPhoVertexReco_step1;
 TH1D                   *fMatchV0_0;
 TH1D                   *fMatchV0_1;
 TH1D                   *fMatchMC_0;
 TH1D                   *fMatchMC_1;
 TH2D                   *fpair;
 TH2D                   *fFakeRejection0;
 TH2D                   *fFakeRejection1;
 TH2D                   *fFakeRejection2;
 TH2D                   *EopFake;
 TH2D                   *EopTrue;
 TH2D                   *MatchFake;
 TH2D                   *MatchTrue;
 TH2D                   *MatchTrCheck;
 TH2D                   *MatchTrEop;

 //<----- correction
 TGraphErrors           *fnSigEtaCorr[7];


  AliAnalysisTaskHFECal(const AliAnalysisTaskHFECal&); // not implemented
  AliAnalysisTaskHFECal& operator=(const AliAnalysisTaskHFECal&); // not implemented
  
  ClassDef(AliAnalysisTaskHFECal, 1); //!example of analysis
};

#endif


