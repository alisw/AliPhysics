#ifndef AliAnalysisTaskEMCalHFEpA_cxx
#define AliAnalysisTaskEMCalHFEpA_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

	////////////////////////////////////////////////////////////////////////
	//                                                                    //
	//      Task for Heavy-flavour electron analysis in pPb collisions    //
	//      (+ Electron-Hadron Jetlike Azimuthal Correlation)             //
	//																	  //
	//		version: March 18, 2014.								      //
	//                                                                    //
	//	    Authors 							                          //
	//		Elienos Pereira de Oliveira Filho (epereira@cern.ch)	      //
	//		Cristiane Jahnke		(cristiane.jahnke@cern.ch)		      //
	//                                                                    //
	////////////////////////////////////////////////////////////////////////

class TH1F;
class TH2F;
class AliESDEvent;
class AliESDtrackCuts;
class AliESDtrack;
class AliHFEcontainer;
class AliHFEcuts;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliCFManager;
class AliPIDResponse;
class AliCentrality;
class AliAODEvent;
class AliVEvent;
class AliAODMCHeader;
class AliSelectNonHFE;
class AliEventPoolManager;
class AliEventPool;
class TObjArray;
	//Lucile
class AliCaloTrackAODReader;
class AliCaloTrackReader;


	//______________________________________________________________________
	//Library
#include "AliAnalysisTaskSE.h"
#include "AliHFEpid.h"
#include "AliLog.h"
	//______________________________________________________________________

	//______________________________________________________________________
class AliAnalysisTaskEMCalHFEpA : public AliAnalysisTaskSE 
{
		//______________________________________________________________________
public:
	AliAnalysisTaskEMCalHFEpA();
	AliAnalysisTaskEMCalHFEpA(const char *name);
	virtual ~AliAnalysisTaskEMCalHFEpA();
	
	virtual void   UserCreateOutputObjects();
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t *);
	
		//Setters
	void SetAssHadronPtRange(Double_t AssHadronPtMin, Double_t AssHadronPtMax) {fAssHadronPtMin = AssHadronPtMin; fAssHadronPtMax = AssHadronPtMax; };
	void SetHFECuts(AliHFEcuts * const cuts) {fCuts = cuts;};
	void SetRejectKinkMother(Bool_t rejectKinkMother = kFALSE) {fRejectKinkMother = rejectKinkMother;};
	void SetCorrelationAnalysis(Bool_t CorrelationFlag=kTRUE) {fCorrelationFlag = CorrelationFlag;};
	void SetMCanalysis() {fIsMC = kTRUE;};
	void SetCentrality(Double_t CentralityMin, Double_t CentralityMax) { fCentralityMin = CentralityMin; fCentralityMax = CentralityMax; fHasCentralitySelection = kTRUE; };
	void SetAODanalysis(Bool_t IsAOD) {fIsAOD = IsAOD;};
	void SetEventMixing(Bool_t EventMixingFlag) { fEventMixingFlag = EventMixingFlag;};
	void SetNonHFEmassCut(Double_t MassCut) { fMassCut = MassCut; fMassCutFlag = kTRUE;};
	void SetEtaCut(Double_t EtaCutMin,Double_t EtaCutMax ) { fEtaCutMin = EtaCutMin; fEtaCutMax = EtaCutMax; };
	
	void SetdPhidEtaCut(Double_t dPhiCut, Double_t dEtaCut ) { fdPhiCut = dPhiCut;fdEtaCut = dEtaCut ;};
	
	void SetEoverPCut(Double_t EoverPCutMin,Double_t EoverPCutMax ) { fEoverPCutMin = EoverPCutMin; fEoverPCutMax = EoverPCutMax; };
	
	void SetM02Cut(Double_t M02CutMin,Double_t M02CutMax ) { fM02CutMin = M02CutMin; fM02CutMax = M02CutMax; };
	void SetM20Cut(Double_t M20CutMin,Double_t M20CutMax ) { fM20CutMin = M20CutMin; fM20CutMax = M20CutMax; };
	
	
	void SetNonHFEangleCut(Double_t AngleCut) { fAngleCut = AngleCut; fAngleCutFlag = kTRUE;};
	void SetNonHFEchi2Cut(Double_t Chi2Cut) { fChi2Cut = Chi2Cut; fChi2CutFlag = kTRUE;};
	void SetNonHFEdcaCut(Double_t DCAcut) { fDCAcut = DCAcut; fDCAcutFlag = kTRUE;};
	void SetUseEMCal() { fUseEMCal=kTRUE;};
	void SetUseTrigger() { fUseTrigger=kTRUE;};
	void SetUseShowerShapeCut(Bool_t UseShowerShapeCut=kFALSE) { fUseShowerShapeCut=UseShowerShapeCut;};
	void SetBackground(Bool_t FillBackground=kFALSE) { fFillBackground=FillBackground;};
	void SetEMCalTriggerEG1() { fEMCEG1=kTRUE; };
	void SetEMCalTriggerEG2() { fEMCEG2=kTRUE; };
	void SetCentralityEstimator(Int_t Estimator) { fEstimator=Estimator; }; //0 = V0A, 1 = Other
	void SetAdditionalCuts(Double_t PtMinAsso, Int_t TpcNclsAsso) {fPtMinAsso = PtMinAsso; fTpcNclsAsso = TpcNclsAsso;};
	void SetSPDCutForHadrons() {fAssocWithSPD = kTRUE;};
	
		//Getters
	AliHFEpid *GetPID() const {return fPID;};
		//______________________________________________________________________
	
		//______________________________________________________________________
private:
	
		//Function to process track cuts
	Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);
		//Function to process eh analysis
	void ElectronHadronCorrelation(AliVTrack *track, Int_t trackIndex, AliVParticle *vtrack);
		//Function to find non-HFE and fill histos
	void Background(AliVTrack *track, Int_t trackIndex, AliVParticle *vtrack, Bool_t IsTPConly);
		//Selected Hadrons, for mixed event analysis
	TObjArray* SelectedHadrons();
		//DiHadron Correlation Background
	void DiHadronCorrelation(AliVTrack *track, Int_t trackIndex);
		//Find Mothers (Finde HFE and NonHFE from MC information)
	Bool_t FindMother(Int_t mcIndex);
	
		//Flags for specifics analysis
	Bool_t 				fCorrelationFlag;
	Bool_t				fIsMC;
	Bool_t				fUseEMCal;
	Bool_t				fUseTrigger;
	Bool_t				fUseShowerShapeCut;
	Bool_t				fFillBackground;
	Bool_t				fAssocWithSPD;
	
	
	Bool_t				fEMCEG1;
	Bool_t				fEMCEG2;
	
		//Used in the function FindMother
	Bool_t				fIsHFE1;
	Bool_t				fIsHFE2;
	Bool_t				fIsNonHFE;
	Bool_t				fIsFromD;
	Bool_t				fIsFromB;
	Bool_t				fIsFromPi0;
	Bool_t				fIsFromEta;
	Bool_t				fIsFromGamma;
	
		//General variables
	AliESDEvent 			*fESD;
	AliAODEvent 		   	*fAOD;				/// new
	AliVEvent 		      	*fVevent;			/// new
	AliESDtrackCuts         *fPartnerCuts;
	TList       			*fOutputList;
	AliPIDResponse 			*fPidResponse;
	AliSelectNonHFE 		*fNonHFE;
	
		//For the case of AOD analysis
	Bool_t					fIsAOD;					//flag for AOD analysis
	
		//For Centrality Selection
	AliCentrality			*fCentrality;
	Double_t				fCentralityMin;
	Double_t				fCentralityMax;
	Bool_t					fHasCentralitySelection;
	TH1F					*fCentralityHist;
	TH1F					*fCentralityHistPass;
	Float_t					fZvtx;	
	Int_t					fEstimator;
	
		//EMCal
	
	AliVCluster				*fClus;
		//AliESDCaloCluster 		*fClusESD;
	
		//Histograms
	TH1F				*fNevent;
	TH1F				*fPtElec_Inc;
	
	TH1F				*fPtPrim;
	TH1F				*fPtSec;
	TH1F				*fPtPrim2;
	TH1F				*fPtSec2;
	
	
	TH1F				*fCharge_n;
	TH1F				*fCharge_p;
	
	TH2D				*fTime;
	TH2D				*fTime2;
	TH2D				*ftimingEle;
	TH2D				*ftimingEle2;	
	
	TH1F				*fPtElec_ULS;
	TH1F				*fPtElec_LS;
	
	TH1F				*fPtElec_ULS_NoPid;
	TH1F				*fPtElec_LS_NoPid;
	
	TH1F				*fPtElec_ULS_MC;
	TH1F				*fPtElec_ULS_MC_weight;
	
	TH1F				*fPtElec_ULS2;
	TH1F				*fPtElec_LS2;
	
	TH1F				*fPtElec_ULS_weight;
	TH1F				*fPtElec_LS_weight;
	TH1F				*fPtElec_ULS2_weight;
	TH1F				*fPtElec_LS2_weight;
	
		//PID Histograms
	
	TH2F				*fTOF01;
	TH2F				*fTOF02;
	TH2F				*fTOF03;

	TH1F				*fpid;		
	
	TH2F				**fEoverP_pt;
	TH2F				**fEoverP_tpc;
	
	TH1F				**fTPC_pt;
	TH2F				**fTPC_p;
	
	TH1F				**fTPCnsigma_pt;
	TH2F				**fTPCnsigma_p;
	TH2F				*fTPCnsigma_pt_2D;
	TH2F				*fShowerShapeCut;
	TH2F				*fShowerShapeM02_EoverP;
	TH2F				*fShowerShapeM20_EoverP;
	TH2F				*fShowerShape_ha;
	TH2F				*fShowerShape_ele;
	
	TH2F				*fTPCnsigma_eta;
	TH2F				*fTPCnsigma_phi;
	
	
	TH1F				**fECluster;
	TH1F				*fECluster_pure;
	TH2F				**fEtaPhi;
	TH2F				*fEtaPhi_num;
	TH2F				*fEtaPhi_den;
	
	TH2F				*fpt_reco_pt_MC_num;
	TH2F				*fpt_reco_pt_MC_den;
	
	TH1F				**fVtxZ;
	
	TH1F				*fVtxZ_new1;
	TH1F				*fVtxZ_new2;
	TH1F				*fVtxZ_new3;
	TH1F				*fVtxZ_new4;
	
	TH1F				**fEtad;
	TH1F				**fNTracks;
	
	TH2F				**fNTracks_pt;
	TH2F				**fNTracks_eta;
	TH2F				**fNTracks_phi;
	
	TH1F				**fNClusters;
	TH2F				**fTPCNcls_EoverP;
	TH2F				**fTPCNcls_pid;
	
	TH1F				**fEta;
	TH1F				**fPhi;
	TH1F				**fR;
	TH2F				**fR_EoverP;
	TH1F				**fNcells;
	TH2F				**fNcells_EoverP;
	TH1F				**fNcells_electrons;
	TH1F				**fNcells_hadrons;
	TH1F				**fECluster_ptbins;
	TH1F				**fEoverP_ptbins;
	TH1F				**fEoverP_wSSCut;
	TH2F				**fM02_EoverP;
	TH2F				**fM20_EoverP;
	TH2F				**fTPCnsigma_eta_electrons;
	TH2F				**fTPCnsigma_eta_hadrons;
	
	TH2F				*fEoverP_pt_pions;
	
	TH2F				*ftpc_p_EoverPcut;
	TH2F				*fnsigma_p_EoverPcut;
	
	TH2F				*fEoverP_pt_pions2;
	TH2F				*fNcells_pt;
	TH2F				*fEoverP_pt_hadrons;
	
		//Electron-Hadron Correlation Histograms
	TH2F				**fCEtaPhi_Inc;
	
	TH2F				**fCEtaPhi_ULS;
	TH2F				**fCEtaPhi_LS;
	TH2F				**fCEtaPhi_ULS_NoP;
	TH2F				**fCEtaPhi_LS_NoP;
	
	TH2F				**fCEtaPhi_ULS_Weight;
	TH2F				**fCEtaPhi_LS_Weight;
	TH2F				**fCEtaPhi_ULS_NoP_Weight;
	TH2F				**fCEtaPhi_LS_NoP_Weight;
	
	TH1F				*fInvMass;
	TH1F				*fInvMassBack;
	TH1F				*fDCA;
	TH1F				*fDCABack;
	TH1F				*fOpAngle;
	TH1F				*fOpAngleBack;
	
	TH1F				*fInvMass2;
	TH1F				*fInvMassBack2;
	TH1F				*fDCA2;
	TH1F				*fDCABack2;
	TH1F				*fOpAngle2;
	TH1F				*fOpAngleBack2;
	
	Double_t			fMassCut;
	Double_t			fEtaCutMin;
	Double_t			fEtaCutMax;
	
	Double_t			fdPhiCut;
	Double_t			fdEtaCut;
	
	Double_t			fEoverPCutMin;
	Double_t			fEoverPCutMax;
	Double_t			fM20CutMin;
	Double_t			fM20CutMax;
	Double_t			fM02CutMin;
	Double_t			fM02CutMax;
	
	Double_t			fAngleCut;
	Double_t			fChi2Cut;
	Double_t			fDCAcut;
	Bool_t				fMassCutFlag;
	Bool_t				fAngleCutFlag;
	Bool_t				fChi2CutFlag;
	Bool_t				fDCAcutFlag;
	
	//Correlation Function
	Double_t			fAssHadronPtMin;
	Double_t			fAssHadronPtMax;
	
		//Non-HFE reconstruction efficiency
	TH1F				*fPtBackgroundBeforeReco;
	TH1F				*fPtBackgroundBeforeReco2;
	TH1F				*fPtBackgroundBeforeReco_weight;
	TH1F				*fPtBackgroundBeforeReco2_weight;
	
	TH2F				*fpT_m_electron;
	TH2F				*fpT_gm_electron;
	
	TH1F				*fPtBackgroundAfterReco;
	
	Double_t			fPtMinAsso;
	Int_t			fTpcNclsAsso;
	
		//Tracking Efficiency
	TH1F				*fPtMCparticleAll;
	TH1F				*fPtMCparticleAll_nonPrimary;
	TH1F				*fPtMCparticleAlle_nonPrimary;
	TH1F				*fPtMCparticleAlle_Primary;

	TH1F				*fPtMCparticleReco;
	TH1F				*fPtMCparticleReco_nonPrimary;
	TH1F				*fPtMCparticleAllHfe1;
	TH1F				*fPtMCparticleRecoHfe1;
	TH1F				*fPtMCparticleAllHfe2;
	TH1F				*fPtMCparticleRecoHfe2;
	TH1F				*fPtMCelectronAfterAll;
	TH1F				*fPtMCelectronAfterAll_unfolding;
	TH1F				*fPtMCelectronAfterAll_nonPrimary;
	TH1F				*fPtMCelectronAfterAll_Primary;
	
	TH1F				*fPtMCpi0;
	TH1F				*fPtMCeta;
	TH1F				*fPtMCpi02;
	TH1F				*fPtMCeta2;
	
	TH1F				*fPtMC_EMCal_All;
	TH1F				*fPtMC_EMCal_Selected;
	TH1F				*fPtMC_TPC_All;
	TH1F				*fPtMC_TPC_Selected;
	TH1F				*fPt_track_match_den;
	TH1F				*fPt_track_match_num;
	
		
	TH1F				*fPtMCWithLabel;
	TH1F				*fPtMCWithoutLabel;
	TH1F				*fPtIsPhysicaPrimary;
	
		//For the HFE package
	AliHFEcuts 			*fCuts;                 		// Cut Collection for HFE
														//Lucile
														//AliCaloTrackAODReader 			*reader; 
	AliCFManager 		*fCFM;                  		// Correction Framework Manager
	AliHFEpid 			*fPID;                  		// PID
	AliHFEpidQAmanager 	*fPIDqa;						// PID QA manager
	
		//Others
	AliStack 			*fMCstack;						//
	Bool_t              fRejectKinkMother;				//
	TParticle 			*fMCtrack;
	TParticle 			*fMCtrackMother;
	TParticle 			*fMCtrackGMother;
	TParticle 			*fMCtrackGGMother;
	TParticle 			*fMCtrackGGGMother;
	TClonesArray 		*fMCarray;
	AliAODMCHeader 		*fMCheader;
	AliAODMCParticle 	*fMCparticle;
	AliAODMCParticle 	*fMCparticle2;
	AliAODMCParticle 	*fMCparticleMother;
	AliAODMCParticle 	*fMCparticleGMother;
	AliAODMCParticle 	*fMCparticleGGMother;
	AliAODMCParticle 	*fMCparticleGGGMother;
	AliMCEventHandler	*fEventHandler;
	AliMCEvent			*fMCevent;
	
		//______________________________________________________________________
		//Mixed event analysis
	AliEventPoolManager *fPoolMgr;
	AliEventPool		*fPool;
	TObjArray			*fTracksClone;
	TObjArray			*fTracks;
	
	TH2F				**fCEtaPhi_Inc_EM;
	
	TH2F				**fCEtaPhi_ULS_EM;
	TH2F				**fCEtaPhi_LS_EM;
	
	TH2F				**fCEtaPhi_ULS_Weight_EM;
	TH2F				**fCEtaPhi_LS_Weight_EM;
	
	TH1F				*fPoolNevents;
	
	Bool_t				fEventMixingFlag;
		//______________________________________________________________________
	
		//______________________________________________________________________
		//Di-hadron correlation
	TH2F				**fCEtaPhi_Inc_DiHadron;
	TH1F				*fPtTrigger_Inc;
		//______________________________________________________________________
	
	AliAnalysisTaskEMCalHFEpA(const AliAnalysisTaskEMCalHFEpA&); 			// not implemented
	AliAnalysisTaskEMCalHFEpA& operator=(const AliAnalysisTaskEMCalHFEpA&); 		// not implemented
	
	ClassDef(AliAnalysisTaskEMCalHFEpA, 1); 								// example of analysis
																			//______________________________________________________________________
};

	///_________________________________________________________________________________________________
	///Class copied from : $ALICE_ROOT/PWGCF/Correlations/DPhi/AliAnalysisTaskLongRangeCorrelations.h
	///Author: Christoph Mayer
class AliEHCParticle : public TObject {
public:
	AliEHCParticle(Double_t eta=0, Double_t phi=0, Double_t pt=0)
    : fEta(eta), fPhi(phi), fPt(pt) {}
	virtual ~AliEHCParticle() {}
	
	Double_t Eta() const { return fEta; }
	Double_t Phi() const { return fPhi; }
	Double_t Pt() const { return fPt; }
	
protected:
private:
	AliEHCParticle(const AliEHCParticle&);
	AliEHCParticle& operator=(const AliEHCParticle&);
	
	Double_t fEta;
	Double_t fPhi;
	Double_t fPt;
	
	ClassDef(AliEHCParticle, 1);
} ;
	///_________________________________________________________________________________________________

#endif
