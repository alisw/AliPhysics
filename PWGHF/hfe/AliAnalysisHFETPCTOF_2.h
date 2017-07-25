#ifndef AliAnalysisHFETPCTOF_2_cxx
#define AliAnalysisHFETPCTOF_2_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//      Task for Heavy-flavour electron analysis in pp collisions     //
//      (and also Pb-Pb)             								  //
//																	  //
//		v1.0														  //
//                                                                    //
//	    Authors 							                          //
//		Camila de Conti (camila.de.conti@cern.ch)				      //
//                                                                    //
////////////////////////////////////////////////////////////////////////

class AliAnalysisUtils;
class TH1F;
class TH2F;
//class TF1;
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


//______________________________________________________________________
//Library
#include "AliAnalysisTaskSE.h"
#include "AliHFEpid.h"
#include "AliLog.h"
//#include "AliMultSelection.h"
//______________________________________________________________________

//______________________________________________________________________
class AliAnalysisHFETPCTOF_2 : public AliAnalysisTaskSE 
{
//______________________________________________________________________
	public:
	AliAnalysisHFETPCTOF_2();
	AliAnalysisHFETPCTOF_2(const char *name);
	virtual ~AliAnalysisHFETPCTOF_2();
  
	virtual void   UserCreateOutputObjects();
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t *);

	//Setters
	void SetHFECuts(AliHFEcuts * const cuts) {fCuts = cuts;};
	void SetRejectKinkMother(Bool_t rejectKinkMother = kFALSE) {fRejectKinkMother = rejectKinkMother;};
	
	void SetMCanalysis() {fIsMC = kTRUE;};
    void SetAODanalysis(Bool_t IsAOD) {fIsAOD = IsAOD;};
    void SetPPanalysis(Bool_t IsPP) {fIsPP = IsPP;};
 
    //Setter for the Hadron contamination function
	void SetHadronFunction(TF1* HadronF) {fHadrons = HadronF;};

	//Setter for the pi0 weight
	void SetPi0Weight(TF1* Pi0F) {fPi0w = Pi0F;};
	void SetPi0Weight2(TF1* Pi0F2) {fPi0w2 = Pi0F2;};
	void SetPi0Weight3(TF1* Pi0F3) {fPi0w3 = Pi0F3;};

	//Setter for the eta weight
	void SetEtaWeight(TF1* EtaF) {fEtaw = EtaF;};
	void SetEtaWeight2(TF1* EtaF2) {fEtaw2 = EtaF2;};
	void SetEtaWeight3(TF1* EtaF3) {fEtaw3 = EtaF3;};
	
	//Setter for the Partner cuts
	void SetPartnerCuts(Float_t Mass, Float_t MinPt, Float_t TpcNclus);
	
	//Setter for the PID cuts (TPC and TOF)
	void SetPIDCuts(Float_t tpcPIDmincut, Float_t tpcPIDmaxcut, Float_t tofPIDmincut, Float_t tofPIDmaxcut);
	
	//Setter for the Eta cut
	void SetEtaCut(Float_t EtaMin, Float_t EtaMax);
			
	//Getters
	AliHFEpid *GetPID() const {return fPID;};
//______________________________________________________________________
  
//______________________________________________________________________
	private:
	
//Function to process track cuts
	Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);

//Find Mothers (Find HFE and NonHFE from MC information)
	Bool_t FindMother(Int_t mcIndex);

//Select HFE for the reconstruction efficiency calculation
	Bool_t IsHFelectronsMC(AliVTrack *track);

//Select pi0 and eta for weight calculation
	Int_t GetPi0EtaType(AliAODMCParticle *pi0eta, TClonesArray *mcArray);
	
// ----- weights for tagging efficiency -----
enum pi0etaType {kNoMother, kNoFeedDown, kNoIsPrimary, kLightMesons, kBeauty, kCharm};

//Select photonic electrons
	Int_t GetElecSourceType(AliAODMCParticle *electron, TClonesArray *mcArray,Double_t &ptm);
	
enum ESourceType {kNoMotherE, kPi0NoFeedDown, kEtaNoFeedDown, kGPi0NoFeedDown, kGEtaNoFeedDown, kDirectGamma, kOthersE};	

// ------------------------------------------
	

//Check if is heavy-flavor electron 
	//Bool_t IsHFelectrons(Int_t mcIndex);

//Flags for specifcs analysis
    Bool_t				fIsMC;  
    Bool_t				fIsPP;  
    
//Weight to normalize the amount of pi and eta    
    //Double_t CalculateWeight(Int_t pdg_particle, Double_t x);

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
	TList       			*fOutputList;
	AliPIDResponse 			*fPidResponse;
	AliSelectNonHFE 		*fNonHFE;
	
	Bool_t				fMassCutFlag;
	Bool_t				fAngleCutFlag;
	Bool_t				fChi2CutFlag;
	Bool_t				fDCAcutFlag;
	
    Double_t		    	fMassCut;
    Double_t			    fAngleCut;
    Double_t			    fChi2Cut;
    Double_t			    fDCAcut;
    
    Double_t			    fPtMinAsso;
    Int_t			        fTpcNclsAsso;
    
    AliESDtrackCuts         *fPartnerCuts;
    
    
	
//For the case of AOD analysis
	Bool_t					fIsAOD;					//flag for AOD analysis
//	
	
//Vertex selection
	Float_t					fZvtx;
	
//EMCal
	//AliESDCaloCluster 		*fClus;
	AliVCluster				*fClus;
	
	
//Histograms for the analysis
	TH1F				*fVertex1;//!
	TH1F				*fNevent; //!
	TH1F				*fNeventT0; //!
	TH1F				*fNevent_3;	//!
	TH1F				*fNevent_T0b; //!
	TH1F				*fNevent_corrcut;//!
	TH1F				*fNevent_no_vertex; //!
	TH1F				*fNevent_no_vertex_2; //!
	TH1F				*fCent;	//!
	TH1F				*fCent2;	//!
	TH2F				*fTPC_p1;//!
	TH2F				*fTPC_p2;//!
	TH2F				*fTPC_p3;//!
	TH1F				*fPt_1;//!
	TH1F				*fPt_2;//!
	TH1F				*fITSnClus_1;//!
	TH1F				*fITSnClus_2;//!
	TH1F				*fTPCnClus_1;//!
	TH1F				*fTPCnClus_2;//!
	TH2F				*fTPCnsigma_p1;//!
	TH2F				*fTPCnsigma_p2;//!
	TH2F				*fTPCnsigma_p3;//!
	TH2F				*fTPCnsigma_pt1;//!
	TH2F				*fTPCnsigma_pt2;//!
	TH2F				*fTPCnsigma_pt3;//!
	TH2F                *fTOFnsigma_p1;//!
	TH2F                *fTOFnsigma_p2;//!
	TH2F                *fTOFnsigma_p3;//!
	TH2F                *fTOFnsigma_pt1;//!
	TH2F                *fTOFnsigma_pt2;//!
	TH2F                *fTOFnsigma_pt3;//!
	TH2F                *fITSnsigma_p1;//!
	TH2F                *fITSnsigma_p2;//!
	TH2F                *fITSnsigma_p3;//!
	TH2F                *fITSnsigma_pt1;//!
	TH2F                *fITSnsigma_pt2;//!
	TH2F                *fITSnsigma_pt3;//!
	TH2F				*fTPCnSigmaPElec_ULS;//!
	TH2F				*fTPCnSigmaPElec_LS;//!
	//TH2F				**fTPCnsigma_EoverP;
	TH2F                *fTPCnsigma_TOFnsigma1;//!
	TH2F                *fTPCnsigma_TOFnsigma2;//!
	TH2F                *fTPCnsigma_TOFnsigma3;//!
	//TH1F				**fVtxZ;
    //TH1F				**fNClusters;
    //TH2F                **fEtaPhiTracks;
    //TH2F                **fEtaPhiClusters;
	TH1F				*fInvMass;//!
	TH1F				*fInvMassBack;//!
	TH1F				*fInvMass1;//!
	TH1F				*fInvMassBack1;//!
	TH1F				*fInvMass2;//!
	TH1F				*fInvMassBack2;//!
	TH1F				*fInvMass3;//!
	TH1F				*fInvMassBack3;//!
	TH1F				*fInvMass4;//!
	TH1F				*fInvMassBack4;//!
	TH2F                *fTPCnsigma_p_after_tof;//!
	TH2F                *fTPCnsigma_pt_after_tof;//!
	TH2F                *fTPCnsigma_p_after_tof_its;//!
	TH2F                *fTPCnsigma_pt_after_tof_its;//!
	
	TH1F				*fDCA;//!
	TH1F				*fDCABack; //!
	TH1F				*fAngle;//!
	TH1F				*fAngleBack;//!

    TH1F                *fPtElec_ULS;//!
    TH1F                *fPtElec_LS;//!
    TH1F                *fPtElec_ULS_MC;//!
    TH1F                *fPtElec_LS_MC;//!
    TH1F                *fPtElec_ULS_MC2;//!
    TH1F                *fPtElec_LS_MC2;//!
    TH1F                *fPtElec_ULS_MC_from_pi0;//!
    TH1F                *fPtElec_LS_MC_from_pi0;//!
    
    TH1F                *fPtElec_ULS_MC_from_pi02;//!
    TH1F                *fPtElec_LS_MC_from_pi02;//!
    
    TH1F                *fPtElec_ULS_MC_from_eta2;//!
    TH1F                *fPtElec_LS_MC_from_eta2;//!
    
    TH1F                *fPtElec_LS_MC_from_gamma;//!
    TH1F                *fPtElec_ULS_MC_from_gamma2;//!
    TH1F                *fPtElec_LS_MC_from_gamma2;//!
    
    TH1F                *fPtElec;//!
    TH1F                *fPElec;//!
    TH1F				*fPtHad_f;//!
    TH1F				*fPHad_f;//!
    TH2F                *fDCAz_pt;//!
    TH2F                *fDCAxy_pt;//!
    TH2F                *fDCAz_pt_gamma;//!
    TH2F                *fDCAxy_pt_gamma;//!
    TH2F                *fDCAz_pt_light;//!
    TH2F                *fDCAxy_pt_light;//!
    TH2F                *fDCAz_pt_B;//!
    TH2F                *fDCAxy_pt_B;//!
    TH2F                *fDCAz_pt_C;//!
    TH2F                *fDCAxy_pt_C;//!
    TH1F                *fPt_elec_phot;//!
    TH1F                *fPt_elec_phot2;//!
    
    TH1F 				*fPt_elec_MB;//!
    TH1F 				*fPt_elec_MB_ULS;//!
    TH1F 				*fPt_elec_MB_LS;//!
	
	TH1F 				*fPt_elec_pi0_MB2;//!
    TH1F 				*fPt_elec_pi0_MB_ULS2;//!
    TH1F 				*fPt_elec_pi0_MB_LS2;//!
	
	TH1F				*fPt_elec_eta_MB2;//!
	TH1F				*fPt_elec_eta_MB_ULS2;//!
	TH1F				*fPt_elec_eta_MB_LS2;//!
	
	TH1F	 			*fPt_elec_gamma_MB;//!
	TH1F	 			*fPt_elec_gamma_MB_ULS;//!
	TH1F	 			*fPt_elec_gamma_MB_LS;//!
    
    TH1F	 			*fPt_elec_gamma_MB2;//!
	TH1F	 			*fPt_elec_gamma_MB_ULS2;//!
	TH1F	 			*fPt_elec_gamma_MB_LS2;//!

    TH1F                *fPt_elec_from_pi02;//!
    TH1F                *fPt_elec_from_eta2;//!
    TH1F                *fPt_elec_from_gamma2;//!
    TH1F                *fPt_HFelec_passTOFandTPC_2;//!	
    TH1F                *fPt_HFelec_passTPC_2;	//!
    TH1F				*fPt_HFelec_pass_track;//!
    TH1F				*fPt_HFelec_pass_trackDCA;//!
    TH1F				*fPt_HFelec_pass_trackITSlayer;//!
    TH1F				*fPt_HFelec_pass_trackFB;//!
    TH1F				*fPt_HFelec_pass_trackTPCITS;//!
    TH1F				*fPt_HFelec_pass_trackKink;//!
    TH1F				*fPt_HFelec_pass_trackandTOF;//!
    TH1F                *fPtMCeta;//!
    TH1F                *fPtMCeta_MB_afterNorm;//!
    TH1F                *fPtMCpi03;//!
    TH1F                *fPtMCeta3;//!
    TH1F                *fPtHFEMC_2;//!

		
//For the HFE package
	AliHFEcuts 			*fCuts;                 		// Cut Collection for HFE
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
	AliAODMCParticle 	*fMCparticleMother;
	AliAODMCParticle 	*fMCparticleGMother;
	AliAODMCParticle 	*fMCparticleGGMother;
	AliAODMCParticle 	*fMCparticleGGGMother;
	AliMCEventHandler	*fEventHandler;
	AliMCEvent			*fMCevent;
	TF1 				*fHadrons;
	TF1					*fPi0w;
	TF1					*fPi0w2;
	TF1					*fPi0w3;
	TF1					*fEtaw;
	TF1					*fEtaw2;
	TF1					*fEtaw3;
	Float_t				fMass;
	Float_t				fMinPt;
	Float_t				fTpcNclusAsso;
	Float_t				ftpcPIDmincut;
	Float_t				ftpcPIDmaxcut;
	Float_t				ftofPIDmincut;
	Float_t				ftofPIDmaxcut;
	Float_t				fEtaMin;
	Float_t				fEtaMax;
	
	
	// ----- weights for tagging efficiency -----
    THnSparseF           *fPi0EtaSpectra;// ! QA for weights
    // ------------------------------------------
	
//______________________________________________________________________

	AliAnalysisHFETPCTOF_2(const AliAnalysisHFETPCTOF_2&); 			// not implemented
	AliAnalysisHFETPCTOF_2& operator=(const AliAnalysisHFETPCTOF_2&); 		// not implemented
  
	ClassDef(AliAnalysisHFETPCTOF_2, 1); 								// example of analysis
//______________________________________________________________________
};

///_________________________________________________________________________________________________

#endif
