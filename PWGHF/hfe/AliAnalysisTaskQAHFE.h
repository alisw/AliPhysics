

#ifndef AliAnalysisTaskQAHFE_cxx
#define AliAnalysisTaskQAHFE_cxx

class TH1F;
class TH2F;
class TH2D;
class TH3F;
class AliAODEvent;
class TList;
class AliESDEvent;
class AliHFEpid;
class AliHFEpidQAmanager;
class AliAnalysisDataSlot;
class AliAnalysisDataContainer;
class AliAnalysisTaskQAHFE;
class AliCFManager;
class AliHFEcuts;
class AliHFEpidQAmanager;
class AliMultSelection;
class AliSelectNonHFE;
class TClonesArray;
class AliAODMCParticle;
class TClonesArray;
class AliAODParticle;
class AliVertexingHFUtils;

#include <TProfile.h>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTask.h"
#include "AliAnalysisManager.h"
#include "AliVEvent.h"
#include "AliAODEvent.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisDataSlot.h"
#include "AliCFManager.h"
//#include "AliMultSelection.h"
#include "AliSelectNonHFE.h"
#include "TClonesArray.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisUtils.h"
#include "AliESDUtils.h"
//#include "AliAODHandler.h"
//#include "AliAnalysisTaskQAHFE.h"

class AliAnalysisTaskQAHFE : public AliAnalysisTaskSE {
 public:

	enum pi0etaType {kNoIsPrimary,kNoMother,kLightMesons,kBeauty, kCharm, kNoFeedDown}; //0,1,2,3,4,5
   enum ESourceType {kNoMotherE, kPi0NoFeedDown, kEtaNoFeedDown, kGPi0NoFeedDown, kGEtaNoFeedDown, kDirectGamma, kOthersE};//0,1,2,3,4,5,6
    
    
	AliAnalysisTaskQAHFE() ;
	AliAnalysisTaskQAHFE(const char *name);
	virtual ~AliAnalysisTaskQAHFE() ;
	  
	virtual void   Init();
	virtual void   UserCreateOutputObjects();
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t *);
	Int_t ClassifyTrack(AliAODTrack* track,const AliVVertex* pVtx);
	Int_t GetNcharged();
	
	//Find Mothers (Find HFE and NonHFE from MC information)
	Int_t GetHFE(AliAODMCParticle *, TClonesArray *);
  
	Int_t GetElecSourceType(AliAODMCParticle *, TClonesArray *,Double_t &ptm);
	Int_t GetPi0EtaType(AliAODMCParticle *, TClonesArray *);
	void SelectPhotonicElectronMC(Int_t ,AliAODTrack *, AliPIDResponse *,Int_t,Int_t);
	void SelectPhotonicElectronR(Int_t itrack, AliAODTrack *track, Int_t motherindex, Int_t pdg,Int_t source,Double_t V0Mmult1 , Double_t SPDntr1);
	Double_t Beta(AliAODTrack *track);
	AliHFEpid *GetPID() const { return fPID; }
	  
	//Setters
	void SetMCAnalysis(Bool_t isMC){fIsMC=isMC;}
	void SetReferenceMultiplicity(Double_t multi){fRefMult=multi;}
		
	void SetMultiplVsZProfile_NoEtaCut(TProfile* hprof){
		if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
		fMultEstimatorAvg[0]=new TProfile(*hprof);
  	}
	void SetMultiplVsZProfile(TProfile* hprof){
		if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
		fMultEstimatorAvg[1]=new TProfile(*hprof);
	}
	//----------Setter for Track and PID cuts
	void SetTrigger(AliVEvent::EOfflineTriggerTypes trigger){ftrigger =trigger;}
	//void SetTrigger(Int_t trigger){ftrigger =trigger;}
	void SetEtaRange(Double_t Etarange){fEtarange=Etarange;}
	void SetMinTPCCluster(Int_t TPCNclus){fTPCNclus=TPCNclus;}
	void SetMinITSCluster(Int_t ITSNclus){fITSNclus=ITSNclus;}
	void SetMinTPCClusterPID(Int_t TPCNclusPID){fTPCNclusPID=TPCNclusPID;}
	void SetHitsOnSPDLayers(Bool_t SPDBoth,Bool_t SPDAny, Bool_t SPDFirst)	{
			fSPDBoth=SPDBoth;
			if(!SPDBoth)fSPDAny=SPDAny;
			if(!SPDBoth && !SPDAny)fSPDFirst=SPDFirst;
	}
	void SetDCACut(Double_t DCAxyCut,Double_t DCAzCut){
		fDCAxyCut=DCAxyCut;
		fDCAzCut=DCAzCut;
	}
	void SetTPCnsigma(Double_t TPCnsigmin,Double_t TPCnsigmax){
			fTPCnsigmin=TPCnsigmin;
			fTPCnsigmax=TPCnsigmax;
	}
	void SetTOFnsigma(Double_t TOFnsig){fTOFnsig=TOFnsig;}
	
	//------------Setters for Photonic Electron Selection Cuts
	void SetInvMassCut(Double_t InvmassCut){fInvmassCut=InvmassCut;}
	void SetAssoTPCclus(Int_t AssoTPCCluster){fAssoTPCCluster=AssoTPCCluster;}
	void SetAssoITSrefit(Bool_t AssoITSRefit){fAssoITSRefit= AssoITSRefit;}
	void SetAssopTMin(Double_t AssopTMin){fAssopTMin = AssopTMin;}
	void SetAssoEtarange(Double_t AssoEtarange){fAssoEtarange=AssoEtarange;}
	void SetAssoTPCnsig(Double_t AssoTPCnsig){fAssoTPCnsig=AssoTPCnsig;}
	void SetFilterBit(Double_t Filterbit){fFilterbit=Filterbit;}
	
	private:
	
	//------------Track and PID cut variables--------------
	AliVEvent::EOfflineTriggerTypes ftrigger;
	//Int_t ftrigger;
	Int_t fTPCNclus;  
	Int_t fITSNclus;  
	Int_t fTPCNclusPID;  
	Bool_t fSPDBoth;  
	Bool_t fSPDAny;  
	Bool_t fSPDFirst;  
	Double_t fDCAxyCut;  
	Double_t fDCAzCut;  
	Double_t fpTMin;  
	Double_t fEtarange; 
	Double_t    fFilterbit;
	Double_t fTPCnsigmin;  
	Double_t fTPCnsigmax;  
	Double_t fTOFnsig;  
	
	Double_t fInvmassCut;	//	  invariant mass cut value
	Int_t fAssoTPCCluster;  
	Bool_t fAssoITSRefit;  
	Double_t fAssopTMin;  
	Double_t fAssoEtarange;  
	Double_t fAssoTPCnsig;  
	
	
	Double_t fDCAxy;  
	Double_t fDCAz;  
  //---------------------------------------------------------------
	AliAODEvent *fAOD;    //! AOD object
	AliHFEpid   *fPID;                  //!PID
	AliPIDResponse   *fpidResponse; //!pid response
	AliAODVertex   *fAODVertex;  //!fAODVertex

	AliSelectNonHFE *fNonHFE; //!
	 
	Bool_t fIsMC;
	TList       *fOutputList; //! Output list

	TH1D        *fcount;//!
	TH1F        *fHistPt; //! Pt spectrum
	TH1F        *fHistMult; //! Multiplicity Distribution
	TH1F        *fTPCSignal; //! TPCSignal vs entries

	TH1F *fVtxZ;//!
	TH1F *fVtxZ_corr;//!
	TH1F *fVtxY;//!
	TH1F *fVtxX;//!
	TH2F *fDCAZvsPt;//!
	TH2F *fDCAXYvsPt;//!

	TH2F	      *fdEdxVsP_TPC;//!
	TH2F	      *fdEdxVsP_TPC_cut ;//!
	TH2F	      *fnSigmaVsP_TPC;//!
	TH2F	      *fnSigmaVsP_TPC_cut;//!
	TH2F	      *fnSigmaVsP_TOF;//!
	TH2F	      *fnBetaVsP_TOF;//!
	
	TH2F	      *fdEdxVsPt_TPC;//!
	TH2F	      *fnSigmaVsPt_TPC;//!
	TH2F	      *fnSigmaVsPt_TOF;//! 
	TH2F	      *fnBetaVsPt_TOF;//!

	TH2F	      *fnSigmaVsEta_TPC;//!
	TH2F	      *fnSigmaVsEta_TPC_cut;//!
  	  
	TH1F        *fNentries;         	//!histogram with number of events on output slot 2 
	TH1F        *fNentries2;         	//!histogram with number of events on output slot 3 
	TF1					*fLandau;//!
	TF1					*fErr;//!
	  
	TH1D *fHadCot_Landau;//!
	TH1D *fHadCot_Err;	//!
	TH1D *fPt_incl_e_Landau;//!
	TH1D *fPt_incl_e_Err;	//!	

  	Double_t fMultiPercentile;
  	Double_t fMultiPercentileSPD;
  	
  	TH1F *fV0Mult; //!
  	TH1F *fV0MultCorr; //!
  	TH1F *fV0AMult; //!
  	TH1F *fV0CMult; //!
  	TH1F *fSPD_tracklet_NoEtaCut; 		  //!
  	TH1F *fSPD_tracklet; 		  //!
  	TH1F *fCent;//!
  	TH1F *fCentSPD;//!

  	TH2F *fPteV0M_Lan;//!
  	TH2F *fPteV0M_Err;//!
  	
  	TH2F *fPteSPD_Lan;//!
  	TH2F *fPteSPD_Err;//!
  	
  	TH2F *fPteV0M;//!
  	TH2F *fPteSPD;//!
  					
  	TH2F *fHadV0MLan;//!
  	TH2F *fHadSPDLan;//!
  	TH2F *fHadV0MErr;//!
  	TH2F *fHadSPDErr;//!
  	
  	/*TH3F *fnSigmaVsP_TPC_cut_multV0M;//!
  	TH3F *fnSigmaVsP_TPC_cut_multSPD;//!
  	TH3F *fnSigmaVsPt_TPC_cut_multV0M;//!
  	TH3F *fnSigmaVsPt_TPC_cut_multSPD;//!
	*/	
	//===============NonHFE========================================
	TH1F        *fPtHFEMCtot;//!
	TH2F        *fPtHFEMCtot_SPD;//!
	TH2F        *fPtHFEMCtot_V0M;//!
	TH1F        *fPtHFEMC;//! HFE pt before track cut	
	TH2F        *fPtHFEMC_SPD;//! HFE pt before track cut	
	TH2F        *fPtHFEMC_V0M;//! HFE pt before track cut	
	
	TH1F *fPtHFEMC_afterfilterbit;//!
	TH1F        *fPtHFEMC_aftertrackcuts;//!  
	TH2F        *fPtHFEMC_aftertrackcuts_SPD;//! 
	TH2F        *fPtHFEMC_aftertrackcuts_V0M;//! 
	TH1F        *fPtHFEMC_aftertracktofcuts;//!
	TH2F        *fPtHFEMC_aftertracktofcuts_SPD;//!
	TH2F        *fPtHFEMC_aftertracktofcuts_V0M;//!

	
	TH1F        *fPte_ULS; //ULS elec Pt
   TH1F        *fPte_LS;// LS elec pt
   TH1F        *fPte_ULS_KF; //ULS elec Pt
   TH1F        *fPte_LS_KF;// LS elec pt  
   TH1F        *fInvmassLS1; //LS Invmass 
   TH1F        *fInvmassULS1;//ULS Invmass
   
   TH2F        *fPte_ULS_multV0M; //ULS elec Pt
   TH2F        *fPte_LS_multV0M;// LS elec pt 
   
   TH2F        *fPte_ULS_multSPD; //ULS elec Pt
   TH2F        *fPte_LS_multSPD;// LS elec pt 
   
   TH1F        *fPt_elec_phot;//
   TH1F        *fPt_elec_phot1;//
   TH1F        *fpt_e_nonphot_MC;//!
   TH1F        *fElecNos;//!
   TH1F        *fPT_elec_MCtrue;//!
   TH1F        *fpt_incl_elec_MC;//!
   TH1F        *fPi0_Pt_MC;//!
   TH1F        *fEta_Pt_MC;//!
   TH1F        *fGamma_Pt_MC;//!
   TH1F        *fPte_LS_invMss_MC;//!
   TH1F        *fPte_ULS_invMss_MC;//!
   TH1F        *fPte_LS_MC;//!
   TH1F        *fPte_ULS_MC;//!
   TH2F			*fPte_LS_MC_multV0M;//!
   TH2F			*fPte_ULS_MC_multV0M;//!
	TH2F			*fPte_LS_MC_multSPD;//!
	TH2F			*fPte_ULS_MC_multSPD;//!
   TH2F			*fInvMULS;//!
   
   TClonesArray 		*fMCarray;//
   AliAODMCParticle 	*fMCparticle;//
   AliAODMCParticle 	*fMCmother;//
   
   //TH3F				*fInvMULS_multV0M;//!
   //TH3F				*fInvMULS_multSPD;//!
   TH2F        *fPt_elec_phot_multV0M;//!
   TH2F        *fPt_elec_phot_multSPD;//!
   TH2F        *fPt_elec_phot1_multV0M;//!
   TH2F        *fPt_elec_phot1_multSPD;//!
   
  
   
   Int_t pdg; 
   
   TH1F        *fPt_elec_pi0_MB;//!
	TH1F        *fPt_elec_pi0_MB_ULS;//!
	TH1F        *fPt_elec_pi0_MB_LS;//!
	TH1F        *fPt_elec_eta_MB;//!
	TH1F        *fPt_elec_eta_MB_ULS;//!
	TH1F        *fPt_elec_eta_MB_LS;//!
	TH1F        *fPt_elec_gamma_MB;//!
	TH1F        *fPt_elec_gamma_MB_ULS;//!
	TH1F        *fPt_elec_gamma_MB_LS;//!
   
   THnSparseF  *fInvMULSnSp;//!
   THnSparseF  *fPi0EtaSpectra;//!
	THnSparseF  *fPi0EtaSpectraSp;//!
	TH2F        *fInvMULSpi0;//!
	TH2F        *fInvMULSeta;//!
	TH2F        *fInvMULSgamma;//!
   
   TH1F        *fPtHFEMC_reco;//! HFE pt after track cut	
   TH1F        *fPtHFEMC_reco_true;//! HFE pt after track cut	
	TH2F        *fPtHFEMC_reco_SPD;//! HFE pt after track cut	
	TH2F        *fPtHFEMC_reco_V0M;//! HFE pt after track cut
		
   TH1F        *pi0MC;//!
	TH1F        *etaMC;//!
	TH1F        *gammaMC;//!
   
    TH1F        *pi0MC_1;//!
	TH1F        *etaMC_1;//!
	TH1F        *gammaMC_1;//!

   
   //-------------------------SPD V0M Multiplicity -------------------------------
   
		
  	TFile *fileEstimator;//!
  	TH1F *fSPDCorrMultDistNoEtaCut_min;//!
  	TH1F *fSPDCorrMultDistNoEtaCut_max;//!
  	TH1F *fSPDCorrMultDist_min;//!
  	TH1F *fSPDCorrMultDist_max;//!
  
  	TH2F *fHistNtrVsZvtxNoEtaCut;//!
  	TH2F *fHistNtrCorrVsZvtxNoEtaCut_min;//!
  	TH2F *fHistNtrCorrVsZvtxNoEtaCut_max;//!
  
  	TH2F *fHistNtrVsZvtx;//!
  	TH2F *fHistNtrCorrVsZvtx_min;//!
  	TH2F *fHistNtrCorrVsZvtx_max;//!
  	TProfile *fMultEstimatorAvg[2]; /// TProfile with mult vs. Z per period
  	TH2F *fSPDCorrMultDist_min_vs_AliSPDMult;//!
  	
  	Double_t 		fRefMult;
  	TH2F *fV0MultVsZvtx;//!
  	TH2F *fV0MultCorrVsZvtx;//!
  	TH2F *fMultSPD_vs_alimult;//!
  	TH2F *fMultV0M_vs_alimult;//!
  	TH2F *fNchVsZvtx;//!
  	TH2F *fSPDNtrCorrVsNch;//!
		TH2F *fV0MCorrVsNch;//!
		TH2F *fSPDtr_V0Mmult;//!
	//	TH3F *fSPDNtrCorrVsV0MCorrVsNch;//!
  
  	
	//-------------------------------------------------------------------------------------
	  
	AliAnalysisTaskQAHFE(const AliAnalysisTaskQAHFE&); // not implemented
	AliAnalysisTaskQAHFE& operator=(const AliAnalysisTaskQAHFE&); // not implemented
	  
	ClassDef(AliAnalysisTaskQAHFE, 1); // ClassDef(ClassName,ClassVersionID)  
};
	
#endif	
