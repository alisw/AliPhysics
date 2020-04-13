

#ifndef AliAnalysisTask_QA_EMCALElectrons_cxx
#define AliAnalysisTask_QA_EMCALElectrons_cxx

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
class AliAnalysisTask_QA_EMCALElectrons;
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
//#include  "AliHelperPID.h"
#include "AliCFManager.h"
//#include "AliMultSelection.h"
#include "AliSelectNonHFE.h"
#include "TClonesArray.h"
#include "AliAODMCParticle.h"
#include "AliAnalysisUtils.h"
#include "AliESDUtils.h"
//#include "AliAODHandler.h"
//#include "AliAnalysisTask_QA_EMCALElectrons.h"

class AliAnalysisTask_QA_EMCALElectrons : public AliAnalysisTaskSE {
 public:
    enum EnhanceSigOrNot {kMB,kEnhance};
	enum pi0etaType {kNoIsPrimary,kNoMother,kLightMesons,kBeauty, kCharm, kNoFeedDown}; //0,1,2,3,4,5
   enum ESourceType {kNoMotherE, kPi0NoFeedDown, kEtaNoFeedDown, kGPi0NoFeedDown, kGEtaNoFeedDown, kDirectGamma, kOthersE};//0,1,2,3,4,5,6
    
    
	AliAnalysisTask_QA_EMCALElectrons() ;
	AliAnalysisTask_QA_EMCALElectrons(const char *name);
	virtual ~AliAnalysisTask_QA_EMCALElectrons() ;
	  
	virtual void   Init();
	virtual void   UserCreateOutputObjects();
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t *);
	void    GetTrkClsEtaPhiDiff(AliVTrack *t, AliVCluster *v, Double_t &phidiff, Double_t &etadiff);
	 Bool_t 		PassEIDCuts(AliAODTrack *track, AliVCluster *clust);// electron identification cuts
	Int_t ClassifyTrack(AliAODTrack* track,const AliVVertex* pVtx);
	//Int_t GetNcharged();
	
	//Find Mothers (Find HFE and NonHFE from MC information)
	Int_t GetHFE(AliAODMCParticle *, TClonesArray *);
  
	Int_t GetElecSourceType(AliAODMCParticle *,Double_t &ptm);
	Int_t GetPi0EtaType(AliAODMCParticle *);

	Bool_t GetNMCPartProduced();
	void    GetPi0EtaWeight(THnSparse *SparseWeight);
	Bool_t  IsNonHFE(AliAODMCParticle *MCPart, Bool_t &fFromMB, Int_t &type, Int_t &iMom, Int_t &MomPDG, Double_t &MomPt);	
	  
	 
	//  void  SetTagEffiFromEnhanceMC(Bool_t EffiFromEnhMC) {fEffiFromEnhMC = EffiFromEnhMC;};
    /*Bool_t  GetNonHFEEffiDenomGenPurMC(AliVTrack *track);
     Bool_t  GetNonHFEEffiRecoTagGenPurMC(AliVTrack *track);
    Bool_t  GetNonHFEEffiULSLSGenPurMC(AliVTrack *track, AliVTrack *Assotrack, Bool_t fFlagLS, Bool_t fFlagULS, Double_t mass);
  */

    Bool_t  GetNonHFEEffiDenom(AliVTrack *track);
    Bool_t  GetNonHFEEffiRecoTag(AliVTrack *track);
    Bool_t  GetNonHFEEffiULSLS(AliVTrack *track, AliVTrack *Assotrack, Bool_t fFlagLS, Bool_t fFlagULS, Double_t mass);
	 void SelectPhotonicElectron(Int_t itrack, AliVTrack *track, Bool_t &fFlagPhotonicElec, Int_t iMC);
	 
	
	void SelectPhotonicElectronR(Int_t itrack, AliAODTrack *track, Int_t motherindex, Int_t pdg,Int_t source , Double_t SPDntr1,Double_t ptmotherwg);
	/*Double_t WeightMCncCorr(Int_t multSPDtr);
	Double_t Beta(AliAODTrack *track);*/
	Double_t GetCorrectedNtracklets(TProfile* estimatorAvg, Double_t uncorrectedNacc, Double_t vtxZ, Double_t refMult);
	TProfile* GetEstimatorHistogram(const AliAODEvent* fAOD);
	AliHFEpid *GetPID() const { return fPID; }
	  
	//Setters
	void SetMCAnalysis(Bool_t isMC){fIsMC=isMC;}
	void SetReferenceMultiplicity(Double_t multi){fRefMult=multi;}
		
	void SetMultiplVsZProfile_16l(TProfile* hprof){
		if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
		fMultEstimatorAvg[0]=new TProfile(*hprof);
  	}
	void SetMultiplVsZProfile_17d20a2_extra(TProfile* hprof){
		if(fMultEstimatorAvg[1]) delete fMultEstimatorAvg[1];
		fMultEstimatorAvg[1]=new TProfile(*hprof);
	}
	void SetMultiplVsZProfile_16k(TProfile* hprof){
		if(fMultEstimatorAvg[2]) delete fMultEstimatorAvg[2];
		fMultEstimatorAvg[2]=new TProfile(*hprof);
  	}
	void SetMultiplVsZProfile_17d20a1_extra(TProfile* hprof){
		if(fMultEstimatorAvg[3]) delete fMultEstimatorAvg[3];
		fMultEstimatorAvg[3]=new TProfile(*hprof);
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
	void SetEopE(Double_t EopEMin,Double_t EopEMax){
			fCutEopEMin=EopEMin;
			fCutEopEMax=EopEMax;
	}

	void SetShowerShapeEM20(Double_t M20Min,Double_t M20Max){
			fCutM20Min=M20Min;
			fCutM20Max=M20Max;
	}
	//-----------------Setter For EMCal--------------------
	void                  SetEMCalTriggerEG1(Bool_t flagTrEG1) { fEMCEG1=flagTrEG1;};
  	void                  SetEMCalTriggerEG2(Bool_t flagTrEG2) { fEMCEG2=flagTrEG2;};
  	void                  SetEMCalTriggerDG1(Bool_t flagTrDG1) { fDCalDG1=flagTrDG1;};
  	void                  SetEMCalTriggerDG2(Bool_t flagTrDG2) { fDCalDG2=flagTrDG2;};
  	//void    		SetTenderSwitch(Bool_t usetender){fUseTender = usetender;};

	
	//------------Setters for Photonic Electron Selection Cuts
	void SetInvMassCut(Double_t InvmassCut){fInvmassCut=InvmassCut;}
	void SetAssoTPCclus(Int_t AssoTPCCluster){fAssoTPCCluster=AssoTPCCluster;}
	void SetAssoITSrefit(Bool_t AssoITSRefit){fAssoITSRefit= AssoITSRefit;}
	void SetAssopTMin(Double_t AssopTMin){fAssopTMin = AssopTMin;}
	void SetAssoEtarange(Double_t AssoEtarange){fAssoEtarange=AssoEtarange;}
	void SetAssoTPCnsig(Double_t AssoTPCnsig){fAssoTPCnsig=AssoTPCnsig;}
	
	
	private:
	
	//------------Track and PID cut variables--------------
	AliVEvent::EOfflineTriggerTypes ftrigger;
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
	Double_t fTPCnsigmin;  
	Double_t fTPCnsigmax;  
	Double_t fCutEopEMin;
  	Double_t fCutEopEMax;
  	Bool_t IsM20;
  	Double_t fCutM20Min;
  	Double_t fCutM20Max;
	Double_t fInvmassCut;	//	  invariant mass cut value
	Int_t fAssoTPCCluster;  
	Bool_t fAssoITSRefit;  
	Double_t fAssopTMin;  
	Double_t fAssoEtarange;  
	Double_t fAssoTPCnsig;  
	
	
	Double_t fDCAxy;  
	Double_t fDCAz;  
	   
   //====EMCal==================
   Bool_t              fFlagClsTypeEMC;//switch to select EMC clusters
   Bool_t              fFlagClsTypeDCAL;//switch to select DCAL clusters

	Bool_t                fEMCEG1;//EMcal Threshold EG1
	Bool_t                fEMCEG2;//EMcal Threshold SetReferenceMultiplicityEG2
	Bool_t                fDCalDG1;//DCal Threshold DG1
  	Bool_t                fDCalDG2;//DCal Threshold DG2
  	
  	TString		fTenderClusterName;//
  	TString	        fTenderTrackName;//
  
  	TClonesArray*         fTracks_tender;//Tender tracks     
  	TClonesArray*         fCaloClusters_tender;//Tender cluster     
   	Bool_t                fUseTender;// switch to add tender
 

 
	
  //---------------------------------------------------------------
	AliAODEvent *fAOD;    //! AOD object
	AliHFEpid   *fPID;                  //!PID
	AliPIDResponse   *fpidResponse; //!pid response
	AliAODVertex   *fAODVertex;  //!fAODVertex
   AliAODMCHeader *fMCHeader;
   
	AliSelectNonHFE *fNonHFE; //!
	 
	Bool_t fIsMC;
	TList       *fOutputList; //! Output list

	TH1D        *fcount;//!
	TH1F        *fHistPt; //! Pt spectrum
	TH1F        *fHistMult; //! Multiplicity Distribution
	TH1F			*feta;//!
	TH1F        *fTPCSignal; //! TPCSignal vs entries

	TH1F *fVtxZ;//!
	TH1F *fVtxZ_corr;//!
	
	TH2F	      *fdEdxVsP_TPC;//!
	TH2F	      *fnSigmaVsP_TPC;//!
	
	TH2F	      *fnSigmaVsPt_TPC;//!
	  	  
	TH1F        *fNentries;         	//!histogram with number of events on output slot 2 
	TH1F        *fNentries2;         	//!histogram with number of events on output slot 3 
	
	   TH2F *fClusEtaPhi;//
   TH1F *fClusT;//
   TH1F *fNCells ;//
   TH1F *fClusE;//
   TH2F *fClusEvsnTracklets;//
	
	TH1F                *fHistPtMatch;//!
    TH2F                *fEMCTrkMatch;//!
    TH1F                *fEMCTrkPt;//!
    TH1F                *fEMCTrketa;//!
    TH1F                *fEMCTrkphi;//!
    TH2F                *fEMCTPCnsig;//!
    TH1F                *fClsEAftMatch;//!
     TH2F                *fClsEAftMatch_SPD;//!
    TH1F                *fClsEopAftMatch;//!
    TH2F                *fClsEtaPhiAftMatch;//!
    THnSparse* 	        fSparseElectron; //! Electron information
    Double_t*             fvalueElectron; //!	
  	//===============NonHFE========================================
	
   	
   	
   TH1F        *fPte_ULS; //! ULS elec Pt
   TH1F        *fPte_LS;//! LS elec pt  
   TH1F        *fInvmassLS1; //! LS Invmass 
   TH1F        *fInvmassULS1;//! ULS Invmass
   
   TH2F        *fPte_ULS_multV0M; //! ULS elec Pt
   TH2F        *fPte_LS_multV0M;//! LS elec pt 
   
   TH2F        *fPte_ULS_multSPD; //! ULS elec Pt
   TH2F        *fPte_LS_multSPD;//! LS elec pt 
   
   

//====================MC================
 	TClonesArray 		*fMCArray;//
   AliAODMCParticle 	*fMCparticle;//
   AliAODMCParticle 	*fMCmother;//
   Int_t pdg; 
	TH1F        *fPtHFEMC;//! HFE pt before track cut	
	TH2F        *fPtHFEMC_SPD;//! HFE pt before track cut	
	TH1F        *fPtHFEMC_reco;//! HFE pt after track cut	
	TH2F        *fPtHFEMC_reco_SPD;//! HFE pt after track cut
	TH1F        *fPtHFEMC_trackcutreco;//!
   TH2F        *fPtHFEMC_trackcutreco_SPD;//!
   TH1F        *fPtHFEMC_trackmatchreco;//!
   TH2F        *fPtHFEMC_trackmatchreco_SPD;//!
   TH1F *fPtHFEMC_TPCEMCreco;//!
	TH2F *fPtHFEMC_TPCEMCreco_SPD;//!
	TH1F *fPtHFEMC_SScutreco;//!
 	TH2F *fPtHFEMC_SScutreco_SPD;//!
	TH1F *fPtHFEMC_TPCreco;//!
	TH2F *fPtHFEMC_TPCreco_SPD;//!
	TH1F        *fPT_elec_MC;//!
	TH1F        *fPT_PIDCaloCut_MC;//!
	TH1F        *fPT_PIDCaloCut_realistic;//!
	TH1F        *fPT_elec_realistic;//!
	TH1F        *fPt_elec_phot1;//!
	TH2F        *fPt_elec_phot1_multSPD;//!
	TH2F			*fInvMULS;//!
	THnSparseF  *fInvMULSnSp;//!
	THnSparseF  *fPi0EtaSpectraSp;//!
	TH1F        *pi0MC;//!
	TH1F        *etaMC;//!
	TH1F        *gammaMC;//!

	
	//++++++++++Weights++++++++++++++++++++
	    Double_t fRefMult;
		TH1F *fSPD_tracklet;//
   TH1F *fSPDCorrMultDist_max;//
   TH1F *fSPDWeightedCorrMultDist_max;//
   	TProfile *fMultEstimatorAvg[4]; /// TProfile with mult vs. Z per period
   	TProfile *Profile_Mean; //
   	TProfile *Profile_MeanCorr; //
   	Double_t SPDntr;//
   	
   	
    Int_t               fNTotMCpart; //! N of total MC particles produced by generator
    Int_t               fNpureMC;//! N of particles from main generator (MB/Enhanced)
    Int_t               fNembMCpi0; //! N > fNembMCpi0 = particles from pi0 generator
    Int_t               fNembMCeta; //! N > fNembMCeta = particles from eta generator
   	Bool_t              fIsFrmEmbPi0;//!
    Bool_t              fIsFrmEmbEta;//!
    Bool_t              fIsFrmPi0;//!
    Bool_t              fIsFrmEta;//!
    Int_t               ftype;//!
    Double_t            fWeight;//!
    Double_t            fWeightPi0;//!
    Double_t            fWeightEta;//!



   	THnSparse           *fSprsPi0EtaWeightCal;//!
   	Bool_t              fCalculateWeight;//
   	
   	 TF1                 *fPi0Weight;//!
    TF1                 *fEtaWeight;//!
   	
   	TH1F                *fRealInclsElecPt;//!
   	 TH1F                *fNonHFeTrkPt;//!
    TH1F                *fNonHFeEmbTrkPt;//!
    TH1F                *fNonHFeEmbWeightTrkPt;//!
    TH1F                *fPi0eEmbWeightTrkPt;//!
    TH1F                *fEtaeEmbWeightTrkPt;//!
    
        TH1F                *fRecoNonHFeTrkPt;//!
    TH1F                *fRecoNonHFeEmbTrkPt;//!
    TH1F                *fRecoNonHFeEmbWeightTrkPt;//!
    TH1F                *fRecoPi0eEmbWeightTrkPt;//!
    TH1F                *fRecoEtaeEmbWeightTrkPt;//!
   
   TH1F                *fULSElecPt;//!
    TH1F                *fLSElecPt;//!
     TH2F        *fInvmassULSPt;//!Invmass of ULS
    TH2F        *fInvmassLSPt;//!Invmass of LS
  
   	

	TH1F                *fNonHFePairInvmassLS;//!
    TH1F                *fNonHFePairInvmassULS;//!
    TH1F                *fNonHFeEmbInvmassLS;//!
    TH1F                *fNonHFeEmbInvmassULS;//!
    TH1F                *fNonHFeEmbWeightInvmassLS;//!
    TH1F                *fNonHFeEmbWeightInvmassULS;//!
    TH1F                *fPi0EmbInvmassLS;//!
    TH1F                *fPi0EmbInvmassULS;//!
    TH1F                *fPi0EmbWeightInvmassLS;//!
    TH1F                *fPi0EmbWeightInvmassULS;//!
    TH1F                *fEtaEmbInvmassLS;//!
    TH1F                *fEtaEmbInvmassULS;//!
    TH1F                *fEtaEmbWeightInvmassLS;//!
    TH1F                *fEtaEmbWeightInvmassULS;//!
    TH1F                *fRecoLSeEmbTrkPt;//!
    TH1F                *fRecoLSeEmbWeightTrkPt;//!
    TH1F                *fRecoPi0LSeEmbWeightTrkPt;//!
    TH1F                *fRecoEtaLSeEmbWeightTrkPt;//!
    TH1F                *fRecoULSeEmbTrkPt;//!
    TH1F                *fRecoULSeEmbWeightTrkPt;//!
    TH1F                *fRecoPi0ULSeEmbWeightTrkPt;//!
    TH1F                *fRecoEtaULSeEmbWeightTrkPt;//!
    
    TH2F                *fNonHFeEmbTrkPt_SPD;//!
	TH2F                *fNonHFeEmbWeightTrkPt_SPD;//!
	TH2F                *fRecoNonHFeTrkPt_SPD;//!
	TH2F                *fRecoNonHFeEmbTrkPt_SPD;//!
	TH2F                *fRecoNonHFeEmbWeightTrkPt_SPD;//!

 	//-------------------------------------------------------------------------------------
	  
	AliAnalysisTask_QA_EMCALElectrons(const AliAnalysisTask_QA_EMCALElectrons&); // not implemented
	AliAnalysisTask_QA_EMCALElectrons& operator=(const AliAnalysisTask_QA_EMCALElectrons&); // not implemented
	  
	ClassDef(AliAnalysisTask_QA_EMCALElectrons, 1); // ClassDef(ClassName,ClassVersionID)  
};
	
#endif	
