#ifndef AliAnalysisTask_JPsi_EMCal_cxx
#define AliAnalysisTask_JPsi_EMCal_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//      Task for J/psi analysis using EMCal and                       //
//      EMCal correction framework                                    //
//                                                                    //
//        v1.0                                                        //
//                                                                    //
//        Authors                                                     //
//                                                                    //
//        Cristiane Jahnke        (cristiane.jahnke@cern.ch)          //
//        22 January, 2021 -> TPC calibrations for 2017 and 2018 data //
//                                                                    //
////////////////////////////////////////////////////////////////////////

class TH1F;
class TH2F;
class TLorentzVector;
class THnSparse;
class TRandom3;
class AliESDEvent;
class AliESDtrackCuts;
class AliESDtrack;
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
#include "AliLog.h"
//______________________________________________________________________

//______________________________________________________________________
class AliAnalysisTask_JPsi_EMCal : public AliAnalysisTaskSE
{
//______________________________________________________________________
	public:
	AliAnalysisTask_JPsi_EMCal();
	AliAnalysisTask_JPsi_EMCal(const char *name);
	virtual ~AliAnalysisTask_JPsi_EMCal();
    virtual void         Init();
	virtual void   UserCreateOutputObjects();
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t *);

	//Setters
	//void SetHFECuts(AliHFEcuts * const cuts) {fCuts = cuts;};
	
	void SetMCanalysis() {fIsMC = kTRUE;};
	
    void SetAODanalysis(Bool_t IsAOD) {fIsAOD = IsAOD;};
	
    //trigger selection
	void SetEMCalTriggerEG1() { fEMCEG1=kTRUE; };
	void SetEMCalTriggerEG2() { fEMCEG2=kTRUE; };
	void SetEMCalTriggerDG1() { fEMCDG1=kTRUE; };
	void SetEMCalTriggerDG2() { fEMCDG2=kTRUE; };
    //trigger selection for EMCal+DCal together
    void SetEMCalTriggerEG1DG1() { fEMCEG1DG1=kTRUE; };
    void SetEMCalTriggerEG2DG2() { fEMCEG2DG2=kTRUE; };
	
	void SetUseTender() { fUseTender=kTRUE;};
    void SetMultiAnalysis() {fMultiAnalysis=kTRUE;};
    
    void Set_Fill_ESparse() {fFill_ESparse=kTRUE;};
    void Set_Fill_ESparseTPC() {fFill_ESparseTPC=kTRUE;};
    void Set_Fill_MSparse() {fFill_MSparse=kTRUE;};
    
    void Set_TPCCalibration(){fIs_TPC_calibration=kTRUE;};
    
    //to select events with high energy cluster (to mimic the trigger)
    void Set_Select_trigger_events2() {fSelect_trigger_events2=kTRUE;};
    void Set_Select_trigger_events1() {fSelect_trigger_events1=kTRUE;};
	//Setters analysis cuts
    
    //event cut
    void SetVertexCut(Double_t VertexCut ) {fVertexCut = VertexCut;};
    
    //track cuts
    void SetEtaCut(Double_t EtaCutMin,Double_t EtaCutMax ) {fEtaCutMin = EtaCutMin; fEtaCutMax = EtaCutMax;};
    
    void SetPtCutMainEle(Double_t PtCutMainEle ) {fPtCutMainEle = PtCutMainEle;};
    void SetPtCutPartner(Double_t PtCutPartner) {fPtCutPartner = PtCutPartner;};
    void SetRejectKinkMother(Bool_t rejectKinkMother = kFALSE) {fRejectKinkMother = rejectKinkMother;};
    void SetAODGlobalTracks(Bool_t AODGlobalTracks = kFALSE) {fAODGlobalTracks = AODGlobalTracks;};
    void SetTPCandITSrefit(Bool_t TPCandITSrefit = kFALSE) {fTPCandITSrefit = TPCandITSrefit;};
    
    void SetITSncls(Int_t ITSncls) {fITSncls = ITSncls;};
    void SetITSpixel(Int_t ITSpixel) {fITSpixel = ITSpixel;};
    void SetTPCncls(Int_t TPCncls) {fTPCncls = TPCncls;};
    void SetTPCnCrossedRows(Int_t TPCnCrossedRows) {fTPCnCrossedRows = TPCnCrossedRows;};
    void SetTPCnclsPID(Int_t TPCnclsPID) {fTPCnclsPID = TPCnclsPID;};
    void SetTPCchi2(Double_t TPCchi2) {fTPCchi2 = TPCchi2;};
    void SetITSchi2(Double_t ITSchi2) {fITSchi2 = ITSchi2;};
    void SetDCACut(Double_t DCAxyCut,Double_t DCAzCut ) {fDCAxyCut = DCAxyCut; fDCAzCut = DCAzCut;};
    
    //SPD corrections
    void SetMultiProfileSPD(TProfile2D * hprof){
        //for(Int_t i=0;i<1;i++){
            if(fMultEstimatorAvg[0]) delete fMultEstimatorAvg[0];
            fMultEstimatorAvg[0]=new TProfile2D(*hprof);
        //}
    }
    Double_t         GetTrackletsMeanCorrection(TProfile2D* estimatorAvg, Double_t uncorrectedNacc, Double_t vtxZ, Double_t refMult, Int_t run_number); /*const*/
    
    Double_t GetTPCCalibration(Int_t runNo, Double_t TPCnsigma0);/*const*/
    
    //V0 correction
    void SetMultiProfileV0(TProfile2D * hprofV0){
        if(fMultEstimatorV0[0]) delete fMultEstimatorV0[0];
        fMultEstimatorV0[0]=new TProfile2D(*hprofV0);
    }
    Double_t  GetV0MeanCorrection(TProfile2D* estimatorV0, Double_t uncorrectedV0, Double_t vtxZ, Double_t refMult_V0, Int_t run_number); /*const*/
    

    
    //TPC PID cuts
    void SetTPCnsigmaCut(Double_t TPCnsigmaCutMin,Double_t TPCnsigmaCutMax ) {fTPCnsigmaCutMin = TPCnsigmaCutMin; fTPCnsigmaCutMax = TPCnsigmaCutMax;};
    
    //EMCal cuts
	void SetEnergyCut(Double_t EnergyCut) {fEnergyCut= EnergyCut;};
	void SetEoverPCut(Double_t EoverPCutMin,Double_t EoverPCutMax ) { fEoverPCutMin = EoverPCutMin; fEoverPCutMax = EoverPCutMax; };
    
    //Mass cut
    void SetMassCut(Double_t MassCutMin,Double_t MassCutMax ) { fMassCutMin = MassCutMin; fMassCutMax = MassCutMax; };
   
    //To apply weights due to J/psi enhancement
    Double_t CalculateWeight(Double_t x);
	//Getters
	//AliHFEpid *GetPID() const {return fPID;};
//______________________________________________________________________
  
//______________________________________________________________________
	private:
	
//Function to process track cuts
	Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);

//Find Mothers (Finde HFE and NonHFE from MC information)
	Bool_t FindMother(Int_t mcIndex);
		//Bool_t TrackCuts(AliVTrack *track, AliAODTrack *atrack, AliESDtrack *etrack, Bool_t fIsAOD);


    Bool_t				fIsMC;
	Bool_t				fUseTender;
    Bool_t              fMultiAnalysis;
    Bool_t              fFill_ESparse;
    Bool_t              fFill_ESparseTPC;
    Bool_t              fFill_MSparse;
    Bool_t              fIs_TPC_calibration;
    Bool_t              fSelect_trigger_events1;
    Bool_t              fSelect_trigger_events2;
    
    //new organization of tender using global variables
    TString        fTenderClusterName;//
    TString            fTenderTrackName;//
    
    TClonesArray*         fTracks_tender;//Tender tracks
    TClonesArray*         fCaloClusters_tender;//Tender cluster


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
    TList                   *fListProfiles; // SPD vtx correction
	AliPIDResponse 			*fPidResponse;
	//AliSelectNonHFE 		*fNonHFE;
	
//For the case of AOD analysis
	Bool_t					fIsAOD;					//flag for AOD analysis
	
//EMCal threshold separation
	Bool_t				fEMCEG1;
	Bool_t				fEMCEG2;
    
    
    Bool_t               fEMCEG1DG1;//to run both EMCal and DCal triggers together
    Bool_t               fEMCEG2DG2;
	
//DCal threshold separation
	Bool_t				fEMCDG1;
	Bool_t				fEMCDG2;
    
    Bool_t                fIsTrack1Emcal;
    Bool_t                fIsTrack1Dcal;
    Bool_t                fIsTrack2Emcal;
    Bool_t                fIsTrack2Dcal;
    
    Bool_t              fIsEMCalCls;
    Bool_t              fIsDCalCls;
    
    //SPD corrections
    TProfile2D*        fMultEstimatorAvg[1];
    Double_t         fRefMult;
    TProfile2D*        GetEstimatorHistogram(const AliAODEvent *fAOD);
    TRandom3*        gRandom;//!< random number generator
	
    
    //V0 corrections
    TProfile2D*        fMultEstimatorV0[1];
    Double_t         fRefMult_V0;
    TProfile2D*        GetEstimatorHistogram_V0(const AliAODEvent *fAOD);
    TRandom3*        gRandom_V0;//!< random number generator
    
//AnalysisCuts
    
    //Event cut
    Double_t            fVertexCut;
    
    //Track cuts
    
    Double_t            fEtaCutMin;
    Double_t            fEtaCutMax;
    Double_t            fPtCutMainEle;
    Double_t            fPtCutPartner;
    Bool_t              fRejectKinkMother;                //
    Bool_t              fAODGlobalTracks;
    Bool_t              fTPCandITSrefit;
    Int_t               fITSncls;
    Int_t               fITSpixel;
    Int_t               fTPCncls;
    Int_t               fTPCnCrossedRows;
    Int_t               fTPCnclsPID;
    Double_t            fTPCchi2;
    Double_t            fITSchi2;
    Double_t            fDCAxyCut;
    Double_t            fDCAzCut;
    
    
    
    Double_t            fTPCnsigmaCutMin;
    Double_t            fTPCnsigmaCutMax;
    
    
    
	Double_t			fEnergyCut;

	Double_t			fEoverPCutMin;
	Double_t			fEoverPCutMax;
    
    Double_t            fMassCutMin;
    Double_t            fMassCutMax;
	
	
	

//Vertex selection
	Float_t					fZvtx;
    
//global multiplicity values
    Double_t                    fV0Mult;
    Double_t                    fSPDMult;
    Double_t                    fV0Mult_corr;
    Double_t                    fV0Mult_corr2;
    Double_t                    fSPDMult_corr;
    
    
	
//EMCal
	//AliESDCaloCluster 	*fClus;
	AliVCluster				*fClus;
	AliVCluster				*fClus2;
	AliAODCaloCluster	    *fClusAOD;
	
//Histograms for the analysis
	TH1F				*fNevent;
    TH1F                *fNevent2;
    TH2F                **fTPC_vs_ITScls;
    TH1F                *fPDG_values;
    TH1F                *fNevent_SPD_multi;
    TH1F                *fNevent_V0_multi;
    
	TH2F				**fEoverP_pt;
	TH2F				**fTPC_p;
    
	TH2F				**fTPCnsigma_p;
    TH2F                *fTPCnsigma_p_beforeCalibration;
    TH2F                *fTPCnsigma_p_afterCalibration;
    
    TH2F                **fTOF_p;
    TH2F                **fTOFnsigma_p;
	
   
    
	
	TH2F				**fTPCnsigma_EoverP;
	TH1F				**fECluster;
	TH1F				**fECluster_emcal;
	TH1F				**fECluster_dcal;
	
	TH1F				*fECluster_pure;
	TH1F				*fECluster_pure_emcal;
    
    TH1F                *fECluster_pure_emcal_SPD1;
    TH1F                *fECluster_pure_emcal_SPD2;
    TH1F                *fECluster_pure_emcal_SPD3;
    TH1F                *fECluster_pure_emcal_SPD4;
    TH1F                *fECluster_pure_emcal_SPD5;
    
    TH1F                *fECluster_pure_emcal_V01;
    TH1F                *fECluster_pure_emcal_V02;
    TH1F                *fECluster_pure_emcal_V03;
    TH1F                *fECluster_pure_emcal_V04;
    TH1F                *fECluster_pure_emcal_V05;
  
	TH1F				*fECluster_pure_dcal;
	
	TH2F				*fEtaPhi_both;
	TH2F				*fEtaPhi_emcal;
	TH2F				*fEtaPhi_dcal;

	
	TH1F				**fTracksPt;
	TH1F				**fTracksQAPt;
    TH1F                **fTracksMCPt;
	
	TH1F				**fVtxZ;
    
    //histos for SPD and V0 multiplicity
    TH2F                 *fVtxZ_V0;
    TH2F                 *fVtxZ_SPD;
    TH2F                 *fV0_SPD;
    TH2F                 *fV0_nch;
    TH2F                 *fSPD_nch;
    
    TH1F                 *fV0;
    TH1F                 *fV01;
    TH1F                 *fV02;
    
    TH1F                 *fSPD;
    TH1F                 *fSPD1;
    TH1F                 *fSPD2;
    
    
    TH1F				**fNClusters;
    
    TH1F                *fNClusters_pure;
    TH2F                *fEoverP_ntracks_matched;
    TH2F                *fEoverP_ncells;

	

	
//Others
	
	TClonesArray 		*fMCarray;
	AliAODMCHeader 		*fMCheader;
	AliAODMCParticle 	*fMCparticle;
	AliAODMCParticle 	*fMCparticleMother;

	
	AliAODMCParticle 	*fMCparticle2;
	AliAODMCParticle 	*fMCparticleMother2;
	
	
	AliAODMCParticle 	*fMCparticleGMother;
	AliAODMCParticle 	*fMCparticleGGMother;
	AliAODMCParticle 	*fMCparticleGGGMother;
	AliMCEventHandler	*fEventHandler;
	AliMCEvent			*fMCevent;
	
	//JPsi histos

	
	//KF
	TH2F				*fHist_InvMass_pt_ULS_KF;
	TH2F				*fHist_InvMass_pt_LS_KF;
    
    TH2F                *fHist_Correlation_leg1_emcal_leg2_not;
    TH2F                *fHist_Correlation_leg1_not_leg2_emcal;
    TH2F                *fHist_Correlation_leg1_emcal_leg2_emcal;
    
    TH2F                *fHist_InvMass_pt_ULS_KF_weight;
    
    //multiplicity histos
    TH2F                *fHist_InvMass_pt_ULS_KF_SPDmulti_1;
    TH2F                *fHist_InvMass_pt_ULS_KF_SPDmulti_2;
    TH2F                *fHist_InvMass_pt_ULS_KF_SPDmulti_3;
    TH2F                *fHist_InvMass_pt_ULS_KF_SPDmulti_4;
    TH2F                *fHist_InvMass_pt_ULS_KF_SPDmulti_5;
    
    TH2F                *fHist_InvMass_pt_ULS_KF_V0multi_1;
    TH2F                *fHist_InvMass_pt_ULS_KF_V0multi_2;
    TH2F                *fHist_InvMass_pt_ULS_KF_V0multi_3;
    TH2F                *fHist_InvMass_pt_ULS_KF_V0multi_4;
    TH2F                *fHist_InvMass_pt_ULS_KF_V0multi_5;
    
    //with weight
    //KF
    
    
    
    //multiplicity histos
    TH2F                *fHist_InvMass_pt_ULS_KF_SPDmulti_1_weight;
    TH2F                *fHist_InvMass_pt_ULS_KF_SPDmulti_2_weight;
    TH2F                *fHist_InvMass_pt_ULS_KF_SPDmulti_3_weight;
    TH2F                *fHist_InvMass_pt_ULS_KF_SPDmulti_4_weight;
    TH2F                *fHist_InvMass_pt_ULS_KF_SPDmulti_5_weight;
    
    TH2F                *fHist_InvMass_pt_ULS_KF_V0multi_1_weight;
    TH2F                *fHist_InvMass_pt_ULS_KF_V0multi_2_weight;
    TH2F                *fHist_InvMass_pt_ULS_KF_V0multi_3_weight;
    TH2F                *fHist_InvMass_pt_ULS_KF_V0multi_4_weight;
    TH2F                *fHist_InvMass_pt_ULS_KF_V0multi_5_weight;
	
    
	//generators
	//BB
	TH2F				*fHist_InvMass_pt_ULS_KF_BB;
	TH2F				*fHist_InvMass_pt_LS_KF_BB;
	//CC
	TH2F				*fHist_InvMass_pt_ULS_KF_CC;
	TH2F				*fHist_InvMass_pt_LS_KF_CC;
	//B
	TH2F				*fHist_InvMass_pt_ULS_KF_B;
	TH2F				*fHist_InvMass_pt_LS_KF_B;
	//JPsi
	TH2F				*fHist_InvMass_pt_ULS_KF_Jpsi;
	TH2F				*fHist_InvMass_pt_LS_KF_Jpsi;
	//BJpsi
	TH2F				*fHist_InvMass_pt_ULS_KF_BJpsi;
	TH2F				*fHist_InvMass_pt_LS_KF_BJpsi;
	
	//leg 1 on EMCal
	TH2F				*fHist_InvMass_pt_ULS1;
	TH2F				*fHist_InvMass_pt_LS1;
	//leg 2 on EMCal
	TH2F				*fHist_InvMass_pt_ULS2;
	TH2F				*fHist_InvMass_pt_LS2;
	//both legs on EMCal
	TH2F				*fHist_InvMass_pt_ULSboth;
	TH2F				*fHist_InvMass_pt_LSboth;
	
	
	TH2F				*fHist_InvMass_pt_ULStpc;
	TH2F				*fHist_InvMass_pt_LStpc;
    
    TH2F                *fHist_InvMass_pt_ULStpc_wMatching;
    TH2F                *fHist_InvMass_pt_LStpc_wMatching;
	
		//new histos
	TH2F                **fdEta_dPhi;
	
	THnSparse  *fSparseElectron;//!Electron info
    THnSparse  *fSparseElectronTPC;//!Electron info
    
	Double_t   *fvalueElectron;//!Electron info
    Double_t   *fvalueElectronTPC;//!Electron info
    
    THnSparse  *fSparseMulti;//!Multiplicity info
    Double_t   *fvalueMulti;//!Multiplicity info
	
	
	
	//MC efficiencies
	TH1F				*fPtMCparticleRecoHfe1;
	TH1F				*fPtMCparticleAllHfe1;
	TH1F				*fPtMCparticleAll_e_from_JPsi;
    TH1F                *fPtMCparticleAll_JPsi_pT;
    
    TH1F                *fPtMCparticleAll_e_from_JPsi_electron;
    TH1F                *fPtMCparticleAll_JPsi_pT_electron;
    TH1F                *fPtMCparticleAll_e_from_JPsi_positron;
    TH1F                *fPtMCparticleAll_JPsi_pT_positron;
    
    TH1F                *fPtMCparticleAll_electrons;
    TH1F                *fPtMCparticleAll_particles;
    
    //tracking efficiency
    TH1F                *fPtMCparticleReco_electrons;
    TH1F                *fPtMCparticleReco_electrons_no_gamma;
    TH1F                *fPtMCparticleReco_particles;
    //TPC PID efficiency
    TH1F                *fPtMCparticle_TPCpid_e_from_JPsi;
    TH1F                *fPtMCparticle_TPCpid_electrons;
    TH1F                *fPtMCparticle_TPCpid_e_from_JPsi_num;
    TH1F                *fPtMCparticle_TPCpid_electrons_num;
    //EMCal PID efficiency
    
    TH1F                *fPtMCparticle_EMCalpid_leg1;
    TH1F                *fPtMCparticle_EMCalpid_leg2;
    
    TH1F                *fPtMCparticle_EMCal_TM_e_from_JPsi;
    TH1F                *fPtMCparticle_EMCal_TM_electrons;
    TH1F                *fPtMCparticle_EMCalpid_leg1_e_from_JPsi;
    TH1F                *fPtMCparticle_EMCalpid_leg2_e_from_JPsi;
    TH1F                *fPtMCparticle_EMCalpid_both_leg1_e_from_JPsi;
    TH1F                *fPtMCparticle_EMCalpid_both_leg2_e_from_JPsi;
    TH1F                *fPtMCparticle_Total_JPsi_pT;
    //J/Psi reco
    TH1F                *fPtMCparticle_JPsi;
    TH1F                *fPtMCparticle_JPsi_num;
    //J/Psi mass cut
    TH1F                *fPtMCparticle_JPsi_mass;
    TH1F                *fPtMCparticle_JPsi_mass_num;
    
    
    
    TH1F                *fPtMCparticleAll_trueJPsi_pT;
    TH1F                *fPtMCparticleAll_trueJPsi_pT_weight;
    TH1F                *fPtMCparticleAll_trueJPsi_pT_weight_prompt;
    
	TH1F				*fPtMCparticleReco_e_from_JPsi;
    
 
    
	TH1F				*fPtMCparticle_Total_e_from_JPsi;
    TH1F                *fPtMCparticle_Total_e_from_JPsi_sameMother;
	TH1F				*fPtMCparticle_TotalplusMass_e_from_JPsi;
    TH1F                *fPtMCparticle_TotalplusMass_e_from_JPsi_sameMother;
    TH1F                *fPtMCparticle_TotalplusMass_JPsi_pT;
    TH1F                *fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother;
    TH1F                *fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight;
    TH1F                *fPtMCparticle_TotalplusMass_JPsi_pT_eSameMother_weight_prompt;
	

//______________________________________________________________________

	AliAnalysisTask_JPsi_EMCal(const AliAnalysisTask_JPsi_EMCal&); 			// not implemented
	AliAnalysisTask_JPsi_EMCal& operator=(const AliAnalysisTask_JPsi_EMCal&); 		// not implemented
  
	ClassDef(AliAnalysisTask_JPsi_EMCal, 1); 								// example of analysis
//______________________________________________________________________
};


///_________________________________________________________________________________________________

#endif
