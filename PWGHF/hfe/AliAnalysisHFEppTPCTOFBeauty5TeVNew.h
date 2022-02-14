#ifndef AliAnalysisHFEppTPCTOFBeauty5TeVNew_cxx
#define AliAnalysisHFEppTPCTOFBeauty5TeVNew_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//      Task for Beauty analysis in p-p collisions  TPC-TOF           //
//                						      //
//								      //
//	           						      //
//                                                                    //
//	    Authors 						      //
//  Sudhir Pandurang Rode (sudhir.pandurang.rode@cern.ch              //
//                                                                    //
////////////////////////////////////////////////////////////////////////

class AliAnalysisUtils;
class TH1F;
class TH2F;
class TGraphErrors;
//class TF1;
class AliESDEvent;
class AliESDtrackCuts;
class AliESDtrack;
class AliHFEcontainer;
class AliHFEcuts;
class AliHFEextraCuts;
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
class AliGenEventHeader;

//______________________________________________________________________
//Library
#include "AliAnalysisTaskSE.h"
#include "AliHFEpid.h"
#include "AliLog.h"
//______________________________________________________________________

//______________________________________________________________________
class AliAnalysisHFEppTPCTOFBeauty5TeVNew : public AliAnalysisTaskSE
{
    //______________________________________________________________________
public:
    
    enum HijingOr {kHijing,kPhytia,kpi0,keta};
    enum ESourceType {kNoMotherE, kPi0NoFeedDown, kEtaNoFeedDown, kGPi0NoFeedDown, kGEtaNoFeedDown, kDirectGamma, kOthersE};
    enum pi0etaType {kNoMother, kNoFeedDown, kNoIsPrimary, kLightMesons, kKaonFromNonHF, kBeauty, kCharm, kKaonFromHF};

    AliAnalysisHFEppTPCTOFBeauty5TeVNew();
    AliAnalysisHFEppTPCTOFBeauty5TeVNew(const char *name);
    virtual ~AliAnalysisHFEppTPCTOFBeauty5TeVNew();
    
    virtual void   UserCreateOutputObjects();
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *);
    
    //Setters
    void SetHFECuts(AliHFEcuts * const cuts) {fCuts = cuts;};
    void SetRejectKinkMother(Bool_t rejectKinkMother = kFALSE) {fRejectKinkMother = rejectKinkMother;};
    
    void SetMCanalysis() {fIsMC = kTRUE;};
    void SetAODanalysis(Bool_t IsAOD) {fIsAOD = IsAOD;};
    void SetPPanalysis(Bool_t IsPP) {fIsPP = IsPP;};
       
    //Setter for the PID cuts (TPC and TOF)
    void SetPIDCuts(Float_t tpcPIDmincut, Float_t tpcPIDmaxcut, Float_t tofPIDmincut, Float_t tofPIDmaxcut);
    
   // void SetAssoPIDCuts(Float_t AssotpcPIDmincut, Float_t AssotpcPIDmaxcut, Float_t AssotofPIDmincut, Float_t AssotofPIDmaxcut);
    
   
    void SetAssoMinNClustersTPC(Float_t kassTPCcl){kAssoTPCcl = kassTPCcl;};//Minimum number of clusters on TPC
    void SetAssoMinNClustersTPCPID(Float_t kassTPCPIDcl){kAssoTPCPIDcl = kassTPCPIDcl;};
    void SetAssoMaxImpactParam(Float_t kassoDCAr, Float_t kassoDCAzd){
    kassoDCAxy = kassoDCAr;
    kassoDCAz = kassoDCAzd;
    }; //DCA to vertex
    
    //Setter for the Eta cut
    void SetEtaCut(Float_t EtaMin, Float_t EtaMax);
  //  void SetAssoEtaCut(Float_t AssoEtaMin, Float_t AssoEtaMax);
    //Setter for the B correction function
    void SetBcorrFunction(TGraphErrors* BcorrF) {fBcorr = BcorrF;};
	
	
    //Setter for the D correction function
    void SetDcorrFunction1(TF1* DcorrF1) {fDcorr1 = DcorrF1;};
    void SetDcorrFunction2(TF1* DcorrF2) {fDcorr2 = DcorrF2;};
    void SetDcorrFunction3(TF1* DcorrF3) {fDcorr3 = DcorrF3;};
    void SetDcorrFunction4(TF1* DcorrF4) {fDcorr4 = DcorrF4;};
    void SetDcorrFunction5(TF1* DcorrF5) {fDcorr5 = DcorrF5;};
    void SetDcorrFunction6(TF1* DcorrF6) {fDcorr6 = DcorrF6;};
    void SetDcorrFunction7(TF1* DcorrF7) {fDcorr7 = DcorrF7;};
    void SetDcorrFunction8(TF1* DcorrF8) {fDcorr8 = DcorrF8;};
    void SetDcorrFunction9(TF1* DcorrF9) {fDcorr9 = DcorrF9;};
    void SetDcorrFunction10(TF1* DcorrF10) {fDcorr10 = DcorrF10;};
    void SetDcorrFunction11(TF1* DcorrF11) {fDcorr11 = DcorrF11;};
    void SetDcorrFunction12(TF1* DcorrF12) {fDcorr12 = DcorrF12;};
    void SetDcorrFunction13(TF1* DcorrF13) {fDcorr13 = DcorrF13;};
    void SetDcorrFunction14(TF1* DcorrF14) {fDcorr14 = DcorrF14;};
    void SetDcorrFunction15(TF1* DcorrF15) {fDcorr15 = DcorrF15;};
    void SetDcorrFunction16(TF1* DcorrF16) {fDcorr16 = DcorrF16;};
    void SetDcorrFunction17(TF1* DcorrF17) {fDcorr17 = DcorrF17;};
    void SetDcorrFunction18(TF1* DcorrF18) {fDcorr18 = DcorrF18;};
    void SetDcorrFunction19(TF1* DcorrF19) {fDcorr19 = DcorrF19;};
    void SetDcorrFunction20(TF1* DcorrF20) {fDcorr20 = DcorrF20;};
    void SetDcorrFunction21(TF1* DcorrF21) {fDcorr21 = DcorrF21;};
    void SetDcorrFunction22(TF1* DcorrF22) {fDcorr22 = DcorrF22;};
   // void SetDcorrFunction(TF1* DcorrF) {fDcorr = DcorrF;};
    
    void SetDcorrFunction(TF1* DcorrF) {fDcorr = DcorrF;};
    
    void SetHCFunction(TF1* HC) {fHC = HC;};
    void SetPi0WeightMB(TGraphErrors* hWeightPi0MB) {hMCWeightPi0MB = hWeightPi0MB;};
    void SetEtaWeightMB(TGraphErrors* hWeightEtaMB) {hMCWeightEtaMB = hWeightEtaMB;};
    void SetPi0WeightEnh(TGraphErrors* hWeightPi0Enh) {hMCWeightPi0Enh = hWeightPi0Enh;};
    void SetEtaWeightEnh(TGraphErrors* hWeightEtaEnh) {hMCWeightEtaEnh = hWeightEtaEnh;};
    
    void SetTaggEffi(Bool_t fSwitch) {fCalculateTaggingEff = fSwitch;};
    void SetElecRecoEffi(Bool_t fSwitch) {fCalculateBeautyElectronTrackEff = fSwitch;};
    void SetBDMesonpTWeightCalc(Bool_t fSwitch) {fCalculateBDMesonpTWeights = fSwitch;};

    void InvMassCheckData(int itrack, AliVTrack *track, Double_t *d0z0, Int_t MagSign);
    void InvMassCheckMC(int itrack, AliVTrack *track, Double_t *d0z0, Int_t MagSign);
    void InvMassCheckMCDenom(AliVTrack *track);
    
    void GetMCTemplateWeight(TClonesArray *mcArray);
  //  void InvMassCheckMCDenomNew(AliVTrack *track);
  //  void InvMassCheckMCNew(int itrack, AliVTrack *track, Double_t *d0z0, Int_t MagSign);

    //Getters
    AliHFEpid *GetPID() const {return fPID;};
    //______________________________________________________________________
    
    //______________________________________________________________________
private:
    
    //Function to process track cuts
    Bool_t ProcessCutStep(Int_t cutStep, AliVParticle *track);
    
    //Find Mothers (Find HFE and NonHFE from MC information)
    Bool_t FindMother(Int_t mcIndex);
    
   
    // Returns the resolution Gaus which is correction to the resolution in different phi regions
    Float_t GetDCAResolMC_phi1(Float_t x);
    Float_t GetDCAResolMC_phi2(Float_t x);
    Float_t GetDCAResolMC_phi3(Float_t x);
    Float_t GetDCAResolMC_phi4(Float_t x);
    
    
    // Returns the mean Gaus which is correction to the mean in different phi regions
    Float_t GetDCAMeanMC_phi1(Float_t x);
    Float_t GetDCAMeanMC_phi2(Float_t x);
    Float_t GetDCAMeanMC_phi3(Float_t x);
    Float_t GetDCAMeanMC_phi4(Float_t x);
    
    
    
    
    //Select HFE for the reconstruction efficiency calculation
    Bool_t IsHFelectronsMC(AliVTrack *track);
    
    //Select pi0 and eta for weight calculation
    Int_t GetPi0EtaType(AliAODMCParticle *pi0eta, TClonesArray *mcArray);
    
    // ----- weights for tagging efficiency -----
    
    //Select photonic electrons
    Int_t GetElecSourceType(AliAODMCParticle *electron, TClonesArray *mcArray,Double_t &ptm);
    
    //Get first mother 
    Int_t GetPrimary(Int_t id, TClonesArray *mcArray);
    
    //Get the number of particles produced by each generator
    Bool_t GetNMCPartProduced();
    
    //Correlation cuts between TPC and SPD vertexes
    Bool_t PassCorrCuts(AliAODEvent *fAOD);
    
    // ------------------------------------------
    
    
    //Check if is heavy-flavor electron
    //Bool_t IsHFelectrons(Int_t mcIndex);
    
    //Flags for specifcs analysis
    Bool_t				fIsMC;
    Bool_t				fIsPP;
    Bool_t				fCalculateBDMesonpTWeights;
    Bool_t				fCalculateTaggingEff;
    Bool_t				fCalculateBeautyElectronTrackEff;
    
    //Weight to normalize the amount of pi and eta
    //Double_t CalculateWeight(Int_t pdg_particle, Double_t x);
    
    //Used in the function FindMother
    Bool_t				fIsHFE1;
    Bool_t				fIsHFE2;
    Bool_t				fIsNonHFE;
    Bool_t				fIsFromD;
    Bool_t				fIsFromBarionB;
	Bool_t				fIsFromMesonB;
    Bool_t				fIsFromBarionBD;
	Bool_t				fIsFromMesonBD;
    Bool_t				fIsFromPi0;
    Bool_t				fIsFromEta;
    Bool_t				fIsFromGamma;
    
    //General variables
    AliESDEvent 			*fESD;
    AliAODEvent 		   	*fAOD;				
    AliVEvent 		      	*fVevent;			
    TList       			*fOutputList;
    AliPIDResponse 			*fPidResponse;
    AliSelectNonHFE 		*fNonHFE;
    AliHFEextraCuts 		*fExtraCuts;    
    
    //For the case of AOD analysis
    Bool_t					fIsAOD;					//flag for AOD analysis
    //
    
    //Vertex selection
    Float_t					fZvtx;//!
    Float_t					kassoDCAxy;//!
    Float_t					kassoDCAz;//!
    Float_t					kAssoTPCcl;//!
    Float_t					kAssoTPCPIDcl;//!
    
    //Histograms for the analysis
    TH1F				*fVertex1;//!
    TH1F				*fNevent; //!
    TH1F				*fNeventT0; //!
    TH1F				*fNevent_3;	//!
    TH1F				*fNevent_T0b; //!
    TH1F				*fNevent_corrcut;//!
    TH1F				*fNevent_no_vertex; //!
    TH1F				*fNevent_passvertex; //!
    TH1F				*fNevent_no_vertex_2; //!
    TH1F				*fNeventAnalized;//!
    TH1F				*fCent;	//!
    TH1F				*fCent2;//!
    TGraphErrors				*hMCWeightPi0MB;
    TGraphErrors				*hMCWeightEtaMB;
    TGraphErrors				*hMCWeightPi0Enh;
    TGraphErrors				*hMCWeightEtaEnh;
    TH2F				*fTPC_p1;//!
    TH2F				*fTPC_p2;//!
    TH2F				*fTPC_p3;//!
    TH1F				*fPt_1;//!
    TH1F				*fPt_2;//!
    TH1F				*fNAnalizedTracks;//!
    TH1F				*fNAnalizedTracksHijing;//!
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
    TH2F                *fTPCnsigma_TOFnsigma1;//!
    TH2F                *fTPCnsigma_TOFnsigma2;//!
    TH2F                *fTPCnsigma_TOFnsigma3;//!
    TH2F                *fTPCnsigma_p_after_tof;//!
    TH2F                *fTPCnsigma_proton_p_after_tof;//!
    TH2F                *fTPCnsigma_p_after_tof_p;//!
    TH2F                *fTPCnsigma_p_after_tof_pion;//!
    TH2F                *fTPCnsigma_p_after_tof_k;//!
    TH2F                *fTPCnsigma_pt_after_tof;//!
    TH2F                *fTPCnsigma_p_after_tof_its;//!
    TH2F                *fTPCnsigma_pt_after_tof_its;//!
    TH1F                *fPtElec;//!
    TH1F                *fPElec;//!
    TH1F                *hPtD0;//!
    TH1F                *hPtLambdaC;//!
    TH1F				*fPtHad_f;//!
    TH1F				*fPHad_f;//!
    TH2F                *fDCAxy_pt_BeforeAllCuts; //!
    TH2F                *fDCAxy_pt_AfterTrkFiltBit;//!
    TH2F                *fDCAxy_pt_AfterMoreCuts;//!
    TH2F                *fDCAxy_pt_AfterSPDLayer;//!
    TH2F		*fDCAxy_pt_AfterTrackDCACuts;//!
    TH2F                *fDCAxy_pt_AfterAllTrackCuts;//! 
    TH2F                *fDCAxy_pt_AfterPIDCuts;//!
    TH2F                *fDCAxy_pt_AfterStepRecPrimCut;//!
    TH2F                *fDCAz_pt_had;//!
    TH2F                *fDCAxy_pt_had;//!
    TH2F                *fDCAz_pt_had_WoPID;//!
    TH2F                *fDCAxy_pt_had_WoPID;//!
    TH2F                *fDCAxy_pt_charmbef;//!
    TH2F                *fDCAxy_pt_charmaft;//!
    TH2F                *fDCAxy_pt_beautybef;//!
    TH2F                *fDCAxy_pt_beautyaft;//!
    TH2F                *fDCAxy_pt_beautybaryons;//!
    TH2F                *fDCAxy_pt_beautybaryons_corr;//!
    TH2F		*fDCAxy_pt_DstarDplusbef;//!
    TH2F		*fDCAxy_pt_Dplusbef;//!
    TH2F		*fDCAxy_pt_DstarDzerobef;//!
    TH2F		*fDCAxy_pt_Dzerobef;//!
    TH2F		*fDCAxy_pt_DstarDsbef;//!
    TH2F		*fDCAxy_pt_Dsbef;//!
    TH2F		*fDCAxy_pt_charmmesonsbef;//!
    TH2F		*fDCAxy_pt_DstarDplusAft;//!
    TH2F		*fDCAxy_pt_DplusAft;//!
    TH2F		*fDCAxy_pt_DstarDzeroAft;//!
    TH2F		*fDCAxy_pt_DzeroAft;//!
    TH2F		*fDCAxy_pt_DstarDsAft;//!
    TH2F		*fDCAxy_pt_DsAft;//!
    TH2F		*fDCAxy_pt_charmmesonsAft;//!
    TH2F		*fDCAxy_pt_Lc;//!
    TH2F		*fDCAxy_pt_charmbaryons;//!
     TH2F		*fDCAxy_pt_Lc_corr;//!
    TH2F		*fDCAxy_pt_charmbaryons_corr;//!
    TH2F                *fDCAxy_pt_MesonB_beautybef;//!
    TH2F                *fDCAxy_pt_MesonB_beautyaft;//!
    TH2F                *fDCAxy_pt_MesonBD_beautybef;//!
    TH2F                *fDCAxy_pt_MesonBD_beautyaft;//!
    TH2F                *fDCAxy_pt_BaryonB_beautybef;//!
    TH2F                *fDCAxy_pt_BaryonBD_beautybef;//!
    TH2F				*fDCAxy_pt_had_onlyDCA;//!
  /*  TH2F				*fDCAxy_pt_had_onlyDCA_phi1;//!
    TH2F				*fDCAxy_pt_had_onlyDCA_phi2;//!
    TH2F				*fDCAxy_pt_had_onlyDCA_phi3;//!
    TH2F				*fDCAxy_pt_had_onlyDCA_phi4;//!
    TH2F				*fDCAxy_pt_had_phi1_ChB;//!
    TH2F				*fDCAxy_pt_had_phi1_B;//!
    TH2F				*fDCAxy_pt_had_phi2_ChB;//!
    TH2F				*fDCAxy_pt_had_phi2_B;//!
    TH2F				*fDCAxy_pt_had_phi3_ChB;//!
    TH2F				*fDCAxy_pt_had_phi3_B;//!
    TH2F				*fDCAxy_pt_had_phi4_ChB;//!
    TH2F				*fDCAxy_pt_had_phi4_B;//!
    TH2F				*fDCAxy_pt_had_ResCorr_phi1;//!
    TH2F				*fDCAxy_pt_had_ResCorr_phi2;//!
    TH2F				*fDCAxy_pt_had_ResCorr_phi3;//!
    TH2F				*fDCAxy_pt_had_ResCorr_phi4;//!
    TH1F				*fResGausCorr_phi1;//!
    TH1F				*fResGausCorr_phi2;//!
    TH1F				*fResGausCorr_phi3;//!
    TH1F				*fResGausCorr_phi4;//!*/
    TH2F				*fDCAxy_pt_had_onlyDCA_WoPID;//!
    TH2F				*fDCAxy_pt_had_onlyDCA_Hijing;//!
    TH2F				*fDCAxy_pt_had_onlyDCA_Phytia;//!
    TH2F                *fDCAz_pt_ele;//!
    TH2F                *fDCAxy_pt_ele;//!
    TH1F                *fPtMCeta;//!
    TH1F                *hCharmMotherPt;//! pt of mothers of eletrons from mesons D
    TH1F                *hCharmMotherPt_corr;//! pt of mothers of eletrons from mesons D corrected statistically
    TH1F                *hCharmMotherPt_corr2;//! pt of mothers of eletrons from mesons D weighted
    
    TH1F		*fInvmassLS;//!
    TH1F		*fInvmassULS;//!
    TH2F		*fULSdcaBelow;//!
    TH2F		*fLSdcaBelow;//!
    TH1F		*fMCInvmassLS;//!
    TH1F		*fMCInvmassULS;//!
    TH2F		*fMCULSdcaBelow;//!
    TH2F		*fMCLSdcaBelow;//!
    
	TH2F                *hBeautyMotherPt;//!
	TH2F 		    *hDCAPtProtons;//!
	TH2F 		    *hDCAPtProtons2;//!
	TH2F 		    *hDCAPtProtons3;//!
        TH1F                *hPtElectrons;//!
        TH1F                *hPtHadrons;//!
        
        TH1F                *fBHadpT;//!
        TH1F                *fBMesonpT;//!
        TH1F                *fBMesonpT_Corr;//!
     //   TH1F                *fBMesonpTG;//!
       // TH1F                *fBMesonpTGG;//!
        TH1F                *fDHadpT;//!
        TH1F                *fDMesonpT;//!
        TH1F                *fBDHadpT;//!
        TH1F                *fD0pT;//!
        TH1F                *fLambdaCpT;//!
        
        
        TH1F                *hDCAElectrons;//!
        TH1F                *hDCAHadrons;//!
	TH2F                *hBeautyMotherPt2Daft;//!
	TH1F                *hBeautyMotherPtbef;//!
	TH1F                *hBeautyMotherPtaft;//!
	TH1F				*fPtBeautyGenerated;//!
	TH1F				*fPtBeautyReconstructedAll;//!
	TH1F				*fPtGeneratedBmesons;//!
	TH1F				*fPtBeautyReconstructedTracks;//!
	TH1F				*fPtBeautyReconstructedTracksPID;//!
	TH1F				*fPtBeautyReconstructedTracksPIDTPC;//!
	TH1F				*fPtBeautyReconstructedTracksPIDTOF;//!
        TH1F				*fMCPtElePi0FoundBef;//! for tagging efficiency
        TH1F				*fMCPtEleEtaFoundBef;//!
        TH1F				*fMCPtEleGammaFound;//!
        TH1F				*fMCPtElePi0GeneratedBef;//!
        TH1F				*fMCPtEleEtaGeneratedBef;//!
        TH1F				*fMCPtEleGammaGenerated;//!
        TH1F				*fMCMBPtElePi0FoundAft;//! for tagging efficiency
        TH1F				*fMCMBPtEleEtaFoundAft;//!
        TH1F				*fMCMBPtElePi0GeneratedAft;//!
        TH1F				*fMCMBPtEleEtaGeneratedAft;//!
        TH1F				*fMCEnhPtElePi0FoundAft;//! for tagging efficiency
        TH1F				*fMCEnhPtEleEtaFoundAft;//!
        TH1F				*fMCEnhPtElePi0GeneratedAft;//!
        TH1F				*fMCEnhPtEleEtaGeneratedAft;//!

    TH1F                *fRealInclsElecPt;//!
/*    TH1F                *fPtElePi0GeneratedAft_NewNoWeight;//!
    TH1F                *fPtElePi0EtaGeneratedAft_NewNoWeight;//!
    TH1F                *fMCMBPtElePi0EtaGeneratedAft_New;//!
    TH1F                *fMCEnhPtElePi0EtaGeneratedAft_New;//!
    TH1F                *fMCMBPtElePi0GeneratedAft_New;//!
    TH1F                *fMCEnhPtElePi0GeneratedAft_New;//!
    TH1F                *fPtEleEtaGeneratedAft_NewNoWeight;//!
    TH1F                *fMCMBPtEleEtaGeneratedAft_New;//!
    TH1F                *fMCEnhPtEleEtaGeneratedAft_New;//!
    
    TH1F                *fPtElePi0FoundAft_NewNoWeight;//!
    TH1F                *fPtElePi0EtaFoundAft_NewNoWeight;//!
    TH1F                *fMCMBPtElePi0EtaFoundAft_New;//!
    TH1F                *fMCEnhPtElePi0EtaFoundAft_New;//!
    TH1F                *fMCMBPtElePi0FoundAft_New;//!
    TH1F                *fMCEnhPtElePi0FoundAft_New;//!
    TH1F                *fPtEleEtaFoundAft_NewNoWeight;//!
    TH1F                *fMCMBPtEleEtaFoundAft_New;//!
    TH1F                *fMCEnhPtEleEtaFoundAft_New;//!
        */
    TH2F				*hCharmMotherPt_vsElecPt;//!
    TH2F				*hElecPt_vsCharmMotherPt;//!
    
    TH2F				*hCharmMotherPt_vsElecPt_corr;//!
    TH2F				*hElecPt_vsCharmMotherPt_corr;//!
    
    TGraphErrors			*fBcorr;
    TF1					*fDcorr;
    TF1					*fDcorr1;
    TF1					*fDcorr2;
    TF1					*fDcorr3;
    TF1					*fDcorr4;
    TF1					*fDcorr5;
    TF1					*fDcorr6;
    TF1					*fDcorr7;
    TF1					*fDcorr8;
    TF1					*fDcorr9;
    TF1					*fDcorr10;
    TF1					*fDcorr11;
    TF1					*fDcorr12;
    TF1					*fDcorr13;
    TF1					*fDcorr14;
    TF1					*fDcorr15;
    TF1					*fDcorr16;
    TF1					*fDcorr17;
    TF1					*fDcorr18;
    TF1					*fDcorr19;
    TF1					*fDcorr20;
    TF1					*fDcorr21;
    TF1					*fDcorr22;
    TF1					*fHC;
    //For the HFE package
    AliHFEcuts 			*fCuts;            		// Cut Collection for HFE
    AliCFManager 		*fCFM;                  // Correction Framework Manager
    AliHFEpid 			*fPID;                  // PID
    AliHFEpidQAmanager 	*fPIDqa;				// PID QA manager
    
    //Others
//    AliStack 			*fMCstack;	//!						
    Bool_t              fRejectKinkMother;	//!			
//    TParticle 			*fMCtrack;//!
//    TParticle 			*fMCtrackMother;//!
//    TParticle 			*fMCtrackGMother;//!
//    TParticle 			*fMCtrackGGMother;//!
//    TParticle 			*fMCtrackGGGMother;//!
    TClonesArray 		*fMCarray;//!
    AliAODMCHeader 		*fMCheader;  //!
    AliAODMCParticle 	*fMCparticle; //!
    AliAODMCParticle 	*fMCparticleMother;//!
    AliAODMCParticle 	*fMCparticleGMother;//!
    AliAODMCParticle 	*fMCparticleGGMother;//!
    AliAODMCParticle 	*fMCparticleGGGMother;//!
    AliMCEventHandler	*fEventHandler;//!
    AliMCEvent			*fMCevent;//!
    Float_t				ftpcPIDmincut;
    Float_t				ftpcPIDmaxcut;
    Float_t				ftofPIDmincut;
    Float_t				ftofPIDmaxcut;
    Float_t				fEtaMin;
    Float_t				fEtaMax;
    Int_t 				TOFcut;
   
	Int_t	            fNTotMCpart; //! N of total MC particles produced by generator
    Int_t               fNpureMC;//! N of particles from main generator (Hijing/Pythia)
    Int_t               fNembMCpi0; //! N > fNembMCpi0 = particles from pi0 generator
    Int_t               fNembMCeta; //! N > fNembMCeta = particles from eta generator
    
    Double_t 		fBWeight;//!
    Double_t 		fDWeight;//!
    THnSparseF           *fD0;//! DCA
    THnSparseF           *fD0_DCorr;//! DCA
    THnSparseF           *fD0_BCorr;//! DCA
    THnSparseF           *fD0Data;//! DCA data
    THnSparseF           *fD0HC;//! DCA HC
    //______________________________________________________________________
    
    AliAnalysisHFEppTPCTOFBeauty5TeVNew(const AliAnalysisHFEppTPCTOFBeauty5TeVNew&); 			// not implemented
    AliAnalysisHFEppTPCTOFBeauty5TeVNew& operator=(const AliAnalysisHFEppTPCTOFBeauty5TeVNew&); 		// not implemented
    
    ClassDef(AliAnalysisHFEppTPCTOFBeauty5TeVNew, 1); 								// example of analysis
    //______________________________________________________________________
};

///_________________________________________________________________________________________________

#endif
