#ifndef AliAnalysisHFETPCTOFBeauty_cxx
#define AliAnalysisHFETPCTOFBeauty_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////////
//                                                                    //
//      Task for Beauty analysis in Pb-Pb collisions   				  //
//      															  //
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
class AliGenEventHeader;

//______________________________________________________________________
//Library
#include "AliAnalysisTaskSE.h"
#include "AliHFEpid.h"
#include "AliLog.h"
//#include "AliMultSelection.h"
//______________________________________________________________________

//______________________________________________________________________
class AliAnalysisHFETPCTOFBeauty : public AliAnalysisTaskSE
{
    //______________________________________________________________________
public:
    
    enum HijingOr {kHijing,kPhytia,kpi0,keta};
    enum ESourceType {kNoMotherE, kPi0NoFeedDown, kEtaNoFeedDown, kGPi0NoFeedDown, kGEtaNoFeedDown, kDirectGamma, kOthersE};
    enum pi0etaType {kNoMother, kNoFeedDown, kNoIsPrimary, kLightMesons, kBeauty, kCharm};

    AliAnalysisHFETPCTOFBeauty();
    AliAnalysisHFETPCTOFBeauty(const char *name);
    virtual ~AliAnalysisHFETPCTOFBeauty();
    
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
    AliAODEvent 		   	*fAOD;				
    AliVEvent 		      	*fVevent;			
    TList       			*fOutputList;
    AliPIDResponse 			*fPidResponse;
    AliSelectNonHFE 		*fNonHFE;
    
    
    //For the case of AOD analysis
    Bool_t					fIsAOD;					//flag for AOD analysis
    //
    
    //Vertex selection
    Float_t					fZvtx;
  
    
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
    TH2F                *fTPCnsigma_TOFnsigma1;//!
    TH2F                *fTPCnsigma_TOFnsigma2;//!
    TH2F                *fTPCnsigma_TOFnsigma3;//!
    TH2F                *fTPCnsigma_p_after_tof;//!
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
    TH2F                *fDCAz_pt_had;//!
    TH2F                *fDCAxy_pt_had;//!
    TH2F                *fDCAz_pt_ele;//!
    TH2F                *fDCAxy_pt_ele;//!
    TH1F                *fPtMCeta;//!
    TH1F                *hCharmMotherPt;//! pt of mothers of eletrons from mesons D
    TH1F                *hCharmMotherPt_corr;//! pt of mothers of eletrons from mesons D corrected statistically
    TH1F                *hCharmMotherPt_corr2;//! pt of mothers of eletrons from mesons D weighted
	TH1F                *hCharmMotherPt_corr3;
	TH1F                *hCharmMotherPt_corr4;
	
	TH1F				*hCharmMotherPt1_corr;
	TH1F				*hCharmMotherPt2_corr;
	TH1F				*hCharmMotherPt3_corr;
	TH1F				*hCharmMotherPt4_corr;
	TH1F				*hCharmMotherPt5_corr;
	TH1F				*hCharmMotherPt6_corr;
	TH1F				*hCharmMotherPt7_corr;
	TH1F				*hCharmMotherPt8_corr;
	TH1F				*hCharmMotherPt9_corr;
	TH1F				*hCharmMotherPt10_corr;
	TH1F				*hCharmMotherPt11_corr;
	TH1F				*hCharmMotherPt12_corr;
	TH1F				*hCharmMotherPt13_corr;
	TH1F				*hCharmMotherPt14_corr;
	TH1F				*hCharmMotherPt15_corr;
	
	TH1F                *hBeautyMotherPt;//!
	TH1F				*fPtBeautyGenerated;
	TH1F				*fPtBeautyReconstructedTracks;
	TH1F				*fPtBeautyReconstructedTracksPID;
    
    TH2F				*hCharmMotherPt_vsElecPt;
    TH2F				*hElecPt_vsCharmMotherPt;
    
    
    //For the HFE package
    AliHFEcuts 			*fCuts;            		// Cut Collection for HFE
    AliCFManager 		*fCFM;                  		// Correction Framework Manager
    AliHFEpid 			*fPID;                  		// PID
    AliHFEpidQAmanager 	*fPIDqa;					// PID QA manager
    
    //Others
    AliStack 			*fMCstack;						//
    Bool_t              fRejectKinkMother;				//
    TParticle 			*fMCtrack;
    TParticle 			*fMCtrackMother;
    TParticle 			*fMCtrackGMother;
    TParticle 			*fMCtrackGGMother;
    TParticle 			*fMCtrackGGGMother;
    TClonesArray 		*fMCarray;
    AliAODMCHeader 		*fMCheader;  //!
    AliAODMCParticle 	*fMCparticle; //!
    AliAODMCParticle 	*fMCparticleMother;
    AliAODMCParticle 	*fMCparticleGMother;
    AliAODMCParticle 	*fMCparticleGGMother;
    AliAODMCParticle 	*fMCparticleGGGMother;
    AliMCEventHandler	*fEventHandler;
    AliMCEvent			*fMCevent;
    Float_t				ftpcPIDmincut;
    Float_t				ftpcPIDmaxcut;
    Float_t				ftofPIDmincut;
    Float_t				ftofPIDmaxcut;
    Float_t				fEtaMin;
    Float_t				fEtaMax;
    
	Int_t	            fNTotMCpart; //! N of total MC particles produced by generator
    Int_t               fNpureMC;//! N of particles from main generator (Hijing/Pythia)
    Int_t               fNembMCpi0; //! N > fNembMCpi0 = particles from pi0 generator
    Int_t               fNembMCeta; //! N > fNembMCeta = particles from eta generator
    
 
    THnSparseF           *fD0;//! DCA
    THnSparseF           *fD0Data;//! DCA data
    //______________________________________________________________________
    
    AliAnalysisHFETPCTOFBeauty(const AliAnalysisHFETPCTOFBeauty&); 			// not implemented
    AliAnalysisHFETPCTOFBeauty& operator=(const AliAnalysisHFETPCTOFBeauty&); 		// not implemented
    
    ClassDef(AliAnalysisHFETPCTOFBeauty, 1); 								// example of analysis
    //______________________________________________________________________
};

///_________________________________________________________________________________________________

#endif
