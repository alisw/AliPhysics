#ifndef AliAnalysisTaskBtoElecPbPbTPCTOF_cxx
#define AliAnalysisTaskBtoElecPbPbTPCTOF_cxx

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////////////
//                                                     		        //
//	 Task for Beauty-hadron decay electrons in Pb-Pb collisions		//
//																	 														  //
//		v1.0								                           						  //
//                                                                //
//	    Authors 						           	                          //
//		Jonghan Park (jonghan@cern.ch)												      //
//                                                                //
////////////////////////////////////////////////////////////////////

class AliAODEvent;
class AliPIDResponse;
class AliAODMCHeader;
class TClonesArray;
class AliHFEextraCuts;
class AliAODMCParticle;
class AliGenEventHeader;
class AliAODv0;
class AliAODv0KineCuts;
class TRandom3;
class THnSparse;
class AliHFEV0taginfo;

//______________________________________________________________________
//Library
#include "AliAnalysisTaskSE.h"
//______________________________________________________________________

//______________________________________________________________________
class AliAnalysisTaskBtoElecPbPbTPCTOF : public AliAnalysisTaskSE 
{
	//______________________________________________________________________
	public:
		enum SourceType {
			kDirectCharm=1,     // electrons from primary charmed hadrons 
			kDirectBeauty=2,    // electrons from primary beauty hadrons
			kBeautyCharm=3,     // electrons from charmed hadrons decaying from the beauty hadrons (B->D->e)
			kGamma=4,           // should be obsolete -> please let me know if you see something! 
			kPi0=5,             // electrons from pi0 Dalitz
			kEta=6,             // electrons from eta Dalitz
			kOmega=7,           // electrons from omega decay (Dalitz and di-electrons)
			kPhi=8,            // electrons from phi decay (di-electron)
			kEtaPrime=9,       // electrons from eta prime decay (Dalitz and 2charged-pions&di-electrons) 
			kRho0=10,           // electrons from rho decay (di-electron)
			kK0s2P=11,          // electrons from secondary pions from K0s
			kK0l2P=12,          // electrons from secondary pions from K0l
			kLamda2P=13,        // electrons from secondary pions from Lamda
			kSigma2P=14,        // electrons from secondary pions from Sigma
			kK2P=15,            // electrons from secondary pions from K  
			kElse=16,            // all the other sources which was not in this enumeration
			kMisID=17,           // not the electrons (hadrons)
			kGammaPi0=18,       // electrons from photon conversion where the photon originated from pi0
			kGammaEta=19,       // electrons from photon conversion where the photon originated from eta
			kGammaOmega=20,     // electrons from photon conversion where the photon originated from omega
			kGammaPhi=21,       // electrons from photon conversion where the photon originated from phi
			kGammaEtaPrime=22,  // electrons from photon conversion where the photon originated from eta prime
			kGammaRho0=23,      // electrons from photon conversion where the photon originated from rho
			kGammaK0s2P=24,     // electrons from photon conversion where the photon originated from secondary pion from K0s
			kGammaK0l2P=25,     // electrons from photon conversion where the photon originated from secondary pion from K0l
			kGammaLamda2P=26,   // electrons from photon conversion where the photon originated from secondary pion from lamda
			kGammaSigma2P=27,   // electrons from photon conversion where the photon originated from secondary pion from sigma
			kGammaK2P=28,        // electrons from photon conversion where the photon originated from secondary pion from K
			kJpsi=29,           // electrons from primary J/psi decay
			kB2Jpsi=30,         // electrons from J/psi decay where the J/psi originated from the beauty hadrons
			kKe3=31,            // Ke3 electrons
			kGammaB2M=32,       // electrons from photon conversion from meson where the meson originated from the beauty hadrons
			kGammaD2M=33,       // electrons from photon conversion from meson where the meson originated from the charm hadrons
			kGammaM2M=34,       // electrons from photon conversion from the light meson decay where the light meson originated from other light meson 
			kB2M=35,            // electrons from the meson where the meson originated from the beauty hadrons 
			kD2M=36,            // electrons from the meson where the meson originated from the charm hadrons
			kM2M=37,            // electrons from the light meson decay where the light meson originated from other light meson 
			kK0L=38,            // K0L -> e +X
			kPromptD0=39,
			kNonPromptD=40,
			kPromptB=41,
			kPromptLc=42,
			kDirectGamma=43,
			kDalitzGamma=44,
			kBaryonGamma=45
		};



	AliAnalysisTaskBtoElecPbPbTPCTOF();
	AliAnalysisTaskBtoElecPbPbTPCTOF(const char *name);
	virtual ~AliAnalysisTaskBtoElecPbPbTPCTOF();
  
	virtual void   UserCreateOutputObjects();
	virtual void   UserExec(Option_t *option);
	virtual void   Terminate(Option_t *);

	//Setters
	
	void SetMCanalysis() { fIsMC = kTRUE; }
	void SetCentrality(double centMin, double centMax) { fCentMin = centMin; fCentMax = centMax; }
	void SetMinTPCNcls(int minTPCncls) { fMinTPCNcls = minTPCncls; }
	void SetMinTPCNclsPID(int minTPCnclsPID) { fMinTPCNclsPID = minTPCnclsPID; }
	void SetMaxTPCchi2(double maxTPCchi2) { fMaxTPCchi2 = maxTPCchi2; }
	void SetMinTPCclsRatio(double minTPCclsratio) { fMinTPCclsRatio = minTPCclsratio; }
	void SetMinITSNcls(int minITSncls) { fMinITSNcls = minITSncls; }
	void SetMaxITSclsFrac(double maxITSclsfrac) { fMaxITSclsFrac = maxITSclsfrac; }
	void SetMaxITSchi2(double maxITSchi2) { fMaxITSchi2 = maxITSchi2; }
	void SetITSlayer(TString itsLayer) { fITSlayer = itsLayer; }
	void SetEtaRange(double eta) { fEta = eta; }
	void SetPtRange(double ptmin, double ptmax) { fMinPt = ptmin; fMaxPt = ptmax; }
	void SetDCARange(double dcaxy, double dcaz) { fDCAxy = dcaxy; fDCAz = dcaz; }
	void SetPIDCuts(double tpcPIDlow, double tpcPIDhigh, double tofPID) { fTPCnsigmaLow = tpcPIDlow; fTPCnsigmaHigh = tpcPIDhigh; fTOFnsigma = tofPID; }
	
	int GetElecSource(const AliAODMCParticle * const mcpart, double &mpt, int &mpdg, double &momTime);
	int GetHeavyFlavours(const AliAODMCParticle * const mcpart, double &mpt, double &meta);
	int GetGammaPt(const AliAODMCParticle * const mcpart, double &mpt);
	void SelectV0Pions(AliAODEvent *evt);

	// Setter for B corr (central)
	void SetBcorrFtn(TF1* BcorrF) { fBmesonCorr = BcorrF; }
	void SetBcorrFtn1(TF1* BcorrF1) { fBmesonCorr1 = BcorrF1; }

	// Setter for D corr (central)
	void SetDcorrFtn(TF1* DcorrF) { fDmesonCorr = DcorrF; }
	void SetDcorrFtn1(TF1* DcorrF1) { fDmesonCorr1 = DcorrF1; }
	void SetDcorrFtn2(TF1* DcorrF2) { fDmesonCorr2 = DcorrF2; }
	void SetDcorrFtn3(TF1* DcorrF3) { fDmesonCorr3 = DcorrF3; }
	void SetDcorrFtn4(TF1* DcorrF4) { fDmesonCorr4 = DcorrF4; }
	void SetDcorrFtn5(TF1* DcorrF5) { fDmesonCorr5 = DcorrF5; }
	void SetDcorrFtn6(TF1* DcorrF6) { fDmesonCorr6 = DcorrF6; }
	void SetDcorrFtn7(TF1* DcorrF7) { fDmesonCorr7 = DcorrF7; }
	void SetDcorrFtn8(TF1* DcorrF8) { fDmesonCorr8 = DcorrF8; }
	void SetDcorrFtn9(TF1* DcorrF9) { fDmesonCorr9 = DcorrF9; }
	void SetDcorrFtn10(TF1* DcorrF10) { fDmesonCorr10 = DcorrF10; }
	void SetDcorrFtn11(TF1* DcorrF11) { fDmesonCorr11 = DcorrF11; }
	void SetDcorrFtn12(TF1* DcorrF12) { fDmesonCorr12 = DcorrF12; }

	// Setter for Lc corr (central)
	void SetLccorrFtn(TF1* LccorrF) {fLcCorr = LccorrF;};

	//______________________________________________________________________
	private:
	
		//Correlation cuts between TPC and SPD vertexes
		bool PassEventCuts(AliAODEvent *event);
		bool PassTrackCuts(AliAODTrack *track);
		bool PassV0PionMinCuts(AliAODTrack *track);
		bool PassV0PionMaxCuts(AliAODTrack *track);

	
		//Flags for analysis
		bool				fIsMC;
		float				fCentrality;
		float				fCentralityCalib;
		double			fCentMin;
		double			fCentMax;
		int					fMinTPCNcls;
		int					fMinTPCNclsPID;
		double			fMaxTPCchi2;
		double			fMinTPCclsRatio;
		int					fMinITSNcls;
		double			fMaxITSclsFrac;
		double			fMaxITSchi2;
		TString			fITSlayer;
		double			fEta;
		double			fMinPt;
		double			fMaxPt;
		double			fDCAxy;
		double			fDCAz;
		double			fTPCnsigmaLow;
		double			fTPCnsigmaHigh;
		double			fTOFnsigma;
    
		//General variables
		AliAODEvent					*fAOD;
		TList								*fOutputList;
		AliPIDResponse			*fPidResponse;
		AliAODMCHeader			*fAODMCHeader;
		TClonesArray				*fAODArrayMCInfo;
		AliHFEextraCuts			*fExtraCuts;
		AliAODMCParticle		*fAODMCParticle;
		AliAODv0						*fAODv0;
		AliAODv0KineCuts		*fAODV0Cuts;
		AliHFEV0taginfo			*fV0Tagger;
		
		// Event
		TH1F				*hCent_nocut;
		TH1F				*hCent_nocut2;
		TH1F				*hCent_cut;
		TH1F				*hNrEvents;

		// Track
		TH1F				*hFilterMask;
		TH1F				*hTPCNcls;
		TH1F				*hTPCclsPID;
		TH1F				*hTPCchi2;
		TH1F				*hTPCclsRatio;
		TH1F				*hITSNcls;
		TH1F				*hITSclsFrac;
		TH1F				*hITSchi2;
		TH1F				*hITSlayer;
		TH1F				*hDCAxy;
		TH1F				*hDCAz;
		TH1F				*hPt;
		TH1F				*hEta;
		TH1F				*hPhi;

		// PID
		TH2F				*hITSnsigma;
		TH2F				*hITSnsigmaTOFcut;
		TH2F				*hITSnsigmaQA;
		TH2F				*hTPCnsigma;
		TH2F				*hTPCnsigmaTOFcut;
		TH2F				*hTPCnsigmaITSTOFcut;
		TH2F				*hTPCnsigmaQA;
		TH2F				*hTPCnsigmaPiQA;
		TH2F				*hTOFnsigma;
		TH2F				*hTOFnsigmaQA;
		
		TH2F				*hV0ElecTOFnsigmaDeno;
		TH2F				*hV0ElecTOFnsigmaNume;
		TH2F				*hV0ElecTPCnsigmaDeno;
		TH2F				*hV0ElecTPCnsigmaNume;
		
		// Efficiency
		TH1F				*hGenBtoElecPt;
		TH1F				*hRecBtoElecPt_nocut;
		TH1F				*hRecBtoElecPt_track;
		TH1F				*hRecBtoElecPt_tof;
		TH1F				*hRecBtoElecPt_tpc;

		// D meson pt correction
		TH1F				*hD0Pt;
		TH1F				*hD0PtCorr;
		
		// B meson pt correction
		TH1F				*hBhadronPt;
		TH1F				*hBhadronPtCorr;

		// Lc pt correction
		TH1F				*hLcPt;
		TH1F				*hLcPtCorr;

		// IP
		TH2F				*dcaTrack;
		TH2F				*dcaDmeson;
		TH2F				*dcaDmesonCorr;
		TH2F				*dcaDzero;
		TH2F				*dcaDplus;
		TH2F				*dcaDsplus;
		TH2F				*dcaLc;
		TH2F				*dcaBeauty;
		TH2F				*dcaBeautyCorr;
		TH2F				*dcaDalitz;
		TH2F				*dcaConv;
		TH2F				*dcaPion1;
		TH2F				*dcaPion2;

		// Conversion multiplicity correction
		TH2F				*hProdR_Pt;
		TH3F				*hV0PionMinCut;
		TH3F				*hV0PionMaxCut;
		TH3F				*hV0PionMult;

		// B corr
		TF1					*fBmesonCorr;
		TF1					*fBmesonCorr1;

		// D corr
		TF1         *fDmesonCorr;
		TF1         *fDmesonCorr1;
		TF1         *fDmesonCorr2;
		TF1         *fDmesonCorr3;
		TF1         *fDmesonCorr4;
		TF1         *fDmesonCorr5;
		TF1         *fDmesonCorr6;
		TF1         *fDmesonCorr7;
		TF1         *fDmesonCorr8;
		TF1         *fDmesonCorr9;
		TF1         *fDmesonCorr10;
		TF1         *fDmesonCorr11;
		TF1         *fDmesonCorr12;
	
		// Lc corr
		TF1         *fLcCorr;

		// tau corr
		TF1					*fD0TauWeight;
		TF1					*fDpTauWeight;
		TF1					*fDsTauWeight;
		TF1					*fB0TauWeight;
		TF1					*fBpTauWeight;
		TF1					*fBsTauWeight;
	
		TRandom3 		*fRnd;

		//______________________________________________________________________

		AliAnalysisTaskBtoElecPbPbTPCTOF(const AliAnalysisTaskBtoElecPbPbTPCTOF&); 			// not implemented
		AliAnalysisTaskBtoElecPbPbTPCTOF& operator=(const AliAnalysisTaskBtoElecPbPbTPCTOF&); 		// not implemented
  
		ClassDef(AliAnalysisTaskBtoElecPbPbTPCTOF, 1); 								// example of analysis
		//______________________________________________________________________
};

///_________________________________________________________________________________________________

#endif
