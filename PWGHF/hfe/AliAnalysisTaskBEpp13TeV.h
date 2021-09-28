/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskBEpp13TeV_H
#define AliAnalysisTaskBEpp13TeV_H

#include "AliAnalysisTaskSE.h"

class AliPIDResponse;
class AliHFEextraCuts;
class AliAODMCHeader;
class TClonesArray;
class AliAODMCParticle;
class TClonesArray;
class AliAODMCParticle;
class TRandom3;

class AliAnalysisTaskBEpp13TeV : public AliAnalysisTaskSE  
{
  public:
                            AliAnalysisTaskBEpp13TeV();
                            AliAnalysisTaskBEpp13TeV(const char *name);
	virtual                 ~AliAnalysisTaskBEpp13TeV();

    virtual void            UserCreateOutputObjects();
    virtual void            UserExec(Option_t* option);
    virtual void            Terminate(Option_t* option);
	
	enum SourceType{
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
	  kK0L=38,            // K0L -> e +X
	  kPromptD0=39,
	  kNonPromptD=40,
	  kPromptB=41,
	  kPromptLc=42,
	  kDirectGamma=43,
	  kDalitzGamma=44,
	  kBaryonGamma=45
	};

	//Setter
	void SetMCanalysis() { fIsMC = true; }
	void SetMinTPCnCrossedRow(int minTPCncls) { fMinTPCnCrossedRow = minTPCncls; }
	void SetMinTPCNclsPID(int minTPCnclsPID) { fMinTPCNclsPID = minTPCnclsPID; }
	void SetMaxTPCchi2(double maxTPCchi2) { fMaxTPCchi2 = maxTPCchi2; }
	void SetMinTPCclsRatio(double minTPCclsratio) { fMinTPCclsRatio = minTPCclsratio; }
	void SetMinITSNcls(int minITSncls) { fMinITSNcls = minITSncls; }
	void SetITSlayer(int itsLayer) { fITSlayer = itsLayer; }
	void SetPIDCuts(double tpcPIDlow, double tpcPIDhigh, double tofPID) { fTPCnsigmaLow = tpcPIDlow;fTPCnsigmaHigh = tpcPIDhigh; fTOFnsigma = tofPID; }


	// Setter for B corr (central)
	void SetBcorrCentLow(TF1* BcorrFcentL) { fBmesonCorrCentLow = BcorrFcentL; }
	void SetBcorrCentHigh(TF1* BcorrFcentH) { fBmesonCorrCentHigh = BcorrFcentH; }
	void SetBcorrMinLow(TF1* BcorrFminL) { fBmesonCorrMinLow = BcorrFminL; }
	void SetBcorrMinHigh(TF1* BcorrFminH) { fBmesonCorrMinHigh = BcorrFminH; }
	void SetBcorrMaxLow(TF1* BcorrFmaxL) { fBmesonCorrMaxLow = BcorrFmaxL; }
	void SetBcorrMaxHigh(TF1* BcorrFmaxH) { fBmesonCorrMaxHigh = BcorrFmaxH; }

	// Setter for D corr (central)
	void SetDcorrFtn(TF1* DcorrF) { fDmesonCorr = DcorrF; }
	void SetDcorrFtnVar1(TF1* DcorrFVar1) { fDmesonCorrVar1 = DcorrFVar1; }
	void SetDcorrFtnVar2(TF1* DcorrFVar2) { fDmesonCorrVar2 = DcorrFVar2; }

	// Setter for Lc corr (central)
	void SetLccorrFtn(TF1* LccorrF) {fLcCorr = LccorrF;};

  private:

	bool PassEventCuts(AliAODEvent *event);
	bool PassPileUpEvent(AliAODEvent *event);
	bool PassTrackCuts(AliAODTrack *track);
	int GetElecSource(const AliAODMCParticle * const mcpart, double &mpt, int &mpdg);
	int GetHeavyFlavours(const AliAODMCParticle * const mcpart, double &hfpt, double &hfeta);

	bool		fIsMC;
	int			fMinTPCnCrossedRow;
	int			fMinTPCNclsPID;
	double		fMaxTPCchi2;
	double		fMinTPCclsRatio;
	int			fMinITSNcls;
	int			fITSlayer;
	double		fTPCnsigmaLow;
	double		fTPCnsigmaHigh;
	double		fTOFnsigma;

    AliAODEvent				*fAOD;           //! input event
    TList					*fOutputList;    //! output list
    AliPIDResponse			*fPIDResponse;
	AliHFEextraCuts			*fExtraCuts;
	//AliAODMCHeader			*fAODMCHeader;
	TClonesArray			*fAODArrayMCInfo;
	AliAODMCParticle		*fAODMCParticle;
	
	TH1F	*hNrEvents;
	
	// track cut QA
	TH1F	*hFilterMask;
	TH1F	*hTPCnCrossedRow;
	TH1F	*hTPCclsPID;
	TH1F	*hTPCchi2;
	TH1F	*hTPCclsRatio;
	TH1F	*hITSNcls;
	TH1F	*hITSlayer;
	TH1F	*hDCAxy;
	TH1F	*hDCAz;
	TH1F	*hPt;
	TH1F	*hEta;
	TH1F	*hPhi;

	// generated level
	TH1F	*hBhadronPt;
	TH1F	*hBhadronPtCorr;
	TH1F	*hD0Pt;
	TH1F	*hD0PtCorr;
	TH1F	*hLcPt;

	TH1F	*hGenBePt;
	TH1F	*hRecBePt_track;
	TH1F	*hRecBePt_tof;
	TH1F	*hRecBePt_tpc;

	// pid cut
	TH2F	*hTPCnsigma;
	TH2F	*hTPCnsigmaTOFcut;
	TH2F	*hTPCnsigmaTOFcutPt;
	TH2F	*hTPCnsigmaQA;
	TH2F	*hTPCnsigmaPiQA;
	TH2F	*hTOFnsigma;
	TH2F	*hTOFnsigmaQA;

	// dca
	TH2F	*dcaTrack;
	TH2F	*dcaPion;
	TH2F	*dcaBeauty;
	TH2F	*dcaBeautyCorr;
	TH2F	*dcaBeautyCorrVar1;
	TH2F	*dcaBeautyCorrVar2;
	TH2F	*dcaBzero;
	TH2F	*dcaBplus;
	TH2F	*dcaBszero;
	TH2F	*dcaLb;
	TH2F	*DelecVsDmother;
	TH2F	*dcaDmeson;
	TH2F	*dcaDmesonCorr;
	TH2F	*dcaDmesonCorrVar1;
	TH2F	*dcaDmesonCorrVar2;
	TH2F	*dcaDzero;
	TH2F	*dcaDplus;
	TH2F	*dcaDsplus;
	TH2F	*dcaLc;
	TH2F	*dcaDalitz;
	TH2F	*dcaConv;
	
	// B corr
    TF1         *fBmesonCorrCentLow;
    TF1         *fBmesonCorrCentHigh;
    TF1         *fBmesonCorrMinLow;
    TF1         *fBmesonCorrMinHigh;
    TF1         *fBmesonCorrMaxLow;
    TF1         *fBmesonCorrMaxHigh;

    // D corr
    TF1         *fDmesonCorr;
    TF1         *fDmesonCorrVar1;
    TF1         *fDmesonCorrVar2;

    // Lc corr
    TF1         *fLcCorr;

    TRandom3    *fRnd;
	
    AliAnalysisTaskBEpp13TeV(const AliAnalysisTaskBEpp13TeV&); // not implemented
    AliAnalysisTaskBEpp13TeV& operator=(const AliAnalysisTaskBEpp13TeV&); // not implemented

    ClassDef(AliAnalysisTaskBEpp13TeV, 1);
};

#endif
