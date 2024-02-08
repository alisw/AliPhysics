/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskWHMult_H
#define AliAnalysisTaskWHMult_H

#include "AliAnalysisTaskSE.h"

class AliPIDResponse;
class AliMultSelection;
class AliVCluster;
class AliAODMCHeader;
class AliAODMCParticle;

class AliAnalysisTaskWHMult : public AliAnalysisTaskSE  
{
    public:
                                AliAnalysisTaskWHMult();
                                AliAnalysisTaskWHMult(const char *name);
        virtual                 ~AliAnalysisTaskWHMult();

        virtual void            UserCreateOutputObjects();
        virtual void            UserExec(Option_t* option);
        virtual void            FindMother(AliAODMCParticle* part, int &label, int &pid, double &ptmom);        //! find mother command
	virtual void		FindFather(AliAODMCParticle* part, int &label, int &pid, double &ptmom);	//! find father command
	virtual Double_t	GetCorrectedNtrackletsD(TTree* tree, Double_t uncorrectedNacc, Int_t BinOfvertexZ, AliAODEvent* fAOD);
	virtual Int_t		CountNch();

	Bool_t                  GetEMCalTriggerEG1() { return fEMCEG1; };
	void                    SetEG1(Bool_t flagEG1) { fEMCEG1= flagEG1;};

	void			SetTrackEff(Int_t ITS, Int_t CrossedRows, Int_t TPC) {NITSClust = ITS, NCrossedRows = CrossedRows, NTPCClust = TPC;};
	void			SetNSigma(Double_t min, Double_t max) {NSigmaMin = min, NSigmaMax = max;};
	void			SetM02(Double_t min, Double_t max) {M02Min = min, M02Max = max;};
	void			SetEoP(Double_t min, Double_t max) {EoPMin = min, EoPMin = max;};
	void			SetEiso(Double_t maxEiso) {EisoMax = maxEiso;};
	void			SetNtrk(Int_t maxNtrk) {NtrkMax = maxNtrk;};
        
        virtual void            Terminate(Option_t* option);

    private:
        AliAODEvent*            fAOD;				//! input event
	TClonesArray*		fTracks_tender;			//! Tender tracks (EMCal Correction framework)
	TClonesArray*		fCaloClusters_tender;		//! Tender cluster (EMCal Correction framework)
	AliVEvent*              fVevent;			//! event object
	AliPIDResponse*		fPIDResponse;			//! pid response object
	AliAODMCParticle*	fMCTrackpart;			//! MC track particle
	TClonesArray*		fMCarray;			//! MC array
	AliAODMCParticle*	fMCpart;			//! MC particle for MC particle loop
	AliAODMCParticle*	fMCparticle;			//! MC particle for Nch
        TList*                  fOutputList;    		//! output list
        TTree*			tree;				//! ref. for mult correction

	Int_t NITSClust, NCrossedRows, NTPCClust;
	Double_t NSigmaMin, NSigmaMax;
	Double_t M02Min, M02Max;
	Double_t EoPMin, EoPMax;
	Double_t EisoMax;
	Int_t NtrkMax;														//Object List Number//

	TH1F*			fNevents;			//! hist of event cut performance						//0
	TH1F*                   pVertex_all;			//! additional histgram in Step3						//1
	TH1F*                   pVertex;			//! after event selection							//2
	TH2F*			EtavsPhi;			//! hist of pseudorapidity and azimuthal distribution with cutVer		//3
	TH2F*			TPCSig;				//! hist of TPC signal (dE/dx) and p/z						//4
	TH2F*			fzvtx_Ntrkl;			//! hist of vertexZ vs tracklets(corresponding to multiplicity)			//5
	TH2F*			fzvtx_Ntrkl_cal;		//! hist of vertexZ vs tracklets after calibration				//6

	TH2F*			fHistNsigmaP;			//! histogram of nsigma vs P							//7
	TH2F*			fPtEoverPE;			//! histogram of Pt vs E/P of electron						//8
        TH2F*			fPtEoverPH;			//! histogram of Pt vs E/P of hadron						//9
	TH1F*			fREisolation;			//! histogram of Eiso								//10
	TH2F*			fdPhi_trkW_Pt;			//! histogram of DeltaPhi w/ e<-W cut vs Pt					//11
	TH2F*			fdPhi_trkHF_Pt;			//! histogram of DeltaPhi w/ e<-HF cut vs Pt					//12
	TH1F*			fHistPt_We;			//! histogram of Pt (W candidate)						//13
	TH1F*			fHistPt_HFe;			//! histogram of Pt (HF candidate)						//14
	TH2F*			fNtrkl_PtOfMaxTrk_W_m;		//! histogram of Ntrkl vs max pT track/e<-W non UE sub.				//15
	TH2F*			fNtrkl_PtOfMaxTrk_W70_m;	//! histogram of Ntrkl vs max pT track/e<-W 70% UE				//16
	TH2F*			fNtrkl_PtOfMaxTrk_W80_m;	//! histogram of Ntrkl vs max pT track/e<-W 80% UE				//17
	TH2F*			fNtrkl_PtOfMaxTrk_W85_m;	//! histogram of Ntrkl vs max pT track/e<-W 85% UE				//18
	TH2F*			fNtrkl_PtOfMaxTrk_W90_m;	//! histogram of Ntrkl vs max pT track/e<-W 90% UE				//19
	TH2F*			fNtrkl_PtOfMaxTrk_W93_m;	//! histogram of Ntrkl vs max pT track/e<-W 93% UE				//20
	TH2F*			fNtrkl_PtOfMaxTrk_W95_m;	//! histogram of Ntrkl vs max pT track/e<-W 95% UE				//21
	TH2F*			fNtrkl_PtOfTrks_W_m;		//! histogram of Ntrkl vs sum of pT track/e<-W					//22
	TH2F*			fNtrkl_PtOfTrks_W70_m;		//! histogram of Ntrkl vs sum of pT track/e<-W 70% UE				//23
	TH2F*			fNtrkl_PtOfTrks_W80_m;		//! histogram of Ntrkl vs sum of pT track/e<-W 80% UE				//24
	TH2F*			fNtrkl_PtOfTrks_W85_m;		//! histogram of Ntrkl vs sum of pT track/e<-W 85% UE				//25
	TH2F*			fNtrkl_PtOfTrks_W90_m;		//! histogram of Ntrkl vs sum of pT track/e<-W 90% UE				//26
	TH2F*			fNtrkl_PtOfTrks_W93_m;		//! histogram of Ntrkl vs sum of pT track/e<-W 93% UE				//27
	TH2F*			fNtrkl_PtOfTrks_W95_m;		//! histogram of Ntrkl vs sum of pT track/e<-W 95% UE				//28
	TH2F*			fNtrkl_PtOfMaxTrk_HF_m;		//! histogram of Ntrkl vs max pT track/e<-HF non UE sub.			//29
	TH2F*			fNtrkl_PtOfMaxTrk_HF70_m;	//! histogram of Ntrkl vs max pT track/e<-HF 70% UE				//30
	TH2F*			fNtrkl_PtOfMaxTrk_HF80_m;	//! histogram of Ntrkl vs max pT track/e<-HF 80% UE				//31
	TH2F*			fNtrkl_PtOfMaxTrk_HF85_m;	//! histogram of Ntrkl vs max pT track/e<-HF 85% UE				//32
	TH2F*			fNtrkl_PtOfMaxTrk_HF90_m;	//! histogram of Ntrkl vs max pT track/e<-HF 90% UE				//33
	TH2F*			fNtrkl_PtOfMaxTrk_HF93_m;	//! histogram of Ntrkl vs max pT track/e<-HF 93% UE				//34
	TH2F*			fNtrkl_PtOfMaxTrk_HF95_m;	//! histogram of Ntrkl vs max pT track/e<-HF 95% UE				//35
	TH2F*			fNtrkl_PtOfTrks_HF_m;		//! histogram of Ntrkl vs sum of pT track/e<-HF					//36
	TH2F*			fNtrkl_PtOfTrks_HF70_m;		//! histogram of Ntrkl vs sum of pT track/e<-HF 70% UE				//37
	TH2F*			fNtrkl_PtOfTrks_HF80_m;		//! histogram of Ntrkl vs sum of pT track/e<-HF 80% UE				//38
	TH2F*			fNtrkl_PtOfTrks_HF85_m;		//! histogram of Ntrkl vs sum of pT track/e<-HF 85% UE				//39
	TH2F*			fNtrkl_PtOfTrks_HF90_m;		//! histogram of Ntrkl vs sum of pT track/e<-HF 90% UE				//40
	TH2F*			fNtrkl_PtOfTrks_HF93_m;		//! histogram of Ntrkl vs sum of pT track/e<-HF 93% UE				//41
	TH2F*			fNtrkl_PtOfTrks_HF95_m;		//! histogram of Ntrkl vs sum of pT track/e<-HF 95% UE				//42
	TH2F*			fNtrkl_PtOfMaxTrk_H_m;		//! histogram of Ntrkl vs max pT track/BG H non UE sub.				//43
	TH2F*			fNtrkl_PtOfMaxTrk_H70_m;	//! histogram of Ntrkl vs max pT track/BG H 70% UE				//44
	TH2F*			fNtrkl_PtOfMaxTrk_H80_m;	//! histogram of Ntrkl vs max pT track/BG H 80% UE				//45
	TH2F*			fNtrkl_PtOfMaxTrk_H85_m;	//! histogram of Ntrkl vs max pT track/BG H 85% UE				//46
	TH2F*			fNtrkl_PtOfMaxTrk_H90_m;	//! histogram of Ntrkl vs max pT track/BG H 90% UE				//47
	TH2F*			fNtrkl_PtOfMaxTrk_H93_m;	//! histogram of Ntrkl vs max pT track/BG H 93% UE				//48
	TH2F*			fNtrkl_PtOfMaxTrk_H95_m;	//! histogram of Ntrkl vs max pT track/BG H 95% UE				//49
	TH2F*			fNtrkl_PtOfTrks_H_m;		//! histogram of Ntrkl vs sum of pT track/BG H					//50
	TH2F*			fNtrkl_PtOfTrks_H70_m;		//! histogram of Ntrkl vs sum of pT track/BG H 70% UE				//51
	TH2F*			fNtrkl_PtOfTrks_H80_m;		//! histogram of Ntrkl vs sum of pT track/BG H 80% UE				//52
	TH2F*			fNtrkl_PtOfTrks_H85_m;		//! histogram of Ntrkl vs sum of pT track/BG H 85% UE				//53
	TH2F*			fNtrkl_PtOfTrks_H90_m;		//! histogram of Ntrkl vs sum of pT track/BG H 90% UE				//54
	TH2F*			fNtrkl_PtOfTrks_H93_m;		//! histogram of Ntrkl vs sum of pT track/BG H 93% UE				//55
	TH2F*			fNtrkl_PtOfTrks_H95_m;		//! histogram of Ntrkl vs sum of pT track/BG H 95% UE				//56
	TH2F*			fNtrkl_PtOfMaxTrk_W_w;		//! histogram of Ntrkl vs max pT track/e<-W non UE sub.				//57
	TH2F*			fNtrkl_PtOfMaxTrk_W70_w;	//! histogram of Ntrkl vs max pT track/e<-W 70% UE				//58
	TH2F*			fNtrkl_PtOfMaxTrk_W80_w;	//! histogram of Ntrkl vs max pT track/e<-W 80% UE				//59
	TH2F*			fNtrkl_PtOfMaxTrk_W85_w;	//! histogram of Ntrkl vs max pT track/e<-W 85% UE				//60
	TH2F*			fNtrkl_PtOfMaxTrk_W90_w;	//! histogram of Ntrkl vs max pT track/e<-W 90% UE				//61
	TH2F*			fNtrkl_PtOfMaxTrk_W93_w;	//! histogram of Ntrkl vs max pT track/e<-W 93% UE				//62
	TH2F*			fNtrkl_PtOfMaxTrk_W95_w;	//! histogram of Ntrkl vs max pT track/e<-W 95% UE				//63
	TH2F*			fNtrkl_PtOfTrks_W_w;		//! histogram of Ntrkl vs sum of pT track/e<-W					//64
	TH2F*			fNtrkl_PtOfTrks_W70_w;		//! histogram of Ntrkl vs sum of pT track/e<-W 70% UE				//65
	TH2F*			fNtrkl_PtOfTrks_W80_w;		//! histogram of Ntrkl vs sum of pT track/e<-W 80% UE				//66
	TH2F*			fNtrkl_PtOfTrks_W85_w;		//! histogram of Ntrkl vs sum of pT track/e<-W 85% UE				//67
	TH2F*			fNtrkl_PtOfTrks_W90_w;		//! histogram of Ntrkl vs sum of pT track/e<-W 90% UE				//68
	TH2F*			fNtrkl_PtOfTrks_W93_w;		//! histogram of Ntrkl vs sum of pT track/e<-W 93% UE				//69
	TH2F*			fNtrkl_PtOfTrks_W95_w;		//! histogram of Ntrkl vs sum of pT track/e<-W 95% UE				//70
	TH2F*			fNtrkl_PtOfMaxTrk_HF_w;		//! histogram of Ntrkl vs max pT track/e<-HF non UE sub.			//71
	TH2F*			fNtrkl_PtOfMaxTrk_HF70_w;	//! histogram of Ntrkl vs max pT track/e<-HF 70% UE				//72
	TH2F*			fNtrkl_PtOfMaxTrk_HF80_w;	//! histogram of Ntrkl vs max pT track/e<-HF 80% UE				//73
	TH2F*			fNtrkl_PtOfMaxTrk_HF85_w;	//! histogram of Ntrkl vs max pT track/e<-HF 85% UE				//74
	TH2F*			fNtrkl_PtOfMaxTrk_HF90_w;	//! histogram of Ntrkl vs max pT track/e<-HF 90% UE				//75
	TH2F*			fNtrkl_PtOfMaxTrk_HF93_w;	//! histogram of Ntrkl vs max pT track/e<-HF 93% UE				//76
	TH2F*			fNtrkl_PtOfMaxTrk_HF95_w;	//! histogram of Ntrkl vs max pT track/e<-HF 95% UE				//77
	TH2F*			fNtrkl_PtOfTrks_HF_w;		//! histogram of Ntrkl vs sum of pT track/e<-HF					//78
	TH2F*			fNtrkl_PtOfTrks_HF70_w;		//! histogram of Ntrkl vs sum of pT track/e<-HF 70% UE				//79
	TH2F*			fNtrkl_PtOfTrks_HF80_w;		//! histogram of Ntrkl vs sum of pT track/e<-HF 80% UE				//80
	TH2F*			fNtrkl_PtOfTrks_HF85_w;		//! histogram of Ntrkl vs sum of pT track/e<-HF 85% UE				//81
	TH2F*			fNtrkl_PtOfTrks_HF90_w;		//! histogram of Ntrkl vs sum of pT track/e<-HF 90% UE				//82
	TH2F*			fNtrkl_PtOfTrks_HF93_w;		//! histogram of Ntrkl vs sum of pT track/e<-HF 93% UE				//83
	TH2F*			fNtrkl_PtOfTrks_HF95_w;		//! histogram of Ntrkl vs sum of pT track/e<-HF 95% UE				//84
	TH2F*			fNtrkl_PtOfMaxTrk_H_w;		//! histogram of Ntrkl vs max pT track/BG H no UE sub.				//85
	TH2F*			fNtrkl_PtOfMaxTrk_H70_w;	//! histogram of Ntrkl vs max pT track/BG H 70% UE				//86
	TH2F*			fNtrkl_PtOfMaxTrk_H80_w;	//! histogram of Ntrkl vs max pT track/BG H 80% UE				//87
	TH2F*			fNtrkl_PtOfMaxTrk_H85_w;	//! histogram of Ntrkl vs max pT track/BG H 85% UE				//88
	TH2F*			fNtrkl_PtOfMaxTrk_H90_w;	//! histogram of Ntrkl vs max pT track/BG H 90% UE				//89
	TH2F*			fNtrkl_PtOfMaxTrk_H93_w;	//! histogram of Ntrkl vs max pT track/BG H 93% UE				//90
	TH2F*			fNtrkl_PtOfMaxTrk_H95_w;	//! histogram of Ntrkl vs max pT track/BG H 95% UE				//91
	TH2F*			fNtrkl_PtOfTrks_H_w;		//! histogram of Ntrkl vs sum of pT track/BG H					//92
	TH2F*			fNtrkl_PtOfTrks_H70_w;		//! histogram of Ntrkl vs sum of pT track/BG H 70% UE				//93
	TH2F*			fNtrkl_PtOfTrks_H80_w;		//! histogram of Ntrkl vs sum of pT track/BG H 80% UE				//94
	TH2F*			fNtrkl_PtOfTrks_H85_w;		//! histogram of Ntrkl vs sum of pT track/BG H 85% UE				//95
	TH2F*			fNtrkl_PtOfTrks_H90_w;		//! histogram of Ntrkl vs sum of pT track/BG H 90% UE				//96
	TH2F*			fNtrkl_PtOfTrks_H93_w;		//! histogram of Ntrkl vs sum of pT track/BG H 93% UE				//97
	TH2F*			fNtrkl_PtOfTrks_H95_w;		//! histogram of Ntrkl vs sum of pT track/BG H 95% UE				//98

	TH2F*			fHistPt_We_Ntrkl;		//! histogram of Pt of e<-W vs tracklets					//99

	TH2F*			fNtrkl_ClustE;			//! histogram for rejection factor						//100
	TH2F*			TPCSigForE;			//! hist of TPC signal (dE/dx) and p						//101
	TH2F*			fNsigmaPtForE;			//! histogram of nsigma vs Pt with e EMCal cut					//102
	TH2F*			fHistUEmult;			//! histogram of Mult vs UE pT							//103

	TH2F*			fNtrkl_LeadHadPt_m;		//! histogram of Ntrkl vs max pT track full in medium				//104
	TH2F*			fNtrkl_LeadHadPt_w;		//! histogram of Ntrkl vs max pT track full in wide				//105

	Bool_t                  fEMCEG1;//EMcal Threshold EG1

        AliAnalysisTaskWHMult(const AliAnalysisTaskWHMult&); // not implemented
        AliAnalysisTaskWHMult& operator=(const AliAnalysisTaskWHMult&); // not implemented

        ClassDef(AliAnalysisTaskWHMult, 1);
};

#endif
