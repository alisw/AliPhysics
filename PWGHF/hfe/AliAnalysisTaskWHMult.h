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
        TList*                  fOutputList;    		//! output list								//Object List Number//
        TTree*			tree;				//! ref. for mult correction
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
	TH1F*			fREisolation[3];		//! histogram of Eiso in R<0.3,0.4,0.5						//10,11,12
	TH2F*			fdPhi_trkW_Pt[3];		//! histogram of DeltaPhi w/ e<-W cut vs Pt R<0.3,0.4,0.5			//13,14,15
	TH2F*			fdPhi_trkHF_Pt[3];		//! histogram of DeltaPhi w/ e<-HF cut vs Pt R<0.3,0.4,0.5			//16,17,18
	TH1F*			fHistPt_We[3];			//! histogram of Pt (W candidate) R<0.3,0.4,0.5					//19,20,21
	TH1F*			fHistPt_HFe[3];			//! histogram of Pt (HF candidate) R<0.3,0.4,0.5				//22,23,24
	TH2F*			fNtrkl_PtOfMaxTrk_W_m[3];	//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5 non UE sub.		//25,26,27
	TH2F*			fNtrkl_PtOfMaxTrk_W70_m[3];	//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5 70% UE		//28,29,30
	TH2F*			fNtrkl_PtOfMaxTrk_W80_m[3];	//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5 80% UE		//31,32,33
	TH2F*			fNtrkl_PtOfMaxTrk_W85_m[3];	//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5 85% UE		//34,35,36
	TH2F*			fNtrkl_PtOfMaxTrk_W90_m[3];	//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5 90% UE		//37,38,39
	TH2F*			fNtrkl_PtOfMaxTrk_W93_m[3];	//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5 93% UE		//40,41,42
	TH2F*			fNtrkl_PtOfMaxTrk_W95_m[3];	//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5 95% UE		//43,44,45
	TH2F*			fNtrkl_PtOfTrks_W_m[3];		//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5			//46,47,48
	TH2F*			fNtrkl_PtOfTrks_W70_m[3];	//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5 70% UE		//49,50,51
	TH2F*			fNtrkl_PtOfTrks_W80_m[3];	//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5 80% UE		//52,53,54
	TH2F*			fNtrkl_PtOfTrks_W85_m[3];	//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5 85% UE		//55,56,57
	TH2F*			fNtrkl_PtOfTrks_W90_m[3];	//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5 90% UE		//58,59,60
	TH2F*			fNtrkl_PtOfTrks_W93_m[3];	//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5 93% UE		//61,62,63
	TH2F*			fNtrkl_PtOfTrks_W95_m[3];	//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5 95% UE		//64,65,66
	TH2F*			fNtrkl_PtOfMaxTrk_HF_m[3];	//! histogram of Ntrkl vs max pT track/e<-HF R<0.3,0.4,0.5 non UE sub.		//67,68,69
	TH2F*			fNtrkl_PtOfMaxTrk_HF70_m[3];	//! histogram of Ntrkl vs max pT track/e<-HF R<0.3,0.4,0.5 70% UE		//70,71,72
	TH2F*			fNtrkl_PtOfMaxTrk_HF80_m[3];	//! histogram of Ntrkl vs max pT track/e<-HF R<0.3,0.4,0.5 80% UE		//73,74,75
	TH2F*			fNtrkl_PtOfMaxTrk_HF85_m[3];	//! histogram of Ntrkl vs max pT track/e<-HF R<0.3,0.4,0.5 85% UE		//76,77,78
	TH2F*			fNtrkl_PtOfMaxTrk_HF90_m[3];	//! histogram of Ntrkl vs max pT track/e<-HF R<0.3,0.4,0.5 90% UE		//79,80,81
	TH2F*			fNtrkl_PtOfMaxTrk_HF93_m[3];	//! histogram of Ntrkl vs max pT track/e<-HF R<0.3,0.4,0.5 93% UE		//82,83,84
	TH2F*			fNtrkl_PtOfMaxTrk_HF95_m[3];	//! histogram of Ntrkl vs max pT track/e<-HF R<0.3,0.4,0.5 95% UE		//85,86,87
	TH2F*			fNtrkl_PtOfTrks_HF_m[3];	//! histogram of Ntrkl vs sum of pT track/e<-HF R<0.3,0.4,0.5			//88,89,90
	TH2F*			fNtrkl_PtOfTrks_HF70_m[3];	//! histogram of Ntrkl vs sum of pT track/e<-HF R<0.3,0.4,0.5 70% UE		//91,92,93
	TH2F*			fNtrkl_PtOfTrks_HF80_m[3];	//! histogram of Ntrkl vs sum of pT track/e<-HF R<0.3,0.4,0.5 80% UE		//94,95,96
	TH2F*			fNtrkl_PtOfTrks_HF85_m[3];	//! histogram of Ntrkl vs sum of pT track/e<-HF R<0.3,0.4,0.5 85% UE		//97,98,99
	TH2F*			fNtrkl_PtOfTrks_HF90_m[3];	//! histogram of Ntrkl vs sum of pT track/e<-HF R<0.3,0.4,0.5 90% UE		//100,101,102
	TH2F*			fNtrkl_PtOfTrks_HF93_m[3];	//! histogram of Ntrkl vs sum of pT track/e<-HF R<0.3,0.4,0.5 93% UE		//103,104,105
	TH2F*			fNtrkl_PtOfTrks_HF95_m[3];	//! histogram of Ntrkl vs sum of pT track/e<-HF R<0.3,0.4,0.5 95% UE		//106,107,108
	TH2F*			fNtrkl_PtOfMaxTrk_H_m[3];	//! histogram of Ntrkl vs max pT track/BG H R<0.3,0.4,0.5 non UE sub.		//109,110,111
	TH2F*			fNtrkl_PtOfMaxTrk_H70_m[3];	//! histogram of Ntrkl vs max pT track/BG H R<0.3,0.4,0.5 70% UE		//112,113,114
	TH2F*			fNtrkl_PtOfMaxTrk_H80_m[3];	//! histogram of Ntrkl vs max pT track/BG H R<0.3,0.4,0.5 80% UE		//115,116,117
	TH2F*			fNtrkl_PtOfMaxTrk_H85_m[3];	//! histogram of Ntrkl vs max pT track/BG H R<0.3,0.4,0.5 85% UE		//118,119,120
	TH2F*			fNtrkl_PtOfMaxTrk_H90_m[3];	//! histogram of Ntrkl vs max pT track/BG H R<0.3,0.4,0.5 90% UE		//121,122,123
	TH2F*			fNtrkl_PtOfMaxTrk_H93_m[3];	//! histogram of Ntrkl vs max pT track/BG H R<0.3,0.4,0.5 93% UE		//124,125,126
	TH2F*			fNtrkl_PtOfMaxTrk_H95_m[3];	//! histogram of Ntrkl vs max pT track/BG H R<0.3,0.4,0.5 95% UE		//127,128,129
	TH2F*			fNtrkl_PtOfTrks_H_m[3];		//! histogram of Ntrkl vs sum of pT track/BG H R<0.3,0.4,0.5			//130,131,132
	TH2F*			fNtrkl_PtOfTrks_H70_m[3];	//! histogram of Ntrkl vs sum of pT track/BG H R<0.3,0.4,0.5 70% UE		//133,134,135
	TH2F*			fNtrkl_PtOfTrks_H80_m[3];	//! histogram of Ntrkl vs sum of pT track/BG H R<0.3,0.4,0.5 80% UE		//136,137,138
	TH2F*			fNtrkl_PtOfTrks_H85_m[3];	//! histogram of Ntrkl vs sum of pT track/BG H R<0.3,0.4,0.5 85% UE		//139,140,141
	TH2F*			fNtrkl_PtOfTrks_H90_m[3];	//! histogram of Ntrkl vs sum of pT track/BG H R<0.3,0.4,0.5 90% UE		//142,143,144
	TH2F*			fNtrkl_PtOfTrks_H93_m[3];	//! histogram of Ntrkl vs sum of pT track/BG H R<0.3,0.4,0.5 93% UE		//145,146,147
	TH2F*			fNtrkl_PtOfTrks_H95_m[3];	//! histogram of Ntrkl vs sum of pT track/BG H R<0.3,0.4,0.5 95% UE		//148,149,150
	TH2F*			fNtrkl_PtOfMaxTrk_W_w[3];	//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5 non UE sub.		//151,152,153
	TH2F*			fNtrkl_PtOfMaxTrk_W70_w[3];	//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5 70% UE		//154,155,156
	TH2F*			fNtrkl_PtOfMaxTrk_W80_w[3];	//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5 80% UE		//157,158,159
	TH2F*			fNtrkl_PtOfMaxTrk_W85_w[3];	//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5 85% UE		//160,161,162
	TH2F*			fNtrkl_PtOfMaxTrk_W90_w[3];	//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5 90% UE		//163,164,165
	TH2F*			fNtrkl_PtOfMaxTrk_W93_w[3];	//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5 93% UE		//166,167,168
	TH2F*			fNtrkl_PtOfMaxTrk_W95_w[3];	//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5 95% UE		//169,170,171
	TH2F*			fNtrkl_PtOfTrks_W_w[3];		//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5			//172,173,174
	TH2F*			fNtrkl_PtOfTrks_W70_w[3];	//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5 70% UE		//175,176,177
	TH2F*			fNtrkl_PtOfTrks_W80_w[3];	//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5 80% UE		//178,179,180
	TH2F*			fNtrkl_PtOfTrks_W85_w[3];	//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5 85% UE		//181,182,183
	TH2F*			fNtrkl_PtOfTrks_W90_w[3];	//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5 90% UE		//184,185,186
	TH2F*			fNtrkl_PtOfTrks_W93_w[3];	//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5 93% UE		//187,188,189
	TH2F*			fNtrkl_PtOfTrks_W95_w[3];	//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5 95% UE		//190,191,192
	TH2F*			fNtrkl_PtOfMaxTrk_HF_w[3];	//! histogram of Ntrkl vs max pT track/e<-HF R<0.3,0.4,0.5 non UE sub.		//193,194,195
	TH2F*			fNtrkl_PtOfMaxTrk_HF70_w[3];	//! histogram of Ntrkl vs max pT track/e<-HF R<0.3,0.4,0.5 70% UE		//196,197,198
	TH2F*			fNtrkl_PtOfMaxTrk_HF80_w[3];	//! histogram of Ntrkl vs max pT track/e<-HF R<0.3,0.4,0.5 80% UE		//199,200,201
	TH2F*			fNtrkl_PtOfMaxTrk_HF85_w[3];	//! histogram of Ntrkl vs max pT track/e<-HF R<0.3,0.4,0.5 85% UE		//202,203,204
	TH2F*			fNtrkl_PtOfMaxTrk_HF90_w[3];	//! histogram of Ntrkl vs max pT track/e<-HF R<0.3,0.4,0.5 90% UE		//205,206,207
	TH2F*			fNtrkl_PtOfMaxTrk_HF93_w[3];	//! histogram of Ntrkl vs max pT track/e<-HF R<0.3,0.4,0.5 93% UE		//208,209,210
	TH2F*			fNtrkl_PtOfMaxTrk_HF95_w[3];	//! histogram of Ntrkl vs max pT track/e<-HF R<0.3,0.4,0.5 95% UE		//211,212,213
	TH2F*			fNtrkl_PtOfTrks_HF_w[3];	//! histogram of Ntrkl vs sum of pT track/e<-HF R<0.3,0.4,0.5			//214,215,216
	TH2F*			fNtrkl_PtOfTrks_HF70_w[3];	//! histogram of Ntrkl vs sum of pT track/e<-HF R<0.3,0.4,0.5 70% UE		//217,218,219
	TH2F*			fNtrkl_PtOfTrks_HF80_w[3];	//! histogram of Ntrkl vs sum of pT track/e<-HF R<0.3,0.4,0.5 80% UE		//220,221,222
	TH2F*			fNtrkl_PtOfTrks_HF85_w[3];	//! histogram of Ntrkl vs sum of pT track/e<-HF R<0.3,0.4,0.5 85% UE		//223,224,225
	TH2F*			fNtrkl_PtOfTrks_HF90_w[3];	//! histogram of Ntrkl vs sum of pT track/e<-HF R<0.3,0.4,0.5 90% UE		//226,227,228
	TH2F*			fNtrkl_PtOfTrks_HF93_w[3];	//! histogram of Ntrkl vs sum of pT track/e<-HF R<0.3,0.4,0.5 93% UE		//229,230,231
	TH2F*			fNtrkl_PtOfTrks_HF95_w[3];	//! histogram of Ntrkl vs sum of pT track/e<-HF R<0.3,0.4,0.5 95% UE		//232,233,234
	TH2F*			fNtrkl_PtOfMaxTrk_H_w[3];	//! histogram of Ntrkl vs max pT track/BG H R<0.3,0.4,0.5 non UE sub.		//235,236,237
	TH2F*			fNtrkl_PtOfMaxTrk_H70_w[3];	//! histogram of Ntrkl vs max pT track/BG H R<0.3,0.4,0.5 70% UE		//238,239,240
	TH2F*			fNtrkl_PtOfMaxTrk_H80_w[3];	//! histogram of Ntrkl vs max pT track/BG H R<0.3,0.4,0.5 80% UE		//241,242,243
	TH2F*			fNtrkl_PtOfMaxTrk_H85_w[3];	//! histogram of Ntrkl vs max pT track/BG H R<0.3,0.4,0.5 85% UE		//244,245,246
	TH2F*			fNtrkl_PtOfMaxTrk_H90_w[3];	//! histogram of Ntrkl vs max pT track/BG H R<0.3,0.4,0.5 90% UE		//247,248,249
	TH2F*			fNtrkl_PtOfMaxTrk_H93_w[3];	//! histogram of Ntrkl vs max pT track/BG H R<0.3,0.4,0.5 93% UE		//250,251,252
	TH2F*			fNtrkl_PtOfMaxTrk_H95_w[3];	//! histogram of Ntrkl vs max pT track/BG H R<0.3,0.4,0.5 95% UE		//253,254,255
	TH2F*			fNtrkl_PtOfTrks_H_w[3];		//! histogram of Ntrkl vs sum of pT track/BG H R<0.3,0.4,0.5			//256,257,258
	TH2F*			fNtrkl_PtOfTrks_H70_w[3];	//! histogram of Ntrkl vs sum of pT track/BG H R<0.3,0.4,0.5 70% UE		//259,260,261
	TH2F*			fNtrkl_PtOfTrks_H80_w[3];	//! histogram of Ntrkl vs sum of pT track/BG H R<0.3,0.4,0.5 80% UE		//262,263,264
	TH2F*			fNtrkl_PtOfTrks_H85_w[3];	//! histogram of Ntrkl vs sum of pT track/BG H R<0.3,0.4,0.5 85% UE		//265,266,267
	TH2F*			fNtrkl_PtOfTrks_H90_w[3];	//! histogram of Ntrkl vs sum of pT track/BG H R<0.3,0.4,0.5 90% UE		//268,269,270
	TH2F*			fNtrkl_PtOfTrks_H93_w[3];	//! histogram of Ntrkl vs sum of pT track/BG H R<0.3,0.4,0.5 93% UE		//271,272,273
	TH2F*			fNtrkl_PtOfTrks_H95_w[3];	//! histogram of Ntrkl vs sum of pT track/BG H R<0.3,0.4,0.5 95% UE		//274,275,276

	TH2F*			fHistPt_We_Ntrkl[3];		//! histogram of Pt of e<-W vs tracklets R<0.3,0.4,0.5				//277,278,279

	TH2F*			fNtrkl_ClustE;			//! histogram for rejection factor						//280
	TH2F*			TPCSigForE;			//! hist of TPC signal (dE/dx) and p						//281
	TH2F*			fNsigmaPtForE;			//! histogram of nsigma vs Pt with e EMCal cut					//282
	TH2F*			fHistUEmult[3];			//! histogram of Mult vs UE pT in R<0.3,0.4,0.5					//283,284,285

	Bool_t                  fEMCEG1;//EMcal Threshold EG1

        AliAnalysisTaskWHMult(const AliAnalysisTaskWHMult&); // not implemented
        AliAnalysisTaskWHMult& operator=(const AliAnalysisTaskWHMult&); // not implemented

        ClassDef(AliAnalysisTaskWHMult, 1);
};

#endif
