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
        TH1F*                   fHistPt;        		//! dummy histogram								//1
	TH1F*			fPt_MCgenWe;			//! MC generated e<-W Pt							//2
	TH1F*			fPt_TrackingMCWe;		//! reconstructed track Pt (MC e<-W)						//3
	TH1F*			fPt_TPCPIDMCWe;			//! dE/dx cut track Pt (MC e<-W)						//4
	TH1F*			fPt_TrackMatchingMCWe;		//! matched track Pt (MC e<-W)							//5
	TH1F*			fPt_EMCalPIDMCWe;		//! M02,E/P,Eiso cut track Pt (MC e<-W)						//6
	TH1F*                   pVertex_all;			//! additional histgram in Step3						//7
	TH1F*                   pVertex;			//! after event selection							//8
	TH2F*			EtavsPhi;			//! hist of pseudorapidity and azimuthal distribution with cutVer		//9
	TH2F*			TPCSig;				//! hist of TPC signal (dE/dx) and p/z						//10
	TH1F*			Cent;				//! hist of the centrality							//11
	TH2F*			fzvtx_Ntrkl;			//! hist of vertexZ vs tracklets(corresponding to multiplicity)			//12
	TH2F*			fzvtx_Ntrkl_cal;		//! hist of vertexZ vs tracklets after calibration				//13
	TH2F*			fNtrklNch;			//! hist of vertexZ calibrated tracklets vs primary charged particles		//14
	TH1F*			fPDG;				//! histogram of pdg code							//15
	TH1F*			fMPDG;				//! histogram of mother's pdg code						//16
	TH1F*			fFPDG;				//! histogram of father's pdg code						//17

	TH1F*			fHistClustEMatch;		//! histogram of Energy in EMCal Cluster					//18
	TH2F*			fHistClustMCM02H;		//! histogram of MC non-electron M02						//19
	TH2F*			fHistClustMCM02E;		//! histogram of MC electron M02						//20

	TH2F*			fHistNsigmaP;			//! histogram of nsigma vs P							//21
	TH2F*			fHistMCNsigmaP;			//! histogram of nsigma vs P of MC events					//22
	TH2F*			fPtEoverPE;			//! histogram of Pt vs E/P of electron						//23
	TH2F*			fPtEoverPMCE;			//! histogram of Pt vs E/P of MC electron					//24
	TH2F*			fPtEoverPEGeo;			//! histogram of Pt vs E/P of electron in M02					//25
	TH2F*			fHistClustM02E;			//! histogram of MC events on M02						//26
        TH2F*			fPtEoverPH;			//! histogram of Pt vs E/P of hadron						//27
	TH1F*			fREisolation[3];		//! histogram of Eiso in R<0.3,0.4,0.5						//28,29,30
        TH1F*                   fREiso_MCW;			//! histogram of Eiso in R<0.3 of W more 10GeV					//31
        TH1F*                   fREiso_MCHF;			//! histogram of Eiso in R<0.3 of HF more 10GeV					//32
	TH1F*			fREiso_MCWhpt;			//! histogram of Eiso in R<0.3 of W more 30GeV					//33
	TH1F*			fREiso_MCHFhpt;			//! histogram of Eiso in R<0.3 of HF more 30GeV					//34
	TH2F*			fdPhi_trkW_Pt[3];		//! histogram of DeltaPhi w/ e<-W cut vs Pt R<0.3,0.4,0.5			//35,36,37
	TH2F*			fdPhi_trkHF_Pt[3];		//! histogram of DeltaPhi w/ e<-HF cut vs Pt R<0.3,0.4,0.5			//38,39,40
	TH2F*			fdPhi_trkW_ePt[3];		//! histogram of DeltaPhi (trk-W) vs Pt of e<-W R<0.3,0.4,0.5			//41,42,43
	TH2F*			fdPhi_trkHF_ePt[3];		//! histogram of DeltaPhi (trk-HF) vs Pt of e<-HF R<0.3,0.4,0.5			//44,45,46
	TH1F*			fHistPt_We[3];			//! histogram of Pt (W candidate) R<0.3,0.4,0.5					//47,48,49
	TH1F*			fHistPt_HFe[3];			//! histogram of Pt (HF candidate) R<0.3,0.4,0.5				//50,51,52
	TH1F*			fPt_maxtrack_W[3];		//! histogram of max pT track/e<-W R<0.3,0.4,0.5				//53,54,55
	TH2F*			fNtrkl_PtOfMaxTrk_W80[3];	//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5 80% UE		//56,57,58
	TH2F*			fNtrkl_PtOfMaxTrk_W85[3];	//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5 85% UE		//59,60,61
	TH2F*			fNtrkl_PtOfMaxTrk_W90[3];	//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5 90% UE		//62,63,64
	TH2F*			fNtrkl_PtOfMaxTrk_W95[3];	//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5 95% UE		//65,66,67
	TH2F*			fNtrkl_PtOfTrks_W[3];		//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5			//68,69,70
	TH2F*			fNtrkl_PtOfTrks_W80[3];		//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5 80% UE		//71,72,73
	TH2F*			fNtrkl_PtOfTrks_W85[3];		//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5 85% UE		//74,75,76
	TH2F*			fNtrkl_PtOfTrks_W90[3];		//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5 90% UE		//77,78,79
	TH2F*			fNtrkl_PtOfTrks_W95[3];		//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5 95% UE		//80,81,82
	TH2F*			fNtrkl_PtOfMaxTrk_HF80[3];	//! histogram of Ntrkl vs max pT track/e<-HF R<0.3,0.4,0.5 80% UE		//83,84,85
	TH2F*			fNtrkl_PtOfMaxTrk_HF85[3];	//! histogram of Ntrkl vs max pT track/e<-HF R<0.3,0.4,0.5 85% UE		//86,87,88
	TH2F*			fNtrkl_PtOfMaxTrk_HF90[3];	//! histogram of Ntrkl vs max pT track/e<-HF R<0.3,0.4,0.5 90% UE		//89,90,91
	TH2F*			fNtrkl_PtOfMaxTrk_HF95[3];	//! histogram of Ntrkl vs max pT track/e<-HF R<0.3,0.4,0.5 95% UE		//92,93,94
	TH2F*			fNtrkl_PtOfTrks_HF[3];		//! histogram of Ntrkl vs sum of pT track/e<-HF R<0.3,0.4,0.5			//95,96,97
	TH2F*			fNtrkl_PtOfTrks_HF80[3];	//! histogram of Ntrkl vs sum of pT track/e<-HF R<0.3,0.4,0.5 80% UE		//98,99,100
	TH2F*			fNtrkl_PtOfTrks_HF85[3];	//! histogram of Ntrkl vs sum of pT track/e<-HF R<0.3,0.4,0.5 85% UE		//101,102,103
	TH2F*			fNtrkl_PtOfTrks_HF90[3];	//! histogram of Ntrkl vs sum of pT track/e<-HF R<0.3,0.4,0.5 90% UE		//104,105,106
	TH2F*			fNtrkl_PtOfTrks_HF95[3];	//! histogram of Ntrkl vs sum of pT track/e<-HF R<0.3,0.4,0.5 95% UE		//107,108,109
	TH2F*			fHistPt_We_Ntrkl[3];		//! histogram of Pt of e<-W vs tracklets R<0.3,0.4,0.5				//110,111,112

	TH2F*			fdPhi_trkW_full[3];		//! tryal hist R<0.3,0.4,0.5							//113,114,115
	TH2F*			fdPhi_trkHF_full[3];		//! Eiso>0.1 tryal hist R<0.3,0.4,0.5						//116,117,118

	TH2F*			fNtrkl_ClustE;			//! histogram for rejection factor						//119
	TH2F*			TPCSigForE;			//! hist of TPC signal (dE/dx) and p						//120
	TH2F*			fNsigmaPtForE;			//! histogram of nsigma vs Pt with e EMCal cut					//121
	TH2F*			fHistNtrk_W[3];			//! histogram of Ntrack in R<0.3,0.4,0.5 of W					//122,123,124
	TH2F*			fHistNtrk_HF[3];		//! histogram of Ntrack in R<0.3,0.4,0.5 of HF					//125,126,127
	TH2F*			fHistEiso_Ntrk[3];		//! histogram of Eiso vs Ntrk in R<0.3,0.4,0.5					//128,129,130
	TH2F*			fHistUEmult[3];			//! histogram of Mult vs UE pT in R<0.3,0.4,0.5					//131,132,133

	Bool_t                  fEMCEG1;//EMcal Threshold EG1

        AliAnalysisTaskWHMult(const AliAnalysisTaskWHMult&); // not implemented
        AliAnalysisTaskWHMult& operator=(const AliAnalysisTaskWHMult&); // not implemented

        ClassDef(AliAnalysisTaskWHMult, 1);
};

#endif
