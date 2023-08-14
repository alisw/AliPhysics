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
	AliVEvent*              fVevent;			//!event object
	AliPIDResponse*		fPIDResponse;			//! pid response object
	AliAODMCParticle*	fMCTrackpart;			//! MC track particle
	TClonesArray*		fMCarray;			//! MC array
	AliAODMCParticle*	fMCparticle;			//! MC particle
        TList*                  fOutputList;    		//! output list								//Object List Number//
        TTree*			tree;				//! ref. for mult correction
	TH1F*			fNevents;			//! hist of event cut performance						//0
        TH1F*                   fHistPt;        		//! dummy histogram								//1
	TH1F*                   pVertex_all;    	    	//! additional histgram in Step3						//2
	TH1F*                   pVertex;        	 	//! after event selection							//3
	TH2F*			EtavsPhi;			//! hist of pseudorapidity and azimuthal distribution with cutVer		//4
	TH2F*			TPCSig;				//! hist of TPC signal (dE/dx) and p/z						//5
	TH1F*			Cent;				//! hist of the centrality							//6
	TH2F*			fzvtx_Ntrkl;			//! hist of vertexZ vs tracklets(corresponding to multiplicity)			//7
	TH2F*			fzvtx_Ntrkl_cal;		//! hist of vertexZ vs tracklets after calibration				//8
	TH2F*			fNtrklNch;			//! hist of vertexZ calibrated tracklets vs primary charged particles		//9
	TH1F*			fPDG;				//! histogram of pdg code							//10
	TH1F*			fMPDG;				//! histogram of mother's pdg code						//11
	TH1F*			fFPDG;				//! histogram of father's pdg code						//12

	TH1F*			fHistClustE;			//! histogram of Energy in EMCal Cluster					//13
	TH1F*			fHistClustEMatch;		//! histogram of Energy in EMCal Cluster					//14
	TH1F*			fHistClustLongAxis;		//! histogram of events on M02							//15
	TH1F*			fHistClustLongAxisE;		//! histogram of events on M02 about Electron					//16

	TH2F*			fHistNsigmaP;			//! histogram of nsigma vs P							//17
	TH2F*			fHistMCNsigmaP;			//! histogram of nsigma vs P of MC events					//18
	TH2F*			fPtEoverPE;			//! histogram of Pt vs E/P of electron						//19
	TH2F*			fPtEoverPMCE;			//! histogram of Pt vs E/P of MC electron					//20
	TH2F*			fPtEoverPEGeo;			//! histogram of Pt vs E/P of electron in M02					//21
	TH1F*			fHistMCClsLAE;			//! histogram of MC events on M02						//22
        TH2F*			fPtEoverPH;			//! histogram of Pt vs E/P of hadron						//23
        TH2F*                   fPtEoverPHGeo;         	 	//! histogram of Pt vs E/P of hadron in M02					//24
	TH1F*			fREisolation[3];		//! histogram of Eiso in R<0.3,0.4,0.5						//25,26,27
        TH1F*                   fREiso_MCW;			//! histogram of Eiso in R<0.3 of W more 10GeV					//28
        TH1F*                   fREiso_MCHF;			//! histogram of Eiso in R<0.3 of HF more 10GeV					//29
	TH1F*			fREiso_MCWhpt;			//! histogram of Eiso in R<0.3 of W more 20GeV					//30
	TH1F*			fREiso_MCHFhpt;			//! histogram of Eiso in R<0.3 of HF more 20GeV					//31
	TH1F*			fREiso_MCWvhpt;			//! histogram of Eiso in R<0.3 of W more 30GeV					//32
	TH1F*			fREiso_MCHFvhpt;		//! histogram of Eiso in R<0.3 of HF more 30GeV					//33
	TH2F*			fdPhi_trkW_Pt[3];		//! histogram of DeltaPhi w/ e<-W cut vs Pt R<0.3,0.4,0.5			//34,35,36
	TH2F*			fdPhi_trkHF_Pt[3];		//! histogram of DeltaPhi w/ e<-HF cut vs Pt R<0.3,0.4,0.5			//37,38,39
	TH2F*			fdPhi_trkW_Pt_hpt[3];		//! histogram of DeltaPhi w/ e<-W cut vs Pt R<0.3,0.4,0.5 over 30GeV/c		//40,41,42
	TH2F*			fdPhi_trkHF_Pt_hpt[3];		//! histogram of DeltaPhi w/ e<-HF cut vs Pt R<0.3,0.4,0.5 over 30GeV/c		//43,44,45
	TH2F*			fdPhi_trkW_ePt[3];		//! histogram of DeltaPhi (trk-W) vs Pt of e<-W R<0.3,0.4,0.5			//46,47,48
	TH2F*			fdPhi_trkHF_ePt[3];		//! histogram of DeltaPhi (trk-HF) vs Pt of e<-HF R<0.3,0.4,0.5			//49,50,51
	TH1F*			fHistPt_We[3];			//! histogram of Pt (e<-W) R<0.3,0.4,0.5					//52,53,54
	TH1F*			fHistPt_HFe[3];			//! histogram of Pt (e<-HF) R<0.3,0.4,0.5					//55,56,57
	TH1F*			fPt_maxtrack_W[3];		//! histogram of max pT track/e<-W R<0.3,0.4,0.5				//58,59,60
	TH1F*			fPt_maxtrack_W_hpt[3];		//! histogram of max pT track/e<-W R<0.3,0.4,0.5 over 30GeV/c			//61,62,63
	TH2F*			fNtrkl_PtOfMaxTrk_W[3];		//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5			//64,65,66
	TH2F*			fNtrkl_PtOfMaxTrk_W_hpt[3];	//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5 over 30GeV/c		//67,68,69
	TH2F*			fHistPt_We_Ntrkl[3];		//! histogram of Pt of e<-W vs tracklets R<0.3,0.4,0.5				//70,71,72

	TH2F*			fdPhi_trkW_full[3];		//! tryal hist R<0.3,0.4,0.5							//73,74,75
	TH2F*			fdPhi_trkHF_full[3];		//! Eiso>0.1 tryal hist R<0.3,0.4,0.5						//76,77,78
	TH2F*			fdPhi_trkW_full_hpt[3];		//! tryal hist R<0.3,0.4,0.5 over 30GeV/c					//79,80,81
	TH2F*			fdPhi_trkHF_full_hpt[3];	//! Eiso>0.1 tryal hist R<0.3,0.4,0.5 over 30GeV/c				//82,83,84

	TH2F*			fNtrkl_ClustE;			//! histogram for rejection factor						//85

	Bool_t                  fEMCEG1;//EMcal Threshold EG1

        AliAnalysisTaskWHMult(const AliAnalysisTaskWHMult&); // not implemented
        AliAnalysisTaskWHMult& operator=(const AliAnalysisTaskWHMult&); // not implemented

        ClassDef(AliAnalysisTaskWHMult, 1);
};

#endif
