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
	//virtual	void		ProcessMCParticles();		//! new function for MC simulation
        virtual void            FindMother(AliAODMCParticle* part, int &label, int &pid, double &ptmom);        //! find mother command
	virtual void		FindFather(AliAODMCParticle* part, int &label, int &pid, double &ptmom);	//! find father command
	virtual Double_t	GetCorrectedNtrackletsD(TTree* tree, Double_t uncorrectedNacc, Int_t BinOfvertexZ);
	virtual Int_t		CountNch();

	Bool_t                  GetEMCalTriggerEG1() { return fEMCEG1; };
	void                    SetEG1(Bool_t flagEG1) { fEMCEG1= flagEG1;};
        
        virtual void            Terminate(Option_t* option);

    private:
        AliAODEvent*            fAOD;           	//! input event
	//AliMCEvent*             fMCEvent;       	//! corresponding MC event
	AliVEvent*              fVevent;        //!event object
	AliPIDResponse*		fPIDResponse;		//! pid response object
	AliAODMCParticle*	fMCTrackpart;		//! MC track particle
	TClonesArray*		fMCarray;		//! MC array
	AliAODMCParticle*	fMCparticle;		//! MC particle
        TList*                  fOutputList;    	//! output list								//Object List Number//
        TTree*                  tree;                   // ref. for mult correction
        TH1F*                   fHistPt;        	//! dummy histogram								//0
	TH1F*                   pVertex;        	//! additional histgram in Step3						//1
	TH1F*                   cutVer;         	//! accepting within 10 cm							//2
	TH2F*			EtavsPhi;		//! hist of pseudorapidity and azimuthal distribution with cutVer		//3
	TH2F*			TPCSig;			//! hist of TPC signal (dE/dx) and p/z						//4
	TH1F*			Cent;			//! hist of the centrality							//5
	TH2F*			fzvtx_Ntrkl;		//! hist of vertexZ vs tracklets(corresponding to multiplicity)			//6
	TH2F*			fzvtx_Ntrkl_cal;	//! hist of vertexZ vs tracklets after calibration				//7
	TH2F*			fNtrklNch;		//! hist of vertexZ calibrated tracklets vs primary charged particles		//8
	TH1F*			fPDG;			//! histogram of pdg code							//9
	TH1F*			fMPDG;			//! histogram of mother's pdg code						//10
	TH1F*			fFPDG;			//! histogram of father's pdg code						//11

	TH1F*			fHistClustE;		//! histogram of Energy in EMCal Cluster					//12
	TH1F*			fHistClustEMatch;	//! histogram of Energy in EMCal Cluster					//13
	TH1F*			fHistClustLongAxis;	//! histogram of events on M02							//14
	TH1F*			fHistClustLongAxisE;	//! histogram of events on M02 about Electron					//15

	TH2F*			fHistNsigmaP;		//! histogram of nsigma vs P							//16
	TH2F*			fHistMCNsigmaP;		//! histogram of nsigma vs P of MC events					//17
	TH2F*			fPtEoverPE;		//! histogram of Pt vs E/P of electron						//18
	TH2F*			fPtEoverPMCE;		//! histogram of Pt vs E/P of MC electron					//19
	TH2F*			fPtEoverPEGeo;		//! histogram of Pt vs E/P of electron in M02					//20
	TH1F*			fHistMCClsLAE;		//! histogram of MC events on M02						//21
        TH2F*			fPtEoverPH;		//! histogram of Pt vs E/P of hadron						//22
        TH2F*                   fPtEoverPHGeo;          //! histogram of Pt vs E/P of hadron in M02					//23
	TH1F*			fREisolation;		//! histogram of Eiso in R<0.3							//24
        TH1F*                   fREisoW;		//! histogram of Eiso in R<0.3 of W more 10GeV					//25
        TH1F*                   fREisoHF;		//! histogram of Eiso in R<0.3 of HF more 10GeV					//26
	TH1F*			fREisoWhpt;		//! histogram of Eiso in R<0.3 of W more 20GeV					//27
	TH1F*			fREisoHFhpt;		//! histogram of Eiso in R<0.3 of HF more 20GeV					//28
	TH1F*			fREisoWvhpt;		//! histogram of Eiso in R<0.3 of W more 30GeV					//29
	TH1F*			fREisoHFvhpt;		//! histogram of Eiso in R<0.3 of HF more 30GeV					//30
	TH2F*			EleTraDelPhi;		//! histogram of DeltaPhi w/ full e<-W cut vs Pt				//31
	TH1F*			fHistPt_We;		//! histogram of Pt (e<-W)							//32
	TH1F*			fHistPt_HFe;		//! histogram of Pt (e<-HF)							//33
	TH1F*			fPt_maxtrack_W;		//! histogram of Rate of Pt(highest Pt track/e<-W)				//34
	TH2F*			fNtrkl_PtOfMaxTrk_W;	//! histogram of Ntrkl vs Rate of Pt(Highest Pt track/e<-W)			//35
	TH2F*			fHistPt_We_Ntrkl;	//! histogram of Pt of e<-W vs tracklets					//36
	TH1F*			fW_true;		//! histogram of Pt of true e<-W						//37
	TH1F*			fW_false;		//! histogram of Pt of false e<-W						//38

	TH1F*			EleTraDelPhi_fullrange;	//! tryal hist

	Bool_t                  fEMCEG1;//EMcal Threshold EG1

        AliAnalysisTaskWHMult(const AliAnalysisTaskWHMult&); // not implemented
        AliAnalysisTaskWHMult& operator=(const AliAnalysisTaskWHMult&); // not implemented

        ClassDef(AliAnalysisTaskWHMult, 1);
};

#endif
