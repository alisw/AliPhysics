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

	TH1F*			fHistClustEMatch;		//! histogram of Energy in EMCal Cluster					//19
	TH2F*			fHistClustMCM02H;		//! histogram of MC non-electron M02						//20
	TH2F*			fHistClustMCM02E;		//! histogram of MC electron M02						//21

	TH2F*			fHistNsigmaP;			//! histogram of nsigma vs P							//22
	TH2F*			fHistMCNsigmaP;			//! histogram of nsigma vs P of MC events					//23
	TH2F*			fPtEoverPE;			//! histogram of Pt vs E/P of electron						//24
	TH2F*			fPtEoverPMCE;			//! histogram of Pt vs E/P of MC electron					//25
	TH2F*			fPtEoverPEGeo;			//! histogram of Pt vs E/P of electron in M02					//26
	TH2F*			fHistClustM02E;			//! histogram of MC events on M02						//27
        TH2F*			fPtEoverPH;			//! histogram of Pt vs E/P of hadron						//28
	TH1F*			fREisolation[3];		//! histogram of Eiso in R<0.3,0.4,0.5						//29,30,31
        TH1F*                   fREiso_MCW;			//! histogram of Eiso in R<0.3 of W more 10GeV					//32
        TH1F*                   fREiso_MCHF;			//! histogram of Eiso in R<0.3 of HF more 10GeV					//33
	TH1F*			fREiso_MCWhpt;			//! histogram of Eiso in R<0.3 of W more 30GeV					//34
	TH1F*			fREiso_MCHFhpt;			//! histogram of Eiso in R<0.3 of HF more 30GeV					//35
	TH2F*			fdPhi_trkW_Pt[3];		//! histogram of DeltaPhi w/ e<-W cut vs Pt R<0.3,0.4,0.5			//36,37,38
	TH2F*			fdPhi_trkHF_Pt[3];		//! histogram of DeltaPhi w/ e<-HF cut vs Pt R<0.3,0.4,0.5			//39,40,41
	TH2F*			fdPhi_trkW_ePt[3];		//! histogram of DeltaPhi (trk-W) vs Pt of e<-W R<0.3,0.4,0.5			//42,43,44
	TH2F*			fdPhi_trkHF_ePt[3];		//! histogram of DeltaPhi (trk-HF) vs Pt of e<-HF R<0.3,0.4,0.5			//45,46,47
	TH1F*			fHistPt_We[3];			//! histogram of Pt (W candidate) R<0.3,0.4,0.5					//48,49,50
	TH1F*			fHistPt_HFe[3];			//! histogram of Pt (HF candidate) R<0.3,0.4,0.5				//51,52,53
	TH1F*			fPt_maxtrack_W[3];		//! histogram of max pT track/e<-W R<0.3,0.4,0.5				//54,55,56
	TH2F*			fNtrkl_PtOfMaxTrk_W[3];		//! histogram of Ntrkl vs max pT track/e<-W R<0.3,0.4,0.5			//57,58,59
	TH2F*			fNtrkl_PtOfTrks_W[3];		//! histogram of Ntrkl vs sum of pT track/e<-W R<0.3,0.4,0.5			//60,61,62
	TH2F*			fHistPt_We_Ntrkl[3];		//! histogram of Pt of e<-W vs tracklets R<0.3,0.4,0.5				//63,64,65

	TH2F*			fdPhi_trkW_full[3];		//! tryal hist R<0.3,0.4,0.5							//66,67,68
	TH2F*			fdPhi_trkHF_full[3];		//! Eiso>0.1 tryal hist R<0.3,0.4,0.5						//69,70,71

	TH2F*			fNtrkl_ClustE;			//! histogram for rejection factor						//72
	TH2F*			TPCSigForE;			//! hist of TPC signal (dE/dx) and p						//73
	TH2F*			fNsigmaPtForE;			//! histogram of nsigma vs Pt with e EMCal cut					//74

	Bool_t                  fEMCEG1;//EMcal Threshold EG1

        AliAnalysisTaskWHMult(const AliAnalysisTaskWHMult&); // not implemented
        AliAnalysisTaskWHMult& operator=(const AliAnalysisTaskWHMult&); // not implemented

        ClassDef(AliAnalysisTaskWHMult, 1);
};

#endif
