/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. */
/* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef AliAnalysisTaskNanoMUON_H
#define AliAnalysisTaskNanoMUON_H

#include "AliAnalysisTaskSE.h"

class AliMuonTrackCuts; 	// Include class for standard muon tack cuts


class AliAnalysisTaskNanoMUON : public AliAnalysisTaskSE  
{
public:
                            AliAnalysisTaskNanoMUON();
                            AliAnalysisTaskNanoMUON(const char *name);
    virtual                 ~AliAnalysisTaskNanoMUON();

    virtual void            UserCreateOutputObjects();
    virtual void            UserExec(Option_t* option);
    virtual void            Terminate(Option_t* option);
    virtual void   			NotifyRun();								// Implement the Notify run to search for the new parameters at each new runs
	void 					TwoMuonAna(Int_t *pos, Int_t *neg);			// Analyses two muons and extracs dimuon information
	void 					TwoMCMuonAna(Int_t *MCpos, Int_t *MCneg);	// Analyses two MC muons and extracs MC dimuon information
	void 					SetPeriod(TString period){fPeriod = period;}  // 0: 2018 q, 1: 2018 r, 2: 2015 o, 90: LHC18l7, 91: LHC16b2
	void 					SetMC(Bool_t flag){fIsMC = flag;}	
	void 					PostAllData();	

    AliMuonTrackCuts* 		fMuonTrackCuts; 					// Use the class as a data member

private:
	TString 				fPeriod;
	Bool_t 					fIsMC;

    AliAODEvent*            fAOD;       		//! input event
    AliMCEvent*				fMC;				//! input MC event

    TList*                  fOutputList; 		//! output list
    TH1F*                   fCounterH; 			//! counter for events passing each cut	
    TH2F*                   fNumberMuonsH; 		//! count good muons per event
    TH2F*                   fNumberMCMuonsH;	//! count MC muons per event
	// TH2F*                   fRAbsMuonH; 		//! distribution of RAbsMuon for selected events
	// TH2F*                   fMuMuMassPtH; 		//! kinematics of dimouns	

	// TH1F*					fZNAEnergyTimingH;	//! ZNA Energy histogram with at least one hit in the timing window
	// TH1F*					fZNCEnergyTimingH;	//! ZNC Energy histogram with at least one hit in the timing window
	// TH1F*					fZNATDCTimingH;	//! ZNA timing of all hits in events with at least one hit in the timing window
	// TH1F*					fZNCTDCTimingH;	//! ZNC timing of all hits in events with at least one hit in the timing window

	// TH1F*					fZNAEnergyTimingAllH;	//! ZNA Energy histogram with all hits in the timing window
	// TH1F*					fZNCEnergyTimingAllH;	//! ZNC Energy histogram with all hits in the timing window
	// TH1F*					fZNATDCTimingAllH;	//! ZNA timing of all hits in events with all hits in the timing window
	// TH1F*					fZNCTDCTimingAllH;	//! ZNC timing of all hits in events with all hits in the timing window

	// TH1F*					fZNAEnergy0NH;	//! ZNA Energy histogram with no hits
	// TH1F*					fZNCEnergy0NH;	//! ZNC Energy histogram with no hits

	TTree *fRecTree; 		//! analysis tree
	Int_t fRunNum;
	// Int_t fTracklets;
	Float_t fZNCEnergy; 
	Float_t fZNAEnergy;
	// Double_t fZPCEnergy; 
	// Double_t fZPAEnergy;
	Float_t fZNATDC[4];
	Float_t fZNCTDC[4];
	// Double_t fZPATDC[4];
	// Double_t fZPCTDC[4];
	Int_t fV0ADecision; 
	Int_t fV0CDecision;
	Int_t fV0AFiredCells; 
	Int_t fV0CFiredCells; 
	Int_t fADADecision; 
	Int_t fADCDecision;
  	// TBits fIR1Map;
  	// TBits fIR2Map;
	Float_t fMuMuPt; 
	// Double_t fMuMuPhi;
	Float_t fMuMuY; 
	Float_t fMuMuM;
	// Double_t fMuPt1; 
	// Double_t fMuPt2;
	// Double_t fMuEta1; 
	// Double_t fMuEta2;
	// Double_t fMuPhi1; 
	// Double_t fMuPhi2;
	// Double_t fMuQ1; 
	// Double_t fMuQ2;

	TClonesArray *fGenPart; 	//! MC particle object
	TTree *fGenTree; 		//! MC tree
	Int_t fMCRunNum;
	Float_t fMCMuMuPt; 
	// Double_t fMCMuMuPhi;
	Float_t fMCMuMuY; 
	Float_t fMCMuMuM;
	// Double_t fMCMuPt1; 
	// Double_t fMCMuPt2;
	// Double_t fMCMuEta1; 
	// Double_t fMCMuEta2;
	// Double_t fMCMuPhi1; 
	// Double_t fMCMuPhi2;
	// Double_t fMCMuPDG1; 
	// Double_t fMCMuPDG2;



    AliAnalysisTaskNanoMUON(const AliAnalysisTaskNanoMUON&); // not implemented
    AliAnalysisTaskNanoMUON& operator=(const AliAnalysisTaskNanoMUON&); // not implemented

    ClassDef(AliAnalysisTaskNanoMUON, 1);
};

#endif
