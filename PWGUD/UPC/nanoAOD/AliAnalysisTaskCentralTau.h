/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice */
/* $Id$ */

#ifndef ALIANALYSISTASKCentralTau_H
#define ALIANALYSISTASKCentralTau_H

class TH1;
class TH2;
class TTree;
class TList;
class TFile;
class TBits;

#include "AliPIDResponse.h"
#include "AliESDtrackCuts.h"
#include "AliTimeRangeCut.h"
#include "AliAnalysisTaskSE.h"
#include "TLorentzVector.h"

class AliAnalysisTaskCentralTau : public AliAnalysisTaskSE {
public:
		AliAnalysisTaskCentralTau();
		AliAnalysisTaskCentralTau(const char *name);
		virtual ~AliAnalysisTaskCentralTau();

		virtual void UserCreateOutputObjects();
		virtual void UserExec(Option_t *option);
		virtual void Terminate(Option_t *);

		Int_t TestPIDhypothesis(AliESDtrack *trk);
		void TPCandTOFsignalInfo(AliESDtrack *trk, Int_t trkID);
		static void SetCrossed(Int_t spd[4], TBits &crossed);
		static Int_t GetChipId(Int_t index, Int_t &chipId2, Bool_t debug=false);
		static Bool_t IsSTGFired(const TBits& bits, Int_t dphiMin=4, Int_t dphiMax=10, Bool_t tolerance = true);
private:

		AliPIDResponse *fPIDResponse;
		AliESDtrackCuts *fTrackCutsBit0;
		AliESDtrackCuts *fTrackCutsBit1;
		AliESDtrackCuts *fTrackCutsBit4;

		TList *fOutputList;		//!
		TList *fOutputPID;   //!
		TH2I *hTriggerCounter;	//!
		TH1I *hParticleTypeCounter; //!
		TTree *tTwoTracks;		//!
		TTree *tPID;    //!

		Double_t  fZNAenergy, fZNCenergy, fZNAtime[4], fZNCtime[4];
		TClonesArray *fESDtracks;
		Float_t fPIDTPC[5][2], fPIDTOF[5][2];
		Int_t fRunNumber, fADAdecision, fADCdecision, fV0Adecision, fV0Cdecision, fSign;
		Int_t fTrackPIDid[2];
		Bool_t fTriggers[10], fTriggerClass[3];

		// PID analysis
		Double_t fPIDpt[2],fPIDmomentum[2],fTPCsignal[2],fTOFsignal[2];
		Int_t fTPCmostProbableTrackType[2], fTOFmostProbableTrackType[2];

		AliAnalysisTaskCentralTau(const AliAnalysisTaskCentralTau&); //not implemented
		AliAnalysisTaskCentralTau& operator =(const AliAnalysisTaskCentralTau&); //not implemented

ClassDef(AliAnalysisTaskCentralTau, 31);
};

#endif
