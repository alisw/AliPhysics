#ifndef ALIANALYSISTASKHFFINDJETS_H
#define ALIANALYSISTASKHFFINDJETS_H

#include "AliAnalysisTaskSE.h"

#include <TMath.h>
#include <TDatabasePDG.h>
#include <TFile.h>
#include <TSystem.h>
#include <TH1F.h>
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliAODVertex.h"
#include "AliESDVertex.h"
#include "AliESDEvent.h"
#include "AliAODRecoDecayHF2Prong.h"
#include "AliAODRecoDecayHF3Prong.h"
#include "AliVertexerTracks.h"
#include "AliFJWrapper.h"
#include "FJ_includes.h"

/* Copyright(c) 1998-2022, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//*************************************************************************
// Class AliAnalysisTaskHFSimpleVertices
// AliAnalysisTaskSE to extract D meson candidates from ESDs
//          
//*************************************************************************


class AliAnalysisTaskHFFindJets : public AliAnalysisTaskSE
{

 public:
	AliAnalysisTaskHFFindJets(); // class constructor
	AliAnalysisTaskHFFindJets(const char *name); // class constructor
	virtual ~AliAnalysisTaskHFFindJets(); // class destructor
	
	virtual void UserCreateOutputObjects(); // called once at beginning of runtime
	virtual void UserExec(Option_t* option); // called for each event
	virtual void Terminate(Option_t* option); // called at end of analysis
	
	void SetReadMC(Bool_t read){fReadMC=read;}
	void InitFromJson(TString esdFile);
	TTree* tree;
	

 private:                                
    double ptmintrack = 0.;
    int do3Prongs = 0;
    bool doJetFinding = true;
    TString triggerstring = ""; 
    
	int minncluTPC;
	float dcatoprimxymin;
	Double_t candpTMin,candpTMax, d_maxr;
	Int_t selectD0, selectD0bar;
	
	AliESDtrackCuts* esdTrackCuts = new AliESDtrackCuts("AliESDtrackCuts", "default");
    
 	Int_t kbitDplus;
	Int_t kbitDs;
	Int_t kbitLc;
	Double_t fMassDzero;
	Double_t fMassDplus;
	Double_t fMassDs;
	Double_t fMassLambdaC;
	
	static const Int_t npTBins = 25;
	static const Int_t nCutVars = 11;
	
	Double_t fCuts[npTBins][nCutVars];
                                                        
	char* GetJsonString(const char* jsonFileName, const char* key);
	int GetJsonInteger(const char* jsonFileName, const char* key);
	bool GetJsonBool(const char* jsonFileName, const char* key);
	float GetJsonFloat(const char* jsonFileName, const char* key);
	void ReadJson();
	void MakeJetFinding(AliESDEvent *esd);
	
	void InitDefault();
	
	TList*  fOutput;    //!<!  list of output histos
	
	TH1F* hpt_nocuts;
	TH1F* htgl_nocuts;
	TH1F* hpt_cuts;
	TH1F* hdcatoprimxy_cuts;
	TH1F* htgl_cuts;
	TH1F* hvx;
	TH1F* hvy;
	TH1F* hvz;
	TH1F* hvx3;
	TH1F* hvy3;
	TH1F* hvz3;
	TH1F* hitsmap;

	TH1F* hvertexx;
	TH1F* hvertexy;
	TH1F* hvertexz;

	TH1F* hdecayxyz;
	TH1F* hdecayxy;
	TH1F* hmass0;
	TH1F* hmassP;
	TH1F* hptD0;
	TH1F* hptprong0;
	TH1F* hptprong1;
	TH1F* hd0;
	TH1F* hd0d0;
	TH1F* hImpParErr;
	TH1F* hDecLenErr;
	TH1F* hDecLenXYErr;
	TH1F* hCovPVXX;
	TH1F* hCovSVXX;

	TH1F* hjetpt;
	TH1F* hjetE; // returns the energy component
	TH1F* hjetpx; // returns the x momentum component
	TH1F* hjetpy; // returns the y momentum component
	TH1F* hjetpz; // returns the z momentum component
	TH1F* hjetphi; // returns the azimuthal angle in range 0 . . . 2π
	TH1F* hjetrap; // returns the rapidity
	TH1F* hjetconstituents;
	TH1F* hjetzg;
	TH1F* hjetrg;
	TH1F* hjetnsd;
	
	Bool_t  fReadMC;             // flag for access to MC
	Bool_t  fUsePhysSel;         // flag use/not use phys sel
	Int_t   fTriggerMask;        // mask used in physics selection
	TH1F* fHistNEvents;                //!<!  histo with N of events

	Int_t GetpTBin(Double_t candpT);
	Bool_t GetTrackMomentumAtSecVert(AliESDtrack* tr, AliAODVertex* secVert, Double_t momentum[3], float fBzkG);
	Bool_t SingleTrkCuts(AliESDtrack* trk, AliESDtrackCuts* esdTrackCuts, AliESDVertex* fV1, Double_t fBzkG);
	Bool_t SingleTrkCutsSimple(AliESDtrack* trk, Int_t minclutpc, int ptmintrack, double dcatoprimxymin, AliESDVertex* fV1, Double_t fBzkG);
	Int_t TwoProngSelectionCuts(AliAODRecoDecayHF2Prong* cand, Double_t candpTMin, Double_t candpTMax);
	AliESDVertex* ReconstructSecondaryVertex(AliVertexerTracks* vt, TObjArray* trkArray, AliESDVertex* primvtx, double rmax);
	AliAODVertex* ConvertToAODVertex(AliESDVertex* trkv);
	Int_t SelectInvMassAndPt3prong(TObjArray* trkArray, AliAODRecoDecay* rd4massCalc3);
	AliAODRecoDecayHF2Prong* Make2Prong(TObjArray* twoTrackArray, AliAODVertex* secVert, Double_t fBzkG);
	AliAODRecoDecayHF3Prong* Make3Prong(TObjArray* threeTrackArray, AliAODVertex* secVert, Double_t fBzkG);
  
	/// \cond CLASSIMP
	ClassDef(AliAnalysisTaskHFFindJets, 1);
	/// \endcond
};
#endif
