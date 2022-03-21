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

class AliAnalysisTaskHFFindJets : public AliAnalysisTaskSE
{

 public:
	AliAnalysisTaskHFFindJets(); // class constructor
	AliAnalysisTaskHFFindJets(const char *name); // class constructor
	virtual ~AliAnalysisTaskHFFindJets(); // class destructor
	
	virtual void UserCreateOutputObjects(); // called once at beginning of runtime
	virtual void UserExec(Option_t* option); // called for each event
	virtual void Terminate(Option_t* option); // called at end of analysis
	
	void InitFromJson(TString esdFile);
	TTree* tree;

 private:                                
    double ptmintrack = 0.;
    int do3Prongs = 0;
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
    //m     dca   cost* ptk  ptpi  d0k            d0pi         d0d0     cosp cosxy normdxy
	Double_t fCuts[npTBins][nCutVars] = {{0.400, 350. * 1E-4, 0.8, 0.5, 0.5, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.80, 0., 0.},   // pt<0.5
                                     {0.400, 350. * 1E-4, 0.8, 0.5, 0.5, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.80, 0., 0.},   // 0.5<pt<1
                                     {0.400, 300. * 1E-4, 0.8, 0.4, 0.4, 1000. * 1E-4, 1000. * 1E-4, -25000. * 1E-8, 0.80, 0., 0.},  // 1<pt<1.5 
                                     {0.400, 300. * 1E-4, 0.8, 0.4, 0.4, 1000. * 1E-4, 1000. * 1E-4, -25000. * 1E-8, 0.80, 0., 0.},  // 1.5<pt<2 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -20000. * 1E-8, 0.90, 0., 0.},  // 2<pt<2.5 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -20000. * 1E-8, 0.90, 0., 0.},  // 2.5<pt<3 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -12000. * 1E-8, 0.85, 0., 0.},  // 3<pt<3.5 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -12000. * 1E-8, 0.85, 0., 0.},  // 3.5<pt<4
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   // 4<pt<4.5 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   // 4.5<pt<5 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   // 5<pt<5.5 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   // 5.5<pt<6 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   // 6<pt<6.5 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -8000. * 1E-8, 0.85, 0., 0.},   // 6.5<pt<7 
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -7000. * 1E-8, 0.85, 0., 0.},   // 7<pt<7.5
                                     {0.400, 300. * 1E-4, 0.8, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -7000. * 1E-8, 0.85, 0., 0.},   // 7.5<pt<8 
                                     {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0.},   // 8<pt<9 
                                     {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0.},   // 9<pt<10 
                                     {0.400, 300. * 1E-4, 0.9, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, -5000. * 1E-8, 0.85, 0., 0.},   // 10<pt<12 
                                     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 10000. * 1E-8, 0.85, 0., 0.},   // 12<pt<16 
                                     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  // 16<pt<20 
                                     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  // 20<pt<24 
                                     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  // 24<pt<36 
                                     {0.400, 300. * 1E-4, 1.0, 0.7, 0.7, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.85, 0., 0.},  // 36<pt<50 
                                     {0.400, 300. * 1E-4, 1.0, 0.6, 0.6, 1000. * 1E-4, 1000. * 1E-4, 999999. * 1E-8, 0.80, 0., 0.}}; // pt>50  ;
                                                        
	char* GetJsonString(const char* jsonFileName, const char* key);
	int GetJsonInteger(const char* jsonFileName, const char* key);
	bool GetJsonBool(const char* jsonFileName, const char* key);
	float GetJsonFloat(const char* jsonFileName, const char* key);
	void ReadJson();
	
	
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
	TH1F* hjetzg;
	TH1F* hjetrg;
	TH1F* hjetnsd;

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
