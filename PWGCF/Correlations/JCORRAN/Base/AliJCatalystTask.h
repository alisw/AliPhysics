#ifndef ALIJCATALYSTTASK_H
#define ALIJCATALYSTTASK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis task : Providing basic event information and tracklist
// Period by Peroid event selection is done here
// Track cuts are applied
// User need to get only tracklist and basic event informations
// This code works for AOD and KineOnly files/ 
// author:  D.J. Kim (dong.jo.kim@cern.ch)
// ALICE Group University of Jyvaskyla, Finland
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>

#include "AliAODMCParticle.h"
#include "AliAnalysisTaskSE.h"
#include "AliESDpid.h"
#include "AliAnalysisUtils.h"
#include "AliVVertex.h"
#include "AliStack.h"
#include "AliJCorrectionMapTask.h"
#include "TH1F.h"
#include "TFile.h"
#include "AliEventCuts.h"

//==============================================================

using namespace std;

class TH1D;
class TH2D;
class TH2F;
class TH3D;
class TList;
class TTree;
class AliMCEvent;
class AliStack;
class AliAODEvent;
class AliAODTrack;
class AliAnalysisFilter;
class AliJTrack;
class TParticle;
class TGraphErrors;
class AliJCatalystTask : public AliAnalysisTaskSE {
public:
	AliJCatalystTask();
	AliJCatalystTask(const char *name);
	AliJCatalystTask(const AliJCatalystTask& ap);
	AliJCatalystTask& operator = (const AliJCatalystTask& ap);
	virtual ~AliJCatalystTask();

	// methods to fill from AliAnalysisTaskSE.
	virtual void UserCreateOutputObjects();
	virtual void Init();
	virtual void LocalInit() { Init(); }
	virtual void UserExec(Option_t *option);
	virtual void Terminate(Option_t* opt="");

	inline void DEBUG(int level, TString msg){ if(level < fDebugLevel){ std::cout<< level << "\t" << msg << endl;};}
	// Getters for other analysis tasks
	// Particle list
	TClonesArray * GetInputList() const{return fInputList;}
	TClonesArray * GetInputListALICE() const{return fInputListALICE;}
	// Getters Event Info, centrality, zvertex, runnumber
	inline float GetCentrality() const{return fcent;};
	inline double GetZVertex() const{return fZvert;};
	inline int GetRunNumber() const{return fRunNum;};
	inline AliAODEvent * GetAODEvent() const{return paodEvent;}

	void SetDebugLevel(int debuglevel){
		fDebugLevel = debuglevel; cout <<"setting Debug Level = " << fDebugLevel << endl;}
	float ReadCentrality(AliAODEvent *aod, TString Trig);
	Bool_t IsGoodEvent(AliAODEvent* aod, Int_t thisCent);
	double GetCentralityFromImpactPar(double ip);
	// Read AOD or KineOnly files
	void ReadAODTracks( AliAODEvent* aod, TClonesArray *fInputList, float fCent);
	void ReadKineTracks( AliMCEvent *mcEvent, TClonesArray *TrackList, TClonesArray *TrackListALICE, float fCent);
	void ReadKineTracks( AliStack *stack, TClonesArray *TrackList, TClonesArray *TrackListALICE, float fCent);
	void SetTestFilterBit( Int_t FilterBit){ fFilterBit = FilterBit; cout << "Settting TestFilterBit = " << FilterBit << endl;}
	void SetNumTPCClusters( UInt_t NumTPCClusters){ fNumTPCClusters = NumTPCClusters; }
	void SetEffConfig( int effMode, int FilterBit );
	inline UInt_t GetEffMode() const{return fEffMode;}
	inline UInt_t GetEffFilterBit() const{return fEffFilterBit;}
	void SetEtaRange( double eta_min, double eta_max ){
		fEta_min = eta_min; fEta_max = eta_max; cout << "setting Eta ragne as " << fEta_min << " ~ " <<	fEta_max << endl;}
	inline double GetEtaMin() const{return fEta_min;}
	inline double GetEtaMax() const{return fEta_max;}
	void SetPtRange( double pt_min, double pt_max){
		fPt_min = pt_min; fPt_max = pt_max; cout << "setting Pt range as " << fPt_min << " ~ " << fPt_max << endl;}
	void SetJCatalystTaskName(TString taskname){fTaskName = taskname;}
	inline TString GetJCatalystTaskName() const{return fTaskName;}
	void ReadVertexInfo( AliAODEvent *aod , double* fvertex);
	Bool_t IsThisAWeakDecayingParticle(AliAODMCParticle *thisGuy);
	Bool_t IsThisAWeakDecayingParticle(AliMCParticle *thisGuy);
	void SetZVertexCut( double zvtxCut ){ fzvtxCut = zvtxCut;
		cout << "setting z vertex cut = " << fzvtxCut << endl;}
  void SetRemoveBadArea( Bool_t shallweremove ){ fremovebadarea = shallweremove;
		cout << "setting RemoveBadArea = " << fremovebadarea << endl;}
  void SetRemoveBadArea18q( Bool_t shallweremove ){ fremovebadarea18q = shallweremove;
		cout << "setting RemoveBadArea18q = " << fremovebadarea18q << endl;}
	inline double GetZVertexCut() const{return fzvtxCut;}
	void SetParticleCharge( int charge ){ fPcharge = charge;
		cout << "setting particle charge = " << charge << endl;}
	void SetCentDetName( TString CentName ){ fCentDetName = CentName;
		cout << "setting : Cenetrality determination =" << fCentDetName.Data() << endl;}
	void SetPhiCorrectionIndex(UInt_t id){phiMapIndex = id;} // need for subwagon
	

	enum{
		FLUC_MC = 0x1,
		FLUC_EXCLUDEWDECAY = 0x2,
		FLUC_KINEONLY = 0x4,
		FLUC_KINEONLYEXT = 0x8,
		FLUC_CENT_FLATTENING = 0x100,
		FLUC_CUT_OUTLIERS = 0x200,
		FLUC_ALICE_IPINFO = 0x400,
		FLUC_EXCLUDE_EPOS = 0x800
		//FLUC_PHI_CORRECTION  = 0x800,
	};
	void AddFlags(UInt_t nflags){flags |= nflags;}
	inline Int_t GetJCatalystEntry(){ return fJCatalystEntry; } // in order to sync event id
	inline bool GetIsGoodEvent(){ return fIsGoodEvent; }
	void SetNoCentralityBin( bool nocent) { fnoCentBin = nocent;}
	inline AliJCorrectionMapTask *GetAliJCorrectionMapTask() {return fJCorMapTask;}

// Methods to provide QA output and additional selection cuts.
  inline TList* GetCataList() const {return fMainList;}
  virtual void InitializeArrays();
  virtual void BookControlHistograms();
	virtual void FillControlHistograms(AliAODTrack *thisTrack, Int_t whichHisto, Float_t cent, Double_t *v);
  void SetSaveAllQA(Bool_t SaveQA){ bSaveAllQA = SaveQA; }
  void SetSaveHMOhist (Bool_t SaveHMO) {bSaveHMOhist = SaveHMO;}
  void SetSaveQCNUA(Bool_t SaveQCNUA){ bSaveQCNUA = SaveQCNUA; }
  Int_t GetCentralityBin(Float_t cent);
  void SetCentrality(Float_t cen0, Float_t cen1, Float_t cen2, Float_t cen3, Float_t cen4, Float_t cen5, Float_t cen6, Float_t cen7, Float_t cen8, Float_t cen9, Float_t cen10, Float_t cen11, Float_t cen12, Float_t cen13, Float_t cen14, Float_t cen15, Float_t cen16 ) {fcent_0 = cen0; fcent_1 = cen1; fcent_2 = cen2; fcent_3 = cen3; fcent_4 = cen4; fcent_5 = cen5; fcent_6 = cen6; fcent_7 = cen7; fcent_8 = cen8; fcent_9 = cen9; fcent_10 = cen10; fcent_11 = cen11; fcent_12 = cen12; fcent_13 = cen13; fcent_14 = cen14; fcent_15 = cen15; fcent_16 = cen16;}
  void SetInitializeCentralityArray(); //Set Centrality array inside the task. Must be called in addTask.
  	void SetChi2Cuts(double chiMin, double chiMax) {
		fChi2perNDF_min = chiMin; fChi2perNDF_max = chiMax;
		cout << "setting chi2perNDF cuts, min = " << fChi2perNDF_min << " and max = " << fChi2perNDF_max << endl;
	}
	void SetDCAxyCut(double DCAxyMax) {fDCAxy_max = DCAxyMax;
		cout << "setting DCAxy cut = " << fDCAxy_max << endl;
	}
	void SetDCAzCut(double DCAzMax) {fDCAz_max = DCAzMax;
		cout << "setting DCAz cut = " << fDCAz_max << endl;
	}
	void SetITSCuts(bool UseITSMinClusters, double ITSMinClusters){
		fUseITSMinClusters = UseITSMinClusters;
		fITSMinClusters = ITSMinClusters;
	}

// Methods to apply tighter cuts on Run2.
	void SetTightCuts(bool usePrimary) {fUseTightCuts = usePrimary;}
	void SetDCABaseCuts(bool usePrimary) {bUseDCAbaseCut = usePrimary;}
	void SetESDpileupCuts(bool ESDpileup, double slope, double intercept, bool saveQA) {fAddESDpileupCuts = ESDpileup;
		fESDpileup_slope = slope; fESDpileup_inter = intercept; fSaveESDpileupQA = saveQA;}
	void SetTPCpileupCuts(bool TPCpileup, bool saveQA) {fAddTPCpileupCuts = TPCpileup; fSaveTPCpileupQA = saveQA;}
	void FillEventQA(AliAODEvent *event, int centBin, int stepBin);

// Methods to use alternative correction weights.
	Int_t GetRunIndex10h(Int_t runNumber);
	void SetInputAlternativeNUAWeights10h(bool UseAltWeight, TString fileWeight);

	Int_t GetRunIndex15o(Int_t runNumber);
	void SetInputCentralityWeight15o(bool useCentWeight, TString fileCentWeight);

private:
	TClonesArray * fInputList;  // tracklist
	TClonesArray * fInputListALICE;  // tracklist ALICE acceptance +-0.8 eta
	//TDirectory *fOutput;     // output
	TString fTaskName; //
	TString fCentDetName; //
	AliAODEvent *paodEvent; //
	float fcent; //
	double fZvert; //
	bool fnoCentBin; // no centrality bin => 1
	int fDebugLevel; //
	UInt_t fEvtNum; //
	UInt_t fFilterBit; //
	UInt_t fNumTPCClusters;
	UInt_t fEffMode; //
	UInt_t fEffFilterBit; //
	int fRunNum; //
	int fPcharge; //
	UInt_t GlobTracks; //
	UInt_t TPCTracks; //
	UInt_t FB32Tracks; //
	UInt_t FB32TOFTracks; //
	double fEta_min; //
	double fEta_max; //
	double fPt_min; //
	double fPt_max; //
	double fzvtxCut; //
	Bool_t fremovebadarea; //
	Bool_t fremovebadarea18q; //

	UInt_t flags; //
	Int_t fJCatalystEntry; //
	bool fIsGoodEvent; //
	AliJCorrectionMapTask *fJCorMapTask; // Correction Map task
	TString fJCorMapTaskName; //
	TH3D *pPhiWeights;
	TGraphErrors *grEffCor; // for one cent
	TAxis *fCentBinEff; // for different cent bin for MC eff
	UInt_t phiMapIndex; //
	Bool_t bUseAlternativeWeights; //
	AliEventCuts *fAliEventCuts;	// Instance of AliEventCuts.

// Data members for the QA of the catalyst.
	TList *fMainList;		// Mother list containing all possible output of the catalyst task.
	Bool_t bSaveAllQA;		// if kTRUE: All Standard QA Histograms are saved (default kFALSE).
	Bool_t bSaveHMOhist;	// if kTRUE: Save the TH2D for the HMO in LHC10h (bSaveAllQA must be kTRUE as well).
	Bool_t bSaveQCNUA; 		// if kTRUE: 
	Int_t fCentralityBins;		// Set to 16, for at maximum 16 bins in case of centrality 0 to 80 in 5% steps. Less bins and different steps may be used.
	Float_t fcent_0, fcent_1, fcent_2, fcent_3, fcent_4, fcent_5, fcent_6, fcent_7, fcent_8, fcent_9, fcent_10, fcent_11, fcent_12, fcent_13, fcent_14, fcent_15, fcent_16;
  		// fcent_i holds the edge of a centrality bin.
  Float_t fcentralityArray[17];		// Number of centrality bins for the control histograms.
  double fChi2perNDF_min;	// Minimum requirement for chi2/ndf for TPC
	double fChi2perNDF_max;	// Maximum requirement for chi2/ndf for TPC
	double fDCAxy_max;	// Maximum requirement for the DCA in transverse plane.
	double fDCAz_max;	// Maximum requirement for the DCA along the beam axis.
	bool bDCABaseCut;
	bool fUseITSMinClusters;	// if true use cut for minimum number of ITS clusters
	double fITSMinClusters;		// minimum number of required ITS clusters

// Data members for the use of tighter cuts in Run2.
	bool fUseTightCuts;		// if kTRUE: apply tighter cuts on DCAxy and goldenChi2
	bool bUseDCAbaseCut;	// 
	bool fAddESDpileupCuts;	// if true: apply a cut on the correlations between ESD and TPConly tracks.
	double fESDpileup_slope;	// Slope of the cut M_ESD >= 15000 + 3.38*M_TPC
	double fESDpileup_inter;	// Intercept of the cut.
	bool fSaveESDpileupQA;	// if true: save the TH2D for the QA.
	bool fAddTPCpileupCuts;	// if true: apply a cut on the correlations between ITS and TPC clusters.
	bool fSaveTPCpileupQA;	// if true: save the TH2D for the QA.

  TList *fControlHistogramsList[16];		//! List to hold all control histograms for a specific centrality bin. Up to 16 centraliy bins possible.
  TList *fControl2DNUAList[16]; 	//!
  TList *fControlProfileList;		//!
  TH1F *fPTHistogram[16][3];		//! 0: P_t Before Track Selection, 1: P_t After Track Selection, 2: After correction.
  TH1F *fPhiHistogram[16][3];		//! 0: Phi Before Track Selection, 1: Phi After Track Selection, 2: after correction.
  TH1F *fEtaHistogram[16][2];		//! 0: Eta Before Track Selection, 1: Eta After Track Selection.
  TH1F *fMultHistogram[16][2];		//! 0: Multiplicity Before Track Selection, 1: Mult. After Track Selection.
  TH1F *fTPCClustersHistogram[16][2];		//! 0: TPC Clusters Before Track Selection, 1: TPC Clusters After Track Selection.
  TH1F *fITSClustersHistogram[16][2];		//! 0: ITS Clusters Before Track Selection, 1: ITS Clusters After Track Selection.
  TH1F *fChiSquareTPCHistogram[16][2];		//! 0: ChiSquare TPC Before Track Selection, 1: ChiSquare TPC After Track Selection.
  TH1F *fChiSquareITSHistogram[16][2];		//! 0: ChiSquare ITS Before Track Selection, 1: ChiSquare ITS After Track Selection.  
  TH1F *fDCAzHistogram[16][2];		//! 0: DCAz Before Track Selection, 1: DCAz After Track Selection.
  TH1F *fDCAxyHistogram[16][2];		//! 0: DCAxy Before Track Selection, 1: DCAxy After Track Selection.
  TH1I *fChargeHistogram[16][2];		//! 0: Charge Before Track Selection, 1: Charge After Track Selection.
  TH1F *fCentralityHistogram[16];		//! Centrality After Corresponding Cut.
  TH1F *fVertexXHistogram[16][2];		//! 0: Vertex X Before Corresponding, 1: Vertex X After Corresponding Cut.
  TH1F *fVertexYHistogram[16][2];		//! 0: Vertex Y Before Corresponding, 1: Vertex Y After Corresponding Cut.
  TH1F *fVertexZHistogram[16][2];		//! 0: Vertex Z Before Corresponding, 1: Vertex Z After Corresponding Cut.
  TH2D *fHMOsHistogram[16][2];		//! 0: Correlations between global and TPC tracks before, 1: after HMO cut.
  TH1F *fHistoPhiWeight[16][90];	// Histograms to save the NUA correction weights per centrality and runs.
  TProfile *fProfileWeights[16];	//! Profiles for the weights to apply per phi bins.
  TH2D *fESDpileupHistogram[16][2];		//! 0: Correlations between ESD and TPC tracks before, 1: after cut.
  TH2I *fTPCpileupHistogram[16][2];		//! 0: Correlations between ITS and TPC clusters before, 1: after cut.

  TH1F *fHistoCentWeight[138];		//! Histograms to save the centrality correction for 15o per run.

  TH2F *f2DEtaPhiHistogram[16][2];  //! 0: Eta-Phi Before correction, 1: Eta-Phi after correction
  TProfile *fProfileCosVSCent[7];		//! 0: n=2, 1: n=3, 2: n=4, 3:n=5, 4:n=6, 5:n=7, 6:n=8
  TProfile *fProfileSinVSCent[7];		//! 0: n=2, 1: n=3, 2: n=4, 3:n=5, 4:n=6, 5:n=7, 6:n=8

  ClassDef(AliJCatalystTask, 9);
};
#endif // AliJCatalystTask_H