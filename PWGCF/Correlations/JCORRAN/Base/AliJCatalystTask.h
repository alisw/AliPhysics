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

//==============================================================

using namespace std;

class TH1D;
class TH2D;
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
	float GetCentrality() const{return fcent;};
	double GetZVertex() const{return fZvert;};
	int GetRunNumber() const{return fRunNum;};
	AliAODEvent * GetAODEvent() const{return paodEvent;}

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
	UInt_t GetEffMode() const{return fEffMode;}
	UInt_t GetEffFilterBit() const{return fEffFilterBit;}
	void SetEtaRange( double eta_min, double eta_max ){
		fEta_min = eta_min; fEta_max = eta_max; cout << "setting Eta ragne as " << fEta_min << " ~ " <<	fEta_max << endl;}
	void SetPtRange( double pt_min, double pt_max){
		fPt_min = pt_min; fPt_max = pt_max; cout << "setting Pt range as " << fPt_min << " ~ " << fPt_max << endl;}
	void SetJCatalystTaskName(TString taskname){fTaskName = taskname;}
	TString GetJCatalystTaskName() const{return fTaskName;}
	void ReadVertexInfo( AliAODEvent *aod , double* fvertex);
	Bool_t IsThisAWeakDecayingParticle(AliAODMCParticle *thisGuy);
	Bool_t IsThisAWeakDecayingParticle(AliMCParticle *thisGuy);
	void SetZVertexCut( double zvtxCut ){ fzvtxCut = zvtxCut;
		cout << "setting z vertex cut = " << fzvtxCut << endl;}
  void SetRemoveBadArea( Bool_t shallweremove ){ fremovebadarea = shallweremove;
		cout << "setting RemoveBadArea = " << fremovebadarea << endl;}
	double GetZVertexCut() const{return fzvtxCut;}
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
		FLUC_ALICE_IPINFO = 0x400
		//FLUC_PHI_CORRECTION  = 0x800,
	};
	void AddFlags(UInt_t nflags){flags |= nflags;}
	Int_t GetJCatalystEntry(){ return fJCatalystEntry; } // in order to sync event id
	bool GetIsGoodEvent(){ return fIsGoodEvent; }
	void SetNoCentralityBin( bool nocent) { fnoCentBin = nocent;}
	AliJCorrectionMapTask *GetAliJCorrectionMapTask() {return fJCorMapTask;}

// Methods to provide QA output and additional selection cuts.
  TList* GetCataList() const {return fMainList;}
  virtual void InitializeArrays();
  virtual void BookControlHistograms();
	virtual void FillControlHistograms(AliAODTrack *thisTrack, Int_t whichHisto, Float_t cent, Double_t *v);
  void SetSaveAllQA(Bool_t SaveQA){ bSaveAllQA = SaveQA; }
  void SetSaveHMOhist (Bool_t SaveHMO) {bSaveHMOhist = SaveHMO;}
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

	UInt_t flags; //
	Int_t fJCatalystEntry; //
	bool fIsGoodEvent; //
	AliJCorrectionMapTask *fJCorMapTask; // Correction Map task
	TString fJCorMapTaskName; //
	TH1 *pPhiWeights;
	TGraphErrors *grEffCor; // for one cent
	TAxis *fCentBinEff; // for different cent bin for MC eff
	UInt_t phiMapIndex; //

// Data members for the QA of the catalyst.
	TList *fMainList;		// Mother list containing all possible output of the catalyst task.
	Bool_t bSaveAllQA;		// if kTRUE: All Standard QA Histograms are saved (default kFALSE).
	Bool_t bSaveHMOhist;	// if kTRUE: Save the TH2D for the HMO in LHC10h (bSaveAllQA must be kTRUE as well).
	Int_t fCentralityBins;		// Set to 16, for at maximum 16 bins in case of centrality 0 to 80 in 5% steps. Less bins and different steps may be used.
	Float_t fcent_0, fcent_1, fcent_2, fcent_3, fcent_4, fcent_5, fcent_6, fcent_7, fcent_8, fcent_9, fcent_10, fcent_11, fcent_12, fcent_13, fcent_14, fcent_15, fcent_16;
  		// fcent_i holds the edge of a centrality bin.
  Float_t fcentralityArray[17];		// Number of centrality bins for the control histograms.
  double fChi2perNDF_min;	// Minimum requirement for chi2/ndf for TPC
	double fChi2perNDF_max;	// Maximum requirement for chi2/ndf for TPC
	double fDCAxy_max;	// Maximum requirement for the DCA in transverse plane.
	double fDCAz_max;	// Maximum requirement for the DCA along the beam axis.

	TList *fControlHistogramsList[16];		//! List to hold all control histograms for a specific centrality bin. Up to 16 centraliy bins possible. 
  TH1F *fPTHistogram[16][2];		//! 0: P_t Before Track Selection, 1: P_t After Track Selection.
  TH1F *fPhiHistogram[16][2];		//! 0: Phi Before Track Selection, 1: Phi After Track Selection.
  TH1F *fEtaHistogram[16][2];		//! 0: Eta Before Track Selection, 1: Eta After Track Selection.
  TH1F *fMultHistogram[16][2];		//! 0: Multiplicity Before Track Selection, 1: Mult. After Track Selection.
  TH1F *fTPCClustersHistogram[16][2];		//! 0: TPC Clusters Before Track Selection, 1: TPC Clusters After Track Selection.
  TH1F *fITSClustersHistogram[16][2];		//! 0: ITS Clusters Before Track Selection, 1: ITS Clusters After Track Selection.
  TH1F *fChiSquareTPCHistogram[16][2];		//! 0: ChiSquare TPC Before Track Selection, 1: ChiSquare TPC After Track Selection.
  TH1F *fDCAzHistogram[16][2];		//! 0: DCAz Before Track Selection, 1: DCAz After Track Selection.
  TH1F *fDCAxyHistogram[16][2];		//! 0: DCAxy Before Track Selection, 1: DCAxy After Track Selection.
  TH1I *fChargeHistogram[16][2];		//! 0: Charge Before Track Selection, 1: Charge After Track Selection.
  TH1F *fCentralityHistogram[16];		//! Centrality After Corresponding Cut.
  TH1F *fVertexXHistogram[16][2];		//! 0: Vertex X Before Corresponding, 1: Vertex X After Corresponding Cut.
  TH1F *fVertexYHistogram[16][2];		//! 0: Vertex Y Before Corresponding, 1: Vertex Y After Corresponding Cut.
  TH1F *fVertexZHistogram[16][2];		//! 0: Vertex Z Before Corresponding, 1: Vertex Z After Corresponding Cut.
  TH2D *fHMOsHistogram[16][2];		//! 0: Correlations between global and TPC tracks before, 1: after HMO cut.

	ClassDef(AliJCatalystTask, 3);
};
#endif // AliJCatalystTask_H
