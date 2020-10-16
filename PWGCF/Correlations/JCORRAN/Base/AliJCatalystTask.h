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
class AliJCatalystTask : public AliAnalysisTaskSE {
public:
	AliJCatalystTask();
	AliJCatalystTask(const char *name);
	AliJCatalystTask(const AliJCatalystTask& ap);
	AliJCatalystTask& operator = (const AliJCatalystTask& ap);
	virtual ~AliJCatalystTask();

	// methods to fill from AliAnalysisTaskSE
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
	Bool_t IsGoodEvent(AliAODEvent* aod);
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
	double GetZVertexCut() const{return fzvtxCut;}
	void SetParticleCharge( int charge ){ fPcharge = charge;
		cout << "setting particle charge = " << charge << endl;}
	void SetCentDetName( TString CentName ){ fCentDetName = CentName;
		cout << "setting : Cenetrality determination =" << fCentDetName.Data() << endl;}
	enum{
		FLUC_MC = 0x1,
		FLUC_EXCLUDEWDECAY = 0x2,
		FLUC_KINEONLY = 0x4,
		FLUC_KINEONLYEXT = 0x8,
		FLUC_CENT_FLATTENING = 0x100,
		FLUC_CUT_OUTLIERS = 0x200,
		FLUC_ALICE_IPINFO = 0x400,
	};
	void AddFlags(UInt_t nflags){flags |= nflags;}
	Int_t GetJCatalystEntry(){ return fJCatalystEntry; } // in order to sync event id
	bool GetIsGoodEvent(){ return fIsGoodEvent; }
	void SetNoCentralityBin( bool nocent) { fnoCentBin = nocent;}
	UInt_t ConnectInputContainer(const TString, const TString);
	void EnablePhiCorrection(const TString);
	void EnableCentFlattening(const TString);
	TH1 * GetCorrectionMap(UInt_t, UInt_t);
	TH1 * GetCentCorrection();

private:
	TClonesArray * fInputList;  // tracklist
	TClonesArray * fInputListALICE;  // tracklist ALICE acceptance +-0.8 eta
	TDirectory *fOutput;     // output
	TString fTaskName; //
	TString fCentDetName; //
	AliAODEvent *paodEvent; //
	std::map<UInt_t, TH1 *> PhiWeightMap[96];
	float fcent; //
	double fZvert; //
	bool fnoCentBin; // no centrality bin => 1
	int fDebugLevel; //
	UInt_t fEvtNum; //
	UInt_t fFilterBit; //
	UInt_t fNumTPCClusters;
	UInt_t fEffMode; //
	UInt_t fEffFilterBit; //
	int fPcharge; //
	int fRunNum; //
	UInt_t GlobTracks; //
	UInt_t TPCTracks; //
	UInt_t FB32Tracks; //
	UInt_t FB32TOFTracks; //
	double fEta_min; //
	double fEta_max; //
	double fPt_min; //
	double fPt_max; //
	double fzvtxCut; //

	UInt_t flags; //
	Int_t fJCatalystEntry; //
	bool fIsGoodEvent; //

	UInt_t inputIndex;
	UInt_t phiInputIndex;
	UInt_t centInputIndex;

	ClassDef(AliJCatalystTask, 1);

};
#endif // AliJCatalystTask_H
