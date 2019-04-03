#ifndef ALIJFFLUCTASK_H
#define ALIJFFLUCTASK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis task for high pt particle correlations
// author: R.Diaz, J. Rak,  D.J. Kim
// ALICE Group University of Jyvaskyla
// Finland
//
// Fill the analysis containers for ESD or AOD
// Note: Adapted for AliAnalysisTaskSE
//////////////////////////////////////////////////////////////////////////////

#include <iostream>
//#include <unordered_map>
#include <map>
#include <stdlib.h>
#include <stdio.h>

#include <AliAnalysisTaskSE.h>
#include <AliAODMCParticle.h>
#include "AliJHistManager.h"
#include <AliESDpid.h>
#include <AliPHOSGeoUtils.h>
#include <AliPIDResponse.h>
#include <AliPIDCombined.h>
#include <AliAnalysisUtils.h>
#include <AliVVertex.h>
#include "AliJConst.h"
#include "AliJFFlucAnalysis.h"

//==============================================================

using namespace std;

//class TF1;
class TH1D;
class TH2D;
class TH3D;
class TList;
class TTree;
class AliMCEvent;
class AliAODEvent;
class AliAODTrack;
class AliAnalysisFilter;
class AliJTrack;
class AliJEventHeader;
class TParticle;

class AliJFFlucTask : public AliAnalysisTaskSE {
public:
	AliJFFlucTask();
	AliJFFlucTask(const char *name);
	AliJFFlucTask(const AliJFFlucTask& ap);
	AliJFFlucTask& operator = (const AliJFFlucTask& ap);
	virtual ~AliJFFlucTask();

	// methods to fill from AliAnalysisTaskSE
	virtual void UserCreateOutputObjects();
	virtual void Init();
	virtual void LocalInit() { Init(); }
	virtual void UserExec(Option_t *option);
	virtual void Terminate(Option_t* opt="");

	void ReadAODTracks( AliAODEvent* aod, TClonesArray *fInputList, float fCent);
	void ReadKineTracks( AliMCEvent *mcEvent, TClonesArray *TrackList, float fCent);
	//float ReadAODCentrality( AliAODEvent* aod, TString Trig );
	//float ReadMultSelectionCentrality( AliAODEvent* aod, TString Trig );
	float ReadCentrality(AliAODEvent *aod, TString Trig);
	Bool_t IsGoodEvent(AliAODEvent* aod);
	//void SetIsMC( Bool_t ismc){
		//IsMC = ismc; cout << "Setting IsMC = " << ismc << endl; };
	//void SetIsKineOnly( Bool_t iskine){
		//IsKineOnly = iskine; cout << "Setting IsKineOnly = " << iskine << endl; };
	double GetCentralityFromImpactPar(double ip);
	//void SetIsWeakDeacyExclude( Bool_t WeakDecay){
		//IsExcludeWeakDecay=WeakDecay; cout << "Setting Exclude Weak Decay Particles = " << WeakDecay <<	endl;}
	void SetTestFilterBit( Int_t FilterBit){ fFilterBit = FilterBit; cout << "setting TestFilterBit = " << FilterBit << endl; }
	void SetEtaRange( double eta_min, double eta_max ){
		fEta_min = eta_min; fEta_max = eta_max; cout << "setting Eta range as " << fEta_min << " ~ " <<	fEta_max << endl;}
	void SetPtRange( double pt_min, double pt_max){
		fPt_min = pt_min; fPt_max = pt_max; cout << "setting Pt range as " << fPt_min << " ~ " << fPt_max << endl;}
	void SetFFlucTaskName(TString taskname){fTaskName = taskname;}
	TString GetFFlucTaskName() const{return fTaskName;}
	AliJFFlucAnalysis * GetAnalysis() const{return fFFlucAna;}
	int GetRunNumber() const{return fRunNum;}
	void ReadVertexInfo( AliAODEvent *aod , double* fvertex);
	Bool_t IsThisAWeakDecayingParticle(AliAODMCParticle *thisGuy);
	Bool_t IsThisAWeakDecayingParticle(AliMCParticle *thisGuy);
	void SetEffConfig( int effMode, int FilterBit );
	void EnableCentFlat(const TString);
	void EnablePhiModule(const TString);
	void EnablePhiCorrection(const TString);
	//void SetIsPhiModule( Bool_t isphi){ IsPhiModule = isphi ;
					//cout << "setting phi modulation = " << isphi << endl; }
	void SetZVertexCut( double zvtxCut ){ fzvtxCut = zvtxCut;
					cout << "setting z vertex cut = " << fzvtxCut << endl;}
	double GetZVertexCut() const{return fzvtxCut;}
	//void SetSCptdep( Bool_t isSCptdep){ IsSCptdep = isSCptdep;
					//cout << "setting : SCpt dep = " << isSCptdep << endl;}
	void SetParticleCharge( int charge ){ fPcharge = charge;
					cout << "setting particle charge = " << charge << endl;}
	//void SetSCwithQC( Bool_t isSCwithQC){ IsSCwithQC = isSCwithQC;
					//cout << "setting : SC with QC = " << isSCwithQC << endl;}
	//void SetEbEWeight( Bool_t isEbEWeighted){ IsEbEWeighted = isEbEWeighted;
					//cout << "setting : EbE weight = " << isEbEWeighted << endl;}
	//void SetCutOnOutliers( Bool_t CutOutliers ){ fCutOutliers = CutOutliers;
					//cout << "setting : Cut on Outliers = " << fCutOutliers << endl;}
	//void SetForceToUseALICEIPinfo( Bool_t ALICEIPinfo ){ fALICEIPinfo = ALICEIPinfo;}

	void SetCentDetName( TString CentName ){ fCentDetName = CentName;
					cout << "setting : Cenetrality determination =" << fCentDetName.Data() << endl;}
	void SetQCetaCut( Double_t QC_eta_min, Double_t QC_eta_max){
					fQC_eta_min=QC_eta_min; fQC_eta_max=QC_eta_max;
					cout << "setting : QC eta range " << fQC_eta_min << "~" << fQC_eta_max << endl;}

	enum{
		FLUC_MC = 0x1,
		FLUC_EXCLUDEWDECAY = 0x2,
		FLUC_KINEONLY = 0x4,
		//FLUC_PHI_MODULATION = 0x8,
		//FLUC_PHI_INVERSE = 0x10,
		FLUC_PHI_CORRECTION = 0x10,
		FLUC_PHI_REJECTION = 0x20,
		FLUC_SCPT = 0x40,
		FLUC_EBE_WEIGHTING = 0x80,
		FLUC_CENT_FLATTENING = 0x100,
		FLUC_CUT_OUTLIERS = 0x200,
		FLUC_ALICE_IPINFO = 0x400,
	};
	void AddFlags(UInt_t nflags){
		flags |= nflags;
	}

private:
	 //TF1 *pfOutlierLowCut, *pfOutlierHighCut;
	 TClonesArray *fInputList;  // tracklist
	 TDirectory *fOutput;     // output
	 AliJFFlucAnalysis *fFFlucAna; // analysis code
	 TH1D *h_ratio;
	 TH1D *h_ModuledPhi[CENTN_NAT][2]; // cent7, sub2
	 TFile *pDataFile[2];
	 //std::unordered_map<UInt_t, TH2D *> PhiWeightMap[CENTN];
	 std::map<UInt_t, TH3D *> PhiWeightMap[CENTN_NAT];

	 TString fTaskName;
	 int fEvtNum;
	 int fFilterBit;
	 int fEffMode;
	 int fEffFilterBit;
	 int fPcharge;
	 int fRunNum;
	 unsigned int GlobTracks;
	 unsigned int TPCTracks;
	 unsigned int FB32Tracks;
	 unsigned int FB32TOFTracks;
	 double fEta_min;
	 double fEta_max;
	 double fPt_min;
	 double fPt_max;
	 double fzvtxCut;

	 Double_t fQC_eta_min;
	 Double_t fQC_eta_max;

	 TString fCentDetName;
	 //TString fInFileName;
	 /*Bool_t IsMC;
	 Bool_t IsKineOnly;
	 Bool_t IsExcludeWeakDecay;
	 Bool_t IsCentFlat;
	 Bool_t IsPhiModule;
	 Bool_t IsSCptdep;
	 //Bool_t IsSCwithQC;
	 Bool_t IsEbEWeighted;
	 Bool_t fCutOutliers;
	 Bool_t fALICEIPinfo;*/
	 UInt_t flags;

	 ClassDef(AliJFFlucTask, 1);

};
#endif // AliJFFlucTask_H
