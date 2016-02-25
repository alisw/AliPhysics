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
#include <fstream>
#include <stdlib.h>
#include <stdio.h>
#include <iomanip>

#include "AliAODMCParticle.h"
#include "AliJHistManager.h"
#include "AliAnalysisTaskSE.h"
#include "AliJConst.h"
#include "AliESDpid.h"
#include "AliPHOSGeoUtils.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"
#include "AliJFFlucAnalysis.h"
#include "AliAnalysisUtils.h"
#include "AliVVertex.h" 

//==============================================================

using namespace std;

class TH1D;
class TH2D;
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
  AliJFFlucTask(const char *name, Bool_t IsMC, Bool_t IsExcludeWeakDecay);
  AliJFFlucTask(const AliJFFlucTask& ap);   
  AliJFFlucTask& operator = (const AliJFFlucTask& ap);
  virtual ~AliJFFlucTask();

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init();   
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t* opt=""  );

  //msong added memeber
  //TClonesArray * ReadAODTracks( AliAODEvent* aod );
  void ReadAODTracks( AliAODEvent* aod, TClonesArray *fInputList);
  void ReadKineTracks( AliMCEvent *mcEvent, TClonesArray *TrackList);
  void SetDebugLevel(int debuglevel){fDebugLevel = debuglevel; cout <<"setting Debug Level = " << fDebugLevel << endl;};
  float ReadAODCentrality( AliAODEvent* aod, TString Trig );
  Bool_t IsGoodEvent( AliAODEvent* aod);
  void SetIsMC( Bool_t ismc){ IsMC = ismc; cout << "Setting IsMC = " << ismc << endl; };
  void SetIsKineOnly( Bool_t iskine){ IsKineOnly = iskine; cout << "Setting IsKineOnly = " << iskine << endl; };
  double GetCentralityFromImpactPar(double ip);
  void SetIsWeakDeacyExclude( Bool_t WeakDecay){ IsExcludeWeakDecay=WeakDecay; cout << "Setting Exclude Weak Decay Particles = " << WeakDecay << endl;}
  void SetTestFilterBit( Int_t FilterBit){ fFilterBit = FilterBit; cout << "Settting TestFilterBit = " << FilterBit << endl; };
  void SetEtaRange( double eta_min, double eta_max ){ fEta_min=eta_min; fEta_max=eta_max; cout << "setting Eta ragne as " << fEta_min << " ~ " << fEta_max << endl;};
  void SetPtRange( double pt_min, double pt_max){ fPt_min = pt_min; fPt_max = pt_max; cout << "setting Pt range as " << fPt_min << " ~ " << fPt_max << endl;};
  inline void DEBUG(int level, TString msg){ if(level < fDebugLevel){ std::cout<< level << "\t" << msg << endl;}};
  void SetFFlucTaskName(TString taskname){fTaskName = taskname;};
  TString GetFFlucTaskName(){return fTaskName;};
  void ReadVertexInfo( AliAODEvent *aod , double* fvertex);
  Bool_t IsThisAWeakDecayingParticle(AliAODMCParticle *thisGuy);
  Bool_t IsThisAWeakDecayingParticle(AliMCParticle *thisGuy);
  void SetIsCentFlat( Bool_t isCentFlat);
  void SetEffConfig( int effMode, int FilterBit );
  void SetInFileName( TString inName); 
  void SetIsPhiModule( Bool_t isphi){ IsPhiModule = isphi ; 
					cout << "setting phi modulation = " << isphi << endl; };
  void SetZVertexCut( double zvtxCut ){ fzvtxCut = zvtxCut; 
					cout << "setting z vertex cut = " << fzvtxCut << endl;};
  void SetSCptdep( Bool_t isSCptdep){ IsSCptdep = isSCptdep;
					cout << "setting : SCpt dep = " << isSCptdep << endl;};
  void SetParticleCharge( int charge ){ fPcharge = charge;
					cout << "setting particle charge = " << charge << endl;};
  void SetSCwithQC( Bool_t isSCwithQC){ IsSCwithQC = isSCwithQC;
					cout << "setting : SC with QC = " << isSCwithQC << endl;};

//  void SetMinTPCNcl( int nTPC ){ fminTPCNcl = nTPC ;
//					cout << "setting Number of Cluster in TPC = " << fNclOfTPC << endl;};



 private:
  TString fTaskName;
  int fDebugLevel;
  int fEvtNum;
  int fFilterBit; 
  int fEffMode;
  int fEffFilterBit;
  int fPcharge;
  double fEta_min;
  double fEta_max;
  double fPt_min;
  double fPt_max;
  double fzvtxCut;
  TString fInFileName;
  Bool_t IsMC;
  Bool_t IsKineOnly;
  Bool_t IsExcludeWeakDecay;
  Bool_t IsCentFlat;
  Bool_t IsPhiModule;
  Bool_t IsSCptdep;
  Bool_t IsSCwithQC;
  TClonesArray * fInputList;  // tracklist 
  TH1D *h_ratio; 
  AliJFFlucAnalysis *fFFlucAna; // analysis code
  TDirectory *fOutput;     // output
  ClassDef(AliJFFlucTask, 1);

  TH1D *h_ModuledPhi[7][2]; // cent7, sub2
 
};
#endif // AliJFFlucTask_H
