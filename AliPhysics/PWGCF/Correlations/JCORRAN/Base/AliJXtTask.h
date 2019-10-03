#ifndef ALIJXTTASK_H
#define ALIJXTTASK_H

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

// root classes
#include <TString.h>
#include <TClonesArray.h>

// AliRoot classes
#include <AliAnalysisTaskSE.h>
#include "AliJEfficiency.h"
 
class TDirectory;
class AliJXtAnalysis;
class AliJRunTable;

//==============================================================

using namespace std;

class AliAODEvent; 

class AliJXtTask : public AliAnalysisTaskSE {

 public:
  AliJXtTask();
  AliJXtTask(const char *name, int CollisionCandidates, Bool_t IsMC);
  AliJXtTask(const AliJXtTask& ap);   
  AliJXtTask& operator = (const AliJXtTask& ap);
  
  virtual ~AliJXtTask();

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init();   
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t* opt="" );

 public:
  double GetInclusiveEfficiency(double pT){ return fEfficiency->GetCorrection(pT, fTrackFilterBit, fCent ); };
  double GetIsolatedEfficiency(double pT){ return fEfficiencyIsolated->GetCorrection(pT, fTrackFilterBit, fCent ); };
  
  //msong added member
  void ReadAODTracks( AliAODEvent* aod, TClonesArray *fInputList);
  void SetDebugLevel(int debuglevel){fDebugLevel = debuglevel; cout <<"setting Debug Level = " << fDebugLevel << endl;};
  float ReadAODCentrality( AliAODEvent* aod, TString Trig );
  void SetIsMC( Bool_t ismc){ fIsMC = ismc; cout << "Settint IsMC = " << ismc << endl; };
  void SetTestFilterBit( Int_t FilterBit){ fTrackFilterBit = FilterBit; cout << "Setting TestFilterBit = " << FilterBit << endl; };
  inline void DEBUG(int level, TString msg){ if(level < fDebugLevel){ std::cout<< level << "\t" << msg << endl;}};
  void SetXtTaskName(TString taskname){fTaskName = taskname;};
  TString GetXtTaskName(){return fTaskName;};
  void ReadVertexInfo( AliAODEvent *aod , double* fvertex);
   
  void SetEtaRange( double etaRange ){ fetaRange = etaRange; };
  void SetZVertexRange( double zvertRange ){ fzvertexRange = zvertRange; };
  
  double GetEtaRange(){ return fetaRange; };
  double GetZVertexRange(){ return fzvertexRange; };
  
  void GetEfficiencyFilterBit(int inputTrackCut );
  
  Bool_t IsGoodEvent( AliAODEvent *event );
  
 private:
  TString fTaskName;	// task name
  double fetaRange;	// allowed eta range, no fiducual cut here
  double fzvertexRange;	// vertex z-range
  int fDebugLevel;	// if != 0, debugging info
  int fEvtNum;		// event number
  int fTriggerBit; 	// choose trigger
  int fTrackFilterBit;	// track selection cut
  int fEffFilterBit;	// bit in the JEfficiency, comes from the track cut
  int frunNumber;	// run number, used by AliJRunTable
  int fEffMode;		// efficiency mode
  double fSQRTS;	// center of mass energy
  double fCent;		// centrality
  Bool_t fFirstEvent;	// is the first event
  Bool_t fIsMC;		// is Monte-Carlo or real data
  Bool_t fIsPP;		// is proton-proton? If yes, no centrality
  AliJRunTable * fRunTable;	//! run information
  AliJEfficiency * fEfficiency; //! efficiency
  AliJEfficiency * fEfficiencyIsolated; //! isolated efficiency
  TClonesArray * fInputList;	//! tracklist  
  AliJXtAnalysis *fXtAna;	//! analysis code
  TDirectory *fOutput;		//! output
  
  ClassDef(AliJXtTask, 1);
};
#endif // ALIJXTTASK_H

