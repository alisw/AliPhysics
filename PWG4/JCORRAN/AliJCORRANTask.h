#ifndef ALIJCORRANTASK_H
#define ALIJCORRANTASK_H

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

#include <TTree.h>
#include <TList.h>
#include "AliAnalysisTaskSE.h"
#include "AliAnalysisFilter.h"
#include "AliMCEvent.h"

#include "AliPhJTrackList.h"
#include "AliPhJMCTrackList.h"
#include "AliPhJPhotonList.h"
#include "AliPhJHeaderList.h"
#include "AliJCORRANSetup.h"

#include "AliJRunHeader.h"
#include "JConst.h"

//==============================================================

using namespace std;

const int kMaxDimBuffer = 300;//max length of a line read to a buffe

class TH1D;
class TH2D;
class TNtuple;
class TList;
class TTree;

class AliMCEvent; 
class AliESDEvent; 
class AliAODEvent; 
class AliPHOSGeoUtils; 
class AliEMCALGeoUtils; 
class AliESDtrackCuts;

class AliJRunHeader;
class AliMCEvent;
class AliJCORRANSetup;
class AliAnalysisFilter;
class AliPhJHeaderList;
class AliPhJPhotonList;
class AliPhJMCTrackList;
class AliPhJTrackList;               



class AliJCORRANTask : public AliAnalysisTaskSE {

public:
  AliJCORRANTask() ;
  AliJCORRANTask(const char *name, TString inputformat, AliESDtrackCuts* esdTrackCuts, Int_t downSc, Double_t lowLPmom, Double_t lowCaloE); //FK//
  AliJCORRANTask(const AliJCORRANTask& ap);   
  AliJCORRANTask& operator = (const AliJCORRANTask& ap);
  virtual ~AliJCORRANTask();
  
  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init(); 	
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t * opt = "");
  
  // methods to read data from ESD
  void ReadESDTracks(const AliESDEvent* esd);
  void ReadESDCaloClusters(const AliESDEvent* esd);
  void ReadESDHeader(const AliESDEvent* esd);
  // methods to read data from AOD
  void ReadAODTracks(const AliAODEvent* aod);
  void ReadAODCaloClusters(const AliAODEvent* aod);
  void ReadAODHeader(const AliAODEvent* aod);
  void ReadFilter();
  //void ReadMCTracks(AliMCEvent* fMC);
  Int_t GetSuperModuleNumber(bool isemcal, Int_t absId);

 
private:

  UInt_t ConvertTriggerMask(/*Long64_t alicetriggermask*/);//Converts alice trigger mask to JCorran trigger mask

  bool StoreDownscaledMinBiasEvent();
  bool ContainsESDHighPtTrack(const AliESDEvent* esd);
  bool ContainsESDHighECaloClusters(const AliESDEvent* esd);

  TString fInputFormat; // specify the input data format (ESD or AOD)

  AliESDtrackCuts* fEsdTrackCuts; //FK//

  Int_t fDownscaling; //FK//

  Double_t fLowerCutOnLPMom; //FK//

  Double_t fLowerCutOnLeadingCaloClusterE;//FK//

  TString fActiveTriggers[kRangeTriggerTableAlice];//alice table mapping trigger bit to trigger name

  TString fTriggerTableJCorran[kRangeTriggerTableJCorran];//JCorran trigger table TBit 0 =MinBias
  
  // jcorran output objects
  TTree*         fTree;        // output tree
  AliPhJTrackList*    fTrackList;  // list of charged track objects
  //AliPhJMCTrackList*    fMCTrackList;  // list of charged track objects
  AliPhJPhotonList*   fPhotonList; // list of photons objects
  AliPhJHeaderList*   fHeaderList; // run, event details

  AliJRunHeader* fAliRunHeader; // run, event details

 
  // QA output 
  TList*    fQAList;        // list to hold all the qa objects
  TNtuple*  fTrackQACuts ;  // tree of track quality cuts applied on ESD 
  TNtuple*  fTrackKineCuts; // tree of track kinematic cuts applied on ESD

  AliPHOSGeoUtils  * fPHOSGeom; //phos geometry matrix 
  AliEMCALGeoUtils * fEMCALGeom; //emcal geometry matrix
   
  ClassDef(AliJCORRANTask, 1); // JCORRAN analysis task 
};
#endif // AliJCORRANTask_H
