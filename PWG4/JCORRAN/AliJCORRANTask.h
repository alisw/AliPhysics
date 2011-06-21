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

//#include <TTree.h>
//#include <TList.h>
#include "AliAnalysisTaskSE.h"
//#include "AliAnalysisFilter.h"
//#include "AliMCEvent.h"

//#include "AliPhJTrackList.h"
//#include "AliPhJMCTrackList.h"
//#include "AliPhJPhotonList.h"
//#include "AliPhJHeaderList.h"

//#include "AliJRunHeader.h"
#include "JConst.h"

//==============================================================

using namespace std;

const int kMaxDimBuffer = 300;//max length of a line read to a buffe

class TH1D;
class TH2D;
class TNtuple;
class TList;
class TTree;
class TFormula;

class AliMCEvent; 
class AliESDEvent; 
class AliAODEvent; 
class AliAODTrack; 
class AliPHOSGeoUtils; 
class AliEMCALGeometry; 
class AliESDtrackCuts;

class AliJRunHeader;
class AliMCEvent;
class AliAnalysisFilter;
class AliPhJHeaderList;
class AliPhJPhotonList;
class AliPhJMCTrackList;
class AliPhJTrackList;               



class AliJCORRANTask : public AliAnalysisTaskSE {

public:
  AliJCORRANTask();
  AliJCORRANTask(const char *name, TString inputformat);
  AliJCORRANTask(const AliJCORRANTask& ap);   
  AliJCORRANTask& operator = (const AliJCORRANTask& ap);
  virtual ~AliJCORRANTask();
  
  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init(); 	
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t * opt = "");
  
  void SetESDtrackCuts(AliESDtrackCuts* esdTrackCuts){ fEsdTrackCuts = esdTrackCuts;}
  void SetDownScalingOfMB(Int_t downSc){ fDownscaling = downSc;} 
  void SetLeadingPaticleMomCut(Double_t lowLPmom){ fLowerCutOnLPMom = lowLPmom;}
  void SetLowerCutOnCaloClusterE(Double_t lowE){ fLowerCutOnCaloClusterE=lowE;}
  void SetLowerCutOnLeadingCaloClusterE(Double_t lowCaloE){fLowerCutOnLeadingCaloClusterE=lowCaloE;} 
  void SetRealOrMC(Bool_t realormc){fIsRealOrMC=realormc;} //flags whether the input 
                                                         //are ESDs from real  exp or MonteCarlo 
  void SetOutputAODName(const char* aodname){ fAODName=aodname;}

  AliEMCALGeometry* GetEMCALGeoUtils (bool doDelete=kFALSE);

private:
  // methods to read data from ESD
  void ReadESDTracks(const AliESDEvent* esd);
  void ReadESDCaloClusters(const AliESDEvent* esd);
  void ReadESDHeader(const AliESDEvent* esd);
  // methods to read data from AOD
  void ReadAODTracks(const AliAODEvent* aod);
  void ReadAODCaloClusters(const AliAODEvent* aod);
  void ReadAODHeader(const AliAODEvent* aod);
  void ReadFilter();
  void ReadMCTracks(AliMCEvent* fMC);
  Int_t GetSuperModuleNumber(bool isemcal, Int_t absId);

 

  UInt_t ConvertTriggerMask();//Converts alice trigger mask to JCorran trigger mask
  //functions used for event selction:
  bool StoreDownscaledMinBiasEvent(); 
  bool ContainsESDHighPtTrack();
  bool ContainsESDHighECaloClusters();
  bool AcceptAODTrack(AliAODTrack* aodTrack);

  // d a t a     m e m b e r s
  TString fInputFormat; // specify the input data format (ESD or AOD)

  AliESDtrackCuts* fEsdTrackCuts; //track selection cuts

  Int_t fDownscaling; //downscaling of usual MB events

  Double_t fLowerCutOnLPMom; // store all events where there is a particle with pT above this threshold

  Double_t fLowerCutOnLeadingCaloClusterE;// store all events where there is a EMCAL clustre with E above this threshold
  Double_t fLowerCutOnCaloClusterE;//store only clusters above this energy
 
  Bool_t  fIsRealOrMC; //flags if the input are real (0) ESDs or MonteCarlo ESDs (1)

  TString fAODName; //output delta AOD name

  TString fActiveTriggers[kRangeTriggerTableAlice];//alice table mapping trigger bit to trigger name

  TString fTriggerTableJCorran[kRangeTriggerTableJCorran];//JCorran trigger table TBit 0 =MinBias
 
  TFormula *f1CutMaxDCAToVertexXYPtDep;  // pt-dep track-to-vertex cut in max absolute distance in xy-plane
 
  // jcorran output objects
  AliPhJTrackList*    fTrackList;   // list of charged track objects
  AliPhJMCTrackList*  fMCTrackList; // list of charged track objects
  AliPhJPhotonList*   fPhotonList;  // list of photons objects
  AliPhJHeaderList*   fHeaderList;  // event details
  AliJRunHeader*      fAliRunHeader;//  run details (mg field, trigger mask,etc...)

  AliPHOSGeoUtils  * fPHOSGeom; //phos geometry matrix 
   
  ClassDef(AliJCORRANTask, 2); 
};
#endif // AliJCORRANTask_H
