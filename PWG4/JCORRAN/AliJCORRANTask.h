#ifndef ALIJCORRANTASK_H
#define ALIJCORRANTASK_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */
//______________________________________________________________________________
// Analysis task for high pt particle correlations 
// author: R.Diaz, J. Rak,  D.J. Kim, F.Krizek
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
#include "AliAODTrack.h"

#include "AliPhJTrackList.h"
#include "AliPhJMCTrackList.h"
#include "AliPhJPhotonList.h"
#include "AliPhJHeaderList.h"
#include "AliJCORRANSetup.h"

#include "AliJRunHeader.h"
#include "JConst.h"

//==============================================================

using namespace std;

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


  //Setters to be used in AddAliJCORRANTask
  void SetAliESDtrackCuts(AliESDtrackCuts* cuts){ fEsdTrackCuts=cuts;} 
  void SetDownscaling(Int_t ds){ fDownscaling=ds;} 
  void SetLowerCutOnLPMom(Double_t lcut){ fLowerCutOnLPMom=lcut;} 
  void SetLowerCutOnLeadingCaloClusterE(Double_t lowE){ fLowerCutOnLeadingCaloClusterE=lowE;} //
  void SetLowerCutOnCaloClusterE(Double_t lowE){ fLowerCutOnCaloClusterE=lowE;} //
  void SetRealOrMC(Bool_t realormc){fIsRealOrMC=realormc;} //flags whether the input 
                                                         //are ESDs from real  exp or MonteCarlo 
  void SetOutputAODName(const char* aodname){ fAODName=aodname;}

private:

  UInt_t ConvertTriggerMask(/*Long64_t alicetriggermask*/);//Converts alice trigger mask to JCorran trigger mask


  bool IsSelectedAODTrack(AliAODTrack   *track); //check if the track fulfils quality cuts on AOD tracks

  bool StoreDownscaledMinBiasEvent(); //functions for offline event selection 
  bool ContainsESDHighPtTrack();
  bool ContainsESDHighECaloClusters();


  TString fInputFormat; // specify the input data format (ESD or AOD)

  //To be set in from AddAliJCORRANTask macro  
  AliESDtrackCuts* fEsdTrackCuts;     //standard track quality cuts
  Int_t    fDownscaling;              // downscaling
  Double_t fLowerCutOnLPMom;          // lower cut on leading particle momentum
  Double_t fLowerCutOnLeadingCaloClusterE; // lower cut on leading calo cluster energy
  Double_t fLowerCutOnCaloClusterE;   //minimal E of cluster to be stored
  Bool_t   fIsRealOrMC;               //flags if the input are real (0) ESDs or MonteCarlo ESDs (1)
  TString fAODName; //output name
  
  TString fActiveTriggers[kRangeTriggerTableAlice];//alice table mapping trigg. bit to trigg. name
  TString fTriggerTableJCorran[kRangeTriggerTableJCorran];//JCorran trigger table TBit 0 =MinBias
 
 
  //output objects
  AliPhJTrackList*   fTrackList;    //list of charged track objects
  AliPhJMCTrackList* fMCTrackList;  //list of charged track objects
  AliPhJPhotonList*  fPhotonList;   //list of photons objects
  AliPhJHeaderList*  fHeaderList;   //event details

  AliJRunHeader*     fAliRunHeader; // run details
 
  AliPHOSGeoUtils  * fPHOSGeom;  //! phos geometry matrix 
  AliEMCALGeoUtils * fEMCALGeom; //! emcal geometry matrix
   
  ClassDef(AliJCORRANTask, 2) // JCORRAN analysis task 
};
#endif // AliJCORRANTask_H
