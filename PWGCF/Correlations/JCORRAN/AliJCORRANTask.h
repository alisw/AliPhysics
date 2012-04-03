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
#include <TVectorT.h>
#include "AliAnalysisTaskSE.h"
//#include "AliAnalysisFilter.h"
//#include "AliMCEvent.h"
#include "AliJRunHeader.h"
#include "AliESDVZERO.h"
#include "AliESDTZERO.h"
//#include "AliESDFMD.h"
#include "AliESDZDC.h"
#include "AliJConst.h"
#include "AliESDpid.h"
#include "AliEMCALGeometry.h"
#include "AliPHOSGeoUtils.h"
#include "AliPIDResponse.h"
#include "AliPIDCombined.h"

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
class AliESDtrackCuts;
class AliESDVZERO;
class AliESDCentrality;

class AliJRunHeader;
class AliMCEvent;
class AliAnalysisFilter;
class AliJTrack;
class AliJEventHeader;



class AliJCORRANTask : public AliAnalysisTaskSE {

 public:
  AliJCORRANTask();
  AliJCORRANTask(const char *name,  TString inputformat);
  AliJCORRANTask(const AliJCORRANTask& ap);   
  AliJCORRANTask& operator = (const AliJCORRANTask& ap);
  virtual ~AliJCORRANTask();

  // methods to fill from AliAnalysisTaskSE
  virtual void UserCreateOutputObjects(); 
  virtual void Init(); 	
  virtual void LocalInit() { Init(); }
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t * opt = "");

  TString GetRunType() const { return fRunType;}
  void SetRunType( const TString type ){ fRunType = type; }
  void SetESDtrackCuts(AliESDtrackCuts* esdTrackCuts){ fEsdTrackCuts = esdTrackCuts;}
  void SetESDFilter( AliAnalysisFilter * filter ){ fESDFilter = filter; }
  void SetRealOrMC(Bool_t realormc){fIsRealOrMC[0]=realormc;} //flags whether the input 
  void SetStoreEventPlaneSource(bool dostore ){ fStoreEventPlaneSource = dostore; }
  bool GetStoreEventPlaneSource(){ return fStoreEventPlaneSource; };
  void SetStoreTPCTrack(bool dostore ){ fStoreTPCTrack = dostore; }
  bool GetStoreTPCTrack(){ return fStoreTPCTrack; };
  //are ESDs from real  exp or MonteCarlo 
  void SetOutputAODName(const char* aodname){ fAODName=aodname;}
  //  AliEMCALGeoUtils* GetEMCALGeoUtils (bool doDelete=kFALSE);

 private:
  AliJEventHeader* ReadCommonHeader(AliVEvent *event);
  // methods to read data from ESD
  void ReadESDTracks(AliESDEvent* esd);
  void ReadESDCaloClusters(const AliESDEvent* esd);
  void ReadESDHeader(AliESDEvent* esd);
  void ReadESDPID(AliESDtrack* track, AliJTrack* ctrack);
  // methods to read data from AOD
  void ReadAODTracks(const AliAODEvent* aod);
  void ReadAODCaloClusters(const AliAODEvent* aod);
  void ReadAODHeader(AliAODEvent* aod);
  void ReadFilter();
  void ReadMCTracks(AliMCEvent* fMC);
  Int_t GetSuperModuleNumber(bool isemcal, Int_t absId);

  UInt_t ConvertTriggerMask();//Converts alice trigger mask to JCorran trigger mask
  //functions used for event selction:
  bool AcceptAODTrack(AliAODTrack* aodTrack);
  void SetOADBPath(const char* path) {fOADBPath=path;}
  const char* GetOADBPath() const { return fOADBPath.Data(); }

  // method to fill jcorran
  bool SetAliceTriggerDef(AliJRunHeader *runHeader);
  bool SetAliceFilterMapDef(AliJRunHeader *runHeader); //TODO Check
  void PrintOut();
  
  // UTILS
  void AddListAODBranch(const char* aname, const char* cname, TClonesArray **obj, int nlist);

  // d a t a     m e m b e r s
  TString fRunType;   // ex) LHC10h
  TString fInputFormat; // specify the input data format (ESD or AOD)
  AliESDtrackCuts* fEsdTrackCuts; //track selection cuts
  AliAnalysisFilter * fESDFilter; //filter set of track selection BS
  TVectorT<double>  fIsRealOrMC; //flags if the input are real (0) ESDs or MonteCarlo ESDs (1)
  TString fAODName; //output delta AOD name
  TString fActiveTriggers[kRangeTriggerTableAlice];//alice table mapping trigger bit to trigger name
  TString fTriggerTableJCorran[kRangeTriggerTableJCorran];//JCorran trigger table TBit 0 =MinBias
  bool fStoreEventPlaneSource;
  bool fStoreTPCTrack;
  TString fOADBPath;

  // jcorran output objects

  TClonesArray *    fTrackList;   // list of charged track objects
  TClonesArray *    fMCTrackList; // list of charged track objects
  TClonesArray *    fPhotonList;  // list of photons objects
  TClonesArray *    fHeaderList;  // event details
  TList *  	    fRunInfoList; // run details

  AliESDpid	      *fPIDesd;
  AliPIDResponse  *fPIDResponse; // PID response object
  AliPIDCombined  *fPIDCombined;

  AliESDVZERO*        fVZEROData;
  AliESDTZERO*        fTZEROData;
  //  AliESDFMD*          fFMDData;
  AliESDZDC*          fZDCData;

  AliJRunHeader*      fAliRunHeader;//  run details (mg field, trigger mask,etc...)
  AliEMCALGeometry * fEMCALGeoUtils; // no AliEMCALGeoUtils.h in trunk aliroot (111130)
  AliPHOSGeoUtils  * fPHOSGeom; //phos geometry matrix 


  ClassDef(AliJCORRANTask, 1); 
};
#endif // AliJCORRANTask_H
