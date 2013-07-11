/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


/*__|______________________________________________________________________________|
 |                              -----Info(i)-----                                  |
 |                                                                                 |
 |  .h of single track  efficiecy class                                            |
 |                                                                                 |
 |   ESDs<-->AODs (ON/OFF)                                                         |
 |                                                      Authors:                   |
 |_____________________________________________________________________________|___*/


#ifndef AliCFSINGLETRACKEFFICIENCYTASK_H
#define AliCFSINGLETRACKEFFICIENCYTASK_H

#include "AliAnalysisTaskSE.h"
#include "AliSingleTrackEffCuts.h"

class TH1I;
class TParticle ;
class TFile ;
class AliStack ;
class AliCFManager;
class AliAODTrack;
class AliESDtrack;
class AliESDtrackCuts;
class AliESDVertex;
class AliVVertex;
class AliVParticle;


class AliCFSingleTrackEfficiencyTask : public AliAnalysisTaskSE {
  public:
    
  enum {
    kStepMCGenCut         = 0,
    kStepMCKineCut        = 1,
    kStepMCAccpCut        = 2,
    kStepReconstructed    = 3,
    kStepRecoKineCuts     = 4,
    kStepReconstructedMC  = 5,
    kStepRecoQualityCuts  = 6,
  };

  AliCFSingleTrackEfficiencyTask();
  AliCFSingleTrackEfficiencyTask(const Char_t* name, AliESDtrackCuts *trackcuts, AliSingleTrackEffCuts * mccuts);
  AliCFSingleTrackEfficiencyTask& operator= (const AliCFSingleTrackEfficiencyTask& c);
  AliCFSingleTrackEfficiencyTask(const AliCFSingleTrackEfficiencyTask& c);
  virtual ~AliCFSingleTrackEfficiencyTask();

  // ANALYSIS FRAMEWORK STUFF to loop on data and fill output objects
  void     UserCreateOutputObjects();
  void     UserExec(Option_t *option);
  void     Terminate(Option_t *);
  void     Init();
  
  // CORRECTION FRAMEWORK RELATED FUNCTIONS
  void           SetCFManager(AliCFManager* io) {fCFManager = io;}   // global correction manager
  AliCFManager * GetCFManager() const {return fCFManager;}           // get corr manager
  void           SetQAList(TList* list) {fQAHistList = list;}

    
  // Data types
  Bool_t IsReadTPCTracks() const {return fReadTPCTracks;}
  Bool_t IsReadAODData()   const {return fReadAODData;}

  //Setters
  void   SetReadTPCTracks (Bool_t flag=kTRUE) {fReadTPCTracks=flag;}
  void   SetReadAODData   (Bool_t flag=kTRUE) {fReadAODData=flag;}
  void   SetFilterBit  (Bool_t flag=kTRUE) {fSetFilterBit=flag;}
  void   SetFilterType (Int_t fbittype) { fbit= fbittype;}
  void   SetTriggerMask(ULong64_t mask=0) { fTriggerMask=mask; }
  void   SetMinNClustersTPCPID(Int_t cl=0) { fMinNclsTPCPID=cl;}
  void   SetRequireTOF(Bool_t tof=kTRUE) {fRequireTOF=tof;}
  void   SetMinRatioTPCclusters(Double_t ratio) {fMinRatioTPCcluster=ratio;}

  //Getters
  ULong64_t GetTriggerMask(){ return fTriggerMask; }
  AliESDtrackCuts *GetTrackCuts(){ return (AliESDtrackCuts*)fTrackCuts; }  
  AliSingleTrackEffCuts *GetSingleTrackEffCuts(){ return (AliSingleTrackEffCuts*)fMCCuts; }


 protected:

  void CheckESDParticles();
  void CheckAODParticles();
  AliESDtrack* ConvertTrack(AliAODTrack *track);

  Bool_t          fReadTPCTracks ; // flag to loop on TPC tracks only (ESD mode only)
  Bool_t          fReadAODData ;   // flag for AOD/ESD input files
  AliCFManager   *fCFManager    ;  // pointer to the CF manager
  TList          *fQAHistList   ;  // list of QA histograms

  AliESDtrackCuts *fTrackCuts;      // track cuts (reconstructed level)
  ULong64_t fTriggerMask;           // trigger mask
  AliSingleTrackEffCuts *fMCCuts;   // Cuts used

  Bool_t          fSetFilterBit ; // 
  Int_t           fbit ;   // 

  Int_t           fMinNclsTPCPID;   // Number of clusters for TPC PID
  Bool_t          fRequireTOF;      // whether or not to require TOF matching of the track
  Double_t        fMinRatioTPCcluster; //Min ratio of TPC clusters found/findable


  //Number of events
  TH1I  *fHistEventsProcessed;     //! simple histo for monitoring the number of events processed
  TH1F  *fElectronPt;              //! histo with final selected electron pt
  TH1F  *fElectronPtStart;         //! histo with final selected electron pt

  ClassDef(AliCFSingleTrackEfficiencyTask,1);

};

#endif
