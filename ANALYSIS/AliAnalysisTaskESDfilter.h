#ifndef ALIANALYSISTASKESDFILTER_H
#define ALIANALYSISTASKESDFILTER_H
 
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TList.h> 
#include <TF1.h> 
#include "AliAnalysisTaskSE.h"
#include "AliESDtrack.h"
#include "AliAODTrack.h"
#include "AliAODPid.h"
#include "AliESDpid.h"

class AliAnalysisFilter;
class AliStack;
class AliESDtrack;
class AliMCEventHandler;
class TRefArray;
class AliAODHeader;
class AliESDtrackCuts;

class AliAnalysisTaskESDfilter : public AliAnalysisTaskSE
{
 public:
    AliAnalysisTaskESDfilter();
    AliAnalysisTaskESDfilter(const char* name);
    virtual ~AliAnalysisTaskESDfilter();
    // Implementation of interface methods
    virtual void   UserCreateOutputObjects();
    virtual void   Init();
    virtual void   LocalInit() {Init();}
    virtual void   UserExec(Option_t *option);
    virtual void   Terminate(Option_t *option);

    virtual void ConvertESDtoAOD();
    // Setters
    virtual void SetTrackFilter   (AliAnalysisFilter*   trackF) {fTrackFilter    =   trackF;}
    virtual void SetTPCOnlyFilterMask (UInt_t filterMask)       {fTPCOnlyFilterMask    =  filterMask;}
    virtual void SetKinkFilter    (AliAnalysisFilter*    KinkF) {fKinkFilter     =    KinkF;}
    virtual void SetV0Filter      (AliAnalysisFilter*      V0F) {fV0Filter       =      V0F;}
    virtual void SetCascadeFilter (AliAnalysisFilter* CascadeF) {fCascadeFilter  = CascadeF;}
    virtual void SetPthreshold    (Double_t p)                  {fHighPthreshold =        p;}
    virtual void SetPshape        (TF1 *func)                   {fPtshape        =     func;}
    virtual void SetEnableFillAOD (Bool_t b)                    {fEnableFillAOD  =     b;}

    virtual void SetAODPID(AliESDtrack *esdtrack, AliAODTrack *aodtrack, AliAODPid *detpid);
    void SetDetectorRawSignals(AliAODPid *aodpid, AliESDtrack *track);

  void PrintTask(Option_t *option="all", Int_t indent=0) const;
  
  void DisableVZERO() { fIsVZEROEnabled = kFALSE; }
  void DisableCascades() { fAreCascadesEnabled = kFALSE; }
  void DisableV0s() { fAreV0sEnabled = kFALSE; }
  void DisableKinks() { fAreKinksEnabled = kFALSE; }
  void DisableTracks() { fAreTracksEnabled = kFALSE; }
  void DisablePmdClusters() { fArePmdClustersEnabled = kFALSE; }
  void DisableCaloClusters() { fAreCaloClustersEnabled = kFALSE; }
  void DisableCells() { fAreEMCALCellsEnabled = fArePHOSCellsEnabled = kFALSE; }
  void DisableTracklets() { fAreTrackletsEnabled = kFALSE; }

  virtual void SetTimeZeroType(AliESDpid::EStartTimeType_t tofTimeZeroType) {fTimeZeroType = tofTimeZeroType;}
  
private:
    AliAnalysisTaskESDfilter(const AliAnalysisTaskESDfilter&);
    AliAnalysisTaskESDfilter& operator=(const AliAnalysisTaskESDfilter&);
    void PrintMCInfo(AliStack *pStack,Int_t label); // for debugging
    Double_t Chi2perNDF(AliESDtrack* track);
    
  AliAODHeader* ConvertHeader(const AliESDEvent& esd);
  void ConvertCascades(const AliESDEvent& esd);
  void ConvertV0s(const AliESDEvent& esd);
  void ConvertKinks(const AliESDEvent& esd);
  void ConvertPrimaryVertices(const AliESDEvent& esd);
  void ConvertTracks(const AliESDEvent& esd);
  void ConvertPmdClusters(const AliESDEvent& esd);
  void ConvertCaloClusters(const AliESDEvent& esd);
  void ConvertEMCALCells(const AliESDEvent& esd);
  void ConvertPHOSCells(const AliESDEvent& esd);
  void ConvertTracklets(const AliESDEvent& esd);
  void ConvertTPCOnlyTracks(const AliESDEvent& esd);
  void ConvertVZERO(const AliESDEvent& esd);
  
  TClonesArray& Tracks();
  TClonesArray& V0s();
  TClonesArray& Vertices();
  TClonesArray& Cascades();
  
  // Filtering
  AliAnalysisFilter* fTrackFilter;      //  Track   Filter
  AliAnalysisFilter* fKinkFilter;       //  Kink    Filter
  AliAnalysisFilter* fV0Filter;         //  V0      Filter
  AliAnalysisFilter* fCascadeFilter;    //  Cascade Filter
  // PID
  Double_t     fHighPthreshold;    //  Pt threshold for detector signal setting
  TF1 *        fPtshape;           //  Pt spectrum distribution
  Bool_t       fEnableFillAOD;     //  value that decides if this task activates AOD filling
  Bool_t* fUsedTrack; //! indices of used tracks
  Bool_t* fUsedKink; //! indices of used kinks
  Bool_t* fUsedV0; //! indices of used V0s
  TRefArray* fAODTrackRefs; // array of track references
  TRefArray* fAODV0VtxRefs; // array of v0 vertices references
  TRefArray* fAODV0Refs ; // array of v0s references
  AliMCEventHandler* fMChandler; // pointer to MC handler (if any)
  Int_t fNumberOfTracks; // current number of tracks
  Int_t fNumberOfPositiveTracks; // current number of positive tracks
  Int_t fNumberOfV0s; // current number of v0s
  Int_t fNumberOfVertices; // current number of vertices
  Int_t fNumberOfCascades; // current number of cascades
  Int_t fNumberOfKinks; // current number of kinks
  Bool_t fOldESDformat; // is the ESD in old format ?
  AliAODVertex* fPrimaryVertex; // pointer to primary vertex of the event
  UInt_t fTPCOnlyFilterMask; //  Filter Mask used to select and store refitted TPC only tracks
  Bool_t fIsVZEROEnabled; // whether or not to fill the vzero branch (true by default)
  Bool_t fAreCascadesEnabled; // whether or not to fill the cascades branch (true by default)
  Bool_t fAreV0sEnabled; // whether or not to fill the v0 branch (true by default)
  Bool_t fAreKinksEnabled; // whether or not to fill the kinks (true by default)
  Bool_t fAreTracksEnabled; // whether or not to fill the (central) tracks branch (true by default)
  Bool_t fArePmdClustersEnabled; // whether or not to fill the pmd clusters (true by default)
  Bool_t fAreCaloClustersEnabled; // whether or not to fill the calo clusters (true by default)
  Bool_t fAreEMCALCellsEnabled; // whether or not to fill the emcal cells (true by default)
  Bool_t fArePHOSCellsEnabled; // whether or not to fill the phos cells (true by default)
  Bool_t fAreTrackletsEnabled; // whether or not to fill the tracklets (true by default)
  AliESDpid* fESDpid; // esd pid
  Bool_t fIsPidOwner; // whether we own fESDpid
  Int_t fTimeZeroType;  //  time zero type 
  AliESDtrackCuts* fTPCaloneTrackCuts; // TPC stand-alone track cuts
  
  ClassDef(AliAnalysisTaskESDfilter, 9); // Analysis task for standard ESD filtering
};
 
#endif
