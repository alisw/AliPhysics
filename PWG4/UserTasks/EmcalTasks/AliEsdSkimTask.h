#ifndef ALIESDSKINTASK_H
#define ALIESDSKINTASK_H

// $Id$

#include "AliPhysicsSelectionTask.h"

class TTree;
class AliESDEvent;
class AliESDtrackCuts;

class AliEsdSkimTask : public AliAnalysisTaskSE {
 public:
  AliEsdSkimTask(const char *opt=0);

  void UserExec(Option_t *opt);
  void UserCreateOutputObjects();
  void SetDoZdc(Bool_t b)          { fDoZDC    = b; }
  void SetDoV0(Bool_t b)           { fDoV0     = b; }
  void SetDoT0(Bool_t b)           { fDoT0     = b; }
  void SetDoTpcV(Bool_t b)         { fDoTPCv   = b; }
  void SetDoSpdV(Bool_t b)         { fDoSPDv   = b; }
  void SetDoPriV(Bool_t b)         { fDoPriv   = b; }
  void SetDoEmC(Bool_t b)          { fDoEmCs   = b; }
  void SetDoPhC(Bool_t b)          { fDoPCs    = b; }
  void SetDoEmT(Bool_t b)          { fDoEmT    = b; }
  void SetDoPhT(Bool_t b)          { fDoPT     = b; }
  void SetDoTracks(Bool_t b)       { fDoTracks = b; }
  void SetDoMult(Bool_t b)         { fDoMult   = b; }
  void SetDoTof(Bool_t b)          { fDoTof    = b; }
  void SetDoPileup(Bool_t b)       { fDoPileup = b; }
  void SetDoClus(Bool_t b)         { fDoClus   = b; }
  void SetCuts(AliESDtrackCuts *c) { fCuts     = c; }

 protected:
  AliESDEvent     *fEvent;    //!esd event
  TTree           *fTree;     //!tree
  AliESDtrackCuts *fCuts;     // track cuts
  Bool_t           fDoZDC;    // do zdc
  Bool_t           fDoV0;     // do vzero
  Bool_t           fDoT0;     // do tzero
  Bool_t           fDoTPCv;   // do tpc vertex
  Bool_t           fDoSPDv;   // do spd vertex
  Bool_t           fDoPriv;   // do primary vertex
  Bool_t           fDoEmCs;   // do emcal cells
  Bool_t           fDoPCs;    // do phos cells
  Bool_t           fDoEmT;    // do emcal trigger
  Bool_t           fDoPT;     // do phos trigger
  Bool_t           fDoTracks; // do tracks
  Bool_t           fDoMult;   // do mult
  Bool_t           fDoTof;    // do TOF
  Bool_t           fDoPileup; // do pileup
  Bool_t           fDoClus;   // do clusters
    
 private:
  AliEsdSkimTask(const AliEsdSkimTask&);            // not implemented
  AliEsdSkimTask &operator=(const AliEsdSkimTask&); // not implemented

 ClassDef(AliEsdSkimTask, 1); // Esd trimming and skimming task
};
#endif
