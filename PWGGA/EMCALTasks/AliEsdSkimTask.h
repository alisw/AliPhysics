#ifndef ALIESDSKINTASK_H
#define ALIESDSKINTASK_H

// $Id$

#include "AliPhysicsSelectionTask.h"
#include "AliESDtrack.h"

class TTree;
class AliESDEvent;
class AliESDtrackCuts;

class AliEsdSkimTask : public AliAnalysisTaskSE {
 public:
  AliEsdSkimTask(const char *opt=0);

  void UserExec(Option_t *opt);
  void UserCreateOutputObjects();
  void SetCuts(AliESDtrackCuts *c) { fCuts          = c; }
  void SetDoClus(Bool_t b)         { fDoClus        = b; }
  void SetDoEmC(Bool_t b)          { fDoEmCs        = b; }
  void SetDoEmT(Bool_t b)          { fDoEmT         = b; }
  void SetDoFmd(Bool_t b)          { fDoFmd         = b; }
  void SetDoMiniTracks(Bool_t b)   { fDoMiniTracks  = b; }
  void SetDoMult(Bool_t b)         { fDoMult        = b; }
  void SetDoPhC(Bool_t b)          { fDoPCs         = b; }
  void SetDoPhT(Bool_t b)          { fDoPT          = b; }
  void SetDoPicoTracks(Bool_t b)   { fDoPicoTracks  = b; }
  void SetDoPileup(Bool_t b)       { fDoPileup      = b; }
  void SetDoPriV(Bool_t b)         { fDoPriv        = b; }
  void SetDoSaveBytes(Bool_t b)    { fDoSaveBytes   = b; }
  void SetDoSpdV(Bool_t b)         { fDoSPDv        = b; }
  void SetDoT0(Bool_t b)           { fDoT0          = b; }
  void SetDoTof(Bool_t b)          { fDoTof         = b; }
  void SetDoTpcV(Bool_t b)         { fDoTPCv        = b; }
  void SetDoTracks(Bool_t b)       { fDoTracks      = b; }
  void SetDoV0(Bool_t b)           { fDoV0          = b; }
  void SetDoZdc(Bool_t b)          { fDoZDC         = b; }
  void SetEmcNames(const char *n)  { fEmcNames      = n; }
  void SetPhosClusOnly(Bool_t b)   { fPhosClusOnly  = b; }
  void SetEmcalClusOnly(Bool_t b)  { fEmcalClusOnly = b; }
  void SetRemoveCP(Bool_t b)       { fRemoveCP      = b; }
  void SetResetCov(Bool_t b)       { fResetCov      = b; }
  void SetTracks(const char *n)    { fTracks        = n; }

 protected:
  AliESDEvent     *fEvent;        //!esd event
  TTree           *fTree;         //!tree
  AliESDtrackCuts *fCuts;         // track cuts
  Bool_t           fDoZDC;        // do zdc
  Bool_t           fDoV0;         // do vzero
  Bool_t           fDoT0;         // do tzero
  Bool_t           fDoTPCv;       // do tpc vertex
  Bool_t           fDoSPDv;       // do spd vertex
  Bool_t           fDoPriv;       // do primary vertex
  Bool_t           fDoEmCs;       // do emcal cells
  Bool_t           fDoPCs;        // do phos cells
  Bool_t           fDoEmT;        // do emcal trigger
  Bool_t           fDoPT;         // do phos trigger
  Bool_t           fDoTracks;     // do tracks
  Bool_t           fDoFmd;        // do fmd
  Bool_t           fDoMult;       // do mult
  Bool_t           fDoTof;        // do TOF
  Bool_t           fDoPileup;     // do pileup
  Bool_t           fDoClus;       // do clusters
  TString          fEmcNames;     // name of clusters
  Bool_t           fDoMiniTracks; // strip down tracks
  TString          fTracks;       // name of tracks (e.g. tracks propagated to EMCAL surface)
  Bool_t           fPhosClusOnly; // if true then only store PHOS clusters
  Bool_t           fEmcalClusOnly;// if true then only store EMCAL clusters
  Bool_t           fDoSaveBytes;  // if true then trim down some of the stored objects (mult, fmd)
  Bool_t           fDoCent;       // do centrality
  Bool_t           fDoRP;         // do reaction plane
  Bool_t           fRemoveCP;     // if false then keep constrained parameters (only reset covariance)
  Bool_t           fResetCov;     // if true reset covariance matrix of track
  Bool_t           fDoPicoTracks; // if true then do pico tracks

 private:
  AliEsdSkimTask(const AliEsdSkimTask&);            // not implemented
  AliEsdSkimTask &operator=(const AliEsdSkimTask&); // not implemented

 ClassDef(AliEsdSkimTask, 4); // Esd trimming and skimming task
};
#endif
