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
  void SetCheckCond(Int_t c)       { fCheckCond     = c; }
  void SetCuts(AliESDtrackCuts *c) { fCuts          = c; }
  void SetDoCent(Bool_t b)         { fDoCent        = b; }
  void SetDoClus(Bool_t b)         { fDoClus        = b; }
  void SetDoEmC(Bool_t b)          { fDoEmCs        = b; }
  void SetDoEmT(Bool_t b)          { fDoEmT         = b; }
  void SetDoFmd(Bool_t b)          { fDoFmd         = b; }
  void SetDoMiniTracks(Bool_t b)   { fDoMiniTracks  = b; }
  void SetDoMult(Bool_t b)         { fDoMult        = b; }
  void SetDoMuonTracks(Bool_t b)   { fDoMuonTracks  = b; }
  void SetDoPhC(Bool_t b)          { fDoPCs         = b; }
  void SetDoPhT(Bool_t b)          { fDoPT          = b; }
  void SetDoPicoTracks(Bool_t b)   { fDoPicoTracks  = b; }
  void SetDoPileup(Bool_t b)       { fDoPileup      = b; }
  void SetDoPriV(Bool_t b)         { fDoPriv        = b; }
  void SetDoRP(Bool_t b)           { fDoRP          = b; }
  void SetDoSaveBytes(Bool_t b)    { fDoSaveBytes   = b; }
  void SetDoSpdV(Bool_t b)         { fDoSPDv        = b; }
  void SetDoT0(Bool_t b)           { fDoT0          = b; }
  void SetDoTof(Bool_t b)          { fDoTof         = b; }
  void SetDoTpcV(Bool_t b)         { fDoTPCv        = b; }
  void SetDoTracks(Bool_t b)       { fDoTracks      = b; }
  void SetDoV0(Bool_t b)           { fDoV0          = b; }
  void SetDoZdc(Bool_t b)          { fDoZDC         = b; }
  void SetEmcNames(const char *n)  { fEmcNames      = n; }
  void SetEmcalClusOnly(Bool_t b)  { fEmcalClusOnly = b; }
  void SetPhosClusOnly(Bool_t b)   { fPhosClusOnly  = b; }
  void SetRemoveCP(Bool_t b)       { fRemoveCP      = b; }
  void SetResetCov(Bool_t b)       { fResetCov      = b; }
  void SetTracks(const char *n)    { fTracks        = n; }
  void SetDoAllTracks(Bool_t b)    { fDoAllTracks   = b; }
  void SetDoV0s(Bool_t b)          { fDoV0s         = b; }
  void SetDoCascades(Bool_t b)     { fDoCascades    = b; }
  void SetDoKinks(Bool_t b)        { fDoKinks       = b; }
  void SetDoErrorLogs(Bool_t b)    { fDoErrorLogs   = b; }

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
  Bool_t           fDoMuonTracks; // do muon tracks
  Bool_t           fDoV0s;        // do v0s 
  Bool_t           fDoCascades;   // do cascades
  Bool_t           fDoKinks;      // do kinks
  Bool_t           fDoErrorLogs;  // do error logs
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
  Bool_t           fDoAllTracks;  // if true then keep full tracks
  Bool_t           fDoPicoTracks; // if true then do pico tracks
  Int_t            fCheckCond;    // if !=0 check certain conditions before event is accepted

 private:
  AliEsdSkimTask(const AliEsdSkimTask&);            // not implemented
  AliEsdSkimTask &operator=(const AliEsdSkimTask&); // not implemented

 ClassDef(AliEsdSkimTask, 6); // Esd trimming and skimming task
};
#endif
