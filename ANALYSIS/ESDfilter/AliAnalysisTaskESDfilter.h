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
  virtual Bool_t Notify();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *option);
  virtual void   ConvertESDtoAOD();

  // Setters
  virtual void SetTrackFilter   (AliAnalysisFilter*   trackF)                {fTrackFilter                 = trackF;}
  virtual void SetTPCOnlyFilterMask (UInt_t filterMask)                      {SetTPCConstrainedFilterMask(filterMask);}
  virtual void SetTPCConstrainedFilterMask (UInt_t filterMask)               {fTPCConstrainedFilterMask    = filterMask;}
  virtual void SetHybridFilterMaskTPCConstrainedGlobal(UInt_t filterMask)    {fHybridFilterMaskTPCCG       = filterMask;}
  virtual void SetWriteHybridTPCConstrainedOnly(bool b)                      {fWriteHybridTPCCOnly         = b;}
  virtual void SetGlobalConstrainedFilterMask (UInt_t filterMask)            {fGlobalConstrainedFilterMask = filterMask;}
  virtual void SetHybridFilterMaskGlobalConstrainedGlobal(UInt_t filterMask) {fHybridFilterMaskGCG         = filterMask;}
  virtual void SetWriteHybridGlobalConstrainedOnly(bool b)    {fWriteHybridGCOnly = b;}
  virtual void SetKinkFilter    (AliAnalysisFilter*    KinkF) {fKinkFilter        = KinkF;}
  virtual void SetV0Filter      (AliAnalysisFilter*      V0F) {fV0Filter          = V0F;}
  virtual void SetCascadeFilter (AliAnalysisFilter* CascadeF) {fCascadeFilter     = CascadeF;}
  virtual void SetPthreshold    (Double_t p)                  {fHighPthreshold    = p;}
  virtual void SetPshape        (TF1 *func)                   {fPtshape           = func;}
  virtual void SetEnableFillAOD (Bool_t b)                    {fEnableFillAOD     = b;}
  virtual void SetAODPID(AliESDtrack *esdtrack, AliAODTrack *aodtrack, AliAODPid *detpid);
  void SetDetectorRawSignals(AliAODPid *aodpid, AliESDtrack *track);
  void SetV0Cuts(const Double_t cuts[7]) {for (Int_t icut = 0; icut<7; icut++) fV0Cuts[icut] = cuts[icut];}
  void SetCascadeCuts(const Double_t cuts[8]) {for (Int_t icut = 0; icut<8; icut++) fCascadeCuts[icut] = cuts[icut];} 
  void GetV0Cuts(Double_t cuts[7]) const {for (Int_t icut = 0; icut<7; icut++) cuts[icut] = fV0Cuts[icut];}
  void GetCascadeCuts(Double_t cuts[8]) const {for (Int_t icut = 0; icut<8; icut++) cuts[icut] = fCascadeCuts[icut];}
  Bool_t AddMetadataToUserInfo();
  void PrintTask(Option_t *option="all", Int_t indent=0) const;
  void DisableVZERO()        {fIsVZEROEnabled = kFALSE;}
  void DisableTZERO()        {fIsTZEROEnabled = kFALSE;}
  void DisableZDC()          {fIsZDCEnabled   = kFALSE;}
  void DisableAD()           {fIsADEnabled   = kFALSE;}
  void DisableCascades()     {fAreCascadesEnabled = kFALSE;}
  void DisableV0s()          {fAreV0sEnabled = kFALSE;}
  void DisableKinks()        {fAreKinksEnabled = kFALSE;}
  void DisableTracks()       {fAreTracksEnabled = kFALSE;}
  void DisablePmdClusters()  {fArePmdClustersEnabled = kFALSE;}
  void DisableCaloClusters() {fAreCaloClustersEnabled = kFALSE;}
  void DisableCells()        {fAreEMCALCellsEnabled = fArePHOSCellsEnabled = kFALSE; }
  void DisableCaloTrigger(TString calo = "PHOS") {if (calo.Contains("EMCAL")) fAreEMCALTriggerEnabled = kFALSE; else fArePHOSTriggerEnabled = kFALSE;}
  void DisableTracklets()    {fAreTrackletsEnabled = kFALSE;}
  void DisableHMPID()        {fIsHMPIDEnabled = kFALSE;} 
  void EnableV0CascadeVerticesReco() {fIsV0CascadeRecoEnabled = kTRUE;}
  void SetPropagateTrackToEMCal(Bool_t propagate) {fDoPropagateTrackToEMCal = propagate;}
  void SetEMCalSurfaceDistance(Double_t d)        {fEMCalSurfaceDistance = d;}
  void SetRefitVertexTracks(Int_t algo=6, Double_t* cuts=0);
  
  void SetMuonCaloPass();
  
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
  void ConvertCaloTrigger(TString calo, const AliESDEvent& esd);
  void ConvertTracklets(const AliESDEvent& esd);
  void ConvertTPCOnlyTracks(const AliESDEvent& esd);
  void ConvertGlobalConstrainedTracks(const AliESDEvent& esd);
  void ConvertVZERO(const AliESDEvent& esd);
  void ConvertTZERO(const AliESDEvent& esd);
  void ConvertZDC(const AliESDEvent& esd);
  void ConvertAD(const AliESDEvent& esd);
  Int_t ConvertHMPID(const AliESDEvent& esd);
  void ConvertTRD(const AliESDEvent& esd);
  void CopyCaloProps(AliESDtrack *esdt, AliAODTrack *aodt);

  TClonesArray& Tracks();
  TClonesArray& V0s();
  TClonesArray& Vertices();
  TClonesArray& Cascades();
  
  // Filtering
  AliAnalysisFilter* fTrackFilter;                 // Track   Filter
  AliAnalysisFilter* fKinkFilter;                  // Kink    Filter
  AliAnalysisFilter* fV0Filter;                    // V0      Filter
  AliAnalysisFilter* fCascadeFilter;               // Cascade Filter
  Double_t           fHighPthreshold;              // Pt threshold for detector signal setting
  TF1 *              fPtshape;                     // Pt spectrum distribution
  Bool_t             fEnableFillAOD;               // value that decides if this task activates AOD filling
  Bool_t*            fUsedTrack;                   //! indices of used tracks
  UInt_t*            fUsedTrackCopy;               //! filterbits of tracks for which a copy was added to the AODs
  Bool_t*            fUsedKink;                    //! indices of used kinks
  Bool_t*            fUsedV0;                      //! indices of used V0s
  TRefArray*         fAODTrackRefs;                // array of track references
  TRefArray*         fAODV0VtxRefs;                // array of v0 vertices references
  TRefArray*         fAODV0Refs;                   // array of v0s references
  AliMCEventHandler* fMChandler;                   // pointer to MC handler (if any)
  Int_t              fNumberOfTracks;              // current number of tracks
  Int_t              fNumberOfPositiveTracks;      // current number of positive tracks
  Int_t              fNumberOfV0s;                 // current number of v0s
  Int_t              fNumberOfVertices;            // current number of vertices
  Int_t              fNumberOfCascades;            // current number of cascades
  Int_t              fNumberOfKinks;               // current number of kinks
  Bool_t             fOldESDformat;                // is the ESD in old format ?
  AliAODVertex*      fPrimaryVertex;               // pointer to primary vertex of the event
  UInt_t             fTPCConstrainedFilterMask;    // Filter Mask used to select and store refitted TPC only tracks
  UInt_t             fHybridFilterMaskTPCCG;       // Filter Mask used to mark global tracks as hybrid
  Bool_t             fWriteHybridTPCCOnly;         // write only the complent tracks not all global constrained
  UInt_t             fGlobalConstrainedFilterMask; // Filter Mask used to select and store refitted TPC only tracks
  UInt_t             fHybridFilterMaskGCG;         // Filter Mask used to mark global tracks as hybrid
  Bool_t             fWriteHybridGCOnly;           // write only the complent tracks not all global constrained
  Bool_t             fIsVZEROEnabled;              // whether or not to fill the vzero branch (true by default)
  Bool_t             fIsTZEROEnabled;              // whether or not to fill the tzero branch (true by default)
  Bool_t             fIsZDCEnabled;                // whether or not to fill the zdc branch (true by default)
  Bool_t             fIsADEnabled;                 // whether or not to fill the ad branch (true by default) 
  Bool_t             fIsHMPIDEnabled;              // whether or not to fill the hmpid branch (true by default) 
  Bool_t             fIsV0CascadeRecoEnabled;      // whether or not to reconstruct again V0s and cascades (false by default)
  Bool_t             fAreCascadesEnabled;          // whether or not to fill the cascades branch (true by default)
  Bool_t             fAreV0sEnabled;               // whether or not to fill the v0 branch (true by default)
  Bool_t             fAreKinksEnabled;             // whether or not to fill the kinks (true by default)
  Bool_t             fAreTracksEnabled;            // whether or not to fill the (central) tracks branch (true by default)
  Bool_t             fArePmdClustersEnabled;       // whether or not to fill the pmd clusters (true by default)
  Bool_t             fAreCaloClustersEnabled;      // whether or not to fill the calo clusters (true by default)
  Bool_t             fAreEMCALCellsEnabled;        // whether or not to fill the emcal cells (true by default)
  Bool_t             fArePHOSCellsEnabled;         // whether or not to fill the phos cells (true by default)
  Bool_t             fAreEMCALTriggerEnabled;      // whether or not to fill the emcal trigger (true by default)
  Bool_t             fArePHOSTriggerEnabled;       // whether or not to fill the phos trigger (true by default)
  Bool_t             fAreTrackletsEnabled;         // whether or not to fill the tracklets (true by default)
  Bool_t             fIsTRDEnabled;                // whether or not to fill on-line tracklets and tracks from TRD (true by default)
  AliESDpid*         fESDpid;                      // esd pid
  Bool_t             fIsPidOwner;                  // whether we own fESDpid
  AliESDtrackCuts*   fTPCaloneTrackCuts;           // TPC stand-alone track cuts
  Double_t           fV0Cuts[7];                   // Array to store the values for the different reco selections V0 related
  Double_t           fCascadeCuts[8];              // Array to store the values for the different reco selections cascades related
  Bool_t             fDoPropagateTrackToEMCal;     // whether or not to propagate the tracks to the EMCal surface -- true by default
  Double_t           fEMCalSurfaceDistance;        // EMCal surface distance from the center of the detector (r = 440 by default)
  Int_t              fRefitVertexTracks;           // request to refit the vertex if >=0 (algoID if cuts not supplied, otherwise ncuts)
  Int_t              fRefitVertexTracksNCuts;      // number of cut parameters
  Double_t*          fRefitVertexTracksCuts;       //[fRefitVertexTracksNCuts] optional cuts for vertex refit
  Bool_t fIsMuonCaloPass; /// whether or not this filtering is used on a muon_calo ESD
  
  ClassDef(AliAnalysisTaskESDfilter, 21); // Analysis task for standard ESD filtering
};

#endif
