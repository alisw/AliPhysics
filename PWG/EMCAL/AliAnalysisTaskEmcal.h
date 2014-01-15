#ifndef ALIANALYSISTASKEMCAL_H
#define ALIANALYSISTASKEMCAL_H

// $Id: AliAnalysisTaskEmcalDev.h 64518 2013-10-14 12:44:52Z loizides $

class TClonesArray;
class TString;
class TList;
class AliEmcalParticle;
class AliMCParticle;
class AliVCluster;
class AliVTrack;
class AliVParticle;
class AliVCaloCells;
class TH1;
class TProfile;
class AliEMCALGeometry;
class AliParticleContainer;
class AliClusterContainer;
class AliGenPythiaEventHeader;
class AliVCaloTrigger;
class AliAnalysisUtils;
class AliEmcalTriggerPatchInfo;

#include "Rtypes.h"

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEmcal : public AliAnalysisTaskSE {
 public:

  enum BeamType {
    kNA       = -1,
    kpp       = 0,
    kAA       = 1,
    kpA       = 2
  };

  enum TriggerType {
    kND       = -1,  //not defined
    kJ1       = 0,
    kJ2       = 1
  };

  AliAnalysisTaskEmcal();
  AliAnalysisTaskEmcal(const char *name, Bool_t histo=kFALSE); 
  virtual ~AliAnalysisTaskEmcal();

  void                        UserExec(Option_t *option);
  void                        UserCreateOutputObjects();
  Bool_t                      UserNotify();

  void                        SetNCentBins(Int_t n)                                 { fNcentBins         = n                              ; }  
  void                        SetCentRange(Double_t min, Double_t max)              { fMinCent           = min  ; fMaxCent = max          ; }
  void                        SetHistoBins(Int_t nbins, Double_t min, Double_t max) { fNbins = nbins; fMinBinPt = min; fMaxBinPt = max    ; }
  void                        SetOffTrigger(UInt_t t)                               { fOffTrigger        = t                              ; }
  void                        SetTrigClass(const char *n)                           { fTrigClass         = n                              ; } 
  void                        SetTriggerTypeSel(TriggerType t)                      { fTriggerTypeSel    = t                              ; } 
  void                        SetVzRange(Double_t min, Double_t max)                { fMinVz             = min  ; fMaxVz   = max          ; }
  void                        SetForceBeamType(BeamType f)                          { fForceBeamType     = f                              ; }
  void                        SetMakeGeneralHistograms(Bool_t g)                    { fGeneralHistograms = g                              ; }
  void                        SetMinPtTrackInEmcal(Double_t min)                    { fMinPtTrackInEmcal = min                            ; }
  void                        SetEventPlaneVsEmcal(Double_t ep)                     { fEventPlaneVsEmcal = ep                             ; }
  void                        SetCentralityEstimator(const char *c)                 { fCentEst           = c                              ; }

  void                        SetMinNTrack(Int_t min)                               { fMinNTrack         = min                            ; }
  void                        SetUseAliAnaUtils(Bool_t b)                           { fUseAliAnaUtils    = b                              ; }

  void                        SetMinMCLabel(Int_t s)                                { fMinMCLabel        = s                              ; }
  void                        SetIsEmbedded(Bool_t i)                               { fIsEmbedded        = i                              ; }
  void                        SetIsPythia(Bool_t i)                                 { fIsPythia          = i                              ; }
  void                        SetMCLabelShift(Int_t s)                              { fMCLabelShift      = s                              ; }

  void                        SetCaloCellsName(const char *n)                       { fCaloCellsName     = n                              ; }
  void                        SetCaloTriggersName(const char *n)                    { fCaloTriggersName  = n                              ; }
  void                        SetCaloTriggerPatchInfoName(const char *n)            { fCaloTriggerPatchInfoName = n                       ; }

  void                        SetTracksName(const char *n)                          { AddParticleContainer(n)                             ; }
  void                        SetClusName(const char *n)                            { AddClusterContainer(n)                              ; }

  void                        SetClusPtCut(Double_t cut, Int_t c=0);
  void                        SetClusTimeCut(Double_t min, Double_t max, Int_t c=0);
  void                        SetTrackPtCut(Double_t cut, Int_t c=0);
  void                        SetTrackEtaLimits(Double_t min, Double_t max, Int_t c=0);
  void                        SetTrackPhiLimits(Double_t min, Double_t max, Int_t c=0);

  AliParticleContainer       *AddParticleContainer(const char *n);
  AliClusterContainer        *AddClusterContainer(const char *n);
  void                        RemoveParticleContainer(Int_t i=0)                      { fParticleCollArray.RemoveAt(i);} 
  void                        RemoveClusterContainer(Int_t i=0)                       { fClusterCollArray.RemoveAt(i);} 

  AliParticleContainer       *GetParticleContainer(Int_t i=0)   const;
  AliClusterContainer        *GetClusterContainer(Int_t i=0)    const;
  AliParticleContainer       *GetParticleContainer(const char* name)  const;
  AliClusterContainer        *GetClusterContainer(const char* name)   const;

  AliEmcalTriggerPatchInfo   *GetMainTriggerPatch();

 protected:
  BeamType                    GetBeamType();
  TriggerType                 GetTriggerType();
  TClonesArray               *GetArrayFromEvent(const char *name, const char *clname=0);
  Bool_t                      PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials, Int_t &pthard);

  TClonesArray               *GetParticleArray(Int_t i=0)                  const;
  TClonesArray               *GetClusterArray(Int_t i=0)                   const;

  AliVParticle               *GetAcceptParticleFromArray(Int_t p, Int_t c=0)     const;
  AliVCluster                *GetAcceptClusterFromArray(Int_t cl, Int_t c=0)     const;

  Int_t                       GetNParticles(Int_t i=0)                           const;
  Int_t                       GetNClusters(Int_t i=0)                            const;

  Bool_t                      AcceptCluster(AliVCluster *clus, Int_t c = 0)      const;
  Bool_t                      AcceptTrack(AliVParticle *track, Int_t c = 0)      const;

  // Virtual functions, to be overloaded in derived classes
  virtual void                ExecOnce();
  virtual Bool_t              FillGeneralHistograms();
  virtual Bool_t              IsEventSelected();
  virtual Bool_t              RetrieveEventObjects();
  virtual Bool_t              FillHistograms()                                     { return kTRUE                 ; }
  virtual Bool_t              Run()                                                { return kTRUE                 ; }


  BeamType                    fForceBeamType;              // forced beam type
  Bool_t                      fGeneralHistograms;          // whether or not it should fill some general histograms
  Bool_t                      fInitialized;                // whether or not the task has been already initialized
  Bool_t                      fCreateHisto;                // whether or not create histograms
  TString                     fCaloCellsName;              // name of calo cell collection
  TString                     fCaloTriggersName;           // name of calo triggers collection
  TString                     fCaloTriggerPatchInfoName;   // trigger patch info array name
  Double_t                    fMinCent;                    // min centrality for event selection
  Double_t                    fMaxCent;                    // max centrality for event selection
  Double_t                    fMinVz;                      // min vertex for event selection
  Double_t                    fMaxVz;                      // max vertex for event selection
  Double_t                    fTrackPtCut;                 // cut on track pt in event selection
  Int_t                       fMinNTrack;                  // minimum nr of tracks in event with pT>fTrackPtCut
  Bool_t                      fUseAliAnaUtils;             //  used for LHC13* data
  AliAnalysisUtils           *fAliAnalysisUtils;           //! vertex selection (optional)

  UInt_t                      fOffTrigger;                 // offline trigger for event selection
  TString                     fTrigClass;                  // trigger class name for event selection
  TriggerType                 fTriggerTypeSel;             // trigger type to select based on trigger patches
  Int_t                       fNbins;                      // no. of pt bins
  Double_t                    fMinBinPt;                   // min pt in histograms
  Double_t                    fMaxBinPt;                   // max pt in histograms
  Double_t                    fMinPtTrackInEmcal;          // min pt track in emcal
  Double_t                    fEventPlaneVsEmcal;          // select events which have a certain event plane wrt the emcal
  Double_t                    fMinEventPlane;              // minimum event plane value
  Double_t                    fMaxEventPlane;              // maximum event plane value
  TString                     fCentEst;                    // name of V0 centrality estimator
  Bool_t                      fIsEmbedded;                 // trigger, embedded signal
  Bool_t                      fIsPythia;                   // trigger, if it is a PYTHIA production
  Int_t                       fSelectPtHardBin;            // select one pt hard bin for analysis
  Int_t                       fMinMCLabel;                 // minimum MC label value for the tracks/clusters being considered MC particles
  Int_t                       fMCLabelShift;               // if MC label > fMCLabelShift, MC label -= fMCLabelShift
  Int_t                       fNcentBins;                  // how many centrality bins
  AliEMCALGeometry           *fGeom;                       //!emcal geometry
  TClonesArray               *fTracks;                     //!tracks
  TClonesArray               *fCaloClusters;               //!clusters
  AliVCaloCells              *fCaloCells;                  //!cells
  AliVCaloTrigger            *fCaloTriggers;               //!calo triggers
  TClonesArray               *fTriggerPatchInfo;           //!trigger patch info array
  Double_t                    fCent;                       //!event centrality
  Int_t                       fCentBin;                    //!event centrality bin
  Double_t                    fEPV0;                       //!event plane V0
  Double_t                    fEPV0A;                      //!event plane V0A
  Double_t                    fEPV0C;                      //!event plane V0C
  Double_t                    fVertex[3];                  //!event vertex
  Int_t                       fNVertCont;                  //!event vertex number of contributors
  BeamType                    fBeamType;                   //!event beam type

  // PYTHIA
  AliGenPythiaEventHeader    *fPythiaHeader;               //!event Pythia header
  Double_t                    fPtHard;                     //!event pt hard
  Int_t                       fPtHardBin;                  //!event pt hard bin
  Int_t                       fNTrials;                    //!event trials

  TObjArray                   fParticleCollArray;          // particle/track collection array
  TObjArray                   fClusterCollArray;           // cluster collection array
  AliEmcalTriggerPatchInfo   *fMainTriggerPatch;           // main trigger patch, will be cached after calling GetMainTriggerPatch() first time
  TriggerType                 fTriggerType;                // trigger type J1 or J2

  // Histograms
  TList                      *fOutput;                     //!output list

  // PYTHIA
  TH1                        *fHistTrialsAfterSel;         //!total number of trials per pt hard bin after selection
  TH1                        *fHistEventsAfterSel;         //!total number of events per pt hard bin after selection
  TH1                        *fHistTrials;                 //!trials from pyxsec.root
  TProfile                   *fHistXsection;               //!x section from pyxsec.root
  TH1                        *fHistEvents;                 //!total number of events per pt hard bin
  TH1                        *fHistPtHard;                 //!pt hard distribution

  // General histograms
  TH1                        *fHistCentrality;             //!event centrality distribution
  TH1                        *fHistZVertex;                //!z vertex position
  TH1                        *fHistEventPlane;             //!event plane distribution
  TH1                        *fHistEventRejection;         //!book keep reasons for rejecting event

 private:
  AliAnalysisTaskEmcal(const AliAnalysisTaskEmcal&);            // not implemented
  AliAnalysisTaskEmcal &operator=(const AliAnalysisTaskEmcal&); // not implemented

  ClassDef(AliAnalysisTaskEmcal, 6) // EMCAL base analysis task
};
#endif
