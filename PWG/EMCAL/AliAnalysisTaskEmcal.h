#ifndef ALIANALYSISTASKEMCAL_H
#define ALIANALYSISTASKEMCAL_H

// $Id: AliAnalysisTaskEmcal.h 64518 2013-10-14 12:44:52Z loizides $

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
    kND       = BIT(31),  //not defined := Sign bit of int
    kJ1       = BIT(0),
    kJ2       = BIT(1),
    kG1		    = BIT(2),
    kG2 	    = BIT(3),
    kL0       = BIT(4)
  };
  enum MainTriggers{
    kPosJ1    = 0,
    kPosJ2    = 1,
    kPosG1		= 2,
    kPosG2 	  = 3,
    kPosL0    = 4,
    kNType    = 5
  };

  AliAnalysisTaskEmcal();
  AliAnalysisTaskEmcal(const char *name, Bool_t histo=kFALSE); 
  virtual ~AliAnalysisTaskEmcal();

  AliParticleContainer       *AddParticleContainer(const char *n);
  AliClusterContainer        *AddClusterContainer(const char *n);
  AliParticleContainer       *GetParticleContainer(Int_t i=0)         const;
  AliClusterContainer        *GetClusterContainer(Int_t i=0)          const;
  AliParticleContainer       *GetParticleContainer(const char* name)  const;
  AliClusterContainer        *GetClusterContainer(const char* name)   const;
  void                        RemoveParticleContainer(Int_t i=0)                    { fParticleCollArray.RemoveAt(i)                      ; } 
  void                        RemoveClusterContainer(Int_t i=0)                     { fClusterCollArray.RemoveAt(i)                       ; } 
  void                        SetCaloCellsName(const char *n)                       { fCaloCellsName     = n                              ; }
  void                        SetCaloTriggerPatchInfoName(const char *n)            { fCaloTriggerPatchInfoName = n                       ; }
  void                        SetCaloTriggersName(const char *n)                    { fCaloTriggersName  = n                              ; }
  void                        SetCentRange(Double_t min, Double_t max)              { fMinCent           = min  ; fMaxCent = max          ; }
  void                        SetCentralityEstimator(const char *c)                 { fCentEst           = c                              ; }
  void                        SetClusName(const char *n)                            { AddClusterContainer(n)                              ; }
  void                        SetClusPtCut(Double_t cut, Int_t c=0);
  void                        SetClusTimeCut(Double_t min, Double_t max, Int_t c=0);
  void                        SetEventPlaneVsEmcal(Double_t ep)                     { fEventPlaneVsEmcal = ep                             ; }
  void                        SetForceBeamType(BeamType f)                          { fForceBeamType     = f                              ; }
  void                        SetHistoBins(Int_t nbins, Double_t min, Double_t max) { fNbins = nbins; fMinBinPt = min; fMaxBinPt = max    ; }
  void                        SetIsEmbedded(Bool_t i)                               { fIsEmbedded        = i                              ; }
  void                        SetIsPythia(Bool_t i)                                 { fIsPythia          = i                              ; }
  void                        SetMCLabelShift(Int_t s)                              { fMCLabelShift      = s                              ; }
  void                        SetMinMCLabel(Int_t s)                                { fMinMCLabel        = s                              ; }
  void                        SetMinNTrack(Int_t min)                               { fMinNTrack         = min                            ; }
  void                        SetMinPtTrackInEmcal(Double_t min)                    { fMinPtTrackInEmcal = min                            ; }
  void                        SetNCentBins(Int_t n)                                 { fNcentBins         = n                              ; } 
  void                        SetNeedEmcalGeom(Bool_t n)                            { fNeedEmcalGeom     = n                              ; }
  void                        SetOffTrigger(UInt_t t)                               { fOffTrigger        = t                              ; }
  void                        SetTrackEtaLimits(Double_t min, Double_t max, Int_t c=0);
  void                        SetTrackPhiLimits(Double_t min, Double_t max, Int_t c=0);
  void                        SetTrackPtCut(Double_t cut, Int_t c=0);
  void                        SetTracksName(const char *n)                          { AddParticleContainer(n)                             ; }
  void                        SetTrigClass(const char *n)                           { fTrigClass         = n                              ; } 
  void                        SetTriggerTypeSel(TriggerType t)                      { fTriggerTypeSel    = t                              ; } 
  void                        SetUseAliAnaUtils(Bool_t b, Bool_t bRejPilup = kTRUE) { fUseAliAnaUtils    = b ; fRejectPileup = bRejPilup  ; }
  void                        SetVzRange(Double_t min, Double_t max)                { fMinVz             = min  ; fMaxVz   = max          ; }
  void                        SetUseSPDTrackletVsClusterBG(Bool_t b)                { fTklVsClusSPDCut   = b                              ; }

 protected:
  void                        SetRejectionReasonLabels(TAxis* axis);
  Double_t*                   GenerateFixedBinArray(Int_t n, Double_t min, Double_t max) const;
  void                        GenerateFixedBinArray(Int_t n, Double_t min, Double_t max, Double_t* array) const;
  void                        SetMakeGeneralHistograms(Bool_t g)                    { fGeneralHistograms = g                              ; }
  Bool_t                      AcceptCluster(AliVCluster *clus, Int_t c = 0)      const;
  Bool_t                      AcceptTrack(AliVParticle *track, Int_t c = 0)      const;
  void                        AddObjectToEvent(TObject *obj);
  AliVParticle               *GetAcceptParticleFromArray(Int_t p, Int_t c=0)     const;
  AliVCluster                *GetAcceptClusterFromArray(Int_t cl, Int_t c=0)     const;
  TClonesArray               *GetArrayFromEvent(const char *name, const char *clname=0);
  BeamType                    GetBeamType();
  TClonesArray               *GetParticleArray(Int_t i=0)                        const;
  TClonesArray               *GetClusterArray(Int_t i=0)                         const;
  Int_t                       GetNParticles(Int_t i=0)                           const;
  Int_t                       GetNClusters(Int_t i=0)                            const;
  AliEmcalTriggerPatchInfo   *GetMainTriggerPatch(Int_t triggersel = kJ1 | kJ2);
  AliEmcalTriggerPatchInfo   *ApplyMainTriggerSelection(Int_t selectionBitmap) const;
  Bool_t					            HasTriggerType(Int_t triggersel);
  ULong_t 					          GetTriggerList();
  Bool_t                      PythiaInfoFromFile(const char* currFile, Float_t &fXsec, Float_t &fTrials, Int_t &pthard);
  void                        UserCreateOutputObjects();
  void                        UserExec(Option_t *option);
  Bool_t                      UserNotify();

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
  Bool_t                      fMainTriggerPatchSet;        // Internal variable indicating whether the array of main trigger patches is initialised
  TString                     fCaloCellsName;              // name of calo cell collection
  TString                     fCaloTriggersName;           // name of calo triggers collection
  TString                     fCaloTriggerPatchInfoName;   // trigger patch info array name
  Double_t                    fMinCent;                    // min centrality for event selection
  Double_t                    fMaxCent;                    // max centrality for event selection
  Double_t                    fMinVz;                      // min vertex for event selection
  Double_t                    fMaxVz;                      // max vertex for event selection
  Double_t                    fTrackPtCut;                 // cut on track pt in event selection
  Int_t                       fMinNTrack;                  // minimum nr of tracks in event with pT>fTrackPtCut
  Bool_t                      fUseAliAnaUtils;             // used for LHC13* data: z-vtx, Ncontributors, z-vtx resolution cuts
  Bool_t                      fRejectPileup;               // Reject pilup using function AliAnalysisUtils::IsPileUpEvent()
  Bool_t                      fTklVsClusSPDCut;            // Apply tracklet-vs-cluster SPD cut to reject background events in pp
  AliAnalysisUtils           *fAliAnalysisUtils;           //! vertex selection (optional)
  UInt_t                      fOffTrigger;                 // offline trigger for event selection
  TString                     fTrigClass;                  // trigger class name for event selection
  Int_t                       fTriggerTypeSel;             // trigger type to select based on trigger patches
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
  Bool_t                      fNeedEmcalGeom;              // whether or not the task needs the emcal geometry
  Bool_t                      fIsEsd;                      //!whether it's an ESD analysis
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
  AliGenPythiaEventHeader    *fPythiaHeader;               //!event Pythia header
  Double_t                    fPtHard;                     //!event pt hard
  Int_t                       fPtHardBin;                  //!event pt hard bin
  Int_t                       fNTrials;                    //!event trials
  Float_t                     fXsection;                   //!x-section from pythia header
  TObjArray                   fParticleCollArray;          // particle/track collection array
  TObjArray                   fClusterCollArray;           // cluster collection array
  AliEmcalTriggerPatchInfo   *fMainTriggerPatch[kNType];   // main trigger patch, will be cached after calling GetMainTriggerPatch() first time
  ULong_t                     fTriggers;                   // list of fired triggers

  TList                      *fOutput;                     //!output list
  TH1                        *fHistEventCount;             //!incoming and selected events
  TH1                        *fHistTrialsAfterSel;         //!total number of trials per pt hard bin after selection
  TH1                        *fHistEventsAfterSel;         //!total number of events per pt hard bin after selection
  TProfile                   *fHistXsectionAfterSel;       //!x section from pythia header
  TH1                        *fHistTrials;                 //!trials from pyxsec.root
  TH1                        *fHistEvents;                 //!total number of events per pt hard bin
  TProfile                   *fHistXsection;               //!x section from pyxsec.root
  TH1                        *fHistPtHard;                 //!pt hard distribution
  TH1                        *fHistCentrality;             //!event centrality distribution
  TH1                        *fHistZVertex;                //!z vertex position
  TH1                        *fHistEventPlane;             //!event plane distribution
  TH1                        *fHistEventRejection;         //!book keep reasons for rejecting event

 private:
  AliAnalysisTaskEmcal(const AliAnalysisTaskEmcal&);            // not implemented
  AliAnalysisTaskEmcal &operator=(const AliAnalysisTaskEmcal&); // not implemented

  ClassDef(AliAnalysisTaskEmcal, 11) // EMCAL base analysis task
};
#endif
