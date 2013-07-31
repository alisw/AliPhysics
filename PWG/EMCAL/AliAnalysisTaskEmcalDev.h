#ifndef ALIANALYSISTASKEMCALDEV_H
#define ALIANALYSISTASKEMCALDEV_H

// $Id$

class TClonesArray;
class TString;
class TList;
class AliEmcalParticle;
class AliMCParticle;
class AliVCluster;
class AliVTrack;
class AliVParticle;
class AliVCaloCells;
class TH1F;
class AliEMCALGeometry;
class AliParticleContainer;
class AliClusterContainer;

#include "Rtypes.h"

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEmcalDev : public AliAnalysisTaskSE {
 public:
  
  enum EmcalAnaType {
    kTPC       = 0,     // TPC acceptance
    kEMCAL     = 1,     // EMCal acceptance
    kUser      = 2,     // User defined acceptance
  };

  enum BeamType {
    kNA       = -1,
    kpp       = 0,
    kAA       = 1,
    kpA       = 2
  };

  AliAnalysisTaskEmcalDev();
  AliAnalysisTaskEmcalDev(const char *name, Bool_t histo=kFALSE); 
  virtual ~AliAnalysisTaskEmcalDev();

  void                        UserExec(Option_t *option);
  void                        UserCreateOutputObjects();

  void                        SetAnaType(EmcalAnaType type)                         { fAnaType           = type                           ; }
  void                        SetNCentBins(Int_t n)                                 { fNcentBins         = n                              ; }  
  void                        SetCentRange(Double_t min, Double_t max)              { fMinCent           = min  ; fMaxCent = max          ; }
  void                        SetClusName(const char *n)                            { fCaloName          = n    ; AddClusterContainer(n);   }
  void                        SetCaloCellsName(const char *n)                       { fCaloCellsName     = n                              ; }
  void                        SetClusPtCut(Double_t cut)                            { fClusPtCut         = cut                            ; }
  void                        SetClusTimeCut(Double_t min, Double_t max)            { fClusTimeCutLow    = min  ; fClusTimeCutUp = max    ; }
  void                        SetHistoBins(Int_t nbins, Double_t min, Double_t max) { fNbins = nbins; fMinBinPt = min; fMaxBinPt = max    ; }
  void                        SetOffTrigger(UInt_t t)                               { fOffTrigger        = t                              ; }
  void                        SetPtCut(Double_t cut)                                { SetClusPtCut(cut)         ; SetTrackPtCut(cut)      ; }
  void                        SetTrackPtCut(Double_t cut)                           { fTrackPtCut        = cut                            ; }
  void                        SetTrackEtaLimits(Double_t min, Double_t max)         { fTrackMaxEta       = max  ; fTrackMinEta      = min ; }
  void                        SetTrackPhiLimits(Double_t min, Double_t max)         { fTrackMaxPhi       = max  ; fTrackMinPhi      = min ; }
  void                        SetTracksName(const char *n)                          { fTracksName        = n    ; AddParticleContainer(n);  }
  void                        SetTrigClass(const char *n)                           { fTrigClass         = n                              ; }  
  void                        SetVzRange(Double_t min, Double_t max)                { fMinVz             = min  ; fMaxVz   = max          ; }
  void                        SetForceBeamType(BeamType f)                          { fForceBeamType     = f                              ; }
  void                        SetMakeGeneralHistograms(Bool_t g)                    { fGeneralHistograms = g                              ; }
  void                        SetMinPtTrackInEmcal(Double_t min)                    { fMinPtTrackInEmcal = min                            ; }
  void                        SetEventPlaneVsEmcal(Double_t ep)                     { fEventPlaneVsEmcal = ep                             ; }
  void                        SetCentralityEstimator(const char *c)                 { fCentEst           = c                              ; }
  void                        SetTrackBitMap(UInt_t m)                              { fTrackBitMap       = m                              ; }
  void                        SetClusterBitMap(UInt_t m)                            { fClusterBitMap     = m                              ; }
  void                        SetParticleBitMap(UInt_t m)                           { fClusterBitMap     = m    ; fTrackBitMap       = m  ; }
  void                        SetMCTrackBitMap(UInt_t m)                            { fMCTrackBitMap     = m                              ; }
  void                        SetMCClusterBitMap(UInt_t m)                          { fMCClusterBitMap   = m                              ; }
  void                        SetMCParticleBitMap(UInt_t m)                         { fMCClusterBitMap   = m    ; fMCTrackBitMap     = m  ; }
  void                        SetMinMCLabel(Int_t s)                                { fMinMCLabel        = s                              ; }

  void                        AddParticleContainer(const char *n);
  void                        AddClusterContainer(const char *n);


  AliParticleContainer       *GetParticleContainer(const Int_t i=0)   const;
  AliClusterContainer        *GetClusterContainer(const Int_t i=0)    const;
  void                        RemoveParticleContainer(Int_t i=0)                      { fParticleCollArray.RemoveAt(i);} 
  void                        RemoveClusterContainer(Int_t i=0)                       { fClusterCollArray.RemoveAt(i);} 


  TClonesArray               *GetParticleArray(const Int_t i=0)       const;
  TClonesArray               *GetClusterArray(const Int_t i=0)        const;

  AliVParticle               *GetAcceptParticleFromArray(Int_t p, Int_t c=0) const;
  AliVCluster                *GetAcceptClusterFromArray(Int_t cl, Int_t c=0) const;

  Int_t                       GetNParticles(Int_t i=0) const;
  Int_t                       GetNClusters(Int_t i=0)  const;

 protected:
  Bool_t                      AcceptCluster(AliVCluster        *clus,  Int_t c = 0)  const;
  Bool_t                      AcceptEmcalPart(AliEmcalParticle *part)                const;
  Bool_t                      AcceptTrack(AliVParticle         *track, Int_t c = 0)  const;
  virtual void                ExecOnce();
  virtual Bool_t              FillGeneralHistograms();
  virtual Bool_t              FillHistograms()                                     { return kTRUE                 ; }
  BeamType                    GetBeamType();
  TClonesArray               *GetArrayFromEvent(const char *name, const char *clname=0);
  virtual Bool_t              IsEventSelected();
  virtual Bool_t              RetrieveEventObjects();
  virtual Bool_t              Run()                                                { return kTRUE                 ; }

  EmcalAnaType                fAnaType;                    // analysis type
  BeamType                    fForceBeamType;              // forced beam type
  Bool_t                      fGeneralHistograms;          // whether or not it should fill some general histograms
  Bool_t                      fInitialized;                // whether or not the task has been already initialized
  Bool_t                      fCreateHisto;                // whether or not create histograms
  TString                     fTracksName;                 // name of track collection
  TString                     fCaloName;                   // name of calo cluster collection
  TString                     fCaloCellsName;              // name of calo cell collection
  Double_t                    fMinCent;                    // min centrality for event selection
  Double_t                    fMaxCent;                    // max centrality for event selection
  Double_t                    fMinVz;                      // min vertex for event selection
  Double_t                    fMaxVz;                      // max vertex for event selection
  UInt_t                      fOffTrigger;                 // offline trigger for event selection
  TString                     fTrigClass;                  // trigger class name for event selection
  Int_t                       fNbins;                      // no. of pt bins
  Double_t                    fMinBinPt;                   // min pt in histograms
  Double_t                    fMaxBinPt;                   // max pt in histograms
  Double_t                    fClusPtCut;                  // cut on cluster pt
  Double_t                    fTrackPtCut;                 // cut on track pt
  Double_t                    fTrackMinEta;                // cut on track eta
  Double_t                    fTrackMaxEta;                // cut on track eta
  Double_t                    fTrackMinPhi;                // cut on track phi
  Double_t                    fTrackMaxPhi;                // cut on track phi
  Double_t                    fClusTimeCutLow;             // low time cut for clusters
  Double_t                    fClusTimeCutUp;              // up time cut for clusters
  Double_t                    fMinPtTrackInEmcal;          // min pt track in emcal
  Double_t                    fEventPlaneVsEmcal;          // select events which have a certain event plane wrt the emcal
  Double_t                    fMinEventPlane;              // minimum event plane value
  Double_t                    fMaxEventPlane;              // maximum event plane value
  TString                     fCentEst;                    // name of V0 centrality estimator
  UInt_t                      fTrackBitMap;                // bit map of accepted tracks (non MC)
  UInt_t                      fClusterBitMap;              // bit map of accepted clusters (non MC)
  UInt_t                      fMCTrackBitMap;              // bit map of accepted MC tracks
  UInt_t                      fMCClusterBitMap;            // bit map of accepted MC clusters
  Int_t                       fMinMCLabel;                 // minimum MC label value for the tracks/clusters being considered MC particles
  Int_t                       fNcentBins;                  //!how many centrality bins
  AliEMCALGeometry           *fGeom;                       //!emcal geometry
  TClonesArray               *fTracks;                     //!tracks
  TClonesArray               *fCaloClusters;               //!clusters
  AliVCaloCells              *fCaloCells;                  //!cells
  Double_t                    fCent;                       //!event centrality
  Int_t                       fCentBin;                    //!event centrality bin
  Double_t                    fEPV0;                       //!event plane V0
  Double_t                    fEPV0A;                      //!event plane V0A
  Double_t                    fEPV0C;                      //!event plane V0C
  Double_t                    fVertex[3];                  //!event vertex
  Int_t                       fNVertCont;                  //!event vertex number of contributors
  BeamType                    fBeamType;                   //!event beam type

  TObjArray                   fParticleCollArray;          // particle/track collection array
  TObjArray                   fClusterCollArray;           // cluster collection array

  TList                      *fOutput;                     //!output list

  TH1F                       *fHistCentrality;             //!Event centrality distribution
  TH1F                       *fHistZVertex;                //!Z vertex position
  TH1F                       *fHistEventPlane;             //!Event plane distribution

 private:
  AliAnalysisTaskEmcalDev(const AliAnalysisTaskEmcalDev&);            // not implemented
  AliAnalysisTaskEmcalDev &operator=(const AliAnalysisTaskEmcalDev&); // not implemented

  ClassDef(AliAnalysisTaskEmcalDev, 1) // EMCAL base analysis task
};
#endif
