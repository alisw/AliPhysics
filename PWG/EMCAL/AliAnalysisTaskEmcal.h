#ifndef ALIANALYSISTASKEMCAL_H
#define ALIANALYSISTASKEMCAL_H

// $Id: AliAnalysisTaskEmcal.h 56756 2012-05-30 05:03:02Z loizides $

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

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskEmcal : public AliAnalysisTaskSE {
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

  AliAnalysisTaskEmcal();
  AliAnalysisTaskEmcal(const char *name, Bool_t histo=kFALSE); 
  virtual ~AliAnalysisTaskEmcal();

  void                        UserExec(Option_t *option);
  void                        UserCreateOutputObjects();

  void                        SetAnaType(EmcalAnaType type)                         { fAnaType           = type ;                         ; }
  void                        SetCentRange(Double_t min, Double_t max)              { fMinCent           = min  ; fMaxCent = max          ; }
  void                        SetClusName(const char *n)                            { fCaloName          = n                              ; }
  void                        SetCaloCellsName(const char *n)                       { fCaloCellsName     = n                              ; }
  void                        SetClusPtCut(Double_t cut)                            { fClusPtCut         = cut                            ; }
  void                        SetClusTimeCut(Double_t min, Double_t max)            { fClusTimeCutLow    = min  ; fClusTimeCutUp = max    ; }
  void                        SetHistoBins(Int_t nbins, Double_t min, Double_t max) { fNbins = nbins; fMinBinPt = min; fMaxBinPt = max    ; }
  void                        SetOffTrigger(UInt_t t)                               { fOffTrigger        = t                              ; }
  void                        SetPtCut(Double_t cut)                                { SetClusPtCut(cut)         ; SetTrackPtCut(cut)      ; }
  void                        SetTrackPtCut(Double_t cut)                           { fTrackPtCut        = cut                            ; }
  void                        SetTrackEtaLimits(Double_t min, Double_t max)         { fTrackMaxEta       = max  ; fTrackMinEta      = min ; }
  void                        SetTrackPhiLimits(Double_t min, Double_t max)         { fTrackMaxPhi       = max  ; fTrackMinPhi      = min ; }
  void                        SetTracksName(const char *n)                          { fTracksName        = n                              ; }
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

 protected:
  Bool_t                      AcceptCluster(AliVCluster        *clus)  const;
  Bool_t                      AcceptEmcalPart(AliEmcalParticle *part)  const;
  Bool_t                      AcceptTrack(AliVParticle         *track) const;
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
  TList                      *fOutput;                     //!output list

  TH1F                       *fHistCentrality;             //!Event centrality distribution
  TH1F                       *fHistZVertex;                //!Z vertex position
  TH1F                       *fHistEventPlane;             //!Event plane distribution

 private:
  AliAnalysisTaskEmcal(const AliAnalysisTaskEmcal&);            // not implemented
  AliAnalysisTaskEmcal &operator=(const AliAnalysisTaskEmcal&); // not implemented

  ClassDef(AliAnalysisTaskEmcal, 13) // EMCAL base analysis task
};
#endif
