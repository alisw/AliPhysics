#ifndef ALIANALYSISETCUTS_H
#define ALIANALYSISETCUTS_H
//_________________________________________________________________________
//  Utility Class for transverse energy studies
//  Base class for ESD & MC analysis
//  - cuts for reconstruction and MonteCarlo 
//  
//*-- Authors: Oystein Djuvsland (Bergen), David Silvermyr (ORNL)
//_________________________________________________________________________

#include "TNamed.h"

class AliAnalysisEtCuts : public TNamed
{
 public:
   
  AliAnalysisEtCuts();
  virtual ~AliAnalysisEtCuts();

  // Getters
  // Common
  Double_t GetCommonEtaCut() const { return fCommonEtaCut; }
  Double_t GetCommonClusterEnergyCut() const { return fCommonClusterEnergyCut; }
  Double_t GetCommonTrackPtCut() const { return fCommonTrackPtCut; }
  Int_t GetCommonSingleCell() const { return fCommonSingleCell; }

  // GeometryPhos
  Double_t GetGeometryPhosEtaAccCut() const { return fGeometryPhosEtaAccCut; }
  Double_t GetGeometryPhosPhiAccMinCut() const { return fGeometryPhosPhiAccMinCut; }
  Double_t GetGeometryPhosPhiAccMaxCut() const { return fGeometryPhosPhiAccMaxCut; }
  Double_t GetGeometryPhosDetectorRadius() const { return fGeometryPhosDetectorRadius; }
  // GeometryEmcal
  Double_t GetGeometryEmcalEtaAccCut() const { return fGeometryEmcalEtaAccCut; }
  Double_t GetGeometryEmcalPhiAccMinCut() const { return fGeometryEmcalPhiAccMinCut; }
  Double_t GetGeometryEmcalPhiAccMaxCut() const { return fGeometryEmcalPhiAccMaxCut; }
  Double_t GetGeometryEmcalDetectorRadius() const { return fGeometryEmcalDetectorRadius; }
  // Reconstructed
  Double_t GetReconstructedVertexXCut() const { return fReconstructedVertexXCut; }
  Double_t GetReconstructedVertexYCut() const { return fReconstructedVertexYCut; }
  Double_t GetReconstructedVertexZCut() const { return fReconstructedVertexZCut; }
  Double_t GetReconstructedIPxyCut() const { return fReconstructedIPxyCut; }
  Double_t GetReconstructedIPzCut() const { return fReconstructedIPzCut; }
  Int_t GetReconstructedNTpcClustersCut() const { return fReconstructedNTpcClustersCut; }
  Int_t GetReconstructedNItsClustersCut() const { return fReconstructedNItsClustersCut; }
  Double_t GetReconstructedPidCut() const { return fReconstructedPidCut; }
  // ReconstructedPhos
  Char_t GetReconstructedPhosClusterType() const { return fReconstructedPhosClusterType; }
  Double_t GetReconstructedPhosClusterEnergyCut() const { return fReconstructedPhosClusterEnergyCut; }
  Double_t GetReconstructedPhosSingleCellEnergyCut() const { return fReconstructedPhosSingleCellEnergyCut; }
  Double_t GetReconstructedPhosTrackDistanceCut() const { return fReconstructedPhosTrackDistanceCut; }
  // ReconstructedEmcal
  Char_t GetReconstructedEmcalClusterType() const { return fReconstructedEmcalClusterType; }
  Double_t GetReconstructedEmcalClusterEnergyCut() const { return fReconstructedEmcalClusterEnergyCut; }
  Double_t GetReconstructedEmcalSingleCellEnergyCut() const { return fReconstructedEmcalSingleCellEnergyCut; }
  Double_t GetReconstructedEmcalTrackDistanceCut() const { return fReconstructedEmcalTrackDistanceCut; }
  // MonteCarlo
  Double_t GetMonteCarloSingleChargedParticle() const { return fMonteCarloSingleChargedParticle; }
  Double_t GetMonteCarloNeutralParticle() const { return fMonteCarloNeutralParticle; }
  // Hist: TTree and histogram info
  Bool_t GetHistMakeTree() const { return fHistMakeTree; }

  // Setters
  // Common
  void SetCommonEtaCut(const Double_t val) { fCommonEtaCut = val; }
  void SetCommonClusterEnergyCut(const Double_t val) { fCommonClusterEnergyCut = val; }
  void SetCommonTrackPtCut(const Double_t val) { fCommonTrackPtCut = val; }
  void SetCommonSingleCell(const Int_t val) { fCommonSingleCell = val;}
  // GeometryPhos
  void SetGeometryPhosEtaAccCut(const Double_t val) { fGeometryPhosEtaAccCut = val; }
  void SetGeometryPhosPhiAccMinCut(const Double_t val) { fGeometryPhosPhiAccMinCut = val; }
  void SetGeometryPhosPhiAccMaxCut(const Double_t val) { fGeometryPhosPhiAccMaxCut = val; }
  void SetGeometryPhosDetectorRadius(const Double_t val) { fGeometryPhosDetectorRadius = val; }
  // GeometryEmcal
  void SetGeometryEmcalEtaAccCut(const Double_t val) { fGeometryEmcalEtaAccCut = val; }
  void SetGeometryEmcalPhiAccMinCut(const Double_t val) { fGeometryEmcalPhiAccMinCut = val; }
  void SetGeometryEmcalPhiAccMaxCut(const Double_t val) { fGeometryEmcalPhiAccMaxCut = val; }
  void SetGeometryEmcalDetectorRadius(const Double_t val) { fGeometryEmcalDetectorRadius = val; }
  // Reconstructed
  void SetReconstructedVertexXCut(const Double_t val) { fReconstructedVertexXCut = val; }
  void SetReconstructedVertexYCut(const Double_t val) { fReconstructedVertexYCut = val; }
  void SetReconstructedVertexZCut(const Double_t val) { fReconstructedVertexZCut = val; }
  void SetReconstructedIPxyCut(const Double_t val) { fReconstructedIPxyCut = val; }
  void SetReconstructedIPzCut(const Double_t val) { fReconstructedIPzCut = val; }
  void SetReconstructedNTpcClustersCut(const Int_t val) { fReconstructedNTpcClustersCut = val; }
  void SetReconstructedNItsClustersCut(const Int_t val) { fReconstructedNItsClustersCut = val; }
  void SetReconstrucedPidCut(const Double_t val) { fReconstructedPidCut = val; }
  // ReconstructedPhos
  void SetReconstructedPhosClusterType(const Char_t val) { fReconstructedPhosClusterType = val; }
  void SetReconstructedPhosClusterEnergyCut(const Double_t val) { fReconstructedPhosClusterEnergyCut = val; }
  void SetReconstructedPhosSingleCellEnergyCut(const Double_t val) { fReconstructedPhosSingleCellEnergyCut = val; }
  void SetReconstructedPhosTrackDistanceCut(const Double_t val) { fReconstructedPhosTrackDistanceCut = val; }
  // ReconstructedEmcal
  void SetReconstructedEmcalClusterType(const Char_t val) { fReconstructedEmcalClusterType = val; }
  void SetReconstructedEmcalClusterEnergyCut(const Double_t val) { fReconstructedEmcalClusterEnergyCut = val; }
  void SetReconstructedEmcalSingleCellEnergyCut(const Double_t val) { fReconstructedEmcalSingleCellEnergyCut = val; }
  void SetReconstructedEmcalTrackDistanceCut(const Double_t val) { fReconstructedEmcalTrackDistanceCut = val; }
  // MonteCarlo
  void SetMonteCarloSingleChargedParticle(const Double_t val) { fMonteCarloSingleChargedParticle = val; }
  void SetMonteCarloNeutralParticle(const Double_t val) { fMonteCarloNeutralParticle = val; }
  // Hist: TTree and histogram info
  void SetHistMakeTree(const Bool_t val) { fHistMakeTree = val; }

 protected:

  // Common   
  Double_t fCommonEtaCut; // Eta cut
  Double_t fCommonClusterEnergyCut; // Cluster Energy cut
  Double_t fCommonTrackPtCut; // Track Pt
  Int_t fCommonSingleCell; // Single Cell (1)
  
  // GeometryPhos
  Double_t fGeometryPhosEtaAccCut; // PHOS Eta Acc cut
  Double_t fGeometryPhosPhiAccMinCut; // PHOS Phi Acc Min cut
  Double_t fGeometryPhosPhiAccMaxCut; // PHOS Phi Acc Max cut
  Double_t fGeometryPhosDetectorRadius; // PHOS Detector Radius 

  // GeometryEmcal
  Double_t fGeometryEmcalEtaAccCut; // EMCal Eta Acc cut
  Double_t fGeometryEmcalPhiAccMinCut; // EMCal Phi Acc Min cut
  Double_t fGeometryEmcalPhiAccMaxCut; // EMCal Phi Acc Max cut
  Double_t fGeometryEmcalDetectorRadius; // EMCal Detector Radius

  // Reconstructed
  Double_t fReconstructedVertexXCut; // vertex X cut
  Double_t fReconstructedVertexYCut; // vertex Y cut
  Double_t fReconstructedVertexZCut; // vertex Z cut
  Double_t fReconstructedIPxyCut; // IP xy cut
  Double_t fReconstructedIPzCut; // IP z cut
  Int_t fReconstructedNTpcClustersCut; // # of TPC clusters cut
  Int_t fReconstructedNItsClustersCut; // # of ITS clusters cut
  Double_t fReconstructedPidCut; // cut on pid prob

  // ReconstructedPhos
  Char_t fReconstructedPhosClusterType; // PHOS cluster type
  Double_t fReconstructedPhosClusterEnergyCut; // PHOS cluster energy
  Double_t fReconstructedPhosSingleCellEnergyCut; // PHOS single cell energy
  Double_t fReconstructedPhosTrackDistanceCut; // PHOS track distance

  // ReconstructedEmcal
  Char_t fReconstructedEmcalClusterType; // EMCal cluster type
  Double_t fReconstructedEmcalClusterEnergyCut; // EMCal cluster energy
  Double_t fReconstructedEmcalSingleCellEnergyCut; // EMCal single cell energy
  Double_t fReconstructedEmcalTrackDistanceCut; // EMCal track distance

  // MonteCarlo
  Double_t fMonteCarloSingleChargedParticle; // MC charged
  Double_t fMonteCarloNeutralParticle; // MC neutral

  // Hist: TTree and histogram info
  Bool_t fHistMakeTree; // whether to make a summary tree or not

private:
  //Declare private to avoid compilation warning
  AliAnalysisEtCuts & operator = (const AliAnalysisEtCuts & g) ;//copy assignment
  AliAnalysisEtCuts(const AliAnalysisEtCuts & g) ; // copy ctor

  ClassDef(AliAnalysisEtCuts, 0);
};

#endif // ALIANALYSISETCUTS_H
