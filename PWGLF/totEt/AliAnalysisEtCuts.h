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
#include <iostream>

class AliAnalysisEtCuts : public TNamed
{
 public:
   
  AliAnalysisEtCuts();
  virtual ~AliAnalysisEtCuts();

  virtual void SetPbPbDefaults();

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
  Char_t GetPhosClusterType() const { return fReconstructedPhosClusterType; }
  Double_t GetReconstructedPhosClusterEnergyCut() const { return fReconstructedPhosClusterEnergyCut; }
  Double_t GetReconstructedPhosSingleCellEnergyCut() const { return fReconstructedPhosSingleCellEnergyCut; }
  Double_t GetPhosTrackDistanceCut() const { return fPhosTrackDistanceCut; }
  Double_t GetPhosTrackDxCut() const { return fPhosTrackDxCut; }
  Double_t GetPhosTrackDzCut() const { return fPhosTrackDzCut; }
  Double_t GetPhosTrackRCut() const { return fPhosTrackRCut; }
  
  Double_t GetPhosBadDistanceCut() const { return fPhosBadDistanceCut; }
  
  // ReconstructedEmcal
  Char_t GetEmcalClusterType() const { return fReconstructedEmcalClusterType; }
  Double_t GetReconstructedEmcalClusterEnergyCut() const { return fReconstructedEmcalClusterEnergyCut; }
  Double_t GetReconstructedEmcalSingleCellEnergyCut() const { return fReconstructedEmcalSingleCellEnergyCut; }
  Double_t GetEmcalTrackDistanceCut() const { return fEmcalTrackDistanceCut; }
  Double_t GetEmcalTrackDxCut() const { return fEmcalTrackDxCut; }
  Double_t GetEmcalTrackDzCut() const { return fEmcalTrackDzCut; }
  
  // MonteCarlo
  Double_t GetMonteCarloSingleChargedParticle() const { return fMonteCarloSingleChargedParticle; }
  Double_t GetMonteCarloNeutralParticle() const { return fMonteCarloNeutralParticle; }
  // Hist: TTree and histogram info
  Bool_t GetHistMakeTree() const { return fHistMakeTree; }
  Bool_t GetHistMakeTreeDeposit() const { return fHistMakeTreeDeposit; }
  //
  Int_t GetHistNbinsMult() const { return fHistNbinsMult; }
  Double_t GetHistMinMult() const { return fHistMinMult; }
  Double_t GetHistMaxMult() const { return fHistMaxMult; }
  //
  Int_t GetHistNbinsTotEt() const { return fHistNbinsTotEt; }
  Double_t GetHistMinTotEt() const { return fHistMinTotEt; }
  Double_t GetHistMaxTotEt() const { return fHistMaxTotEt; }
  //
  Int_t GetHistNbinsParticleEt() const { return fHistNbinsParticleEt; }
  Double_t GetHistMinParticleEt() const { return fHistMinParticleEt; }
  Double_t GetHistMaxParticleEt() const { return fHistMaxParticleEt; }
  //
  Int_t GetHistNbinsParticlePt() const { return fHistNbinsParticlePt; }
  Double_t GetHistMinParticlePt() const { return fHistMinParticlePt; }
  Double_t GetHistMaxParticlePt() const { return fHistMaxParticlePt; }
  
  
  
  Short_t GetDetectorPhos() const { return fgkDetectorPhos; }
  Short_t GetDetectorEmcal() const { return fgkDetectorEmcal; }
  
  Double_t GetPrimaryVertexCutXY() const { return fPrimaryVertexCutXY; }
  Double_t GetPrimaryVertexCutZ() const { return fPrimaryVertexCutZ; }
  

  // Setters
  // Common
  void SetCommonEtaCut(Double_t val) { fCommonEtaCut = val; }
  void SetCommonClusterEnergyCut(Double_t val) { fCommonClusterEnergyCut = val; }
  void SetCommonTrackPtCut(Double_t val) { fCommonTrackPtCut = val; }
  void SetCommonSingleCell(Int_t val) { fCommonSingleCell = val;}
  // GeometryPhos
  void SetGeometryPhosEtaAccCut(Double_t val) { fGeometryPhosEtaAccCut = val; }
  void SetGeometryPhosPhiAccMinCut(Double_t val) { fGeometryPhosPhiAccMinCut = val; }
  void SetGeometryPhosPhiAccMaxCut(Double_t val) { fGeometryPhosPhiAccMaxCut = val; }
  void SetGeometryPhosDetectorRadius(Double_t val) { fGeometryPhosDetectorRadius = val; }
  // GeometryEmcal
  void SetGeometryEmcalEtaAccCut(Double_t val) { fGeometryEmcalEtaAccCut = val; }
  void SetGeometryEmcalPhiAccMinCut(Double_t val) { fGeometryEmcalPhiAccMinCut = val; }
  void SetGeometryEmcalPhiAccMaxCut(Double_t val) { fGeometryEmcalPhiAccMaxCut = val; }
  void SetGeometryEmcalDetectorRadius(Double_t val) { fGeometryEmcalDetectorRadius = val; }
  // Reconstructed
  void SetReconstructedVertexXCut(Double_t val) { fReconstructedVertexXCut = val; }
  void SetReconstructedVertexYCut(Double_t val) { fReconstructedVertexYCut = val; }
  void SetReconstructedVertexZCut(Double_t val) { fReconstructedVertexZCut = val; }
  void SetReconstructedIPxyCut(Double_t val) { fReconstructedIPxyCut = val; }
  void SetReconstructedIPzCut(Double_t val) { fReconstructedIPzCut = val; }
  void SetReconstructedNTpcClustersCut(Int_t val) { fReconstructedNTpcClustersCut = val; }
  void SetReconstructedNItsClustersCut(Int_t val) { fReconstructedNItsClustersCut = val; }
  void SetReconstrucedPidCut(Double_t val) { fReconstructedPidCut = val; }
  // ReconstructedPhos
  void SetReconstructedPhosClusterType(Char_t val) { fReconstructedPhosClusterType = val; }
  void SetReconstructedPhosClusterEnergyCut(Double_t val) { fReconstructedPhosClusterEnergyCut = val; }
  void SetReconstructedPhosSingleCellEnergyCut(Double_t val) { fReconstructedPhosSingleCellEnergyCut = val; }
  void SetPhosTrackDistanceCut(Double_t val) { fPhosTrackDistanceCut = val; }
  void SetPhosTrackDxCut(Double_t val) { fPhosTrackDxCut = val; }
  void SetPhosTrackDzCut(Double_t val) { fPhosTrackDzCut = val; }
  void SetPhosTrackRCut(Double_t val) { std::cout << "Setting: " << val << std::endl; fPhosTrackRCut = val; }
  
  void SetPhosBadDistanceCut(Double_t val) { fPhosBadDistanceCut = val; }
  
  // ReconstructedEmcal
  void SetReconstructedEmcalClusterType(Char_t val) { fReconstructedEmcalClusterType = val; }
  void SetReconstructedEmcalClusterEnergyCut(Double_t val) { fReconstructedEmcalClusterEnergyCut = val; }
  void SetReconstructedEmcalSingleCellEnergyCut(Double_t val) { fReconstructedEmcalSingleCellEnergyCut = val; }
  void SetEmcalTrackDistanceCut(Double_t val) { fEmcalTrackDistanceCut = val; }
  // MonteCarlo
  void SetMonteCarloSingleChargedParticle(Double_t val) { fMonteCarloSingleChargedParticle = val; }
  void SetMonteCarloNeutralParticle(Double_t val) { fMonteCarloNeutralParticle = val; }
  // Hist: TTree and histogram info
  void SetHistMakeTree(Bool_t val) { fHistMakeTree = val; }
  void SetHistMakeTreeDeposit(Bool_t val) { fHistMakeTreeDeposit = val; }
  //
  void SetHistNbinsMult(Int_t val) { fHistNbinsMult = val; }
  void SetHistMinMult(Double_t val) { fHistMinMult = val; }
  void SetHistMaxMult(Double_t val) { fHistMaxMult = val; }
  //
  void SetHistNbinsTotEt(Int_t val) { fHistNbinsTotEt = val; }
  void SetHistMinTotEt(Double_t val) { fHistMinTotEt = val; }
  void SetHistMaxTotEt(Double_t val) { fHistMaxTotEt = val; }
  //
  void SetHistNbinsParticleEt(Int_t val) { fHistNbinsParticleEt = val; }
  void SetHistMinParticleEt(Double_t val) { fHistMinParticleEt = val; }
  void SetHistMaxParticleEt(Double_t val) { fHistMaxParticleEt = val; }
  //
  void SetHistNbinsParticlePt(Int_t val) { fHistNbinsParticlePt = val; }
  void SetHistMinParticlePt(Double_t val) { fHistMinParticlePt = val; }
  void SetHistMaxParticlePt(Double_t val) { fHistMaxParticlePt = val; }

  void SetPrimaryVertexCutXY(Double_t val) { fPrimaryVertexCutXY = val; }
  void SetPrimaryVertexCutZ(Double_t val) { fPrimaryVertexCutZ = val; }
  
  
  
 protected:

  // Common   
  Double_t fCommonEtaCut; // Eta cut
  Double_t fCommonClusterEnergyCut; // Cluster Energy cut
  Double_t fCommonTrackPtCut; // Track Pt
  Int_t fCommonSingleCell; // Single Cell (1)
  Double_t fEmcalTrackDistanceCut; // EMCal track distance
  Double_t fEmcalTrackDxCut; // EMCal track distance in x 
  Double_t fEmcalTrackDzCut; // EMCal track distance in z
  
  Double_t fPhosTrackDistanceCut; // PHOS track distance  
  Double_t fPhosTrackDxCut; // PHOS track distance in x
  Double_t fPhosTrackDzCut; // PHOS track distance  in z
  Double_t fPhosTrackRCut; // PHOS track distance  in r (using the parametrized track distance)
 
 Double_t fPhosBadDistanceCut; // PHOS distance to bad channel 
 
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
  Double_t fReconstructedPhosTrackDistanceTightCut; // PHOS track distance
  Double_t fReconstructedPhosTrackDistanceMediumCut; // PHOS track distance
  Double_t fReconstructedPhosTrackDistanceLooseCut; // PHOS track distance

  // ReconstructedEmcal
  Char_t fReconstructedEmcalClusterType; // EMCal cluster type
  Double_t fReconstructedEmcalClusterEnergyCut; // EMCal cluster energy
  Double_t fReconstructedEmcalSingleCellEnergyCut; // EMCal single cell energy
  Double_t fReconstructedEmcalTrackDistanceTightCut; // EMCAL track distance
  Double_t fReconstructedEmcalTrackDistanceMediumCut; // EMCAL track distance
  Double_t fReconstructedEmcalTrackDistanceLooseCut; // EMCAL track distance

  // MonteCarlo
  Double_t fMonteCarloSingleChargedParticle; // MC charged
  Double_t fMonteCarloNeutralParticle; // MC neutral

  // Hist: TTree and histogram info
  Bool_t fHistMakeTree; // whether to make a summary tree or not
  Bool_t fHistMakeTreeDeposit; // whether to make a summary tree of energy deposit or not
  
  Int_t fHistNbinsMult; // number of bins in multiplicity histograms
  Double_t fHistMinMult; // minimum value in multiplicity histograms
  Double_t fHistMaxMult; // maximum value in multiplicity histograms

  Int_t fHistNbinsTotEt; // number of bins in event Et histograms
  Double_t fHistMinTotEt; // minimum value in event Et histograms
  Double_t fHistMaxTotEt; // maximum value in event Et histograms

  Int_t fHistNbinsParticleEt; // number of bins in particle Et histograms
  Double_t fHistMinParticleEt; // minimum value in particle Et histograms
  Double_t fHistMaxParticleEt; // maximum value in particle Et histograms

  Int_t fHistNbinsParticlePt; // number of bins in particle Pt histograms
  Double_t fHistMinParticlePt; // minimum value in particle Pt histograms
  Double_t fHistMaxParticlePt; // maximum value in particle Pt histograms

// Detector definition
  static const Short_t fgkDetectorPhos = -1; // PHOS 
  static const Short_t fgkDetectorEmcal = 1; // EMCAL 
  
  Double_t fPrimaryVertexCutXY; // Cut to decide if particle is from primary vertex
  Double_t fPrimaryVertexCutZ; // Cut to decide if particle is from primary vertex
  

private:
  //Declare private to avoid compilation warning
  AliAnalysisEtCuts & operator = (const AliAnalysisEtCuts & g) ;//copy assignment
  AliAnalysisEtCuts(const AliAnalysisEtCuts & g) ; // copy ctor

  ClassDef(AliAnalysisEtCuts, 1);
};

#endif // ALIANALYSISETCUTS_H
