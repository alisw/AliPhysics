/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// A class for calculating acceptance correction


#ifndef ALIJACCEPTANCECORRECTION_H
#define ALIJACCEPTANCECORRECTION_H

#ifndef ROOT_TObject
#include <TObject.h>
#endif

#include "AliJHistManager.h"
#include "AliJCard.h"
#include "AliJConst.h"

class AliJAcceptanceCorrection{

public:


  AliJAcceptanceCorrection();                     // default constructor
  AliJAcceptanceCorrection(AliJCard *inputCard);  // constructor with card definition
  AliJAcceptanceCorrection(const AliJAcceptanceCorrection& a); // copy constructor
  ~AliJAcceptanceCorrection();    //destructor

  void ReadMixedEventHistograms(const char *fileName);      // Read acceptance histograms from input file
  double GetAcceptanceCorrectionTriangle(double deltaEta);  // Calculate acceptance correction from deltaEta triangle
  double GetAcceptanceCorrectionTriangleInclusive(double deltaEta, int centralityBin, int triggerBin, int assocBin);  // Calculate the acceptance correction from a one dimensional input histogram
  double GetAcceptanceCorrectionTraditionalInclusive(double deltaEta, double deltaPhi, int centralityBin, int triggerBin);  // Calculate the acceptance correction from a two dimensional input histogram
  double GetAcceptanceCorrectionTraditionalInclusive(double deltaEta, double deltaPhi, int centralityBin, int zVertexBin, int triggerBin);  // Calculate the acceptance correction from a two dimensional input histogram
  double GetAcceptanceCorrection3DNearSideCalculation(double deltaEta, double deltaPhi); // Calculate the acceptance correction from geometry for 3D near side
  double GetAcceptanceCorrection3DNearSideInclusive(double deltaEta, double deltaPhi, int centralityBin, int triggerBin);  // Calculate the acceptance correction from a z-vertex summed two dimensional histogram
  double GetAcceptanceCorrection3DNearSideInclusive(double deltaEta, double deltaPhi, int centralityBin, int zVertexBin, int triggerBin);  // Calculate the acceptance correction from a two dimensional histogram

  double GetAcceptanceCorrectionTraditional(int samplingMethod, double deltaEta, double deltaPhi, int centralityBin, int triggerBin);  // Acceptance correction for traditional near side
  double GetAcceptanceCorrectionTraditional(int samplingMethod, double deltaEta, double deltaPhi, int centralityBin, int zVertexBin, int triggerBin);  // Acceptance correction for traditional near side
  double GetAcceptanceCorrection3DNearSide(int samplingMethod, double deltaEta, double deltaPhi, int centralityBin, int triggerBin);  // Acceptance correction for 3D near side
  double GetAcceptanceCorrection3DNearSide(int samplingMethod, double deltaEta, double deltaPhi, int centralityBin, int zVertexBin, int triggerBin);  // Acceptance correction for 3D near side
  
  double GetAcceptanceCorrection(int nearSideDefinition, int samplingMethod, double deltaEta, double deltaPhi, int centralityBin, int triggerBin);  // Acceptance correction using near side and sampling definitions as given in parameters
  double GetAcceptanceCorrection(int nearSideDefinition, int samplingMethod, double deltaEta, double deltaPhi, int centralityBin, int zVertexBin, int triggerBin);  // Acceptance correction using near side and sampling definitions as given in parameters
  
  AliJAcceptanceCorrection& operator=(const AliJAcceptanceCorrection& a); // Equal sign operator
  TH2D *GetAcceptanceHistogram(){return fDEtaDPhi3DNearAcceptanceCalculation;} // Getter for calculated acceptance

  void SetMinCountsPerBinInclusive(int minCounts){ fMinCountsPerBinInclusive = minCounts; } // Setter for fMinCountsPerBinInclusive
  
  void SetLeadingParticle(bool leadingParticle){ fLeadingParticleCorrelation = leadingParticle; } // Setter for fLeadingParticleCorrelation
  void SetTestMode(bool mode){ fTestMode = mode; } // Setter for fTestMode

private:
  void NormalizeAcceptanceTraditional(AliJTH1D &acceptanceHisto, corrType assocType); // Normalize one dimensional histograms to interval [0,1]
  void NormalizeAcceptanceInclusive(AliJTH2D &acceptanceHisto, corrType assocType, double peakValue); // Normalize two dimensional histograms to interval [0,peakValue]
  void NormalizeAcceptanceTraditionalInclusive(AliJTH2D &acceptanceHisto, corrType assocType); // Normalize two dimensional histograms to interval [0,1]
  void NormalizeAcceptance3DNearSideInclusive(AliJTH2D &acceptanceHisto, corrType assocType); // Normalize two dimensional 3D near side histograms according to acceptance limits
  void Generate3DAcceptanceCorrection(); // Calculate 3D near side acceptance correction and store it in 2D histogram
  int GetRebin(double counts, int nBins, int dimension); // Get rebinning factor for histogram
  void RebinAndNormalize(TH2 *histogram, double peakValue); // Rebin and normalize two dimensional histogram
  double GetAcceptanceCorrection3DNearSideInclusiveBin(double deltaEta, double deltaPhi, int centralityBin, int zVertexBin, int triggerBin, int firstBin);  // Common correction getter for z-vertex summed and z-vertex binned histograms
  double GetAcceptanceCorrectionTraditionalInclusiveBin(double deltaEta, double deltaPhi, int centralityBin, int zVertexBin, int triggerBin, int firstBin);  // Common correction getter for z-vertex summed and z-vertex binned histograms
  
  AliJCard *fCard;                    // Card containing analysis details
  AliJTH1D fDEtaNearAcceptance;       // DeltaEta acceptance histogram for the near side
  AliJTH2D fDEtaDPhiNearAcceptance;   // DeltaEta DeltaPhi acceptance histogram for the near side
  AliJTH2D fDEtaDPhi3DNearAcceptance; // DeltaEta DeltaPhi acceptance histogram for the 3D near side
  
  TH2D *fDEtaDPhi3DNearAcceptanceCalculation; // Calculated acceptance correction histogram for 3D near side
  
  int fMinCountsPerBinInclusive; // Minimum number of counts per histogram bin in inclusive deltaEta deltaPhi histograms
  
  bool fDEtaNearLoaded; // Tells whether the one dimensional traditional near side mixed event histogram is loaded or not
  bool fDEtaDPhiNearLoaded; // Tells whether the two dimensional traditional near side mixed event histogram is loaded or not
  bool fDEtaDPhi3DNearLoaded; // Tells whether the two dimensional 3D near side mixed event histogram is loaded or not
  bool fLeadingParticleCorrelation; // Tells whether the acceptance histograms are done with leading particle
  bool fTestMode; // True = mixed event distributions not properly corrected in 3D near side. This setting can be used to test the effect of doing this correctly vs. doing thing wrong
};

#endif
