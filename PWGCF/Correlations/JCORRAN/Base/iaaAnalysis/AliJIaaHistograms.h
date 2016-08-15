/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Container class for histograms needed in the jT analysis.

//===========================================================
// AliJIaaHistograms.h
//
// author: Jussi Viinikainen
//===========================================================

#ifndef ALIJIAAHISTOGRAMS_H
#define ALIJIAAHISTOGRAMS_H

#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <stdio.h>

#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>
#include <TFile.h>
#include <TList.h>
#include <TLorentzVector.h>

#include "../AliJConst.h"
#include "../AliJHistogramInterface.h"

class AliJCard;
class AliJBaseTrack;
class AliJPhoton;
class AliJTrack;

using namespace std;

class AliJIaaHistograms : public AliJHistogramInterface{
  
public:
  AliJIaaHistograms(AliJCard* cardP); //constructor
  AliJIaaHistograms(const AliJIaaHistograms& obj); // copy constructor
	virtual ~AliJIaaHistograms();    //destructor
  AliJIaaHistograms& operator=(const AliJIaaHistograms& obj); // equal sign operator
  
  // create histograms
  void CreateCorrelationHistograms(); // Correlation histogram creation
  void CreateEventTrackHistos(); // Event track histogram creation

  void ReadInclusiveHistos(const char *inclusFileName); // Inclusive histogram reader

  //bool Is2DHistosEnabled(){ return fenable2DHistos; } // Getter for 2D histogram enabling flag
  //void Set2DHistoCreate(bool isenable) { fenable2DHistos = isenable; } // Set 2D histograms for creation
  //void SetAcceptanceCorrectionQA(bool isenable) { fEnableAcceptanceQAHistos = isenable; } // Set flag for QA histogram creation
  
public:

  //==================================================
  // 1D correlation histograms
  //==================================================

  AliJTH1D fhDphiAssoc;
  AliJTH1D fhDEtaNear;
  AliJTH1D fhDEtaNearM;
  // JV // One dimensional deltaEta histograms for acceptance correction
  AliJTH1D fhDetaNearMixAcceptance;   // Mixed event uncorrected deltaEta histogram for acceptance correction

  //==================================================
  // 2D correlation histograms
  //==================================================

  AliJTH2D fhDphiDetaPta;      // 2D histogram (deltaPhi-deltaEta plane) in [fTyp][fCentralityBin][fZBin][fpttBin][fptaBin] bins
  //AliJTH2D fhDphiDetaPtaRgap;    // same with 3D near-side definition
  
  //==================================================
  // Trigger and associated particle specra
  //==================================================
  
  AliJTH1D fhIphiTrigg; // Phi distribution for trigger particle
  AliJTH1D fhIetaTrigg; // Eta distribution for trigger particle
  AliJTH1D fhIphiAssoc; // Phi distribution for associated particle
  AliJTH1D fhIetaAssoc; // Eta distribution for associated particle
  AliJTH1D fhTriggPtBin; // pTt distribution inside a trigger pT bin
  AliJTH1D fhAssocPtBin; // pTa distribution in pTt and pTa bins
  
  //==================================================
  // Inclusive spectra
  //==================================================

  AliJHistManager *fHmgInclusive; // Histogram manager for inclusive histograms
  AliJTH1D fhIetaTriggFromFile;   // Trigger inclusive eta distribution
  AliJTH1D fhIetaAssocFromFile;   // Associated inclusive eta distribution
  AliJTH1D fhIphiTriggFromFile;   // Trigger inclusive phi distribution
  AliJTH1D fhIphiAssocFromFile;   // Associated inclusive phi distribution
  
  //===================================================
  // Event/Track histograms
  //===================================================

  AliJTH1D fhLPpt; // pT distribution of leading particles
  AliJTH1D fhChargedPt, fhChargedPtNoCorr; // Corrected and raw pT distribution of charged particles
  AliJTProfile fhTrackingEfficiency; // Tracking efficiency
  AliJTH1D fhChargedEta; // Charged particle pseudorapidity distribution
  AliJTH1D fhLPeta; // Leading particle eta distribution
  AliJTH1D fhChargedMult; // Charged particle multiplicity distribution
  AliJTH1D fhCentr;
  AliJTH1D fhiCentr;
  AliJTH1D fhZVert; // z-vertex distribution
  

	//===================================================
	// Acceptance correction QA histograms
	// Not needed in the actual analysis, only for check
	//===================================================

	//AliJTH1D fhAcceptanceTraditional;   // Acceptance correction check for traditional 1D triangle
	//AliJTH2D fhAcceptanceTraditional2D; // Acceptance correction check for traditional deltaEta deltaPhi distribution
	//AliJTH2D fhAcceptance3DNearSide; // Acceptance correction check for 3D near side deltaEta deltaPhi distribution
	//AliJTH2D fhAcceptanceTraditional2DZ; // Z-vertex binned acceptance correction check for traditional deltaEta deltaPhi distribution
	//AliJTH2D fhAcceptance3DNearSideZ; // Z-vertex binned acceptance correction check for 3D near side deltaEta deltaPhi distribution

protected:

  double fmaxEtaRange;              // maximum eta range
  
};

#endif






















