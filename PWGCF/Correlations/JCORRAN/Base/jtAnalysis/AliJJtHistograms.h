/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Container class for histograms needed in the jT analysis.

//===========================================================
// AliJJtHistograms.h
//
// author: Jussi Viinikainen
//===========================================================

#ifndef ALIJJTHISTOGRAMS_H
#define ALIJJTHISTOGRAMS_H

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

class AliJJtHistograms : public AliJHistogramInterface{
  
public:
  AliJJtHistograms(AliJCard* cardP); //constructor
  AliJJtHistograms(const AliJJtHistograms& obj); // copy constructor
    virtual ~AliJJtHistograms();    //destructor
  AliJJtHistograms& operator=(const AliJJtHistograms& obj); // equal sign operator
  
  // create histograms
  void CreateCorrelationHistograms(); // Correlation histogram creation
  void CreateEventTrackHistos(); // Event track histogram creation

  void ReadInclusiveHistos(const char *inclusFileName); // Inclusive histogram reader

  bool Is2DHistosEnabled(){ return fenable2DHistos; } // Getter for 2D histogram enabling flag
  void Set2DHistoCreate(bool isenable) { fenable2DHistos = isenable; } // Set 2D histograms for creation
  void SetAcceptanceCorrectionQA(bool isenable) { fEnableAcceptanceQAHistos = isenable; } // Set flag for QA histogram creation
  
public:
  
  // 2D dphi deta histograms for acceptance correction
  AliJTH2D fhDphiDetaKlong;    //! 2D histogram from deltaPhi-deltaEta plane in klong bins
  AliJTH2D fhDphiDetaXlong;    //! 2D histogram from deltaPhi-deltaEta plane in xlong bins
  AliJTH2D fhDphiDetaPta;      //! 2D histogram from deltaPhi-deltaEta plane in pta bins
  
  // JV // One dimensional deltaEta histograms for acceptance correction
  AliJTH1D fhDetaNearMixAcceptance;   //! Mixed event uncorrected deltaEta histogram for acceptance correction
  
  // JV // 2D dphi deta histos to study better the background
  AliJTH2D fhDphiDetaBgKlongEta; //! 2D histogram from deltaPhi-deltaEta plane in klong and eta gap bins
  AliJTH2D fhDphiDetaBgKlongR;   //! 2D histogram from deltaPhi-deltaEta plane in klong and R gap bins
  AliJTH2D fhDphiDetaBgKlongPhi; //! 2D histogram from deltaPhi-deltaEta plane in klong and phi gap bins
  AliJTH2D fhDphiDetaBgXlongEta; //! 2D histogram from deltaPhi-deltaEta plane in xlong and eta gap bins
  AliJTH2D fhDphiDetaBgXlongR;   //! 2D histogram from deltaPhi-deltaEta plane in xlong and R gap bins
  AliJTH2D fhDphiDetaBgXlongPhi; //! 2D histogram from deltaPhi-deltaEta plane in xlong and phi gap bins
  AliJTH2D fhDphiDetaBgPtaEta;   //! 2D histogram from deltaPhi-deltaEta plane in pta and eta gap bins
  AliJTH2D fhDphiDetaBgPtaR;     //! 2D histogram from deltaPhi-deltaEta plane in pta and R gap bins
  AliJTH2D fhDphiDetaBgPtaPhi;   //! 2D histogram from deltaPhi-deltaEta plane in pta and phi gap bins
  
  // JV // pTa distributions in background bins
  AliJTH1D fhBgAssocKlongEta; //! background pta distribution in klong and eta gap bins
  AliJTH1D fhBgAssocKlongR;   //! background pta distribution in klong and R gap bins
  AliJTH1D fhBgAssocKlongPhi; //! background pta distribution in klong and phi gap bins
  AliJTH1D fhBgAssocXlongEta; //! background pta distribution in xlong and eta gap bins
  AliJTH1D fhBgAssocXlongR;   //! background pta distribution in xlong and R gap bins
  AliJTH1D fhBgAssocXlongPhi; //! background pta distribution in xlong and phi gap bins
  AliJTH1D fhBgAssocPtaEta;   //! background pta distribution in pta and eta gap bins
  AliJTH1D fhBgAssocPtaR;     //! background pta distribution in pta and R gap bins
  AliJTH1D fhBgAssocPtaPhi;   //! background pta distribution in pta and phi gap bins
  
  // JV // Invariant mass histograms
  AliJTH1D fhInvariantMassXe;    //! Invariant mass histogram in xlong bins
  AliJTH1D fhInvariantMassKlong; //! Invariant mass histogram in klong bins
  AliJTH1D fhInvariantMassPta;   //! Invariant mass histogram in pta bins

  // Like sign bins for invariant mass
  AliJTH1D fhInvariantMassXeLikeSign;    //! Invariant mass histogram in xlong bins for like sign pairs
  AliJTH1D fhInvariantMassKlongLikeSign; //! Invariant mass histogram in klong bins for like sign pair
  AliJTH1D fhInvariantMassPtaLikeSign;   //! Invariant mass histogram in pta bins for like sign pairs

  // Unlike sign bins for invariant mass
  AliJTH1D fhInvariantMassXeUnlikeSign;    //! Invariant mass histogram in xlong bins for unlike sign pairs
  AliJTH1D fhInvariantMassKlongUnlikeSign; //! Invariant mass histogram in klong bins for unlike sign pair
  AliJTH1D fhInvariantMassPtaUnlikeSign;   //! Invariant mass histogram in pta bins for unlike sign pairs
  
  //==================================================
  // Trigger and associated particle specra
  //==================================================
  
  AliJTH1D fhIphiTrigg; //! Phi distribution for trigger particle
  AliJTH1D fhIetaTrigg; //! Eta distribution for trigger particle
  AliJTH1D fhIphiAssoc; //! Phi distribution for associated particle
  AliJTH1D fhIetaAssoc; //! Eta distribution for associated particle
  AliJTH1D fhTriggPtBin; //! pTt distribution inside a trigger pT bin
  AliJTH1D fhAssocPtBin; //! pTa distribution in pTt and pTa bins
  
  //==================================================
  // jT histograms, disributions and background
  //==================================================
  
  // jT distribution and background in xlong bins
  AliJTH1D     fhJT;      //! jt distribution in xlong bins
  AliJTH1D     fhJTBg;    //! eta gap background in xlong bins
  AliJTH1D     fhJTBgR;   //! R gap background in xlong bins
  AliJTH1D     fhJTBgPhi; //! phi gap background in xlong bins
 
  AliJTH1D     fhJTLikeSign;      //! jt distribution for like sign pairs
  AliJTH1D     fhJTBgLikeSign;    //! eta gap background for like sign pair
  AliJTH1D     fhJTBgRLikeSign;   //! R gap background for like sign pairs
  AliJTH1D     fhJTBgPhiLikeSign; //! phi gap background for like sign pairs

  AliJTH1D     fhJTUnlikeSign;      //! jt distribution for unlike sign pairs
  AliJTH1D     fhJTBgUnlikeSign;    //! eta gap background for unlike sign pair
  AliJTH1D     fhJTBgRUnlikeSign;   //! R gap background for unlike sign pairs
  AliJTH1D     fhJTBgPhiUnlikeSign; //! phi gap background for unlike sign pairs

  // jT distribution and background in klong bins
  AliJTH1D     fhJTKlong;      //! jT distribution in klong bins
  AliJTH1D     fhJTKlongBg;    //! eta gap background in klong bins
  AliJTH1D     fhJTKlongBgR;   //! R gap background in klong bins
  AliJTH1D     fhJTKlongBgPhi; //! phi gap background in klong bins
 
  AliJTH1D     fhJTKlongLikeSign;      //! jt distribution for like sign pairs
  AliJTH1D     fhJTKlongBgLikeSign;    //! eta gap background for like sign pair
  AliJTH1D     fhJTKlongBgRLikeSign;   //! R gap background for like sign pairs
  AliJTH1D     fhJTKlongBgPhiLikeSign; //! R-gap background for like sign pairs

  AliJTH1D     fhJTKlongUnlikeSign;      //! jt distribution for unlike sign pairs
  AliJTH1D     fhJTKlongBgUnlikeSign;    //! eta gap background for unlike sign pair
  AliJTH1D     fhJTKlongBgRUnlikeSign;   //! R gap background for unlike sign pairs
  AliJTH1D     fhJTKlongBgPhiUnlikeSign; //! phi gap background for unlike sign pairs

  // jT distribution and background in pta bins
  AliJTH1D     fhJTPta;      //! jT distribution in pTa bins
  AliJTH1D     fhJTPtaBg;    //! eta gap background in pTa bins
  AliJTH1D     fhJTPtaBgR;   //! R gap background in pTa bins
  AliJTH1D     fhJTPtaBgPhi; //! phi gap background in pTa bins
 
  AliJTH1D     fhJTPtaLikeSign;      //! jt distribution for like sign pairs
  AliJTH1D     fhJTPtaBgLikeSign;    //! eta gap background for like sign pair
  AliJTH1D     fhJTPtaBgRLikeSign;   //! R gap background for like sign pairs
  AliJTH1D     fhJTPtaBgPhiLikeSign; //! Rphi gap background for like sign pairs

  AliJTH1D     fhJTPtaUnlikeSign;      //! jt distribution for unlike sign pairs
  AliJTH1D     fhJTPtaBgUnlikeSign;    //! eta gap background for unlike sign pair
  AliJTH1D     fhJTPtaBgRUnlikeSign;   //! R gap background for unlike sign pairs
  AliJTH1D     fhJTPtaBgPhiUnlikeSign; //! phi gap background for unlike sign pairs
  
  //==================================================
  // Inclusive spectra
  //==================================================

  AliJHistManager *fHmgInclusive; //! Histogram manager for inclusive histograms
  AliJTH1D fhIetaTriggFromFile;   //! Trigger inclusive eta distribution
  AliJTH1D fhIetaAssocFromFile;   //! Associated inclusive eta distribution
  AliJTH1D fhIphiTriggFromFile;   //! Trigger inclusive phi distribution
  AliJTH1D fhIphiAssocFromFile;   //! Associated inclusive phi distribution
  
  //===================================================
  // Event/Track histograms
  //===================================================
  
  AliJTH1D fhLPpt; //! pT distribution of leading particles
  AliJTH1D fhChargedPt, fhChargedPtNoCorr; //! Corrected and raw pT distribution of charged particles
  AliJTProfile fhTrackingEfficiency; //! Tracking efficiency
  AliJTH1D fhChargedEta; //! Charged particle pseudorapidity distribution
  AliJTH1D fhLPeta; //! Leading particle eta distribution
  AliJTH1D fhChargedMult; //! Charged particle multiplicity distribution
  
  AliJTH1D fhZVert; //! z-vertex distribution

  //===================================================
  // Acceptance correction QA histograms
  // Not needed in the actual analysis, only for check
  //===================================================
  
  AliJTH1D fhAcceptanceTraditional;   //! Acceptance correction check for traditional 1D triangle
  AliJTH2D fhAcceptanceTraditional2D; //! Acceptance correction check for traditional deltaEta deltaPhi distribution
  AliJTH2D fhAcceptance3DNearSide; //! Acceptance correction check for 3D near side deltaEta deltaPhi distribution
  AliJTH2D fhAcceptanceTraditional2DZ; //! Z-vertex binned acceptance correction check for traditional deltaEta deltaPhi distribution
  AliJTH2D fhAcceptance3DNearSideZ; //! Z-vertex binned acceptance correction check for 3D near side deltaEta deltaPhi distribution

  
protected:
  double fmaxEtaRange;              // maximum eta range
  bool fenable2DHistos;             // enable the filling of two dimensional histograms
  bool fEnableAcceptanceQAHistos;   // enable the filling of acceptance correction quality control histograms
  
};

#endif






















