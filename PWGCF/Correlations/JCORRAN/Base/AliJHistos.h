/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Interface for all analysis histograms.
// This is the minimum amount of histograms all analysis must have defined

//===========================================================
// AliJHistogramInterface.h
//
//   J
//===========================================================

#ifndef ALIJHISTOS_H
#define ALIJHISTOS_H

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

#include "AliJConst.h"
#include "AliJHistogramInterface.h"
#include "AliJHistManager.h"

class AliJCard;
class AliJBaseTrack;
class AliJPhoton;
class AliJTrack;

using namespace std;

class AliJHistos : public AliJHistogramInterface {
  
public:
  AliJHistos(AliJCard* cardP); //constructor
  virtual ~AliJHistos();    //destructor
  AliJHistos(const AliJHistos& obj);
  AliJHistos& operator=(const AliJHistos& obj);
  
  // create histograms
  void CreateAzimuthCorrHistos();
  void CreateIAAMoons();
  void CreateXEHistos();
  void CreateXtHistos();
  
  void CreateEventTrackHistos();
  void CreateJetHistos();
  
  void CreatePairPtCosThetaStar();
  void CreatePtCorrHistos();
  void CreateRunByRunHistos(int runID, int runcounter) const;
  
  void ReadInclusiveHistos(const char *inclusFileName);

  bool Is2DHistosEnabled(){ return fenable2DHistos; }
  void Set2DHistoCreate(bool isenable) { fenable2DHistos = isenable; }
  
public:
  
  AliJTProfile fhMixStat; //! comment me
  
  //==Pt stat. fcorrelations ===============================================
  AliJTH1D fTestHist; //! test only, remove it not needed
  AliJTH1D fhPtNear; //! comment me
  AliJTH1D fhPtFar ; //! comment me
  
  //assorted
  AliJTH1D fhPhi; //! comment me
  AliJTH1D fhDphiAssoc; //! comment me
  AliJTH1D fhDphiAssocXEbin; //! comment me
  //2D
  AliJTH2D fhDphiAssoc2DIAA; //! comment me
  AliJTH2D fhDphiAssoc2D; //! comment me
  
  // JV // 2D dphi deta histograms for acceptance correction
  AliJTH2D fhDphiDetaXlong;    //! 2D histogram from deltaPhi-deltaEta plane in xlong bins
  AliJTH2D fhDphiDetaPta;      //! 2D histogram from deltaPhi-deltaEta plane in pta bins
  
  // JV // One dimensional deltaEta histograms for acceptance correction
  AliJTH1D fhDetaNearMixAcceptance;   //! Mixed event uncorrected deltaEta histogram for acceptance correction
  AliJTH1D fhDeta3DNearMixAcceptance; //! Mixed event uncorrected deltaEta histogram in 3D near side for acceptance corfection

  AliJTH1D fhDphiAssocIsolTrigg ; //! comment me
  AliJTProfile fhMeanPtAssoc ; //! comment me
  AliJTProfile fhMeanZtAssoc ; //! comment me
  AliJTH1D fhPtAssocUE; //!FK// Eta Gap dependent UE
  AliJTH1D fhPtAssocUEIsolTrigg; //!FK// trigger isolated hadron
  AliJTH1D fhPtAssocN, fhPtAssocF; //! comment needed
  
  //cosThetaStar
  AliJTH1D fhCosThetaStar; //! comment me
  AliJTH2D fhCMSrap;           //! comment me
  AliJTProfile fpCMSrap; //! comment me
  AliJTH1D fhInvMass, fhPairPtMass, fhPairDPhi, fhPairDpT; //! comment me
  AliJTH1D fhPairPtDphi, fhPairPt; //! comment me
  
  AliJTH1D fhDEtaNear; //! comment me
  AliJTH1D fhDEtaNearM; //! comment me
  AliJTH1D fhDEtaNearXEbin; //! comment me
  AliJTH1D fhDEtaNearMXEbin; //! comment me
  
  AliJTH1D fhDRNearPt; //! comment me
  AliJTH1D fhDRFarPt; //! comment me
  AliJTH1D fhDRNearPtMoon; //! comment me
  AliJTH1D fhDRFarPtMoon; //! comment me
  AliJTH1D fhDRNearPtMoonM; //! comment me
  AliJTH1D fhDRFarPtMoonM ; //! comment me
  
  //AliJTH1D *hDEtaUE;
  AliJTH1D fhDEtaFar; //! comment me
  AliJTH1D fhIphiTrigg; //! comment me
  AliJTH1D fhIetaTrigg; //! comment me
  AliJTH1D fhIphiAssoc; //! comment me
  AliJTH1D fhIetaAssoc; //! comment me
  AliJTH1D fhFixPtBin; //! comment me
  AliJTH1D fhTriggPtBin; //! comment me
  AliJTH1D fhTriggPtBinIsolTrigg; //!FK// trigger isolated hadron
  AliJTH1D fhTriggMult; //! comment me
  AliJTH1D fhAssocPtBin; //! comment me
  //TH1D *hInvMassLike[2][kMaxNoCentrBin][kPtDim], *hInvMassUnLike[2][kMaxNoCentrBin][kPtDim];
  
  //==================================================
  // xe slopes - done manually for the pp purpose only
  // you have to determine 7 trigger bins
  //==================================================
  AliJTH1D     fhxEN, fhxEF,  fhxEFIsolTrigg; //! comment me
  AliJTH1D	   fhxEPtBin; //! all=0,near=1,far=2
  
  //==================================================
  //FK//mix2    inclusve spectra
  //==================================================

  AliJHistManager * fHmgInclusive;
  AliJTH1D fhIetaTriggFromFile; //!FK//mix2
  AliJTH1D fhIetaAssocFromFile  ;//!FK//mix2
  AliJTH1D fhIphiTriggFromFile  ;//!FK//mix2
  AliJTH1D fhIphiAssocFromFile  ;//!FK//mix2
  AliJTH1D fhDphiAssocMixFromFile  ;//!FK//mix2
  
  AliJTH1D fhDEtaNearMixFromFile; //!
  AliJTH1D fhDEta3DNearMixFromFile; //!
  
  //TH1D *fhDEtaNearMixFromFile[kMaxNoCentrBin][kPtDim][kPtDim];   // mixed event near side delta eta distribution
  //TH1D *fhDEta3DNearMixFromFile[kMaxNoCentrBin][kPtDim][kPtDim]; // mixed event 3D near side delta eta distribution
  
  //===================================================
  // Event/Track histograms
  //===================================================
  AliJTH1D fhLPpt; //! comment me
  AliJTH1D fhLPpairPt; //! comment me
  AliJTH1D fhChargedPt, fhChargedPtNoCorr; //! comment me
  AliJTH1D fhChargedPtJacek; //! comment me
  AliJTH1D fhChargedPtJacekPos; //! comment me
  AliJTH1D fhChargedPtJacekNeg; //! comment me
  AliJTH1D fhChargedPtJacekEta; //! 3 bins in eta
  AliJTH1D fhChargedPtFiete; //! comment me
  AliJTH1D fhVdelta2,fhVdelta3, fhVN; //! comment me
  AliJTProfile fhTrackingEfficiency; //! comment needed
  AliJTProfile fpV2, fpV3, fpVdeltaNorm; //! comment me
  AliJTH1D fhChargedEta; //! comment me
  AliJTH1D fhLPeta; //! comment me
  AliJTH1D fhAssocMult; //! comment me
  AliJTH1D fhChargedMult, fhChargedMultCut; //! comment me
  AliJTH2D fhChargedMultCent; //! comment me
  
  // XAliJT histogrmas
  AliJTH1D     fhXt;  //! comment needed
  AliJTH1D     fhXtWeighted; //! comment needed
  AliJTH1D     fhXtWeightedHT; //! HT pions
  AliJTH1D 	 fhPtForXt; //! comment me
  AliJTProfile fhConeActivity;          //! pT sum in cone, to be compared to the ALICE UE results
  AliJTProfile fhConeActivityIsolated;  //! activity for isolated triggers
  AliJTProfile fhPerpConeActivity;      //! pT sum in cone perpendicular to the leading particle
  AliJTProfile fhPerpConeActivityIsolated;  //! same as above but for isolated leading particle
  
  AliJTH1D  fhV0AMult; //! comment needed
  
  AliJTH1D fhZVertRawErr; //! comment me
  AliJTH1D fhZVert; //! comment me
  AliJTH1D fhCentr; //! comment me
  AliJTH1D fhiCentr; //! comment me
  AliJTH1D fhEventPerRun; //! comment me
  
  AliJTH1D fhIsolatedLPpt; //! comment me
  AliJTH1D fhBkgActivity; //! comment me
  
  // D.J
  //===================================================
  // Jet with LP Histograms
  //===================================================
  AliJTH1D fhDphiLPJet; //! comment me
  AliJTH1D fhDEtaLPJet; //! comment me
  AliJTH1D fhDPtLPJet; //! comment me
  AliJTH1D fhLPJetPTt; //! comment me
  AliJTH1D fhLPJetPt; //! comment me
  AliJTH1D fhLPJetEtaPTt; //! comment me
  AliJTH1D fhLPJetRapidityPTt; //! comment me
  AliJTH1D fhLPJetMassPTt; //! comment me
  
  AliJTH1D fhLeadingJetWLPPTt; //! comment me
  
  AliJTH1D fhJetPt; //! comment me
  AliJTH1D fhLeadingJetPt; //! comment me
  AliJTH1D fhLeadingJetWLPPt; //! comment me
  
  AliJTH1D fhDiJetAsym; //! dijet asymmetry in cbin
  AliJTH1D fhJetMassPTt; //! comment me
  AliJTH1D fhJetUEPt; //! comment me
  
  //===================================================
  // Jet Histograms
  //===================================================
  AliJTH1D fhJetDphi; //! comment me
  AliJTH1D fhJetDeta; //! comment me
  AliJTH1D fhJetMultPt; //! comment me
  
  //==================================================
  // Background study
  //==================================================
  AliJTH1D fhJetRho; //! comment me
  AliJTH1D fhJetRhoSigma; //! comment me
  
  AliJTH1D fhJetPartMult; //! comment me
  
  AliJTH1D fhRecoDiJetM; //! comment me
  AliJTH1D fhRecoDiJetdPhi; //! comment me
  AliJTH1D fhRecoDiJetkT; //! comment me
  
  //===================================================
  // parton 71 Histogram
  //===================================================
  AliJTH1D fhNParton71 ; //! comment me
  AliJTH1D fhNStringGroup; //! comment me
  AliJTH1D fhNStringGroupFrom; //! comment me
  AliJTH1D fhNTracksInStringGroupFrom; //! comment me
  AliJTH1D fhRapidity71From; //! comment me
  AliJTH1D fhPt71From; //! comment me
  
  //===================================================
  // PHENIX histograms
  //===================================================
  //==Run-by-Run calib ================================
  //TH1D *fTofPbSc[kMaxNoRuns]; TH1D *fTofPbGl[kMaxNoRuns]; // comment me
  //==Pt and FlipSlide Spectra=============================================
  //TH1D *fhCglPt3PC[kMaxNoCentrBin],    *fhCglPtFlip3PC[kMaxNoCentrBin]; // comment me
  //TH1D *fhCglPt32PC[kMaxNoCentrBin],   *fhCglPtFlip32PC[kMaxNoCentrBin]; // comment me
  
protected:
  double fmaxEtaRange;                       // maximum eta range
  double fmaxTriggEtaRange;                  // should be the same as above. Use for GeoAccCorr
  double ftriggFiducCut;                     // fiducial cut for the trigger part in eta. Not in use I think (Jan)
  int fnUE, fnUEfar;                         // logarithmic binning for some pT and UE histos
  double fUEBinsx[101], fUEBinsxFar[101];    // logarithmic bins for the underlaying event
  double fLowRange, fHighRange;              // lower and upper range for dphi histos
  bool fenable2DHistos;                      // enable the filling of two dimensional histograms
  
  // Manual bin definitions
  int fNJacek;        // Number of bins in Jacek binning
  double *fPttJacek;  // Bin borders in Jacek binning
  int fNEta;          // Number of bins in eta binning
  double *fEta;       // Bin borders in eta binning
  int fNJanFiete;     // Number of bins in JanFiete binning
  double *fJanFiete;  // Bin borders in JanFiete binning
  
private:
  void NormalizeAcceptanceHistos(AliJTH1D &acceptanceHisto, corrType assocType);
};

#endif






















