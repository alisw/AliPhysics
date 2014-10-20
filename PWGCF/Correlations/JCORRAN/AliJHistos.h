/* Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice */

// Short comment describing what this class does needed!

//===========================================================
// AliJHistos.h
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
#include "AliJHistManager.h"

class AliJCard;
class AliJBaseTrack;
class AliJPhoton;
class AliJTrack;

using namespace std;

inline ostream &operator << (ostream &out_file, const TLorentzVector &Vec)
{
  out_file<<"Px="<<Vec.Px()<<" Py="<<Vec.Py()<<" Pz="<<Vec.Pz()<<" E="<<Vec.E()<<" M="<<Vec.M()<<endl;
  out_file<<"Theta="<<Vec.Theta()<<" Phi="<<Vec.Phi()<<" p="<<Vec.Rho()<<endl;
  return(out_file);
}

class AliJHistos {
  
public:
  AliJHistos(AliJCard* cardP); //constructor
  virtual ~AliJHistos(){;}    //destructor
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
  
public:
  AliJCard  *fCard;       // card
  
  AliJHistManager * fHMG; // Histogram manager
  AliJBin fCentBin;       // Bin of Centrality
  AliJBin fVtxBin;        // Bin of Z-vertex Bin
  AliJBin fPTtBin;        // Bin of pT Trigged
  AliJBin fPTaBin;        // Bin of pT Associated
  AliJBin fXEBin;         // Bin of XE
  AliJBin fKLongBin;      // Bin of Klong
  AliJBin fRGapBin;       // Bin of R-gap
  AliJBin fEtaGapBin;     // Bin of Eta-gap
  AliJBin fPhiGapBin;     // Bin of Phi-gap
  AliJBin fMassBin;       // Bin of Mass
  AliJBin fTypBin;        // Bin of type ( data, mixed ), but is being used for any 2 dimension
  AliJBin fTypBin3;       // Bin of type3, is being used for any 3 dimension
  AliJBin fPairPtBin;     // Bin of pT pair
  
  AliJTProfile fhMixStat; // comment me
  
  //==Pt stat. fcorrelations ===============================================
  AliJTH1D fTestHist; // test only, remove it not needed
  AliJTH1D fhPtNear; // comment me
  AliJTH1D fhPtFar ; // comment me
  
  //assorted
  AliJTH1D fhPhi; // comment me
  AliJTH1D fhDphiAssoc; // comment me
  AliJTH1D fhDphiAssocXEbin; // comment me
                               //2D
  AliJTH2D fhDphiAssoc2DIAA; // comment me
  AliJTH2D fhDphiAssoc2D; // comment me
  
  // JV // 2D dphi deta histos to study better the background
  AliJTH2D fhDphiDetaKlong;  // 2D histogram from deltaPhi-deltaEta plane in klong and eta gap bins
  AliJTH2D fhDphiDetaKlongR; // 2D histogram from deltaPhi-deltaEta plane in klong and R gap bins
  AliJTH2D fhDphiDetaXe;     // 2D histogram from deltaPhi-deltaEta plane in xlong and eta gap bins
  AliJTH2D fhDphiDetaXeR;    // 2D histogram from deltaPhi-deltaEta plane in xlong and R gap bins
  AliJTH2D fhDphiDetaPta;    // 2D histogram from deltaPhi-deltaEta plane in pta and eta gap bins
  AliJTH2D fhDphiDetaPtaR;   // 2D histogram from deltaPhi-deltaEta plane in pta and R gap bins
  
  // JV // pTa distributions in background bins
  AliJTH1D fhBgAssocKlong;   // background pta distribution in klong and eta gap bins
  AliJTH1D fhBgAssocKlongR;  // background pta distribution in klong and R gap bins
  AliJTH1D fhBgAssocXe;      // background pta distribution in xlong and eta gap bins
  AliJTH1D fhBgAssocXeR;     // background pta distribution in xlong and R gap bins
  AliJTH1D fhBgAssocPta;     // background pta distribution in pta and eta gap bins
  AliJTH1D fhBgAssocPtaR;    // background pta distribution in pta and R gap bins
  
  AliJTH1D fhDphiAssocIsolTrigg ; // comment me
  AliJTProfile fhMeanPtAssoc ; // comment me
  AliJTProfile fhMeanZtAssoc ; // comment me
  AliJTH1D fhPtAssocUE; //FK// Eta Gap dependent UE
  AliJTH1D fhPtaEtaGapN; //FK// pta spectrum with Eta Gap
  AliJTH1D fhPtaRGapN; //FK// pta spectrum with Eta Gap
  AliJTH1D fhPtAssocUEIsolTrigg; //FK// trigger isolated hadron
  AliJTH1D fhPtAssocN, fhPtAssocF; // comment needed
  
  //cosThetaStar
  AliJTH1D fhCosThetaStar; // comment me
  AliJTH2D fhCMSrap;           // comment me
  AliJTProfile fpCMSrap; // comment me
  AliJTH1D fhInvMass, fhPairPtMass, fhPairDPhi, fhPairDpT; // comment me
  AliJTH1D fhPairPtDphi, fhPairPt; // comment me
  
  AliJTH1D fhDEtaNear; // comment me
  AliJTH1D fhDEtaNearM; // comment me
  AliJTH1D fhDEtaNearXEbin; // comment me
  AliJTH1D fhDEtaNearMXEbin; // comment me
  
  AliJTH1D fhDRNearPt; // comment me
  AliJTH1D fhDRFarPt; // comment me
  AliJTH1D fhDRNearPtMoon; // comment me
  AliJTH1D fhDRFarPtMoon; // comment me
  AliJTH1D fhDRNearPtMoonM; // comment me
  AliJTH1D fhDRFarPtMoonM ; // comment me
  
  //AliJTH1D *hDEtaUE;
  AliJTH1D fhDEtaFar; // comment me
  AliJTH1D fhIphiTrigg; // comment me
  AliJTH1D fhIetaTrigg; // comment me
  AliJTH1D fhIphiAssoc; // comment me
  AliJTH1D fhIetaAssoc; // comment me
  AliJTH1D fhFixPtBin; // comment me
  AliJTH1D fhTriggPtBin; // comment me
  AliJTH1D fhTriggPtBinIsolTrigg; //FK// trigger isolated hadron
  AliJTH1D fhTriggMult; // comment me
  AliJTH1D fhAssocPtBin; // comment me
  //TH1D *hInvMassLike[2][kMaxNoCentrBin][kPtDim], *hInvMassUnLike[2][kMaxNoCentrBin][kPtDim];
  
  //==================================================
  // xe slopes - done manually for the pp purpose only
  // you have to determine 7 trigger bins
  //==================================================
  AliJTH1D     fhxEN, fhxEF,  fhxEFIsolTrigg; // comment me
  AliJTH1D     fhPoutF; // comment me
  AliJTH1D	   fhxEPtBin; // all=0,near=1,far=2
                        // xe bins
  AliJTH1D     fhJT;  // comment me
  AliJTH1D     fhJTBg;  // comment me
  AliJTH1D     fhJTBgR;  // comment me
                         // klong bins
  AliJTH1D     fhJTKlong;  // comment me
  AliJTH1D     fhJTKlongBg;  // comment me
  AliJTH1D     fhJTKlongBgR;  // comment me
                              // pta bins
  AliJTH1D     fhJTPta;  // comment me
  AliJTH1D     fhJTPtaBg;  // comment me
  AliJTH1D     fhJTPtaBgR;  // comment me
  
  //FK//mix2    inclusve spectra
  TH1D *fhIetaTriggFromFile  [kMaxNoCentrBin][kPtDim];//FK//mix2
  TH1D *fhIetaAssocFromFile  [kMaxNoCentrBin][kPtDim];//FK//mix2
  TH1D *fhIphiTriggFromFile  [kMaxNoCentrBin][kPtDim];//FK//mix2
  TH1D *fhIphiAssocFromFile  [kMaxNoCentrBin][kPtDim];//FK//mix2
  TH1D *fhDphiAssocMixFromFile  [kMaxNoCentrBin][kPtDim][kPtDim];//FK//mix2
  
  TH1D *fhDEtaNearMixFromFile[kMaxNoCentrBin][kPtDim][kPtDim]; // comment me
  
  //===================================================
  // Event/Track histograms
  //===================================================
  AliJTH1D fhLPpt; // comment me
  AliJTH1D fhLPpairPt; // comment me
  AliJTH1D fhChargedPt, fhChargedPtNoCorr; // comment me
  AliJTH1D fhChargedPtJacek; // comment me
  AliJTH1D fhChargedPtJacekEta; // 3 bins in eta
  AliJTH1D fhChargedPtFiete; // comment me
  AliJTH1D fhVdelta2,fhVdelta3, fhVN; // comment me
  AliJTProfile fhTrackingEfficiency; // comment needed
  AliJTProfile fpV2, fpV3, fpVdeltaNorm; // comment me
  AliJTH1D fhChargedEta; // comment me
  AliJTH1D fhLPeta; // comment me
  AliJTH1D fhAssocMult; // comment me
  AliJTH1D fhChargedMult, fhChargedMultCut; // comment me
  AliJTH2D fhChargedMultCent; // comment me
  
  // XAliJT histogrmas
  AliJTH1D     fhXt;  // comment needed
  AliJTH1D     fhXtWeighted; // comment needed
  AliJTH1D     fhXtWeightedHT; // HT pions
  AliJTH1D 	 fhPtForXt; // comment me
  AliJTProfile fhConeActivity;          // pT sum in cone, to be compared to the ALICE UE results
  AliJTProfile fhConeActivityIsolated;  // activity for isolated triggers
  AliJTProfile fhPerpConeActivity;      // pT sum in cone perpendicular to the leading particle
  AliJTProfile fhPerpConeActivityIsolated;  // same as above but for isolated leading particle
  
  AliJTH1D  fhV0AMult; // comment needed
  
  AliJTH1D fhZVertRaw; // comment me
  AliJTH1D fhZVertRawErr; // comment me
  AliJTH1D fhZVert; // comment me
  AliJTH1D fhCentr; // comment me
  AliJTH1D fhiCentr; // comment me
  AliJTH1D fhEventPerRun; // comment me
  AliJTH2D fhVertexZTriggVtx; // comment me
  
  AliJTH1D fhIsolatedLPpt; // comment me
  AliJTH1D fhBkgActivity; // comment me
  
  // D.J
  //===================================================
  // Jet with LP Histograms
  //===================================================
  AliJTH1D fhDphiLPJet; // comment me
  AliJTH1D fhDEtaLPJet; // comment me
  AliJTH1D fhDPtLPJet; // comment me
  AliJTH1D fhLPJetPTt; // comment me
  AliJTH1D fhLPJetPt; // comment me
  AliJTH1D fhLPJetEtaPTt; // comment me
  AliJTH1D fhLPJetRapidityPTt; // comment me
  AliJTH1D fhLPJetMassPTt; // comment me
  
  AliJTH1D fhLeadingJetWLPPTt; // comment me
  
  AliJTH1D fhJetPt; // comment me
  AliJTH1D fhLeadingJetPt; // comment me
  AliJTH1D fhLeadingJetWLPPt; // comment me
  
  AliJTH1D fhJetAssymPTt; // comment me
  AliJTH1D fhJetMassPTt; // comment me
  AliJTH1D fhJetUEPt; // comment me
  
  //===================================================
  // Jet Histograms
  //===================================================
  AliJTH1D fhJetDphi; // comment me
  AliJTH1D fhJetDeta; // comment me
  AliJTH1D fhJetMultPt; // comment me
  
  //==================================================
  // Background study
  //==================================================
  AliJTH1D fhJetRho; // comment me
  AliJTH1D fhJetRhoSigma; // comment me
  
  AliJTH1D fhJetPartMult; // comment me
  
  AliJTH1D fhRecoDiJetM; // comment me
  AliJTH1D fhRecoDiJetdPhi; // comment me
  AliJTH1D fhRecoDiJetkT; // comment me
  
  //===================================================
  // parton 71 Histogram
  //===================================================
  AliJTH1D fhNParton71 ; // comment me
  AliJTH1D fhNStringGroup; // comment me
  AliJTH1D fhNStringGroupFrom; // comment me
  AliJTH1D fhNTracksInStringGroupFrom; // comment me
  AliJTH1D fhRapidity71From; // comment me
  AliJTH1D fhPt71From; // comment me
  
  //===================================================
  // PHENIX histograms
  //===================================================
  //==Run-by-Run calib ================================
  TH1D *fTofPbSc[kMaxNoRuns]; TH1D *fTofPbGl[kMaxNoRuns]; // comment me
  //==Pt and FlipSlide Spectra=============================================
  TH1D *fhCglPt3PC[kMaxNoCentrBin],    *fhCglPtFlip3PC[kMaxNoCentrBin]; // comment me
  TH1D *fhCglPt32PC[kMaxNoCentrBin],   *fhCglPtFlip32PC[kMaxNoCentrBin]; // comment me


  AliJTH1D fhTrackSelection; // checking bit convention
  
  // Manual bin definitions
  int fNJacek;        // Number of bins in Jacek binning
  double *fPttJacek;  // Bin borders in Jacek binning
  int fNEta;          // Number of bins in eta binning
  double *fEta;       // Bin borders in eta binning
  int fNJanFiete;     // Number of bins in JanFiete binning
  double *fJanFiete;  // Bin borders in JanFiete binning
  
  // additional histos
  AliJTH1D fhEvents; // comment needed
  AliJTH1D fhEventTrigger; // comment needed
  
protected:
  double fmaxEtaRange;                       // maximum eta range
  double fmaxTriggEtaRange;                  // should be the same as above. Use for GeoAccCorr
  double ftriggFiducCut;                     // fiducial cut for the trigger part in eta. Not in use I think (Jan)
  int fnUE, fnUEfar;                           // logarithmic binning for some pT and UE histos
  double fUEBinsx[101], fUEBinsxFar[101];      // logarithmic bins for the underlaying event
  double fLowRange, fHighRange;                // lower and upper range for dphi histos
};

#endif






















