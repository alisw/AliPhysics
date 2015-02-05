#ifndef ALIJETFASTSIMULATION_H
#define ALIJETFASTSIMULATION_H

// $Id$

class TClonesArray;
class TRandom3;
class AliVParticle;
class AliPicoTrack;

#include "AliAnalysisTaskEmcal.h"

class AliJetFastSimulation : public AliAnalysisTaskEmcal {
 public:
  AliJetFastSimulation();
  AliJetFastSimulation(const char *name); 
  virtual ~AliJetFastSimulation();

  virtual void           LocalInit();
  virtual void           UserCreateOutputObjects();

  void                   SetTracksOutName(const char *n)          { fTracksOutName   = n;    }
  void                   SetNTrackClasses(Int_t i)                { fNTrackClasses   = i;    }
  void                   SetFixedTrackEfficiency(Double_t eff)    { fEfficiencyFixed = eff ; }

  void                   SetUseTrResolutionFromOADB(Bool_t b=kTRUE, TString path="$ALICE_PHYSICS/OADB/PWGJE/Resolution/PtResol_LHCh_Cent0-10_v1.root") {fUseTrPtResolutionFromOADB = b; fPathTrPtResolution=path;}
  void                   SetUseTrEfficiencyFromOADB(Bool_t b=kTRUE, TString path="$ALICE_PHYSICS/OADB/PWGJE/Efficiency/Efficiency_LHC11a2aj_Cent0_v1.root") {fUseTrEfficiencyFromOADB = b; fPathTrEfficiency=path;}
  void                   SetSmearResolution(Bool_t b)                               { fUseTrPtResolutionSmearing = b ;}
  void                   SetDiceEfficiency(Int_t b)                                 { fUseDiceEfficiency         = b ;}
  void                   SetDiceEfficiencyMinPt(Double_t pt)                        { fDiceEfficiencyMinPt       = pt;}
  void                   SetUncertEfficiency(Double_t uncerteff)                   { fUncertEfficiency           =uncerteff;}      
 protected:
  void                   ExecOnce();
  Bool_t                 Run();

  void                   SimulateTracks();
  Bool_t                 DiceEfficiency(AliPicoTrack *vp, Double_t eff[3], Double_t rnd);
  AliPicoTrack          *SmearPt(AliPicoTrack *vp, Double_t eff[3], Double_t rnd);
  Double_t               GetMomentumSmearing(Int_t cat, Double_t pt);
  void                   FitMomentumResolution();
  void                   LoadTrEfficiencyRootFileFromOADB();
  void                   LoadTrPtResolutionRootFileFromOADB();
  void                   SetMomentumResolutionHybrid(TProfile *p1, TProfile *p2, TProfile *p3);
  void                   SetEfficiencyHybrid(TH1 *h1, TH1 *h2, TH1 *h3);

  TString                fTracksOutName;       // name of output track collection
  TClonesArray          *fTracksOut;           //!output track collection
  Int_t                  fNTrackClasses;       // number of track classes
  TRandom3 *fRandom;                           //! random number generator
  Double_t  fEfficiencyFixed;                  // fixed efficiency for all pT and all types of tracks
  TProfile *fMomResH1;                         // Momentum resolution from TrackQA Hybrid Category 1
  TProfile *fMomResH2;                         // Momentum resolution from TrackQA Hybrid Category 2
  TProfile *fMomResH3;                         // Momentum resolution from TrackQA Hybrid Category 3
  TF1      *fMomResH1Fit;                      // fit to momentum resolution
  TF1      *fMomResH2Fit;                      // fit to momentum resolution
  TF1      *fMomResH3Fit;                      // fit to momentum resolution
  TH1      *fhEffH1;                           // Efficiency for Spectra Hybrid Category 1
  TH1      *fhEffH2;                           // Efficiency for Spectra Hybrid Category 2
  TH1      *fhEffH3;                           // Efficiency for Spectra Hybrid Category 3
  Bool_t    fUseTrPtResolutionSmearing;        // Apply momentum smearing on track level
  Int_t     fUseDiceEfficiency;                // Flag to apply efficiency on track level by dicing 0: no dicing; 1: dicing wrt to input branch;
  Double_t  fDiceEfficiencyMinPt;              // Only do efficiency dicing for tracks above this pt
  Double_t  fUncertEfficiency;                  //tracking efficiency uncertainty, usually +-0.4%
  Bool_t    fUseTrPtResolutionFromOADB;        // Load track pt resolution root file from OADB path
  Bool_t    fUseTrEfficiencyFromOADB;          // Load tracking efficiency root file from OADB path
  TString   fPathTrPtResolution;               // OADB path to root file
  TString   fPathTrEfficiency;                 // OADB path to root file

  //Output objects
  TH1F     *fHistPtDet;                        //!pT spectrum of detector level particles
  TH2F     *fh2PtGenPtSmeared;                 //! Control histo smeared momentum
  TProfile *fp1Efficiency;                     //! Control profile efficiency
  TProfile *fp1PtResolution;                   //! Control profile for pT resolution

  
 private:
  AliJetFastSimulation(const AliJetFastSimulation&);            // not implemented
  AliJetFastSimulation &operator=(const AliJetFastSimulation&); // not implemented

  ClassDef(AliJetFastSimulation, 1) // Jet fast simulation task
};
#endif
