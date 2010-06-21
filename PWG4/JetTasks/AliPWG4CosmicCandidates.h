#ifndef AliPWG4CosmicCandidates_cxx
#define AliPWG4CosmicCandidates_cxx

// Analysis task looking for cosmic candidates
// Authors: Marta Verweij marta.verweij@cern.ch

class TH1F;
class TH2F;
class TH3F;
class TList;
class AliESDEvent;
class AliESDfriend;
class AliESDfriendTrack;
class AliMCEvent;
class AliVEvent;
class AliESDtrackCuts;
class AliESDtrack;

#include "AliAnalysisTaskSE.h"

class AliPWG4CosmicCandidates : public AliAnalysisTaskSE {
 public:
  AliPWG4CosmicCandidates();
  AliPWG4CosmicCandidates(const char *name);
  AliPWG4CosmicCandidates(const AliPWG4CosmicCandidates &res);
  AliPWG4CosmicCandidates& operator=(const AliPWG4CosmicCandidates& trclass);
  virtual ~AliPWG4CosmicCandidates() {;}

  virtual void   LocalInit();
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  
  //Setters
  void SetCuts(AliESDtrackCuts* trackCuts) {fTrackCuts = trackCuts;}
  void SetPtMin(Double_t ptmin)            {fPtMin = ptmin;}
  void SetMaxCosmicAngle(Double_t angle)   {fMaxCosmicAngle = angle;}

 private:
  AliESDtrackCuts *fTrackCuts;    // Standard trackCuts for global tracks

  Double_t fPtMin;                // Minimal pt for cosmic candidate 
  Double_t fMaxCosmicAngle;       // Max deviation from pi (angle between two tracks) in case of cosmic candidate

  TH1F *fNEventAll;                             //! Event counter
  TH1F *fNEventSel;                             //! Event counter: Selected events for analysis
  
  TH1F *fPtSignedCosmicCandidates;              //! Cosmic Candidates
  TH1F *fDeltaPtCosmicCandidates;               //! Cosmic Candidates Delta Pt
  TH2F *fDeltaPhiSumEta;                        //! Cosmic Candidates Delta Phi vs Sum Eta
  TH2F *fDCAZCosmicCandidates;                  //! Cosmic Candidates DCAZ track1 vs track2
  TH2F *fDCARCosmicCandidates;                  //! Cosmic Candidates DCAR track1 vs track2
  TH1F *fTheta;                                 //! Angle \theta between cosmic candidates in 3D space
  TH1F *fThetaZoom;                             //! Angle between cosmic candidates in 3D space zoomed into back-to-back region
  TH3F *fThetaPt1Pt2;                           //! Angle theta vs Pt1 vs Pt2
  TH3F *fThetaPt1Pt2Signed;                     //! Angle theta vs Pt1 vs Pt2
  TH3F *fDeltaPhiSumEtaPt1;                     //! Delta Phi vs Sum Eta vs Pt1
  TH3F *fDeltaPhiSumEtaPt2;                     //! Delta Phi vs Sum Eta vs Pt2
  TH3F *fThetaDCAZ1DCAZ2;                       //! Angle theta vs DCAZ1 vs DCAZ2
  TH1F *fRisol;                                 //! Isolation R
  TH2F *fRisolTheta;                            //! Isolation R vs Theta

  TList *fHistListCosmics;                      //! List of Histograms for cosmic candidates  

  
  ClassDef(AliPWG4CosmicCandidates, 1);
};

#endif
