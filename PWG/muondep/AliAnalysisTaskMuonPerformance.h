#ifndef ALIANALYSISTASKMUONPERFORMANCE_H
#define ALIANALYSISTASKMUONPERFORMANCE_H

/// \ingroup "PWG3muon"
/// \class AliAnalysisTaskMuonPerformance
/// \brief Analysis task for 
///
//  Author D.Stocco and P.Pillot, Subatech, Nantes


#include "AliAnalysisTaskSE.h"
#include "AliLog.h"

class TObjArray;
class TH1;
class TH2;
class TGraphAsymmErrors;
class TCanvas;
class AliMUONRecoParam;
class AliCFContainer;
class AliMCParticle;
class AliESDMuonTrack;
class AliCFEffGrid;

class AliAnalysisTaskMuonPerformance : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskMuonPerformance();
  AliAnalysisTaskMuonPerformance(const char *name);
  virtual ~AliAnalysisTaskMuonPerformance();
  
  /// Set location of the default OCDB storage (if not set use "raw://")
  void SetDefaultStorage(const char* ocdbPath) { fDefaultStorage = ocdbPath; }
  
  // Set the binning to be used to study the detector resolution versus momentum
  void SetPBins(Int_t nBins, Double_t min, Double_t max);
  
  /// set the flag to add or not the systematic shifts of the residuals to the resolution
  void CorrectClusterResForSystematics(Bool_t flag = kTRUE) { fCorrectForSystematics = flag; }
  
  /// set the flag to fit or not the cluster residuals to extract means and sigmas
  void FitClusterResiduals(Bool_t flag = kTRUE) { fFitResiduals = flag; }
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *option);
  virtual void   NotifyRun();
  
  enum {
    kVarPt,         ///< Transverse momentum
    kVarEta,        ///< Pseudo-rapidity
    kVarPhi,        ///< Azimuthal angle
    kVarThetaZones, ///< Theta at absorber end (4 zones)
    kVarCharge,     ///< Particle charge
    kVarHasTracker, ///< Is tracker track
    kVarTrigger,    ///< Trigger info
    kVarMotherType, ///< Mother type
    kVarMatchMC,    ///< MC matching flag
    kNvars          ///< THnSparse dimensions
  };
  
  enum {
    kStepReconstructed, ///< Reconstructed tracks
    kStepGeneratedMC,   ///< Generated tracks (MC)
    kNsteps             ///< Number of steps
  };
  
 private:
  
  AliAnalysisTaskMuonPerformance(const AliAnalysisTaskMuonPerformance&);
  AliAnalysisTaskMuonPerformance& operator=(const AliAnalysisTaskMuonPerformance&);
  
  Bool_t  GetEfficiency(AliCFEffGrid* efficiency, Double_t& calcEff, Double_t& calcEffErr);
  Int_t   RecoTrackMother(AliMCParticle* mcParticle);
  Float_t GetBinThetaAbsEnd(Float_t RAtAbsEnd, Bool_t isTheta = kFALSE);
  void    FillContainerInfo(Double_t* containerInput, AliESDMuonTrack* esdTrack, Int_t mcID);
  
  void    FitLandauGausResVsP(TH2* h, const char* fitting, TGraphAsymmErrors* gMean, TGraphAsymmErrors* gMostProb, TGraphAsymmErrors* gSigma);
  void    FitGausResVsMom(TH2* h, const Double_t mean0, const Double_t sigma0, const char* fitting, TGraphAsymmErrors* gMean, TGraphAsymmErrors* gSigma);
  void    FitPDCAVsMom(TH2* h, const char* fitting, TGraphAsymmErrors* gMean, TGraphAsymmErrors* gSigma);
  void    FitClusterResidual(TH1* h, Int_t i, Double_t& sigma, TGraphErrors* gMean, TGraphErrors* gSigma);
  
  TCanvas* DrawVsAng(const char* name, const char* title, TH1* h1, TH2* h2);
  TCanvas* DrawVsPos(const char* name, const char* title, TH2* h1, TH2* h2, TH2* h3);
  TCanvas* DrawFitLandauGausResPVsP(const char* name, const char* title, TH2* h, const Int_t nBins, const char* fitting);
  TCanvas* DrawResPVsP(const char* name, const char* title, TH2* h, const Int_t nBins);
  
  void Zoom(TH1* h, Double_t fractionCut = 0.01);
  
  enum {
    kNoMatchTrig,  ///< No match with trigger
    kAllPtTrig,    ///< Match All Pt
    kLowPtTrig,    ///< Match Low Pt
    kHighPtTrig,   ///< Match High Pt
    kNtrigCuts     ///< Total number of trigger types
  };
  
  enum {
    kCharmMu,       ///< Mu from charm
    kBeautyMu,      ///< Mu from beauty
    kPrimaryMu,     ///< Primary mu
    kSecondaryMu,   ///< Secondary mu
    kRecoHadron,    ///< Reconstructed hadron
    kUnknownPart,   ///< Particle that fails matching kine
    kNtrackSources  ///< Total number of track sources
  };
  
  enum {
    kNoMatch,     ///< not matched with either reconstructible track or triggerable track
    kTrackerOnly, ///< matched with reconstructible track only
    kMatchedSame, ///< matched with reconstructible track and triggerable track of same ID
    kMatchedDiff, ///< matched with reconstructible track and triggerable track of different ID
    kTriggerOnly, ///< matched with triggerable track only
    kNMatchMC     ///< Total number of MC matching flags
  };
  
  // Histograms for trigger resolution
  enum {
    kResTrigX11,                          ///< Residual of x position in first trigger chamber
    kResTrigY11,                          ///< Residual of y position in first trigger chamber
    kResTrigSlopeY                        ///< Residual of trigger track slope
  };
  
  // Histograms for tracker resolution
  enum {
    kResPAtVtx,                           ///< momentum residual at vertex
    kResPAtVtxVsP,                        ///< momentum residual at vertex versus P
    kResPAtVtxVsPIn23deg,                 ///< momentum residual at vertex versus P for tracks in ]2,3] degrees at absorber end
    kResPAtVtxVsPIn310deg,                ///< momentum residual at vertex versus P for tracks in ]3,10[ degrees at absorber end
    kResPAtVtxVsPIn02degMC,               ///< momentum residual at vertex versus P for tracks with MC angle < 2 degrees
    kResPAtVtxVsPosAbsEndIn02degMC,       ///< momentum residual at vertex versus position at absorber end for tracks with MC angle <= 2 degrees
    kResPAtVtxVsPosAbsEndIn23degMC,       ///< momentum residual at vertex versus position at absorber end for tracks with MC angle in ]2,3] degrees
    kResPAtVtxVsPosAbsEndIn310degMC,      ///< momentum residual at vertex versus position at absorber end for tracks with MC angle in ]3,10[ degrees
    kResPAtVtxVsAngleAtAbsEnd,            ///< momentum residual at vertex versus angle at absorber end
    kResPAtVtxVsMCAngle,                  ///< momentum residual at vertex versus MC angle
    kResPAtVtxVsAngleAtAbsEndVsP,         ///< momentum residual at vertex versus angle at absorber end versus momentum
    kResPtAtVtxVsPt,                      ///< transverse momentum residual at vertex versus pT
    
    kResPAt1stCl,                         ///< momentum residual at first cluster
    kResPAt1stClVsP,                      ///< momentum residual at first cluster versus P
    kResPtAt1stClVsPt,                    ///< transverse momentum residual at first cluster versus pT
    
    kResSlopeXAtVtx,                      ///< slope-X residual at vertex
    kResSlopeYAtVtx,                      ///< slope-Y residual at vertex
    kResSlopeXAtVtxVsP,                   ///< slope-X residual at vertex versus P
    kResSlopeYAtVtxVsP,                   ///< slope-Y residual at vertex versus P
    kResSlopeXAtVtxVsPosAbsEndIn02degMC,  ///< slope-X residual at vertex versus position at absorber end for tracks with MC angle <= 2 degrees
    kResSlopeYAtVtxVsPosAbsEndIn02degMC,  ///< slope-Y residual at vertex versus position at absorber end for tracks with MC angle <= 2 degrees
    kResSlopeXAtVtxVsPosAbsEndIn23degMC,  ///< slope-X residual at vertex versus position at absorber end for tracks with MC angle in ]2,3] degrees
    kResSlopeYAtVtxVsPosAbsEndIn23degMC,  ///< slope-Y residual at vertex versus position at absorber end for tracks with MC angle in ]2,3] degrees
    kResSlopeXAtVtxVsPosAbsEndIn310degMC, ///< slope-X residual at vertex versus position at absorber end for tracks with MC angle in ]3,10[ degrees
    kResSlopeYAtVtxVsPosAbsEndIn310degMC, ///< slope-Y residual at vertex versus position at absorber end for tracks with MC angle in ]3,10[ degrees
    kResSlopeXAtVtxVsAngleAtAbsEnd,       ///< slope-X residual at vertex versus angle at absorber end
    kResSlopeYAtVtxVsAngleAtAbsEnd,       ///< slope-Y residual at vertex versus angle at absorber end
    kResSlopeXAtVtxVsMCAngle,             ///< slope-X residual at vertex versus MC angle
    kResSlopeYAtVtxVsMCAngle,             ///< slope-Y residual at vertex versus MC angle
    
    kResSlopeXAt1stCl,                    ///< slope-X residual at first cluster
    kResSlopeYAt1stCl,                    ///< slope-Y residual at first cluster
    kResSlopeXAt1stClVsP,                 ///< slope-X residual at first cluster versus P
    kResSlopeYAt1stClVsP,                 ///< slope-Y residual at first cluster versus P
    
    kResEtaAtVtx,                         ///< eta residual at vertex
    kResEtaAtVtxVsP,                      ///< eta residual at vertex versus P
    kResEtaAtVtxVsPosAbsEndIn02degMC,     ///< eta residual at vertex versus position at absorber end for tracks with MC angle <= 2 degrees
    kResEtaAtVtxVsPosAbsEndIn23degMC,     ///< eta residual at vertex versus position at absorber end for tracks with MC angle in ]2,3] degrees
    kResEtaAtVtxVsPosAbsEndIn310degMC,    ///< eta residual at vertex versus position at absorber end for tracks with MC angle in ]3,10[ degrees
    kResEtaAtVtxVsAngleAtAbsEnd,          ///< eta residual at vertex versus angle at absorber end
    kResEtaAtVtxVsMCAngle,                ///< eta residual at vertex versus MC angle
    
    kResPhiAtVtx,                         ///< phi residual at vertex
    kResPhiAtVtxVsP,                      ///< phi residual at vertex versus P
    kResPhiAtVtxVsPosAbsEndIn02degMC,     ///< phi residual at vertex versus position at absorber end for tracks with MC angle <= 2 degrees
    kResPhiAtVtxVsPosAbsEndIn23degMC,     ///< phi residual at vertex versus position at absorber end for tracks with MC angle in ]2,3] degrees
    kResPhiAtVtxVsPosAbsEndIn310degMC,    ///< phi residual at vertex versus position at absorber end for tracks with MC angle in ]3,10[ degrees
    kResPhiAtVtxVsAngleAtAbsEnd,          ///< phi residual at vertex versus angle at absorber end
    kResPhiAtVtxVsMCAngle,                ///< phi residual at vertex versus MC angle
    
    kPDCA,                                ///< P*DCA distribution
    kPDCAVsPIn23deg,                      ///< P * DCA distribution versus P for tracks in ]2,3] degrees at absorber end
    kPDCAVsPIn310deg,                     ///< P * DCA distribution versus P for tracks in ]3,10[ degrees at absorber end
    kPDCAVsPosAbsEndIn02degMC,            ///< P * DCA distribution versus position at absorber end for tracks with MC angle <= 2 degrees
    kPDCAVsPosAbsEndIn23degMC,            ///< P * DCA distribution versus position at absorber end for tracks with MC angle in ]2,3] degrees
    kPDCAVsPosAbsEndIn310degMC,           ///< P * DCA distribution versus position at absorber end for tracks with MC angle in ]3,10[ degrees
    kPDCAVsAngleAtAbsEnd,                 ///< P * DCA distribution versus angle at absorber end
    kPDCAVsMCAngle,                       ///< P * DCA distribution versus MC angle
    
    kPMCSAngVsPIn23deg,                   ///< P * MCS deviation angle distribution versus P for tracks in ]2,3] degrees at absorber end
    kPMCSAngVsPIn310deg,                  ///< P * MCS deviation angle distribution versus P for tracks in ]3,10[ degrees at absorber end
    
    kResClXVsCh,                          ///< cluster residual-X versus chamber
    kResClYVsCh,                          ///< cluster residual-Y versus chamber
    kResClXVsDE,                          ///< cluster residual-X versus DE
    kResClYVsDE                           ///< cluster residual-Y versus DE
  };
  
  // Graphs and canvases for momentum resolution at vertex
  enum {
    kMeanResPAtVtxVsP,                    ///< mean momentum residual at vertex versus P
    kMostProbResPAtVtxVsP,                ///< most probable momentum residual at vertex versus P
    kSigmaResPAtVtxVsP,                   ///< relative momentum resolution at vertex versus P
    kcResPAtVtx,                          ///< momentum residual at vertex in 3 angular regions
    kcResPAtVtxMC,                        ///< momentum residual at vertex in 3 MC angular regions
    kcResPAtVtxVsPosAbsEndMC,             ///< momentum residual at vertex versus position at absorber end in 3 MC angular regions
    kcResPAtVtxVsPIn23deg,                ///< momentum residual for tracks between 2 and 3 degrees
    kcResPAtVtxVsPIn310deg,               ///< momentum residual for tracks between 3 and 10 degrees
    kcResPAtVtxVsPIn02degMC               ///< momentum residuals for tracks with MC angle < 2 degrees
  };
  
  // Graphs for momentum resolution at first cluster
  enum {
    kMeanResPAt1stClVsP,                  ///< mean momentum residual at first cluster versus P
    kSigmaResPAt1stClVsP                  ///< relative momentum resolution at first cluster versus P
  };
  
  // Graphs and canvases for slope resolution at vertex
  enum {
    kMeanResSlopeXAtVtxVsP,               ///< mean slope-X residual at vertex versus P
    kMeanResSlopeYAtVtxVsP,               ///< mean slope-Y residual at vertex versus P
    kSigmaResSlopeXAtVtxVsP,              ///< slope-X resolution at vertex versus P
    kSigmaResSlopeYAtVtxVsP,              ///< slope-Y resolution at vertex versus P
    kcResSlopeXAtVtx,                     ///< slope-X residual at vertex in 3 angular regions
    kcResSlopeYAtVtx,                     ///< slope-Y residual at vertex in 3 angular regions
    kcResSlopeXAtVtxMC,                   ///< slope-X residual at vertex in 3 MC angular regions
    kcResSlopeYAtVtxMC,                   ///< slope-Y residual at vertex in 3 MC angular regions
    kcResSlopeXAtVtxVsPosAbsEndMC,        ///< slope-X residual at vertex versus position at absorber end in 3 MC angular regions
    kcResSlopeYAtVtxVsPosAbsEndMC         ///< slope-Y residual at vertex versus position at absorber end in 3 MC angular regions
  };
  
  // Graphs for slope resolution at first cluster
  enum {
    kMeanResSlopeXAt1stClVsP,             ///< mean slope-X residual at first cluster versus P
    kMeanResSlopeYAt1stClVsP,             ///< mean slope-Y residual at first cluster versus P
    kSigmaResSlopeXAt1stClVsP,            ///< slope-X resolution at first cluster versus P
    kSigmaResSlopeYAt1stClVsP             ///< slope-Y resolution at first cluster versus P
  };
  
  // Graphs and canvases for eta resolution at vertex
  enum {
    kMeanResEtaAtVtxVsP,                  ///< mean eta residual at vertex versus P
    kSigmaResEtaAtVtxVsP,                 ///< eta resolution at vertex
    kcResEtaAtVtx,                        ///< eta residual at vertex in 3 angular regions
    kcResEtaAtVtxMC,                      ///< eta residual at vertex in 3 MC angular regions
    kcResEtaAtVtxVsPosAbsEndMC            ///< eta residual at vertex versus position at absorber end in 3 MC angular regions
  };
  
  // Graphs and canvases for phi resolution at vertex
  enum {
    kMeanResPhiAtVtxVsP,                  ///< mean phi residual at vertex versus P
    kSigmaResPhiAtVtxVsP,                 ///< phi resolution at vertex
    kcResPhiAtVtx,                        ///< phi residual at vertex in 3 angular regions
    kcResPhiAtVtxMC,                      ///< phi residual at vertex in 3 MC angular regions
    kcResPhiAtVtxVsPosAbsEndMC            ///< phi residual at vertex versus position at absorber end in 3 MC angular regions
  };
  
  // Graphs and cavases for DCA resolution and MCS dispersion
  enum {
    kMeanPDCAVsPIn23deg,                  ///< mean P * DCA versus P for tracks in ]2,3] degrees at absorber end
    kSigmaPDCAVsPIn23deg,                 ///< P * DCA resolution versus P for tracks in ]2,3] degrees at absorber end
    kMeanPDCAVsPIn310deg,                 ///< mean P * DCA versus P for tracks in ]3,10[ degrees at absorber end
    kSigmaPDCAVsPIn310deg,                ///< P * DCA resolution versus P for tracks in ]3,10[ degrees at absorber end
    kMeanPMCSAngVsPIn23deg,               ///< mean P * MCS deviation angle versus P for tracks in ]2,3] degrees at absorber end
    kSigmaPMCSAngVsPIn23deg,              ///< P * MCS deviation angle dispersion versus P for tracks in ]2,3] degrees at absorber end
    kMeanPMCSAngVsPIn310deg,              ///< mean P * MCS deviation angle versus P for tracks in ]3,10[ degrees at absorber end
    kSigmaPMCSAngVsPIn310deg,             ///< P * MCS deviation angle dispersion versus P for tracks in ]3,10[ degrees at absorber end
    kcPDCA,                               ///< P * DCA in 3 angular regions
    kcPDCAMC,                             ///< P * DCA in 3 MC angular regions
    kcPDCAVsPosAbsEndMC                   ///< P * DCA versus position at absorber end in 3 MC angular regions
  };
  
  // Graphs for cluster resolution
  enum {
    kMeanResClXVsCh,                      ///< mean cluster residual-X versus chamber
    kMeanResClYVsCh,                      ///< mean cluster residual-Y versus chamber
    kSigmaResClXVsCh,                     ///< cluster resolution-X versus chamber
    kSigmaResClYVsCh,                     ///< cluster resolution-Y versus chamber
    kMeanResClXVsDE,                      ///< mean cluster residual-X versus DE
    kMeanResClYVsDE,                      ///< mean cluster residual-Y versus DE
    kSigmaResClXVsDE,                     ///< cluster resolution-X versus DE
    kSigmaResClYVsDE                      ///< cluster resolution-Y versus DE
  };
  
  TString  fDefaultStorage;        ///< location of the default OCDB storage
  Int_t    fNPBins;                ///< number of momentum bins
  Double_t fPRange[2];             ///< momentum range
  Bool_t   fCorrectForSystematics; ///< add or not the systematic shifts of the residuals to the resolution
  Bool_t   fFitResiduals;          ///< fit or not the cluster residuals to extract means and sigmas
  UInt_t   fRequestedStationMask;  //!< mask of requested stations
  Bool_t   fRequest2ChInSameSt45;  //!< 2 fired chambers requested in the same station (4 or 5) or not
  Double_t fSigmaCut;              //!< sigma cut to associate clusters with TrackRefs
  Double_t fSigmaCutTrig;          //!< sigma cut to associate trigger track to triggerable track
  Double_t fClusterMaxRes[2];      //!< highest chamber resolution in both directions
  Int_t    fNDE;                   //!< total number of DE
  Int_t    fDEIndices[1100];       //!< index of DE in histograms refered by ID
  Int_t    fDEIds[200];            //!< ID of DE refered by index in histograms
  
  AliCFContainer* fCFContainer; //!< Pointer to the CF container
  TObjArray* fEfficiencyList;   //!< List of histograms for tracker/trigger efficiencies
  TObjArray* fTriggerList;      //!< List of histograms for trigger resolution
  TObjArray* fTrackerList;      //!< List of histograms for tracker resolution
  TObjArray* fPAtVtxList;       //!< List of graph and canvas about momentum resolution at vertex
  TObjArray* fSlopeAtVtxList;   //!< List of graph and canvas about slope resolution at vertex
  TObjArray* fEtaAtVtxList;     //!< List of graph and canvas about eta resolution at vertex
  TObjArray* fPhiAtVtxList;     //!< List of graph and canvas about phi resolution at vertex
  TObjArray* fPAt1stClList;     //!< List of graph and canvas about momentum resolution at first cluster
  TObjArray* fSlopeAt1stClList; //!< List of graph and canvas about slope resolution at first cluster
  TObjArray* fDCAList;          //!< List of graph and canvas about DCA
  TObjArray* fClusterList;      //!< List of graph and canvas about cluster resolution
  
  ClassDef(AliAnalysisTaskMuonPerformance, 1); // Muon performance analysis
};

inline void AliAnalysisTaskMuonPerformance::SetPBins(Int_t nBins, Double_t pMin, Double_t pMax)
{
  /// Set the binning to be used to study the detector resolution versus momentum
  if (nBins > 0) fNPBins = nBins;
  else AliError("Incorrect number of momentum bins");
  if (pMin >= 0. && pMax > pMin) {
    fPRange[0] = pMin;
    fPRange[1] = pMax;
  } else AliError("Incorrect momentum range");
}


#endif
