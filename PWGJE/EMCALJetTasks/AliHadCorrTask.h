#ifndef ALIHADCORRTASK_H
#define ALIHADCORRTASK_H

// $Id$

class TClonesArray;
class TList;
class TH1;
class TH2;
class AliEmcalParticle;
class TString;

#include "AliAnalysisTaskEmcal.h"

class AliHadCorrTask : public AliAnalysisTaskEmcal {

 public:
  AliHadCorrTask();
  AliHadCorrTask(const char *name, Bool_t histo=kFALSE); 
  virtual ~AliHadCorrTask();

  void                   UserCreateOutputObjects();

  void                   SetEexcl(Double_t Emin)                 { fEexclCell      = Emin ; }
  void                   SetEtaMatch(Double_t eta)               { fEtaMatch       = eta  ; }
  void                   SetHadCorr(Double_t c)                  { fHadCorr        = c    ; }
  void                   SetOutClusName(const char *n)           { fOutCaloName    = n    ; }
  void                   SetPhiMatch(Double_t phi)               { fPhiMatch       = phi  ; }
  void                   SetTrackClus(Int_t c)                   { fDoTrackClus    = c    ; }
  void                   SetDoExact(Bool_t d)                    { fDoExact        = d    ; }

 protected:
  Double_t               ApplyHadCorrOneTrack(AliEmcalParticle *emccluster, Double_t hadCorr);
  Double_t               ApplyHadCorrAllTracks(AliEmcalParticle *emccluster, Double_t hadCorr);
  void                   DoMatchedTracksLoop(AliEmcalParticle *emccluster, Double_t &totalTrkP, Int_t &Nmatches, Double_t &trkPMCfrac, Int_t &NMCmatches);
  void                   DoTrackLoop();
  Double_t               GetEtaSigma(Int_t pbin)                   const;
  Int_t                  GetMomBin(Double_t pt)                    const;
  Double_t               GetPhiMean(Int_t pbin, Int_t centbin)     const;
  Double_t               GetPhiSigma(Int_t pbin, Int_t centbin)    const;
  Bool_t                 Run()                                          ;
  void                   ExecOnce()                                     ;

  // Task configuration
  TString                fOutCaloName;               // name of output clusters
  Double_t               fPhiMatch;                  // phi match value (pp=0.050)
  Double_t               fEtaMatch;                  // eta match value (pp=0.025)
  Int_t                  fDoTrackClus;               // loop over tracks first
  Double_t               fHadCorr;                   // hadronic correction (fraction)
  Double_t               fEexclCell;                 // energy/cell that we cannot subtract from the clusters
  Bool_t                 fDoExact;                   // do exact correction (embedding only) 

  // Service fields (non-streamed)
  Bool_t                 fEsdMode;                   //!ESD/AOD mode
  TClonesArray          *fOutClusters;               //!output cluster collection

  // QA plots
  TH2                   *fHistMatchEtaPhi[8][9][2];  //!deta vs. dphi of matched cluster-track pairs
  TH2                   *fHistMatchEvsP[4];          //!cluster energy vs. track momentum of matched pairs
  TH2                   *fHistNMatchEnergy[4];       //!n matches vs. cluster energy
  TH2                   *fHistNCellsEnergy[4][4];    //!n cells vs. cluster energy
  TH2                   *fHistMatchdRvsEP[4];        //!matching distance vs. E/P
  TH1                   *fHistNclusvsCent;           //!n clusters vs. centrality
  TH1                   *fHistNclusMatchvsCent;      //!n clusters matched to some track vs. centrality
  TH1                   *fHistEbefore;               //!average energy of clusters before correction vs. centrality
  TH1                   *fHistEafter;                //!average energy of clusters after correction vs. centrality
  TH2                   *fHistEoPCent;               //!E/P vs. centrality
  TH2                   *fHistNMatchCent;            //!n matches vs. centraity
  TH2                   *fHistNClusMatchCent;        //!n clusters macthed to some track (tracks allowed to match more than one cluster)
  TH1                   *fHistEsubPch[8];            //!Esub vs. total momentum of matched tracks (only 1 match)
  TH2                   *fHistEsubPchRat[8];         //!Esub/momentum of matched tracks vs. total momentum of matched tracks (only 1 match)
  TH2                   *fHistEsubPchRatAll[8];      //!Esub/momentum of matched tracks vs. total momentum of matched tracks (all number of matches)

  // Embedding QA plots
  TH2                   *fHistEmbTrackMatchesOversub[4];    //!Over-subtracted energy with embedded track matches (non-embedded matches < 5%)
  TH2                   *fHistNonEmbTrackMatchesOversub[4]; //!Over-subtracted energy with non-embedded track matches (embedded matches < 5%)
  TH2                   *fHistOversub[4];                   //!Over-subtracted energy
  TH2                   *fHistOversubOverClusE[4];          //!Over-subtracted energy / cluster energy

 private:
  AliHadCorrTask(const AliHadCorrTask&);            // not implemented
  AliHadCorrTask &operator=(const AliHadCorrTask&); // not implemented

  ClassDef(AliHadCorrTask, 12) // Hadronic correction task
};
#endif
