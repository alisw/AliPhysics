#ifndef ALIHADCORRTASK_H
#define ALIHADCORRTASK_H

// $Id$

class TClonesArray;
class TList;
class TH1F;
class TH2F;
class AliEmcalParticle;
class TString;

#include "AliAnalysisTaskEmcal.h"

class AliHadCorrTask : public AliAnalysisTaskEmcal {

 public:
  AliHadCorrTask();
  AliHadCorrTask(const char *name); 
  AliHadCorrTask(const char *name, Bool_t histo); 
  virtual ~AliHadCorrTask();

  void         UserCreateOutputObjects();
  void         Terminate(Option_t *);

  void         SetOutClusName(const char *n)           { fOutCaloName    = n    ; }
  void         SetEtaMatch(Double_t eta)               { fEtaMatch       = eta  ; }
  void         SetPhiMatch(Double_t phi)               { fPhiMatch       = phi  ; }
  void         SetTrackClus(Int_t c)                   { fDoTrackClus    = c    ; }
  void         SetHadCorr(Double_t c)                  { fHadCorr        = c    ; }
  void         SetEexcl(Double_t Emin)                 { fEexclCell      = Emin ; }

 protected:

  virtual Bool_t       Run()                                          ;
  virtual Bool_t       FillHistograms()                { return kTRUE ; }
  Int_t                GetMomBin(Double_t pt)                    const;
  Double_t             GetEtaSigma(Int_t pbin)                   const;
  Double_t             GetPhiMean(Int_t pbin, Int_t centbin)     const;
  Double_t             GetPhiSigma(Int_t pbin, Int_t centbin)    const;
  void                 DoTrackClusLoop()                              ;
  void                 DoMatchedTracksLoop(AliEmcalParticle *emccluster, Double_t &totalTrkP, Int_t &Nmatches);
  Double_t             ApplyHadCorrOneTrack(AliEmcalParticle *emccluster, Double_t hadCorr)                   ;
  Double_t             ApplyHadCorrAllTracks(AliEmcalParticle *emccluster, Double_t hadCorr)                  ;


  TString                fOutCaloName;            // name of output clusters
  Double_t               fPhiMatch;               // phi match value (pp=0.050)
  Double_t               fEtaMatch;               // eta match value (pp=0.025)
  Int_t                  fDoTrackClus;            // loop over tracks first
  Double_t               fHadCorr;                // hadronic correction (fraction)
  Double_t               fEexclCell;              // Energy/cell that we cannot subtract from the clusters

  TClonesArray          *fOutClusters;            //!output cluster collection

  TH2F                  *fHistMatchEtaPhi[8][9][2];  //!output histograms
  TH2F                  *fHistMatchEvsP[4];          //!output histograms
  TH2F                  *fHistNMatchEnergy[4];       //!output histograms
  TH2F                  *fHistNCellsEnergy[8][4];    //!output histograms
  TH2F                  *fHistMatchdRvsEP[4];        //!output histograms
  TH1F                  *fHistNclusvsCent;           //!output histograms
  TH1F                  *fHistNclusMatchvsCent;      //!output histograms
  TH1F                  *fHistEbefore;               //!output histograms
  TH1F                  *fHistEafter;                //!output histograms
  TH2F                  *fHistEoPCent;               //!output histograms
  TH2F                  *fHistNMatchCent;            //!output histograms
  TH2F                  *fHistNMatchCent_trk;        //!output histograms
  TH1F                  *fHistEsubPch[8];            //!output histograms
  TH2F                  *fHistEsubPchRat[8];         //!output histograms
  TH1F                  *fHistCentrality;            //!output histograms
  TH2F                  *fHistNoMatchEtaPhi;         //!output histograms

 private:
  AliHadCorrTask(const AliHadCorrTask&);            // not implemented
  AliHadCorrTask &operator=(const AliHadCorrTask&); // not implemented

  ClassDef(AliHadCorrTask, 9) // Hadronic correction task
};
#endif
