#ifndef ALIHADCORRTASK_H
#define ALIHADCORRTASK_H

// $Id$

class TClonesArray;
class TList;
class TH1;
class TH2;

#include "AliAnalysisTaskSE.h"

class AliHadCorrTask : public AliAnalysisTaskSE {
 public:
  AliHadCorrTask();
  AliHadCorrTask(const char *name); 
  virtual ~AliHadCorrTask();

  void         UserCreateOutputObjects();
  void         UserExec(Option_t *option);
  void         Terminate(Option_t *option);

  void         SetClusName(const char *n)            { fCaloName    = n;   }
  void         SetEtaMatch(Double_t eta)             { fEtaMatch    = eta; }
  void         SetHadCorr(Double_t c)                { fHadCorr     = c;   }
  void         SetMinPt(Double_t min)                { fMinPt       = min; }
  void         SetOutClusName(const char *n)         { fOutCaloName = n;   }
  void         SetPhiMatch(Double_t phi)             { fPhiMatch    = phi; }
  void         SetTracksName(const char *n)          { fTracksName  = n;   }

 protected:
  Int_t        GetCentBin(Double_t cent) const;
  Int_t        GetMomBin(Double_t pt)    const;

  TString                fTracksName;             // name of track collection
  TString                fCaloName;               // name of calo cluster collection
  TString                fOutCaloName;            // name of output clusters
  Double_t               fPhiMatch;               // phi match value (pp=0.050)
  Double_t               fEtaMatch;               // eta match value (pp=0.025)
  Double_t               fHadCorr;                // hadronic correction (fraction)
  Double_t               fMinPt;                  // minimum pt (on tracks and clusters)
  TClonesArray          *fOutClusters;            //!output cluster collection
  TList                 *fOutputList;             //!output list
  TH2                   *fHistMatchEtaPhi[8][5];  //!output histograms
  TH2                   *fHistMatchEvsP[4];       //!output histograms
  TH2                   *fHistMatchdRvsEP[4];     //!output histograms
  TH1                   *fHistNclusvsCent;        //!output histograms
  TH1                   *fHistNclusMatchvsCent;   //!output histograms
  TH1                   *fHistEbefore;            //!output histograms
  TH1                   *fHistEafter;             //!output histograms
  TH2                   *fHistEoPCent;            //!output histograms
  TH2                   *fHistNMatchCent;         //!output histograms
  TH2                   *fHistEsubPch[4];         //!output histograms

 private:
  AliHadCorrTask(const AliHadCorrTask&);            // not implemented
  AliHadCorrTask &operator=(const AliHadCorrTask&); // not implemented

  ClassDef(AliHadCorrTask, 3) // Hadronic correction task
};
#endif
