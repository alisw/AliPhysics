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

  void         SetClusName(const char *n)            { fCaloName   = n; }
  void         SetHadCorr(Double_t c)                { fHadCorr    = c; }
  void         SetMinPt(Double_t min)                { fMinPt = min;  }
  void         SetTracksName(const char *n)          { fTracksName = n; }
  void         SetOutClusName(const char *n)         { fOutCaloName = n;}

 protected:
  Int_t        GetCentBin(Double_t cent) const;
  Int_t        GetMomBin(Double_t pt)    const;

  TString                fTracksName;             // name of track collection
  TString                fCaloName;               // name of calo cluster collection
  TString                fOutCaloName;            // name of output clusters
  Double_t               fHadCorr;                // hadronic correction (fraction)
  Double_t               fMinPt;                  // minimum pt (on tracks and clusters)
  TClonesArray          *fOutClusters;            //!output cluster collection
  TList                 *fOutputList;             //!output list
  TH2                   *fHistMatchEtaPhi[4][5];  //!output histograms
  TH2                   *fHistMatchEvsP[4];       //!output histograms
  TH2                   *fHistMatchdRvsEP[4];     //!output histograms
  TH1                   *fHistNclusvsCent;        //!output histograms
  TH1                   *fHistNclusMatchvsCent;   //!output histograms
  TH1                   *fHistEbefore;            //!output histograms
  TH1                   *fHistEafter;             //!output histograms
  TH2                   *fHistEoPCent;            //!output histograms
  TH2                   *fHistNMatchCent;         //!output histograms

 private:
  AliHadCorrTask(const AliHadCorrTask&);            // not implemented
  AliHadCorrTask &operator=(const AliHadCorrTask&); // not implemented

  ClassDef(AliHadCorrTask, 2) // Hadronic correction task
};
#endif
