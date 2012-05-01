#ifndef ALIHADCORRTASK_H
#define ALIHADCORRTASK_H

// $Id$

class TClonesArray;
class AliESDtrackCuts;
class TList;
class TH1;
class TH2;

#include "AliAnalysisTaskSE.h"

class AliHadCorrTask : public AliAnalysisTaskSE {
 public:
  AliHadCorrTask(const char *name); 
  AliHadCorrTask();
  virtual ~AliHadCorrTask();

  void         UserCreateOutputObjects();
  void         UserExec(Option_t *option);
  void         Terminate(Option_t *option);

  void         SetClusName(const char *n)            { fCaloName   = n; }
  void         SetHadCorr(Double_t c)                { fHadCorr    = c; }
  void         SetTracksName(const char *n)          { fTracksName = n; }
  void         SetMinPt(Double_t min)                { fMinPt = min;  }
  void         SetOutClusName(const char *n)         { fOutCaloName = n;}

 protected:
  void FindJets(TObjArray *tracks, TObjArray *clus, Int_t algo, Double_t radius, Float_t fCent);

  TString                fTracksName;             // name of track collection
  TString                fCaloName;               // name of calo cluster collection

  Double_t               fHadCorr;                // hadronic correction
  TClonesArray          *fJets;                   //!jet collection
  Double_t               fMinPt; 
  TString                fOutCaloName;            //name of output clusters
  TClonesArray          *fOutClusters;           //output cluster collection

  TH2                   *fHistMatchEtaPhi[4][5];
  TList                 *fOutputList;
  TH2                   *fHistMatchEvsP[4];
  TH2                   *fHistMatchdRvsEP[4];
  TH1                   *fHistNclusvsCent;
  TH1                   *fHistNclusMatchvsCent;
  TH1                   *fHistEbefore;
  TH1                   *fHistEafter;
  TH2                   *fHistEoPCent;
  TH2                   *fHistNMatchCent;

 private:
  AliHadCorrTask(const AliHadCorrTask&);            // not implemented
  AliHadCorrTask &operator=(const AliHadCorrTask&); // not implemented

  ClassDef(AliHadCorrTask, 1) 
};
#endif
