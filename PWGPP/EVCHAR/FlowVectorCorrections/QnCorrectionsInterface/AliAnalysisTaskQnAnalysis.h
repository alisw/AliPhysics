/*
***********************************************************
  event plane corrections framework
  contact: jaap onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
  2014/12/10
  *********************************************************
*/

//#include "AliSysInfo.h"

#ifndef ALIANALYSISTASKQNANALYSIS_H
#define ALIANALYSISTASKQNANALYSIS_H

#include "TFile.h"
#include "TTree.h"
#include "AliAnalysisTaskSE.h"

class AliAnalysis;
class AliQnCorrectionsCuts;
class AliQnCorrectionsFillEvent;
class AliQnCorrectionsHistos;
class AliQnCorrectionsQnVector;
class TList;
class TProfile;

//_________________________________________________________
class AliAnalysisTaskQnAnalysis : public AliAnalysisTaskSE {


  enum Qvectors{
    kTPC=0,
    kVZEROA,
    kVZEROC,
    kFMDA,
    kFMDC,
    kNqvectors
  };

  enum Constants{
    kNcorrelation=6,
    kNharmonics=4
  };

 public:
  AliAnalysisTaskQnAnalysis();
  AliAnalysisTaskQnAnalysis(const char *name);
  ~AliAnalysisTaskQnAnalysis();


  virtual void UserExec(Option_t *);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();

  AliQnCorrectionsHistos* GetHistograms() {return fEventPlaneHistos;}
  AliQnCorrectionsFillEvent* GetFillEvent() {return fFillEvent;}
  AliQnCorrectionsCuts* EventCuts()  const {return fEventCuts;}
  Bool_t IsEventSelected(Float_t* values);
  AliQnCorrectionsQnVector* GetQvector(TList* l, TString det, Int_t maxcorrection);

  void SetEventCuts(AliQnCorrectionsCuts* cuts)  {fEventCuts = cuts;}

 private:
  TList* fEventQAList;
  AliQnCorrectionsCuts * fEventCuts;
  AliQnCorrectionsHistos* fEventPlaneHistos;
  AliQnCorrectionsFillEvent* fFillEvent;
  
  AliAnalysisTaskQnAnalysis(const AliAnalysisTaskQnAnalysis &c);
  AliAnalysisTaskQnAnalysis& operator= (const AliAnalysisTaskQnAnalysis &c);

  TProfile* fVn[kNcorrelation][kNharmonics][4];//!
  TProfile* fCorrelation[kNcorrelation*3][kNharmonics][4];//!

  Int_t fEPflow[kNcorrelation][2];
  Int_t fEPcorrelation[kNcorrelation*3][2];

  ClassDef(AliAnalysisTaskQnAnalysis, 1);
};

#endif


