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
#include "AliQnCorrectionsFillEventTask.h"

class AliAnalysis;
class AliQnCorrectionsManager;
class AliQnCorrectionsCutsSet;
class AliQnCorrectionsHistos;
class TList;
class TProfile;
class TGraphErrors;


//_________________________________________________________
class AliAnalysisTaskQnVectorAnalysis : public AliQnCorrectionsFillEventTask {

public:

  enum enumTrackDetectors{
    kTPC=0,
    kSPD,
    nTrackDetectors
  };

  enum enumEPdetectors{
    kVZEROA=0,
    kVZEROC,
    kTZEROA,
    kTZEROC,
    kFMDA,
    kFMDC,
//    kRawFMDA,
//    kRawFMDC,
    nEPDetectors
  };

  enum Constants{
    kNharmonics=4,
    kNresolutionComponents=3, /* valid for 3-(sub-event)detector method */
    kNxy=2
  };

  enum CorrelationConstants{
    kXX = 0,
    kXY,
    kYX,
    kYY,
    kNcorrelationComponents,
  };

  /* this is only valid if we are using 3-(sub-event)detector method
   * for evaluating the detector resolution
   */
  enum SubEventConstants {
    kABXX = 0,
    kABYY,
    kACXX,
    kACYY,
    kBCXX,
    kBCYY,
    nCorrelationPerDetector
  };

  AliAnalysisTaskQnVectorAnalysis();
  AliAnalysisTaskQnVectorAnalysis(const char *name);
  virtual ~AliAnalysisTaskQnVectorAnalysis();


  virtual void UserExec(Option_t *);
  virtual void UserCreateOutputObjects();
  virtual void FinishTaskOutput();

  AliQnCorrectionsHistos* GetHistograms() {return fEventPlaneHistos;}
  AliQnCorrectionsCutsSet* EventCuts()  const {return fEventCuts;}
  Bool_t IsEventSelected(Float_t* values);

  void SetEventCuts(AliQnCorrectionsCutsSet* cuts)  {fEventCuts = cuts;}
  void SetCentralityVariable(Int_t var) { fCentralityVariable = var; }
  void SetExpectedCorrectionPass(const char *pass) { fExpectedCorrectionPass = pass; }
  void SetAlternativeCorrectionPass(const char *pass) { fAlternativeCorrectionPass = pass; }

 private:
  TList* fEventQAList;
  AliQnCorrectionsCutsSet *fEventCuts;
  AliQnCorrectionsHistos* fEventPlaneHistos;

  AliAnalysisTaskQnVectorAnalysis(const AliAnalysisTaskQnVectorAnalysis &c);
  AliAnalysisTaskQnVectorAnalysis& operator= (const AliAnalysisTaskQnVectorAnalysis &c);

  TProfile* fVn[nTrackDetectors*nEPDetectors][kNharmonics][kNcorrelationComponents];

  Int_t fNDetectorResolutions;
  TGraphErrors* fDetectorResolution[nTrackDetectors*nEPDetectors][kNharmonics];
  Int_t fDetectorResolutionContributors[nTrackDetectors*nEPDetectors][kNresolutionComponents];
  TProfile *fDetectorResolutionCorrelations[nTrackDetectors*nEPDetectors][nCorrelationPerDetector][kNharmonics];

  TString fTrackDetectorNameInFile[nTrackDetectors];
  TString fEPDetectorNameInFile[nEPDetectors];

  Int_t fCentralityVariable;
  TString fExpectedCorrectionPass;
  TString fAlternativeCorrectionPass;

  ClassDef(AliAnalysisTaskQnVectorAnalysis, 1);
};

#endif

