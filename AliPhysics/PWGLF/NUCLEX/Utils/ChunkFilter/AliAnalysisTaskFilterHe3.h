#ifndef AliAnalysisTaskFilterHe3_H
#define AliAnalysisTaskFilterHe3_H

#include "AliAnalysisTaskSE.h"
#include "TString.h"
#include "TTree.h"

class AliESDEvent;
class AliESDtrackCuts;
class AliMCEvent;
class TList;
class TH1F;
class TH2F;
class TH3F;
class AliPIDResponse;
class AliMultSelection;


class AliAnalysisTaskFilterHe3 : public AliAnalysisTaskSE
{
 public:
  AliAnalysisTaskFilterHe3();
  AliAnalysisTaskFilterHe3(const char* name);
  virtual ~AliAnalysisTaskFilterHe3();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t* option);
  virtual void Terminate(Option_t* option);

  void SetUseMultTaskCentrality(Bool_t useMultTaskCentrality = kTRUE) {fUseMultTaskCentrality = useMultTaskCentrality;};

 private:
  AliESDEvent* fESD;                           //! input event
  TList* fOutputList;                          //! output list
  TTree * fListOfFiles;                        //! list of files 
  AliESDtrackCuts* fESDtrackCuts;              //! input track cuts (secondary & primary)
  AliESDtrackCuts* fESDtrackCutsPrimary;       //! input track cuts (only primary)
  AliPIDResponse*        fPIDResponse;         //! pid response object
  AliMultSelection*      fMultSel;             //! centrality and multiplicity selection
  Bool_t fUseMultTaskCentrality;               // flag: true if centrality should be taken from AliMultSelectionTask, false if traditional method should be used
  //
  Int_t fEventIdFile;                        //! event id in file
  TString fFileName;                           //! chunk file name
  //
  // histograms
  //
  TH1F* fHistZv;                               //! z-Vertex distribution
  TH2F* fHistdEdxData;                         //! PID histogram dEdx all particles
  TH2F* fHistTof;                              //! PID histogram TOF all particles
  TH3F* fHistdEdxDeuteronParam;                //! PID-QA histogram for deuteron Bethe-Bloch parameterisation
  TH3F* fHistdEdxHe3Param;                     //! PID-QA histogram for he3 Bethe-Bloch parameterisation
  TH3F* fHistdEdxTritonParam;                  //! PID-QA histogram for triton Bethe-Bloch parameterisation
  TH1F* fHistCent;                             //! centrality histogram
  //
  //
  // Bethe-Bloch parameterisations (hard-coded for the time being)
  //
  Double_t fParamDeuteron[5];
  Double_t fParamTriton[5];
  Double_t fParamHe3[5];
  Double_t fParamDeuteronMC[5];
  Double_t fParamTritonMC[5];
  Double_t fParamHe3MC[5];
  //
  //
  AliAnalysisTaskFilterHe3(const AliAnalysisTaskFilterHe3&);
  AliAnalysisTaskFilterHe3& operator=(const AliAnalysisTaskFilterHe3&);
  //
  ClassDef(AliAnalysisTaskFilterHe3, 1);
  
};

#endif
  
  
