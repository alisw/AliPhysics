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
  AliAnalysisTaskFilterHe3(const char *name);
  virtual ~AliAnalysisTaskFilterHe3();

  virtual void UserCreateOutputObjects();
  virtual void UserExec(Option_t *option);
  virtual void Terminate(Option_t *option);

  void SetUseMultTaskCentrality(Bool_t useMultTaskCentrality = kTRUE) { fUseMultTaskCentrality = useMultTaskCentrality; };

  // Setters for some cuts
  void SetMinNsig3He(Double_t nsmin = 3) { fMinNSigma3He = nsmin; }
  void SetMaxNsig3He(Double_t nsmax = 10) { fMaxNSigma3He = nsmax; }
  void SetNclsTPC(Double_t ncltpc = 50) { fMinNclsTPC = ncltpc; }
  void SetPtotmin(Double_t ptotmin = 0.5) { fMinPtot = ptotmin; }
  void SetPtotmax(Double_t ptotmax = 20) { fMaxPtot = ptotmax; }
  void SetFillSecifTOF (Bool_t flgsectof = kFALSE) {fillSecifTOF = flgsectof; }

private:
  AliESDEvent *fESD;                     //! input event
  TList *fOutputList;                    //! output list
  TTree *fListOfFiles;                   //! list of files
  AliESDtrackCuts *fESDtrackCuts;        //! input track cuts (secondary & primary)
  AliESDtrackCuts *fESDtrackCutsPrimary; //! input track cuts (only primary)
  AliPIDResponse *fPIDResponse;          //! pid response object
  AliMultSelection *fMultSel;            //! centrality and multiplicity selection
  Bool_t fUseMultTaskCentrality;         // flag: true if centrality should be taken from AliMultSelectionTask, false if traditional method should be used
  //
  Int_t fEventIdFile; //! event id in file
  TString fFileName;  //! chunk file name

  // Cut params
  Double_t fMinNSigma3He;
  Double_t fMaxNSigma3He;
  Double_t fMinNclsTPC;
  Double_t fMinPtot;
  Double_t fMaxPtot;
  Bool_t fillSecifTOF;

  //
  // histograms
  //
  TH1F *fHistZv;                //! z-Vertex distribution
  TH2F *fHistdEdxData;          //! PID histogram dEdx all particles
  TH2F *fHistTof;               //! PID histogram TOF all particles
  TH3F *fHistdEdxDeuteronParam; //! PID-QA histogram for deuteron Bethe-Bloch parameterisation
  TH3F *fHistdEdxHe3Param;      //! PID-QA histogram for he3 Bethe-Bloch parameterisation
  TH3F *fHistdEdxTritonParam;   //! PID-QA histogram for triton Bethe-Bloch parameterisation
  TH1F *fHistCent;              //! centrality histogram
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
  AliAnalysisTaskFilterHe3(const AliAnalysisTaskFilterHe3 &);
  AliAnalysisTaskFilterHe3 &operator=(const AliAnalysisTaskFilterHe3 &);
  //
  ClassDef(AliAnalysisTaskFilterHe3, 1);
};

#endif
