#ifndef AliAnalysisTaskProtons_cxx
#define AliAnalysisTaskProtons_cxx

// Analysis task creating a the 2d y-p_t spectrum of p and antip
// Author: Panos Cristakoglou

class TList;
class AliESDEvent;
#include "TF1.h"

#include "PWG2spectra/SPECTRA/AliProtonAnalysis.h"
#include "AliAnalysisTask.h"

class AliAnalysisTaskProtons : public AliAnalysisTask {
 public:
  AliAnalysisTaskProtons(const char *name = "AliAnalysisTaskProtons");
  virtual ~AliAnalysisTaskProtons() {}
  
  virtual void   ConnectInputData(Option_t *);
  virtual void   CreateOutputObjects();
  virtual void   Exec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetPriorProbabilityFunctions(TF1 *felectrons, 
				    TF1 *fmuons,
				    TF1 *fpions,
				    TF1 *fkaons,
				    TF1 *fprotons) {
    fFunctionUsed = kTRUE;
    fElectronFunction = felectrons;
    fMuonFunction = fmuons;
    fPionFunction = fpions;
    fKaonFunction = fkaons;
    fProtonFunction = fprotons;
  }
  
 private:
  AliESDEvent *fESD;    //ESD object
  TList  *fList; //TList output object
  AliProtonAnalysis *fAnalysis; //analysis object
  TF1 *fElectronFunction; //TF1 for e
  TF1 *fMuonFunction; //TF1 for mu
  TF1 *fPionFunction; //TF1 for pi
  TF1 *fKaonFunction; //TF1 for K
  TF1 *fProtonFunction; //TF1 for p

  Bool_t fFunctionUsed; //kTRUE if Functions are used

  AliAnalysisTaskProtons(const AliAnalysisTaskProtons&); // not implemented
  AliAnalysisTaskProtons& operator=(const AliAnalysisTaskProtons&); // not implemented
  
  ClassDef(AliAnalysisTaskProtons, 1); // example of analysis
};

#endif
