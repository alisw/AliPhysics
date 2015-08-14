#ifndef AliAnalysisTaskHMTFMCMultEst_cxx
#define AliAnalysisTaskHMTFMCMultEst_cxx

class TH1F;
class TH1I;
class TGraphErrors;

#include "AliAnalysisTaskSE.h"
#include "AliMultiplicityEstimators.h"

class AliAnalysisTaskHMTFMCMultEst : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskHMTFMCMultEst();
  AliAnalysisTaskHMTFMCMultEst(const char *name );
  virtual ~AliAnalysisTaskHMTFMCMultEst() {};

  void AddEstimator(const char* n);
  void InitEstimators();
  void SetRequireINELgt0(Bool_t b){fRequireINELgt0 = b;};
  MultiplicityEstimatorBase* MakeEstimator(const TString& name);
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

 private:
  TList *fMyOut;             // Output list
  TList *fRunconditions;     // TString with run conditions
  TList *fEstimatorsList;   // List to get the estimators out in terminate
  TString fEstimatorNames;
  Bool_t fRequireINELgt0;
  TNtuple *fParticleCounter;
  // MultiplicityEstimatorBase* festi;

  std::vector<MultiplicityEstimatorBase*> festimators;

  // Declaring these shuts up warnings from Weffc++
  AliAnalysisTaskHMTFMCMultEst(const AliAnalysisTaskHMTFMCMultEst&); // not implemented
  AliAnalysisTaskHMTFMCMultEst& operator=(const AliAnalysisTaskHMTFMCMultEst&); // not implemented

  ClassDef(AliAnalysisTaskHMTFMCMultEst, 1); // example of analysis
};

#endif
