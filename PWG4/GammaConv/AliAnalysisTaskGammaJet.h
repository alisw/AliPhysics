#ifndef AliAnalysisTaskGammaJet_cxx
#define AliAnalysisTaskGammaJet_cxx

// example of an analysis task creating a p_t spectrum
// Authors: Panos Cristakoglou, Jan Fiete Grosse-Oetringhaus, Christian Klein-Boesing

class TH1F;
class AliESDEvent;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskGammaJet : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskGammaJet(const char *name);
  virtual ~AliAnalysisTaskGammaJet() {}
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);

  void SetDeltaAODFileName(TString string) { fDeltaAODFileName = string;}
  
 private:

  TClonesArray * GetConversionGammas();
  
  AliESDEvent *fESD;    //! ESD object
  TList       *fOutputList; //! Output list
  TH1F        *fHistPt; //! Pt spectrum
  TH1F        *fHistPtPhos; //! Pt spectrum
  TH1F        *fHistPtEmcal; //! Pt spectrum
  TH1F        *fHistPtJets; //! Pt spectrum
   
  TString     fDeltaAODFileName;//! File where Gamma Conv AOD is located, if not in default AOD

  AliAnalysisTaskGammaJet(); //not implemented
  AliAnalysisTaskGammaJet(const AliAnalysisTaskGammaJet&); // not implemented
  AliAnalysisTaskGammaJet& operator=(const AliAnalysisTaskGammaJet&); // not implemented
  
  ClassDef(AliAnalysisTaskGammaJet, 1); // example of analysis
};

#endif
