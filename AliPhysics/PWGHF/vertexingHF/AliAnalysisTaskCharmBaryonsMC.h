#ifndef AliAnalysisTaskCharmBaryonsMC_H
#define AliAnalysisTaskCharmBaryonsMC_H

//
// Task used to analize simulations at generation level (i.e. only
// needs galice.root and Kinematics.root).
// 

class TH1F;
class TH1I;
class TGraphErrors;

#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCharmBaryonsMC : public AliAnalysisTaskSE {
 public:
  AliAnalysisTaskCharmBaryonsMC();
  AliAnalysisTaskCharmBaryonsMC(const char *name );
  virtual ~AliAnalysisTaskCharmBaryonsMC();
  
  virtual void   UserCreateOutputObjects();
  virtual void   UserExec(Option_t *option);
  virtual void   Terminate(Option_t *);
  void BookHistograms();

  void GetHistory(AliMCParticle *part, AliMCEvent *mcevt, Int_t *pdgarray, Int_t *labelarray, Int_t &ngen);
  Bool_t Match(Int_t a, Int_t b);
  Bool_t Match(Float_t a, Float_t b);
  Bool_t FromBottom(Int_t *hitory);
  void Finalize();
  TList * GetList() const { return fMyOut;} 
 private:
  
  AliAnalysisTaskCharmBaryonsMC(const AliAnalysisTaskCharmBaryonsMC&); // not implemented
  AliAnalysisTaskCharmBaryonsMC& operator=(const AliAnalysisTaskCharmBaryonsMC&); // not implemented

  TList * fMyOut; //!<! list of output histos
  TH1F * fHistEvt; //!<! number of events 
  TH1F * fHistMult; //!<! multiplicity
  TH1D *fHistPtPromptD0; //!<! D0 pT distribution in |y|<0.5
  TH1D *fHistPtPromptLc; //!<! Lc pT distribution in |y|<0.5
  TH1D *fHistPtPromptXic0; //!<! Xic0 pT distribution in |y|<0.5
  TH1D *fHistPtFeeddownD0; //!<! D0 pT distribution in |y|<0.5
  TH1D *fHistPtFeeddownLc; //!<! Lc pT distribution in |y|<0.5
  TH1D *fHistPtFeeddownXic0; //!<! Xic0 pT distribution in |y|<0.5
  TH1D *fHistPtInclusiveD0; //!<! D0 pT distribution in |y|<0.5
  TH1D *fHistPtInclusiveLc; //!<! Lc pT distribution in |y|<0.5
  TH1D *fHistPtInclusiveXic0; //!<! Xic0 pT distribution in |y|<0.5
  TH2D *fHistPtvsRapidityPromptD0; //!<! D0 pT distribution in |y|<0.5
  TH2D *fHistPtvsRapidityPromptLc; //!<! Lc pT distribution in |y|<0.5
  TH2D *fHistPtvsRapidityPromptXic0; //!<! Xic0 pT distribution in |y|<0.5
  TH2D *fHistPtvsRapidityFeeddownD0; //!<! D0 pT distribution in |y|<0.5
  TH2D *fHistPtvsRapidityFeeddownLc; //!<! Lc pT distribution in |y|<0.5
  TH2D *fHistPtvsRapidityFeeddownXic0; //!<! Xic0 pT distribution in |y|<0.5
  TH2D *fHistPtvsRapidityInclusiveD0; //!<! D0 pT distribution in |y|<0.5
  TH2D *fHistPtvsRapidityInclusiveLc; //!<! Lc pT distribution in |y|<0.5
  TH2D *fHistPtvsRapidityInclusiveXic0; //!<! Xic0 pT distribution in |y|<0.5

  ClassDef(AliAnalysisTaskCharmBaryonsMC, 1); 
};

#endif
