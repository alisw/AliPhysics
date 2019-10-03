#ifndef ALIMCTRUTHCENT_H
#define ALIMCTRUTHCENT_H

class TList;
class TH1;
class TH2;

class AliESDEvent;
class AliMCEvent;
class AliVParticle;

#include "AliAnalysisTaskSE.h"

class AliMCTruthCent : public AliAnalysisTaskSE {
 public:
  AliMCTruthCent();
  AliMCTruthCent(const char *name);
  virtual ~AliMCTruthCent();

  void SetV0ARange(Double_t mL, Double_t mH)     { fV0ALo = mL; fV0AHi = mH; }
  void SetV0CRange(Double_t mL, Double_t mH)     { fV0CLo = mL; fV0CHi = mH; }
  void SetV0MRange(Double_t mL, Double_t mH)     { fV0MLo = mL; fV0MHi = mH; }
  void SetFillHistos()                           { fFillHistos=kTRUE; DefineOutput(1, TList::Class()); }
  void UserCreateOutputObjects();
  void UserExec(Option_t *option);

 protected:
  AliVParticle      *GetTrack(Int_t i);
  
  TList            * fOutputList;           //! Output list
  Bool_t             fFillHistos;           //! flag to fill the QA histos
  TH1D             * fHMultV0A;             //!
  TH1D             * fHMultV0C;             //!
  TH1D             * fHMultV0M;             //!
  TH2D             * fHMultV0AvsV0C;        //!
  TH1D             * fHCentV0A;             //!
  TH1D             * fHCentV0C;             //!
  TH1D             * fHCentV0M;             //!
  TH2D             * fHCentV0AvsV0C;        //!
  Double_t           fV0ALo;                //! for linear centrality approximation
  Double_t           fV0AHi;                //!
  Double_t           fV0CLo;                //!
  Double_t           fV0CHi;                //!
  Double_t           fV0MLo;                //!
  Double_t           fV0MHi;                //!
  
 private:
  AliMCTruthCent(const AliMCTruthCent&);            // not implemented
  AliMCTruthCent &operator=(const AliMCTruthCent&); // not implemented

  ClassDef(AliMCTruthCent, 1); // Task to select tracks in MC events
};
#endif
