#ifndef AliAnalysisTaskEmcalJetHadCorQA_h
#define AliAnalysisTaskEmcalJetHadCorQA_h

class TH1F;
class TH2F;
class TH3F;
class THnSparse;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalJetHadCorQA : public AliAnalysisTaskEmcalJet {
 public:
  AliAnalysisTaskEmcalJetHadCorQA();
  AliAnalysisTaskEmcalJetHadCorQA(const char *name);
  virtual ~AliAnalysisTaskEmcalJetHadCorQA() {}
  
  
  virtual void           UserCreateOutputObjects();
  void                   SetCalo2Name(const char *n) {fCalo2Name = n;}

 protected:
  Bool_t                 Run();
  void                   ExecOnce();
  virtual Int_t          GetCentBin(Double_t cent) const;
  Float_t                RelativePhi(Double_t mphi,Double_t vphi) const;

 private:
  TString             fCalo2Name;
  TClonesArray       *fCaloClusters2;
  TH2F                  *fHistRhovsCent;           //!rho vs cent
  TH2F                  *fHistNjetvsCent;          //!number of jets versus Centrality
  TH2F               *fHistNEFvsPt[3];      //!
  TH2F               *fHistNTMatchvsPt[3];      //!  
  TH2F               *fHistNCMatchvsPt[3];      //!
  TH2F               *fHistHadCorvsPt[3];      //!
  TH2F               *fHistNEFvsPtBias[3];      //!
  TH2F               *fHistNTMatchvsPtBias[3];      //!  
  TH2F               *fHistNCMatchvsPtBias[3];      //!
  TH2F               *fHistHadCorvsPtBias[3];      //!
  TH3F               *fHistNTMatchvsPtvsNtack0;      //!


  AliAnalysisTaskEmcalJetHadCorQA(const AliAnalysisTaskEmcalJetHadCorQA&); // not implemented
  AliAnalysisTaskEmcalJetHadCorQA& operator=(const AliAnalysisTaskEmcalJetHadCorQA&); // not implemented
  
  ClassDef(AliAnalysisTaskEmcalJetHadCorQA, 1); // Emcal jet spectra task
};
#endif
