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
  void                   SetMCName(const char *n) {fMCParticlesName = n;}

 protected:
  Bool_t                 Run();
  void                   ExecOnce();
  virtual Int_t          GetCentBin(Double_t cent) const;
  Float_t                RelativePhi(Double_t mphi,Double_t vphi) const;

 private:
  TString             fCalo2Name;
  TClonesArray       *fCaloClusters2;
  TString             fMCParticlesName;
  TClonesArray       *fMCParticles;
  TH2F                  *fHistRhovsCent;           //!rho vs cent
  TH2F                  *fHistNjetvsCent;          //!number of jets versus Centrality
  TH2F               *fHistNEFvsPt[3];          //!
  TH2F               *fHistNTMatchvsPt[3];      //!  
  TH2F               *fHistNCMatchvsPt[3];      //!
  TH2F               *fHistHadCorvsPt[3];       //!
  TH2F               *fHistNconvsPt[3];         //!
  TH2F               *fHistNtvsPt[3];         //!
  TH2F               *fHistNcvsPt[3];         //!
  TH2F               *fHistNEFvsPtBias[3];      //!
  TH2F               *fHistNTMatchvsPtBias[3];  //!  
  TH2F               *fHistNCMatchvsPtBias[3];      //!
  TH2F               *fHistHadCorvsPtBias[3];      //!
  TH2F               *fHistNconvsPtBias[3];         //!
  TH2F               *fHistNtvsPtBias[3];         //!
  TH2F               *fHistNcvsPtBias[3];         //!
  TH3F               *fHistNTMatchvsPtvsNtack0;      //!


  AliAnalysisTaskEmcalJetHadCorQA(const AliAnalysisTaskEmcalJetHadCorQA&); // not implemented
  AliAnalysisTaskEmcalJetHadCorQA& operator=(const AliAnalysisTaskEmcalJetHadCorQA&); // not implemented
  
  ClassDef(AliAnalysisTaskEmcalJetHadCorQA, 2); // Hadronic Correction for jet task
};
#endif
