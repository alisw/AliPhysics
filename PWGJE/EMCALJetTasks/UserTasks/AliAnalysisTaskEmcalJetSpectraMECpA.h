#ifndef AliAnalysisTaskEmcalJetSpectraMECpA_h
#define AliAnalysisTaskEmcalJetSpectraMECpA_h

// $Id: AliAnalysisTaskEmcalJetSpectraMECpA.h 3010 2012-06-10 05:40:56Z loizides $


class TH1F;
class TH2F;
class TH3F;
class THnSparse;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalJetSpectraMECpA : public AliAnalysisTaskEmcalJet {
 public:
  AliAnalysisTaskEmcalJetSpectraMECpA();
  AliAnalysisTaskEmcalJetSpectraMECpA(const char *name);
  virtual ~AliAnalysisTaskEmcalJetSpectraMECpA() {}
  
  
  virtual void           UserCreateOutputObjects();

 protected:
  Bool_t                 Run();
  virtual Int_t          GetCentBin(Double_t cent) const;
  Float_t                RelativePhi(Double_t mphi,Double_t vphi) const;

 private:
  TH2F                  *fHistRhovsCent;  //!
  TH2F                  *fHistRhoScvsCent;  //!
  TH2F                  *fHistNjetvsCent; //!number of jets versus Centrality
  TH2F                  *fHistJetPtvsTrackPt[7];//!
  TH2F                  *fHistJetPtScvsTrackPt[7];//!
  TH2F                  *fHistRawJetPtvsTrackPt[7];//!
  TH1F                  *fHistTrackPt[7];//!
  TH1F                  *fHistEP0[7];//!
  TH1F                  *fHistEP0A[7];//!
  TH1F                  *fHistEP0C[7];//!
  TH2F                  *fHistEPAvsC[7];//!
  TH2F                  *fHistJetPtvsdEP[7];//!
  TH2F                  *fHistJetPtvsdEPBias[7];//!
  TH2F                  *fHistJetPtvsEP[7];//!
  TH2F                  *fHistJetPtvsEPBias[7];//!
  TH2F                  *fHistRhovsEP[7]; //!
  TH3F                  *fHistJetPtEtaPhi[7];




  AliAnalysisTaskEmcalJetSpectraMECpA(const AliAnalysisTaskEmcalJetSpectraMECpA&); // not implemented
  AliAnalysisTaskEmcalJetSpectraMECpA& operator=(const AliAnalysisTaskEmcalJetSpectraMECpA&); // not implemented
  
  ClassDef(AliAnalysisTaskEmcalJetSpectraMECpA, 4); // Emcal jet spectra task
};
#endif
