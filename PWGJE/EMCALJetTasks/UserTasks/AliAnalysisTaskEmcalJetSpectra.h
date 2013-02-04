#ifndef AliAnalysisTaskEmcalJetSpectra_h
#define AliAnalysisTaskEmcalJetSpectra_h

// $Id$


class TH1F;
class TH2F;
class THnSparse;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalJetSpectra : public AliAnalysisTaskEmcalJet {
 public:
  AliAnalysisTaskEmcalJetSpectra();
  AliAnalysisTaskEmcalJetSpectra(const char *name);
  virtual ~AliAnalysisTaskEmcalJetSpectra() {}
  
  
  virtual void           UserCreateOutputObjects();

 protected:
  Bool_t                 Run();
  virtual Int_t          GetCentBin(Double_t cent) const;
  Float_t                RelativePhi(Double_t mphi,Double_t vphi) const;

 private:
  TH2F                  *fHistRhovsCent; //!
  TH2F                  *fHistNjetvsCent;          //!number of jets versus Centrality
  TH2F                  *fHistJetPtvsTrackPt[6];//!
  TH2F                  *fHistRawJetPtvsTrackPt[6];//!
  TH1F                  *fHistTrackPt[6];//!
  TH1F                  *fHistEP0[6];//!
  TH1F                  *fHistEP0A[6];//!
  TH1F                  *fHistEP0C[6];//!
  TH2F                  *fHistEPAvsC[6];//!
  TH2F                  *fHistJetPtvsdEP[6];//!
  TH2F                  *fHistJetPtvsdEPBias[6];//!
  TH2F                  *fHistJetPtvsEP[6];//!
  TH2F                  *fHistJetPtvsEPBias[6];//!
  TH2F                  *fHistRhovsEP[6]; //!




  AliAnalysisTaskEmcalJetSpectra(const AliAnalysisTaskEmcalJetSpectra&); // not implemented
  AliAnalysisTaskEmcalJetSpectra& operator=(const AliAnalysisTaskEmcalJetSpectra&); // not implemented
  
  ClassDef(AliAnalysisTaskEmcalJetSpectra, 4); // Emcal jet spectra task
};
#endif
