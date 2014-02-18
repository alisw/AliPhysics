#ifndef ALIANALYSISTASKEMCALJETMASSRESPONSE_H
#define ALIANALYSISTASKEMCALJETMASSRESPONSE_H

class TH1;
class TH2;
class TH3;
class TH3F;
class THnSparse;
class TClonesArray;
class TArrayI;
class AliAnalysisManager;
class AliJetContainer;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalJetMassResponse : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskEmcalJetMassResponse();
  AliAnalysisTaskEmcalJetMassResponse(const char *name);
  virtual ~AliAnalysisTaskEmcalJetMassResponse();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainerBase(Int_t c)                             { fContainerBase     = c   ; }
  void SetMinFractionShared(Double_t f)                         { fMinFractionShared = f   ; }
  void SetJetMassAverage(Double_t m)                            { fJetMassAvg        = m   ; }

 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  Double_t                            GetJetMass(AliEmcalJet *jet);

  Int_t                               fContainerBase;              // jets to be analyzed
  Double_t                            fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  Double_t                            fJetMassAvg;                 // average jet mass

  TH3F            **fh3PtJet1DeltaPtDeltaM;                        //!pt jet1 vs delta-pt vs delta-M
  TH3F            **fh3PtJet2DeltaPtDeltaM;                        //!pt jet2 vs delta-pt vs delta-M
  TH3F            **fh3PtJet1MJet1MJet2;                           //!pt jet1 vs jet mass jet1 vs jet mass jet2
  TH3F            **fh3PtJet2MJet1MJet2;                           //!pt jet2 vs jet mass jet1 vs jet mass jet2

 private:
  AliAnalysisTaskEmcalJetMassResponse(const AliAnalysisTaskEmcalJetMassResponse&);            // not implemented
  AliAnalysisTaskEmcalJetMassResponse &operator=(const AliAnalysisTaskEmcalJetMassResponse&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetMassResponse, 1)
};
#endif

