#ifndef ALIANALYSISTASKEMCALJETMASSRESPONSE_H
#define ALIANALYSISTASKEMCALJETMASSRESPONSE_H

class TH1;
class TH2;
class TH3;
class TH3F;
class THnSparse;
class TF1;
class TLorentzVector;
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
  void SetJetMassAverageFunc(TF1 *f)                            { f1JetMassAvg       = f   ; }

 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  TLorentzVector                      GetSubtractedVector(AliEmcalJet *jet);
  TLorentzVector                      GetSubtractedVectorCheat(AliEmcalJet *jet);
  Double_t                            GetJetMass(AliEmcalJet *jet);

  Int_t                               fContainerBase;              // jets to be analyzed
  Double_t                            fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  TF1                                *f1JetMassAvg;                // parametrization of average jet mass

  TH2F            **fh2PtJet1DeltaMNoSub;                          //!pt jet1 vs delta-pt vs delta-M
  TH2F            **fh2PtJet2DeltaMNoSub;                          //!pt jet2 vs delta-pt vs delta-M

  TH3F            **fh3PtJet1DeltaPtDeltaMCheat;                   //!pt jet1 vs delta-pt vs delta-M
  TH3F            **fh3PtJet2DeltaPtDeltaMCheat;                   //!pt jet2 vs delta-pt vs delta-M

  TH3F            **fh3PtJet1DeltaPtDeltaM;                        //!pt jet1 vs delta-pt vs delta-M
  TH3F            **fh3PtJet2DeltaPtDeltaM;                        //!pt jet2 vs delta-pt vs delta-M
  TH3F            **fh3PtJet1MJet1MJet2;                           //!pt jet1 vs jet mass jet1 vs jet mass jet2
  TH3F            **fh3PtJet2MJet1MJet2;                           //!pt jet2 vs jet mass jet1 vs jet mass jet2

  TH2F            **fh2PtJet1DeltaPtVecSub;                        //!pt jet1 (AA) vs delta pT while using vector subtraction

 private:
  AliAnalysisTaskEmcalJetMassResponse(const AliAnalysisTaskEmcalJetMassResponse&);            // not implemented
  AliAnalysisTaskEmcalJetMassResponse &operator=(const AliAnalysisTaskEmcalJetMassResponse&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetMassResponse, 1)
};
#endif

