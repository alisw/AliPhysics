#ifndef ALIJETRESPONSEMAKER_H
#define ALIJETRESPONSEMAKER_H

// $Id$

class TClonesArray;
class TH2;
class THnSparse;
class AliNamedArrayI;

#include "AliEmcalJet.h"
#include "AliAnalysisTaskEmcalJet.h"

class AliJetResponseMaker : public AliAnalysisTaskEmcalJet {
 public:
  AliJetResponseMaker();
  AliJetResponseMaker(const char *name);
  virtual ~AliJetResponseMaker();

  enum MatchingType{
    kNoMatching = 0,
    kGeometrical = 1,
    kMCLabel = 2,
    kSameCollections = 3
  };

  void                        UserCreateOutputObjects();

  void                        SetMatching(MatchingType t, Double_t p1=1, Double_t p2=1)       { fMatching = t; fMatchingPar1 = p1; fMatchingPar2 = p2; }
  void                        SetPtHardBin(Int_t b)                                           { fSelectPtHardBin   = b         ; }
  void                        SetUseCellsToMatch(Bool_t i)                                    { fUseCellsToMatch   = i         ; }
  void                        SetMinJetMCPt(Float_t pt)                                       { fMinJetMCPt        = pt        ; }
  void                        SetHistoType(Int_t b)                                           { fHistoType         = b         ; }
  void                        SetDeltaPtAxis(Int_t b)                                         { fDeltaPtAxis       = b         ; }
  void                        SetDeltaEtaDeltaPhiAxis(Int_t b)                                { fDeltaEtaDeltaPhiAxis= b       ; }
  void                        SetNEFAxis(Int_t b)                                             { fNEFAxis           = b         ; }
  void                        SetZAxis(Int_t b)                                               { fZAxis             = b         ; }
  void                        SetDoJet2Histogram(Int_t b)                                     { fDoJet2Histogram   = b         ; }

 protected:
  void                        ExecOnce();
  void                        DoJetLoop();
  Bool_t                      FillHistograms();
  Bool_t                      Run();
  Bool_t                      DoJetMatching();
  void                        SetMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, MatchingType matching);
  void                        GetGeometricalMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d) const;
  void                        GetMCLabelMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d1, Double_t &d2) const;
  void                        GetSameCollectionsMatchingLevel(AliEmcalJet *jet1, AliEmcalJet *jet2, Double_t &d1, Double_t &d2) const;
  void                        FillMatchingHistos(Double_t Pt1, Double_t Pt2, Double_t Eta1, Double_t Eta2, Double_t Phi1, Double_t Phi2, 
						 Double_t A1, Double_t A2, Double_t d, Double_t CE1, Double_t CE2, Double_t CorrPt1, Double_t CorrPt2, 
						 Double_t MCPt1, Double_t NEF1, Double_t NEF2, Double_t Z1, Double_t Z2);
  void                        FillJetHisto(Double_t Phi, Double_t Eta, Double_t Pt, Double_t A, Double_t NEF, Double_t Z, Double_t CorrPt, Double_t MCPt, Int_t Set);
  void                        AllocateTH2();
  void                        AllocateTHnSparse();

  MatchingType                fMatching;                               // matching type
  Double_t                    fMatchingPar1;                           // matching parameter for jet1-jet2 matching
  Double_t                    fMatchingPar2;                           // matching parameter for jet2-jet1 matching
  Bool_t                      fUseCellsToMatch;                        // use cells instead of clusters to match jets (slower but sometimes needed)
  Double_t                    fMinJetMCPt;                             // minimum jet MC pt
  Int_t                       fHistoType;                              // histogram type (0=TH2, 1=THnSparse)
  Int_t                       fDeltaPtAxis;                            // add delta pt axis in THnSparse (default=0)
  Int_t                       fDeltaEtaDeltaPhiAxis;                   // add delta eta and delta phi axes in THnSparse (default=0)
  Int_t                       fNEFAxis;                                // add NEF axis in matching THnSparse (default=0)
  Int_t                       fZAxis;                                  // add Z axis in matching THnSparse (default=0)
  Int_t                       fDoJet2Histogram;                        // add unbiased jet2 histogram (potentially memory consuming if on particle level)

  Bool_t                      fIsJet1Rho;                              //!whether the jet1 collection has to be average subtracted
  Bool_t                      fIsJet2Rho;                              //!whether the jet2 collection has to be average subtracted

  // General histos
  TH2                        *fHistLeadingJets1PtArea;                 //!leading jet pt vs. area histogram 1
  TH2                        *fHistLeadingJets1CorrPtArea;             //!leading jet corr pt vs. area histogram 1
  TH2                        *fHistLeadingJets2PtArea;                 //!leading jet pt vs. area histogram 2
  TH2                        *fHistLeadingJets2CorrPtArea;             //!leading jet corr pt vs. area histogram 2
  TH2                        *fHistLeadingJets2PtAreaAcceptance;       //!leading jet pt vs. area histogram 2 using jet 1 cuts (acceptance, leading hadron bias, ...)
  TH2                        *fHistLeadingJets2CorrPtAreaAcceptance;   //!leading jet corr pt vs. area histogram 2 using jet 1 cuts (acceptance, leading hadron bias, ...)

  // THnSparse
  THnSparse                  *fHistJets1;                              //!jet1 THnSparse
  THnSparse                  *fHistJets2;                              //!jet2 THnSparse
  THnSparse                  *fHistJets2Acceptance;                    //!jet2 acceptance THnSparse
  THnSparse                  *fHistMatching;                           //!matching THnSparse

  // Jets 1
  TH2                        *fHistJets1PhiEta;                        //!phi-eta distribution of jets 1
  TH2                        *fHistJets1PtArea;                        //!inclusive jet pt vs. area histogram 1
  TH2                        *fHistJets1CorrPtArea;                    //!inclusive jet pt vs. area histogram 1
  TH2                        *fHistJets1NEFvsPt;                       //!Jet neutral energy fraction vs. jet pt 1
  TH2                        *fHistJets1CEFvsCEFPt;                    //!Jet charged energy fraction vs. charged jet pt 1
  TH2                        *fHistJets1ZvsPt;                         //!Constituent Pt over Jet Pt ratio vs. jet pt 1

  // Jets 2
  TH2                        *fHistJets2PhiEta;                        //!phi-eta distribution of jets 2
  TH2                        *fHistJets2PtArea;                        //!inclusive jet pt vs. area histogram 2
  TH2                        *fHistJets2CorrPtArea;                    //!inclusive jet pt vs. area histogram 2
  TH2                        *fHistJets2PhiEtaAcceptance;              //!phi-eta distribution of jets 2 using jet 1 cuts (acceptance, leading hadron bias, ...)
  TH2                        *fHistJets2PtAreaAcceptance;              //!inclusive jet pt vs. area histogram 2 using jet 1 cuts (acceptance, leading hadron bias, ...)
  TH2                        *fHistJets2CorrPtAreaAcceptance;          //!inclusive jet pt vs. area histogram 2 using jet 1 cuts (acceptance, leading hadron bias, ...)
  TH2                        *fHistJets2NEFvsPt;                       //!Jet neutral energy fraction vs. jet pt 2
  TH2                        *fHistJets2CEFvsCEFPt;                    //!Jet charged energy fraction vs. charged jet pt 2
  TH2                        *fHistJets2ZvsPt;                         //!Constituent Pt over Jet Pt ratio vs. jet pt 2

  // Jet1-Jet2 matching
  TH2                        *fHistCommonEnergy1vsJet1Pt;              //!common energy 1 (%) vs jet 1 pt
  TH2                        *fHistCommonEnergy2vsJet2Pt;              //!common energy 2 (%) vs jet 2 pt
  TH2                        *fHistDistancevsJet1Pt;                   //!distance vs jet 1 pt
  TH2                        *fHistDistancevsJet2Pt;                   //!distance vs jet 2 pt
  TH2                        *fHistDistancevsCommonEnergy1;            //!distance vs common energy 1 (%)
  TH2                        *fHistDistancevsCommonEnergy2;            //!distance vs common energy 2 (%)
  TH2                        *fHistCommonEnergy1vsCommonEnergy2;       //!common energy 1 (%) vs common energy 2 (%)
  TH2                        *fHistDeltaEtaDeltaPhi;                   //!delta eta vs delta phi of matched jets

  TH2                        *fHistJet2PtOverJet1PtvsJet2Pt;           //!jet 2 pt over jet 1 pt vs jet 2 pt
  TH2                        *fHistJet1PtOverJet2PtvsJet1Pt;           //!jet 1 pt over jet 2 pt vs jet 1 pt
  TH2                        *fHistDeltaPtvsJet1Pt;                    //!delta pt between matched jets vs jet 1 pt
  TH2                        *fHistDeltaPtvsJet2Pt;                    //!delta pt between matched jets vs jet 2 pt
  TH2                        *fHistDeltaPtOverJet1PtvsJet1Pt;          //!delta pt / jet 1 pt between matched jets vs jet 1 pt
  TH2                        *fHistDeltaPtOverJet2PtvsJet2Pt;          //!delta pt / jet 2 pt between matched jets vs jet 2 pt
  TH2                        *fHistDeltaPtvsDistance;                  //!delta pt between matched jets vs distance
  TH2                        *fHistDeltaPtvsCommonEnergy1;             //!delta pt between matched jets vs common energy 1 (%)
  TH2                        *fHistDeltaPtvsCommonEnergy2;             //!delta pt between matched jets vs common energy 2 (%)
  TH2                        *fHistDeltaPtvsArea1;                     //!delta pt between matched jets vs jet 1 area
  TH2                        *fHistDeltaPtvsArea2;                     //!delta pt between matched jets vs jet 2 area
  TH2                        *fHistDeltaPtvsDeltaArea;                 //!delta pt between matched jets vs delta area
  TH2                        *fHistJet1PtvsJet2Pt;                     //!correlation jet 1 pt vs jet 2 pt

  TH2                        *fHistDeltaCorrPtOverJet1CorrPtvsJet1CorrPt;//!delta pt corr / jet 1 corr pt between matched jets vs jet 1 corr pt
  TH2                        *fHistDeltaCorrPtOverJet2CorrPtvsJet2CorrPt;//!delta pt corr / jet 2 corr pt between matched jets vs jet 2 corr pt
  TH2                        *fHistDeltaCorrPtvsJet1CorrPt;            //!delta pt corr between matched jets vs jet 1 corr pt 
  TH2                        *fHistDeltaCorrPtvsJet2CorrPt;            //!delta pt corr between matched jets vs jet 2 corr pt
  TH2                        *fHistDeltaCorrPtvsDistance;              //!delta pt corr between matched jets vs distance
  TH2                        *fHistDeltaCorrPtvsCommonEnergy1;         //!delta pt corr between matched jets vs common energy 1 (%)
  TH2                        *fHistDeltaCorrPtvsCommonEnergy2;         //!delta pt corr between matched jets vs common energy 2 (%)
  TH2                        *fHistDeltaCorrPtvsArea1;                 //!delta pt corr between matched jets vs jet 1 area
  TH2                        *fHistDeltaCorrPtvsArea2;                 //!delta pt corr between matched jets vs jet 2 area
  TH2                        *fHistDeltaCorrPtvsDeltaArea;             //!delta pt corr between matched jets vs delta area
  TH2                        *fHistJet1CorrPtvsJet2CorrPt;             //!correlation jet 1 corr pt vs jet 2 corr pt

  TH2                        *fHistDeltaMCPtOverJet1MCPtvsJet1MCPt;    //!jet 1 MC pt - jet2 pt / jet 1 MC pt vs jet 1 pt
  TH2                        *fHistDeltaMCPtOverJet2PtvsJet2Pt;        //!jet 1 MC pt - jet2 pt / jet 2 pt vs jet 2 pt
  TH2                        *fHistDeltaMCPtvsJet1MCPt;                //!jet 1 MC pt - jet2 pt vs jet 1 MC pt
  TH2                        *fHistDeltaMCPtvsJet2Pt;                  //!jet 1 MC pt - jet2 pt vs jet 2 pt
  TH2                        *fHistDeltaMCPtvsDistance;                //!jet 1 MC pt - jet2 pt vs distance
  TH2                        *fHistDeltaMCPtvsCommonEnergy1;           //!jet 1 MC pt - jet2 pt vs common energy 1 (%)
  TH2                        *fHistDeltaMCPtvsCommonEnergy2;           //!jet 1 MC pt - jet2 pt vs common energy 2 (%)
  TH2                        *fHistDeltaMCPtvsArea1;                   //!jet 1 MC pt - jet2 pt vs jet 1 area
  TH2                        *fHistDeltaMCPtvsArea2;                   //!jet 1 MC pt - jet2 pt vs jet 2 area
  TH2                        *fHistDeltaMCPtvsDeltaArea;               //!jet 1 MC pt - jet2 pt vs delta area
  TH2                        *fHistJet1MCPtvsJet2Pt;                   //!correlation jet 1 MC pt vs jet 2 pt

 private:
  AliJetResponseMaker(const AliJetResponseMaker&);            // not implemented
  AliJetResponseMaker &operator=(const AliJetResponseMaker&); // not implemented

  ClassDef(AliJetResponseMaker, 23) // Jet response matrix producing task
};
#endif
