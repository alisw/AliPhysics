#ifndef ALIJETRESPONSEMAKER_H
#define ALIJETRESPONSEMAKER_H
/**************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

//-----------------------------------------------------------------------
// Author : Salvatore Aiola, Yale University, salvatore.aiola@cern.ch
//-----------------------------------------------------------------------

class TClonesArray;
class TH2;
class THnSparse;
class AliNamedArrayI;

#include "AliEmcalJet.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "AliEmcalEmbeddingQA.h"

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
  void                        SetFlavourZAxis(Int_t b)                                        { fFlavourZAxis      = b         ; }
  void                        SetFlavourPtAxis(Int_t b)                                       { fFlavourPtAxis     = b         ; }
  void                        SetZgAxis(Int_t b)                                              { fZgAxis            = b         ; }
  void                        SetdRAxis(Int_t b)                                              { fdRAxis            = b         ; }
  void                        SetPtgAxis(Int_t b)                                             { fPtgAxis           = b         ; }
  void                        SetDBCAxis(Int_t b)                                             { fDBCAxis           = b         ; }
  void                        SetJetRelativeEPAngleAxis(Int_t b)                              { fJetRelativeEPAngle = b        ; }

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
  void                        FillMatchingHistos(AliEmcalJet* jet1, AliEmcalJet* jet2, Double_t d, Double_t CE1, Double_t CE2);
  void                        FillJetHisto(AliEmcalJet* jet, Int_t Set);
  void                        AllocateTH2();
  void                        AllocateTHnSparse();
  Double_t                    GetRelativeEPAngle(Double_t jetAngle, Double_t epAngle) const;

  MatchingType                fMatching;                               // matching type
  Double_t                    fMatchingPar1;                           // matching parameter for jet1-jet2 matching
  Double_t                    fMatchingPar2;                           // matching parameter for jet2-jet1 matching
  Bool_t                      fUseCellsToMatch;                        // use cells instead of clusters to match jets (slower but sometimes needed)
  Double_t                    fMinJetMCPt;                             // minimum jet MC pt
  AliEmcalEmbeddingQA         fEmbeddingQA;                            //!<! Embedding QA hists (will only be added if embedding)
  Int_t                       fHistoType;                              // histogram type (0=TH2, 1=THnSparse)
  Int_t                       fDeltaPtAxis;                            // add delta pt axis in THnSparse (default=0)
  Int_t                       fDeltaEtaDeltaPhiAxis;                   // add delta eta and delta phi axes in THnSparse (default=0)
  Int_t                       fNEFAxis;                                // add NEF axis in matching THnSparse (default=0)
  Int_t                       fZAxis;                                  // add Z axis in matching THnSparse (default=0)
  Int_t                       fFlavourZAxis;                           // add flavour Z axis in matching THnSparse (default=0)
  Int_t                       fFlavourPtAxis;                          // add flavour pt axis in matching THnSparse (default=0)
  Int_t                       fZgAxis;                                 // add Zg axis in matching THnSparse (default=0)
  Int_t                       fdRAxis;                                 // add dR axis in matching THnSparse (default=0)
  Int_t                       fPtgAxis;                                // add Ptg axis in matching THnSparse (default=0)
  Int_t                       fDBCAxis;                                // add DBC (number of soft dropped branches) axis in matching THnSparse (default=0)
  Int_t                       fJetRelativeEPAngle;                     ///< add jet angle relative to the EP in matching THnSparse (default=0)

  Bool_t                      fIsJet1Rho;                              //!whether the jet1 collection has to be average subtracted
  Bool_t                      fIsJet2Rho;                              //!whether the jet2 collection has to be average subtracted

  TH2                        *fHistRejectionReason1;                   //!Rejection reason vs. jet pt
  TH2                        *fHistRejectionReason2;                   //!Rejection reason vs. jet pt

  // THnSparse
  THnSparse                  *fHistJets1;                              //!jet1 THnSparse
  THnSparse                  *fHistJets2;                              //!jet2 THnSparse
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

  ClassDef(AliJetResponseMaker, 29) // Jet response matrix producing task
};
#endif
