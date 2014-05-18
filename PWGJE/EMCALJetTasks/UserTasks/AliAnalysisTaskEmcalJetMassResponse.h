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
class TTree;
class TLorentzVector;
class AliAnalysisManager;
class AliVParticle;
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
  void SetCreateTree(Bool_t b)                                  { fCreateTree        = b   ; }

  void SetJetContainerBase(Int_t c)                             { fContainerBase     = c   ; }
  void SetMinFractionShared(Double_t f)                         { fMinFractionShared = f   ; }
  void SetJetMassAverageFunc(TF1 *f)                            { f1JetMassAvg       = f   ; }
  void SetSingleTrackEmbedding(Bool_t b)                        { fSingleTrackEmb    = b   ; }
  void SetSubtractMasslessParticleJet(Bool_t b)                 { fSubtractMassless  = b   ; }

 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  TLorentzVector                      GetSubtractedVector(AliEmcalJet *jet);
  TLorentzVector                      GetSubtractedVectorCheat(AliEmcalJet *jet);
  TLorentzVector                      GetBkgVector(AliEmcalJet *jet, AliJetContainer *cont);
  TLorentzVector                      GetBkgVectorCheat(AliEmcalJet *jet);
  Double_t                            GetJetMass(AliEmcalJet *jet);
  Double_t                            GetJetMassMasslessConstituents(AliEmcalJet *jet);
  AliVParticle*                       GetEmbeddedConstituent(AliEmcalJet *jet);

  Int_t                               fContainerBase;              // jets to be analyzed
  Double_t                            fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  TF1                                *f1JetMassAvg;                // parametrization of average jet mass
  Bool_t                              fSingleTrackEmb;             // single track embedding
  Bool_t                              fSubtractMassless;           // subtract mass of jet assuming massless particles
  Bool_t                              fCreateTree;                 // create output tree

  TH2F            **fh2PtJet1DeltaMNoSub;                          //!pt jet1 vs delta-pt vs delta-M
  TH2F            **fh2PtJet2DeltaMNoSub;                          //!pt jet2 vs delta-pt vs delta-M

  TH3F            **fh3PtJet1DeltaPtDeltaMCheat;                   //!pt jet1 vs delta-pt vs delta-M
  TH3F            **fh3PtJet2DeltaPtDeltaMCheat;                   //!pt jet2 vs delta-pt vs delta-M

  TH3F            **fh3PtJet1DeltaPtDeltaM;                        //!pt jet1 vs delta-pt vs delta-M
  TH3F            **fh3PtJet2DeltaPtDeltaM;                        //!pt jet2 vs delta-pt vs delta-M
  TH2F            **fh2PtJet1DeltaE;                               //!pt jet1 vs delta-E
  TH2F            **fh2PtJet2DeltaE;                               //!pt jet2 vs delta-E
  TH2F            **fh2PtJet1DeltaP;                               //!pt jet1 vs delta-P
  TH2F            **fh2PtJet2DeltaP;                               //!pt jet2 vs delta-P
  TH2F            **fh2PtJet2DeltaM;                               //!pt jet2 vs delta-M
  TH3F            **fh3PtJet1MJet1MJet2;                           //!pt jet1 vs jet mass jet1 vs jet mass jet2
  TH3F            **fh3PtJet2MJet1MJet2;                           //!pt jet2 vs jet mass jet1 vs jet mass jet2

  TH3F            **fh3PtJet1DeltaPtDeltaMRho;                     //!pt jet1 vs delta-pt vs delta-M
  TH2F            **fh2PtJet1DeltaERho;                            //!pt jet1 vs delta-E
  TH3F            **fh3PtJet2DeltaPtDeltaMRho;                     //!pt jet2 vs delta-pt vs delta-M
  TH2F            **fh2PtJet2DeltaPxRho;                           //!pt jet2 vs delta-px
  TH2F            **fh2PtJet2DeltaPyRho;                           //!pt jet2 vs delta-py
  TH2F            **fh2PtJet2DeltaPzRho;                           //!pt jet2 vs delta-pz
  TH2F            **fh2PtJet2DeltaERho;                            //!pt jet2 vs delta-E
  TH2F            **fh2PtJet2DeltaMRho;                            //!pt jet2 vs delta-M
  TH2F            **fh2PtJet2DeltaPtRho;                           //!pt jet2 vs delta-pT
  TH3F            **fh3PtJet2DeltaEDeltaMRho;                      //!pt jet2 vs delta-E vs delta-M
  TH3F            **fh3PtJet2DeltaPDeltaMRho;                      //!pt jet2 vs delta-P vs delta-M

  TH2F            **fh2PtJet1DeltaPtVecSub;                        //!pt jet1 (AA) vs delta pT while using vector subtraction

  TTree           *fTreeJetBkg;                                    //!tree with jet and bkg variables
  TLorentzVector  *fJet1Vec;                                       // jet1(AA) vector  
  TLorentzVector  *fJet2Vec;                                       // jet2(probe) vector
  TLorentzVector  *fBkgVec;                                        // bkg vector
  Float_t         fArea;                                           // area
  Float_t         fAreaPhi;                                        // area phi
  Float_t         fAreaEta;                                        // area eta
  Float_t         fRho;                                            // rho
  Float_t         fRhoM;                                           // rho_m
  Int_t           fNConst;                                         // N constituents in jet1
  Float_t         fJetMassMassless;                                // jet mass for massless constituents

 private:
  AliAnalysisTaskEmcalJetMassResponse(const AliAnalysisTaskEmcalJetMassResponse&);            // not implemented
  AliAnalysisTaskEmcalJetMassResponse &operator=(const AliAnalysisTaskEmcalJetMassResponse&); // not implemented

  ClassDef(AliAnalysisTaskEmcalJetMassResponse, 1)
};
#endif

