#ifndef ALIANALYSISTASKJETSHAPEDERIV_H
#define ALIANALYSISTASKJETSHAPEDERIV_H

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

namespace fastjet {
  class PseudoJet;
  class GenericSubtractor;
}

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskJetShapeDeriv : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskJetShapeDeriv();
  AliAnalysisTaskJetShapeDeriv(const char *name);
  virtual ~AliAnalysisTaskJetShapeDeriv();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetCreateTree(Bool_t b)                                  { fCreateTree        = b   ; }

  void SetJetContainerBase(Int_t c)                             { fContainerBase     = c   ; }
  void SetMinFractionShared(Double_t f)                         { fMinFractionShared = f   ; }
  void SetSingleTrackEmbedding(Bool_t b)                        { fSingleTrackEmb    = b   ; }

 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  AliVParticle*                       GetEmbeddedConstituent(AliEmcalJet *jet);

  Int_t                               fContainerBase;              // jets to be analyzed
  Double_t                            fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  Bool_t                              fSingleTrackEmb;             // single track embedding
  Bool_t                              fCreateTree;                 // create output tree

  TTree           *fTreeJetBkg;                                    //!tree with jet and bkg variables
  TLorentzVector  *fJet1Vec;                                       // jet1(AA) vector  
  TLorentzVector  *fJet2Vec;                                       // jet2(probe) vector
  Float_t         fArea;                                           // area
  Float_t         fAreaPhi;                                        // area phi
  Float_t         fAreaEta;                                        // area eta
  Float_t         fRho;                                            // rho
  Float_t         fRhoM;                                           // rho_m
  Int_t           fNConst;                                         // N constituents in jet1
  Float_t         fM1st;                                           // 1st order subtracted jet mass
  Float_t         fM2nd;                                           // 2nd order subtracted jet mass
  Float_t         fDeriv1st;                                       // 1st derivative
  Float_t         fDeriv2nd;                                       // 2nd derivative
  Int_t           fMatch;                                          // 1: matched to MC jet; 0: no match

  TH2F          **fh2MSubMatch;                                    //! subtracted jet mass vs match index (0: no match; 1:match)
  TH2F          **fh2MSubPtRawAll;                                 //! subtracted jet mass vs subtracted jet pT
  TH2F          **fh2MSubPtRawMatch;                               //! subtracted jet mass vs subtracted jet pT for matched jets
  TH2F          **fh2MSubPtTrue;                                   //! subtracted jet mass vs true jet pT for matched jets
  TH2F          **fh2MTruePtTrue;                                  //! true jet mass vs true jet pT for matched jets
  TH2F          **fh2PtTrueDeltaM;                                 //! true jet pT vs (Msub - Mtrue)
  THnSparse     **fhnMassResponse;                                 //! Msub vs Mtrue vs PtCorr vs PtTrue

  TH2F          **fh2PtTrueSubFacV1;                               //! true pT vs -(rho+rhom)*V1
  TH2F          **fh2PtRawSubFacV1;                                //! raw pT vs -(rho+rhom)*V1
  TH2F          **fh2PtCorrSubFacV1;                               //! subtracted pT vs -(rho+rhom)*V1
  TH2F          **fh2NConstSubFacV1;                               //! N constituents vs -(rho+rhom)*V1
  TH2F          **fh2PtTrueSubFacV2;                               //! true pT vs 0.5(rho+rhom)^2*V2
  TH2F          **fh2PtRawSubFacV2;                                //! raw pT vs 0.5(rho+rhom)^2*V2
  TH2F          **fh2PtCorrSubFacV2;                               //! subtracted pT vs 0.5(rho+rhom)^2*V2
  TH2F          **fh2NConstSubFacV2;                               //! N constituents vs 0.5(rho+rhom)^2*V2

 private:
  AliAnalysisTaskJetShapeDeriv(const AliAnalysisTaskJetShapeDeriv&);            // not implemented
  AliAnalysisTaskJetShapeDeriv &operator=(const AliAnalysisTaskJetShapeDeriv&); // not implemented

  ClassDef(AliAnalysisTaskJetShapeDeriv, 1)
};
#endif



