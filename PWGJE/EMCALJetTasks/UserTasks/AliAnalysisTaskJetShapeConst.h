#ifndef ALIANALYSISTASKJETSHAPECONST_H
#define ALIANALYSISTASKJETSHAPECONST_H

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

class AliAnalysisTaskJetShapeConst : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskJetShapeConst();
  AliAnalysisTaskJetShapeConst(const char *name);
  virtual ~AliAnalysisTaskJetShapeConst();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetCreateTree(Bool_t b)                                  { fCreateTree        = b   ; }

  void SetJetContainerBase(Int_t c)                             { fContainerBase     = c   ; }
  void SetJetContainerSub(Int_t c)                              { fContainerSub      = c   ; }
  void SetMinFractionShared(Double_t f)                         { fMinFractionShared = f   ; }
  void SetSingleTrackEmbedding(Bool_t b)                        { fSingleTrackEmb    = b   ; }

 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  AliVParticle*                       GetEmbeddedConstituent(AliEmcalJet *jet);

  Int_t                               fContainerBase;              // jets to be analyzed
  Int_t                               fContainerSub;               // subtracted jets to be analyzed
  Double_t                            fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  Bool_t                              fSingleTrackEmb;             // single track embedding
  Bool_t                              fCreateTree;                 // create output tree

  TTree           *fTreeJetBkg;                                    //!tree with jet and bkg variables
  TLorentzVector  *fJet1Vec;                                       // jet1(AA) vector  
  TLorentzVector  *fJet2Vec;                                       // jet2(probe) vector
  TLorentzVector  *fJetSubVec;                                     // subtracted AA jet vector
  Float_t         fArea;                                           // area
  Float_t         fAreaPhi;                                        // area phi
  Float_t         fAreaEta;                                        // area eta
  Float_t         fRho;                                            // rho
  Float_t         fRhoM;                                           // rho_m
  Int_t           fNConst;                                         // N constituents in jet1
  Int_t           fMatch;                                          // 1: matched to MC jet; 0: no match

  TH2F          **fh2MSubMatch;                                    //! subtracted jet mass vs match index (0: no match; 1:match)
  TH2F          **fh2MSubPtRawAll;                                 //! subtracted jet mass vs subtracted jet pT
  TH2F          **fh2MSubPtRawMatch;                               //! subtracted jet mass vs subtracted jet pT for matched jets
  TH2F          **fh2MSubPtTrue;                                   //! subtracted jet mass vs true jet pT for matched jets
  TH2F          **fh2MTruePtTrue;                                  //! true jet mass vs true jet pT for matched jets
  TH2F          **fh2PtTrueDeltaM;                                 //! true jet pT vs (Msub - Mtrue)
  TH2F          **fh2PtTrueDeltaMRel;                              //! true jet pT vs (Msub - Mtrue)/Mtrue
  THnSparse     **fhnMassResponse;                                 //! Msub vs Mtrue vs PtCorr vs PtTrue

 private:
  AliAnalysisTaskJetShapeConst(const AliAnalysisTaskJetShapeConst&);            // not implemented
  AliAnalysisTaskJetShapeConst &operator=(const AliAnalysisTaskJetShapeConst&); // not implemented

  ClassDef(AliAnalysisTaskJetShapeConst, 1)
};
#endif



