#ifndef ALIANALYSISTASKJETSHAPEDERIV_H
#define ALIANALYSISTASKJETSHAPEDERIV_H

class TH1;
class TH2;
class TH3;
class THnSparse;
class TF1;
class TArrayI;
class TTree;
class TLorentzVector;
class AliAnalysisManager;
class AliVParticle;
class AliJetContainer;

#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskJetShapeDeriv : public AliAnalysisTaskEmcalJet {
 public:
  enum JetMassVarType {
    kMass   = 0,   //jet mass
    kRatMPt = 1    //ratio mass/pt jet
  };
  enum ResponseReference {
    kDet    = 0,   //detector level
    kPart   = 1    //particle level
  };

  AliAnalysisTaskJetShapeDeriv();
  AliAnalysisTaskJetShapeDeriv(const char *name);
  virtual ~AliAnalysisTaskJetShapeDeriv();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetCreateTree(Bool_t b)                                  { fCreateTree        = b   ; }

  void SetJetContainerBase(Int_t c)                             { fContainerBase     = c   ; }
  void SetJetContainerNoEmb(Int_t c)                            { fContainerNoEmb    = c   ; }
  void SetMinFractionShared(Double_t f)                         { fMinFractionShared = f   ; }
  void SetSingleTrackEmbedding(Bool_t b, Int_t min = 99999, Int_t max = 999999) { fSingleTrackEmb = b; fMinLabelEmb = min; fMaxLabelEmb = max; }
  void SetJetMassVarType(JetMassVarType t)                      { fJetMassVarType    = t   ; }
  void SetResponseReference(ResponseReference r)                { fResponseReference = r   ; }

  void SetUseSumw2(Bool_t b)                                    { fUseSumw2          = b   ; }
  void SetPartialExclusion(Bool_t b)                            { fPartialExclusion  = b   ; }

  void SetSmallSystRanges(Bool_t small = kTRUE)                 { fSmallSyst = small; }

 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  AliVParticle*                       GetEmbeddedConstituent(AliEmcalJet *jet);

 private:
  Int_t                               fContainerBase;              // jets to be analyzed
  Int_t                               fContainerNoEmb;             // subtracted jets from Pb-Pb only events
  Double_t                            fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  Bool_t                              fSingleTrackEmb;             // single track embedding
  Bool_t                              fCreateTree;                 // create output tree
  JetMassVarType                      fJetMassVarType;             // observable to use
  ResponseReference                   fResponseReference;          // true axis of response matrix
  Bool_t                              fUseSumw2;                   // activate sumw2 for output histograms
  Bool_t                              fPartialExclusion;           // randomly esclude areas according to Ncoll
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
  Int_t           fMinLabelEmb;                                    // min label of embedded particles
  Int_t           fMaxLabelEmb;                                    // max label of embedded particles
  Bool_t          fSmallSyst;                                      // flag for the axes ranges in pPb
  TH2F          **fh2MSubMatch;                                    //! subtracted jet mass vs match index (0: no match; 1:match)
  TH2F          **fh2MSubPtRawAll;                                 //! subtracted jet mass vs subtracted jet pT
  TH3F          **fh3MSubPtRawDRMatch;                             //! subtracted jet mass vs subtracted jet pT vs distance to leading Pb-Pb jet
  TH3F          **fh3MSubPtTrueLeadPt;                             //! subtracted jet mass vs true jet pT vs LeadPt for matched jets for matched jets 
  TH3F          **fh3MTruePtTrueLeadPt;                            //! true jet mass vs true jet pT vs LeadPt for matched jets for matched jets
  TH3F          **fh3PtTrueDeltaMLeadPt;                           //! true jet pT vs (Msub - Mtrue) vs LeadPt for matched jets
  TH3F          **fh3PtTrueDeltaMRelLeadPt;                        //! true jet pT vs (Msub - Mtrue)/Mtrue vs LeadPt for matched jets
  THnSparse     **fhnMassResponse;                                 //! Msub vs Mtrue vs PtCorr vs PtTrue
  THnSparse     **fhnDeltaMass;                                    //! DeltaM, DeltapT 
  THnSparse      *fhnDeltaMassAndBkgInfo;                          //! DeltaM, DeltapT bkg-unsubtracted M and pT, rho and rhom 

  TH2F          **fh2PtTrueSubFacV1;                               //! true pT vs -(rho+rhom)*V1
  TH2F          **fh2PtRawSubFacV1;                                //! raw pT vs -(rho+rhom)*V1
  TH2F          **fh2PtCorrSubFacV1;                               //! subtracted pT vs -(rho+rhom)*V1
  TH2F          **fh2NConstSubFacV1;                               //! N constituents vs -(rho+rhom)*V1
  TH2F          **fh2PtTrueSubFacV2;                               //! true pT vs 0.5(rho+rhom)^2*V2
  TH2F          **fh2PtRawSubFacV2;                                //! raw pT vs 0.5(rho+rhom)^2*V2
  TH2F          **fh2PtCorrSubFacV2;                               //! subtracted pT vs 0.5(rho+rhom)^2*V2
  TH2F          **fh2NConstSubFacV2;                               //! N constituents vs 0.5(rho+rhom)^2*V2

  AliAnalysisTaskJetShapeDeriv(const AliAnalysisTaskJetShapeDeriv&);            // not implemented
  AliAnalysisTaskJetShapeDeriv &operator=(const AliAnalysisTaskJetShapeDeriv&); // not implemented

  ClassDef(AliAnalysisTaskJetShapeDeriv, 9)
};
#endif



