#ifndef ALIANALYSISTASKJETSHAPECONST_H
#define ALIANALYSISTASKJETSHAPECONST_H

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

class AliAnalysisTaskJetShapeConst : public AliAnalysisTaskEmcalJet {
 public:
  enum JetMassVarType {
    kMass   = 0,  //jet mass
    kRatMPt = 1    //ratio mass/pt jet
  };
  enum ResponseReference {
    kDet    = 0,   //detector level
    kPart   = 1    //particle level
  };

  AliAnalysisTaskJetShapeConst();
  AliAnalysisTaskJetShapeConst(const char *name);
  virtual ~AliAnalysisTaskJetShapeConst();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetCreateTree(Bool_t b)                                  { fCreateTree        = b   ; }

  void SetJetContainerBase(Int_t c)                             { fContainerBase     = c   ; }
  void SetJetContainerSub(Int_t c)                              { fContainerSub      = c   ; }
  void SetJetContainerNoEmb(Int_t c)                            { fContainerNoEmb    = c   ; }
  void SetMinFractionShared(Double_t f)                         { fMinFractionShared = f   ; }
  void SetSingleTrackEmbedding(Bool_t b, Int_t min = 99999, Int_t max = 999999) { fSingleTrackEmb = b; fMinLabelEmb = min; fMaxLabelEmb = max; }
  void SetJetMassVarType(JetMassVarType t)                      { fJetMassVarType    = t   ; }
  void SetResponseReference(ResponseReference r)                { fResponseReference = r   ; }
  void SetUseSumw2(Bool_t b)                                    { fUseSumw2          = b   ; }
  void SetSmallSystRanges(Bool_t small = kTRUE)                 { fSmallSyst = small; }
  
 protected:
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  AliVParticle*                       GetEmbeddedConstituent(AliEmcalJet *jet);

 private:
  Int_t                               fContainerBase;              // jets to be analyzed
  Int_t                               fContainerSub;               // subtracted jets to be analyzed
  Int_t                               fContainerNoEmb;             // subtracted jets from Pb-Pb only events
  Double_t                            fMinFractionShared;          // only fill histos for jets if shared fraction larger than X
  Bool_t                              fSingleTrackEmb;             // single track embedding
  Bool_t                              fCreateTree;                 // create output tree
  JetMassVarType                      fJetMassVarType;             // observable to use
  ResponseReference                   fResponseReference;          // true axis of response matrix
  Bool_t                              fUseSumw2;                   // activate sumw2 for output histograms

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
  THnSparse     **fhnMassResponse;                                 //! Msub vs Mtrue vs PtCorr vs PtTrue vs DR
  THnSparse     **fhnDeltaMass;                                    //! deltaM vs deltapT
  THnSparse      *fhnDeltaMassAndBkgInfo;                          //! DeltaM, DeltapT bkg-unsubtracted M and pT, rho and rhom 

  AliAnalysisTaskJetShapeConst(const AliAnalysisTaskJetShapeConst&);            // not implemented
  AliAnalysisTaskJetShapeConst &operator=(const AliAnalysisTaskJetShapeConst&); // not implemented

  ClassDef(AliAnalysisTaskJetShapeConst, 9)
};
#endif



