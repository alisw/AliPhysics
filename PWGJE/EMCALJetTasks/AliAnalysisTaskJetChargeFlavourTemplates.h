#ifndef ALIANALYSISTASKJETCHARGEFLAVOURTEMPLATES_H
#define ALIANALYSISTASKJETCHARGEFLAVOURTEMPLATES_H


class TTree;
class AliAnalysisManager;

#include "AliAnalysisTaskEmcalJet.h"


class AliAnalysisTaskJetChargeFlavourTemplates : public AliAnalysisTaskEmcalJet {
 public:


  AliAnalysisTaskJetChargeFlavourTemplates();
  AliAnalysisTaskJetChargeFlavourTemplates(const char *name);
  virtual ~AliAnalysisTaskJetChargeFlavourTemplates();

  void                                UserCreateOutputObjects();
  void                                Terminate(Option_t *option);

  //Setters
  void SetJetContainer(Int_t c)            { fContainer = c; }
  void SetJetPtMinThreshold(Float_t f)     { fPtThreshold = f; }
  void SetCentralitySelectionOn(Bool_t t)  { fCentSelectOn = t; }
  void SetMinCentrality(Float_t t)         { fCentMin = t; }
  void SetMaxCentrality(Float_t t)         { fCentMax = t; }
  void SetJetRadius(Double_t t)            { fJetRadius = t; }
  void SetK(Double_t t)                    { JetChargeK = t; }
  void SetMotherFraction(Double_t t)       { MotherFraction = t; }
  void SetJetMidPt(Double_t t)             { JetMidPt = t; }
  void SetJetHighPt(Double_t t)            { JetHighPt = t; }

 protected:
  // Obligatory Functions
  Bool_t                              RetrieveEventObjects();
  Bool_t                              Run();
  Bool_t                              FillHistograms();
  // Essential variables
  // Jet container to be analysed: 0 for raw, 1 for subtracted
  Int_t                               fContainer;                           //
  // Jet lower Pt threshold
  Float_t                             fPtThreshold;                         //
  // Switch on/off centrality selection
  Bool_t                              fCentSelectOn;                        //
  // Min centrality value
  Float_t                             fCentMin;                             //
  // Max centrality value
  Float_t                             fCentMax;                             //
  // Jet radius
  Double_t                            fJetRadius;
  // Value of K the scaling factor on jet charge default k = 0.5
  Double_t                            JetChargeK = 0.5;                     //
  // Set Fraction of the Mother particles which need to agree
  Double_t                            MotherFraction = 0.8;                 //
  // Set the Mid Pt jet threshold Defult = 40
  Double_t                            JetMidPt = 40;                        //
  // Set the Highest jet Pt Band threshold Default = 40
  Double_t                            JetHighPt = 80;                       //


  // Splits Jets into too high and Lower Pt Bins
  // User histograms and output tree
  // Histograms first


  TH1F                                *fhJetPt;                             //!
  TH1F                                *fhJetPhi;                            //!
  TH1F                                *fhJetEta;                            //!

  TH1F                                *JC;                                  //!

  TH1F                                *JCUp;                                //!
  TH1F                                *JCDown;                              //!
  TH1F                                *JCGluon;                             //!
  TH1F                                *JCOther;                             //!
  TH1F                                *JCUnmatched;                         //!

  TH1F                                *JCLow;                               //!

  TH1F                                *JCUpLow;                             //!
  TH1F                                *JCDownLow;                           //!
  TH1F                                *JCGluonLow;                          //!
  TH1F                                *JCOtherLow;                          //!
  TH1F                                *JCUnmatchedLow;                      //!

  TH1F                                *JCMid;                               //!

  TH1F                                *JCUpMid;                             //!
  TH1F                                *JCDownMid;                           //!
  TH1F                                *JCGluonMid;                          //!
  TH1F                                *JCOtherMid;                          //!
  TH1F                                *JCUnmatchedMid;                      //!

  TH1F                                *JCHigh;                              //!

  TH1F                                *JCUpHigh;                            //!
  TH1F                                *JCDownHigh;                          //!
  TH1F                                *JCGluonHigh;                         //!
  TH1F                                *JCOtherHigh;                         //!
  TH1F                                *JCUnmatchedHigh;                     //!

  TH1F                                *fhParticleJetPt;                     //!
  TH1F                                *fhParticleJetPhi;                    //!
  TH1F                                *fhParticleJetEta;                    //!

  TH1F                                *JCParticle;                          //!

  TH1F                                *JCParticleUp;                        //!
  TH1F                                *JCParticleDown;                      //!
  TH1F                                *JCParticleGluon;                     //!
  TH1F                                *JCParticleOther;                     //!
  TH1F                                *JCParticleUnmatched;                 //!

  TH1F                                *JCParticleLow;                       //!

  TH1F                                *JCParticleUpLow;                     //!
  TH1F                                *JCParticleDownLow;                   //!
  TH1F                                *JCParticleGluonLow;                  //!
  TH1F                                *JCParticleOtherLow;                  //!
  TH1F                                *JCParticleUnmatchedLow;              //!

  TH1F                                *JCParticleMid;                       //!

  TH1F                                *JCParticleUpMid;                     //!
  TH1F                                *JCParticleDownMid;                   //!
  TH1F                                *JCParticleGluonMid;                  //!
  TH1F                                *JCParticleOtherMid;                  //!
  TH1F                                *JCParticleUnmatchedMid;              //!

  TH1F                                *JCParticleHigh;                      //!

  TH1F                                *JCParticleUpHigh;                    //!
  TH1F                                *JCParticleDownHigh;                  //!
  TH1F                                *JCParticleGluonHigh;                 //!
  TH1F                                *JCParticleOtherHigh;                 //!
  TH1F                                *JCParticleUnmatchedHigh;             //!

  TH1F                                *PtComparison;                        //!
  TH1F                                *JCComparison;                        //!
  TH1F                                *JCComparisonUp;                      //!
  TH1F                                *JCComparisonDown;                    //!
  TH1F                                *JCComparisonGluon;                   //!
  TH1F                                *JCComparisonOther;                   //!
  TH1F                                *JCComparisonUnmatched;               //!

  TH2F                                *Pt2DCompare;                         //!
  TH2F                                *JC2DCompare;                         //!

  TH2F                                *PtComparisonVsJCDiff;                      //!

  TH2F                                *PtComparisonVsJCDiffUp;                    //!
  TH2F                                *PtComparisonVsJCDiffDown;                  //!
  TH2F                                *PtComparisonVsJCDiffGluon;                 //!
  TH2F                                *PtComparisonVsJCDiffOther;                 //!
  TH2F                                *PtComparisonVsJCDiffUnmatched;             //!


  TH2F                                *JCComparisonVsPtDiff;                        //!

  TH2F                                *JCComparisonVsPtDiffUp;                      //!
  TH2F                                *JCComparisonVsPtDiffDown;                    //!
  TH2F                                *JCComparisonVsPtDiffGluon;                   //!
  TH2F                                *JCComparisonVsPtDiffOther;                   //!
  TH2F                                *JCComparisonVsPtDiffUnmatched;               //!

/*
  TH3F                                *ParticlePtAndJC;                       //!

  TH3F                                *ParticlePtAndJCUp;                     //!
  TH3F                                *ParticlePtAndJCDown;                   //!
  TH3F                                *ParticlePtAndJCGluon;                  //!
  TH3F                                *ParticlePtAndJCOther;                  //!
  TH3F                                *ParticlePtAndJCUnmatched;              //!
*/
  // Here is the TTree
  TTree                               *fTreeJets;                                                  //!
  // These are the branch variables; there are nBranches of them
  static const Int_t nBranchesJetChargeFlavourTemplates = 25;                                      //
  Double_t                            fTreeBranch[nBranchesJetChargeFlavourTemplates];             //
  TChain                              *pChain;                                                     //



  ClassDef(AliAnalysisTaskJetChargeFlavourTemplates, 4)
};
#endif
