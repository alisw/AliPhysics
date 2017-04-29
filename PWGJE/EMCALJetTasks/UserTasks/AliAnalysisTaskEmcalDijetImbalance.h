#ifndef AliAnalysisTaskEmcalDijetImbalance_H
#define AliAnalysisTaskEmcalDijetImbalance_H
/**
 * \file AliAnalysisTaskEmcalDijetImbalance.h
 * \brief Di-jet imbalance analysis with full jets.
 *
 * This header file declares the class AliAnalysisTaskEmcalDijetImbalance.
 * Structure of class based on AliAnalysisTaskEmcalJetSample, AliAnalysisTaskEmcalJetQA.
 *
 * \author James Mulligan <james.mulligan@yale.edu>, Yale University
 * \date Jun 27, 2016
 */

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "THistManager.h"

#include "AliEventCuts.h"
#include "AliAnalysisTaskEmcalJet.h"

class AliAnalysisTaskEmcalDijetImbalance : public AliAnalysisTaskEmcalJet {
 public:
  
  // Struct to store all relevant info of a di-jet pair
  struct Dijet_t {
    Dijet_t() : leadingHadronCutType(0), trigJet(0), trigJetPt(0), trigJetPhi(0), trigJetEta(0),
                assJet(0), assJetPt(0), assJetPhi(0), assJetEta(0), isAccepted(kFALSE), deltaPhi(0), deltaEta(0), AJ(0), xJ(0), kTy(0) {}
    
    Int_t leadingHadronCutType;
    
    AliEmcalJet* trigJet;
    Double_t trigJetPt;
    Double_t trigJetPhi;
    Double_t trigJetEta;
    
    AliEmcalJet* assJet;
    Double_t assJetPt;
    Double_t assJetPhi;
    Double_t assJetEta;
    
    Bool_t isAccepted;
    
    Double_t deltaPhi;
    Double_t deltaEta;
    Double_t AJ;
    Double_t xJ;
    Double_t kTy;
    
    void clear() {leadingHadronCutType=0; trigJet=0; trigJetPt=0; trigJetPhi=0; trigJetEta=0;
                    assJet=0; assJetPt=0; assJetPhi=0; assJetEta=0; isAccepted=kFALSE; deltaPhi=0; deltaEta=0; AJ=0; xJ=0; kTy=0;}
  };

  AliAnalysisTaskEmcalDijetImbalance()                                          ;
  AliAnalysisTaskEmcalDijetImbalance(const char *name)                          ;
  virtual ~AliAnalysisTaskEmcalDijetImbalance()                                 ;

  void UserCreateOutputObjects()                                                ;
  void LoadBackgroundScalingHistogram(const char* path = "alien:///alice/cern.ch/user/j/jmulliga/BackgroundScalingWeights.root", const char* name = "hBackgroundScalingWeights");
  
  // Setters
  void SetDeltaPhiCut(Double_t d)                           { fDeltaPhiMin = d; }
  void SetMaxPt(Double_t d)                                 { fMaxPt = d; }
  void SetPlotJetHistograms(Bool_t b)                       { fPlotJetHistograms = b; }
  void SetPlotDijetCandHistograms(Bool_t b)                 { fPlotDijetCandHistograms = b; }
  void SetPlotDijetImbalanceHistograms(Bool_t b)            { fPlotDijetImbalanceHistograms = b; }
  void SetComputeBackground(Bool_t b)                       { fComputeBackground = b; }
  void SetDoMomentumBalance(Bool_t b)                       { fDoMomentumBalance = b; }
  void SetDoGeometricalMatching(Bool_t b, Double_t r, Double_t trackThresh, Double_t clusThresh)
    { fDoGeometricalMatching = b; fMatchingJetR = r; fTrackConstituentThreshold = trackThresh; fClusterConstituentThreshold = clusThresh;}
  void SetLoadBackgroundScalingWeights(Bool_t b)            { fLoadBackgroundScalingWeights = b; }
  void SetComputeMBDownscaling(Bool_t b)                    { fComputeMBDownscaling = b; }
  void SetMinTrigJetPt(Double_t p)                          { fMinTrigJetPt = p; }
  void SetMinAssJetPt(Double_t p)                           { fMinAssJetPt = p; }
  void SetDijetLeadingHadronPt(Double_t pt)                 { fDijetLeadingHadronPt = pt; }
  void SetUseAliEventCuts(Bool_t b)                         { fUseAliEventCuts = b; }
  void SetUseManualEvtCuts(Bool_t input)                    { fUseManualEventCuts = input;}
  void SetNEtaBins(Int_t n)                                 { fNEtaBins = n; }
  void SetNPhiBins(Int_t n)                                 { fNPhiBins = n; }

 protected:
  void                        ExecOnce()                                        ;
  Bool_t                      FillHistograms()                                  ;
  Bool_t                      Run()                                             ;
  Bool_t                      IsEventSelected()                                 ;
  void                        RunChanged(Int_t run)                             ;

  // Analysis and plotting functions
  void                        GenerateHistoBins()                               ;
  void                        AllocateJetHistograms()                           ;
  void                        AllocateDijetCandHistograms()                     ;
  void                        AllocateDijetImbalanceHistograms()                ;
  void                        AllocateMomentumBalanceHistograms()               ;
  void                        AllocateGeometricalMatchingHistograms()           ;
  void                        FindDijet(AliJetContainer* jetCont, Int_t leadingHadronCutBin);
  void                        ComputeBackground(AliJetContainer* jetCont)       ;
  void                        DoMomentumBalance(TString histname)               ;
  void                        DoGeometricalMatching();
  void                        FindMatchingDijet(AliJetContainer* jetCont)       ;
  void                        FillJetHistograms()                               ;
  void                        FillDijetCandHistograms(TString histname);
  void                        FillDijetImbalanceHistograms(TString histname)    ;
  void                        FillMomentumBalanceHistograms(TString histname, Double_t deltaPhi, Double_t trackPt, Double_t balancePt);
  void                        FillGeometricalMatchingHistograms();
  
  // Utility functions
  Double_t                    GetJetPt(AliJetContainer* jetCont, AliEmcalJet* jet);
  AliEmcalJet*                GetLeadingJet(AliJetContainer* jetCont);
  Double_t                    GetDeltaR(AliEmcalJet* jet1, AliEmcalJet* jet2);
  Double_t                    GetDeltaR(AliTLorentzVector* part, Double_t etaRef, Double_t phiRef);
  
  // Analysis parameters
  Double_t                    fDeltaPhiMin;                         ///< minimum delta phi between di-jets
  Double_t                    fMinTrigJetPt;                        ///< leading jet min pT in a dijet pair
  Double_t                    fMinAssJetPt;                         ///< subleading jet min pT in a dijet pair, for it to be accepted
  Double_t                    fDijetLeadingHadronPt;                ///< leading hadron pT threshold for leading jet in dijet
  Double_t                    fMatchingJetR;                        ///< jet R for matching study
  Double_t                    fTrackConstituentThreshold;           ///< constituent threshold for matching study
  Double_t                    fClusterConstituentThreshold;         ///< constituent threshold for matching study
  Dijet_t                     fDijet;                               //!<! dijet candidate (per event)
  Dijet_t                     fMatchingDijet;                       //!<! low-threshold matching dijet, for matching study
  Int_t                       fNEtaBins;                            ///< Number of eta bins in DCal region (for background/correction)
  Int_t                       fNPhiBins;                            ///< Number of phi bins in DCal region (for background/correction)
  TH1D*                       fBackgroundScalingWeights;            ///< Histogram storing eta-phi weights for full-jet background scale factors
  TH1D*                       fGapJetScalingWeights;                ///< Histogram storing eta-phi weights scaling jets near the gap region

  // Analysis configuration and plotting options
  Bool_t                      fPlotJetHistograms;                   ///< Set whether to enable inclusive jet histograms
  Bool_t                      fPlotDijetCandHistograms;             ///< Set whether to enable dijet pair histograms
  Bool_t                      fPlotDijetImbalanceHistograms;        ///< Set whether to enable dijet imbalance histograms
  Bool_t                      fComputeBackground;                   ///< Set whether to enable study of background
  Bool_t                      fDoMomentumBalance;                   ///< Set whether to enable momentum balance study
  Bool_t                      fDoGeometricalMatching;               ///< Set whether to enable constituent study with geometrical matching
  Bool_t                      fLoadBackgroundScalingWeights;        ///< Flag to load eta-phi weights for full-jet background scale factors
  Bool_t                      fComputeMBDownscaling;                ///< Set whether to compute and plot MB downscaling factors

  // Plotting parameters
  Float_t                     fMaxPt;                               ///< Histogram pt limit
  Int_t                       fNCentHistBins;                       //!<! number of cent bins
  Double_t*                   fCentHistBins;                        //!<! cent bins
  
  // Event selection
  Bool_t                      fUseAliEventCuts;                     ///< Flag to use AliEventCuts (otherwise AliAnalysisTaskEmcal will be used)
  AliEventCuts                fEventCuts;                           ///< event selection utility
  TList                      *fEventCutList;                        //!<! Output list for event cut histograms
  Bool_t                      fUseManualEventCuts;                  ///< Flag to use manual event cuts
  
  // Trigger parameters
  Double_t                    fMBUpscaleFactor;                     //!<! inverse of downscale factor, for MB trigger
  
  // Hist manager
  THistManager                fHistManager;                         ///< Histogram manager

 private:
  AliAnalysisTaskEmcalDijetImbalance(const AliAnalysisTaskEmcalDijetImbalance&)           ; // not implemented
  AliAnalysisTaskEmcalDijetImbalance &operator=(const AliAnalysisTaskEmcalDijetImbalance&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalDijetImbalance, 7);
  /// \endcond
};
#endif
