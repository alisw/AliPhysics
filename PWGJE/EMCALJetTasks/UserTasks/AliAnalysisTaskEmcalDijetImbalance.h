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
    Dijet_t() : leadingHadronCutType(0), trigJetMinPtType(0), assJetMinPtType(0), trigJet(0), trigJetPt(0), trigJetPhi(0), trigJetEta(0),
                assJet(0), assJetPt(0), assJetPhi(0), assJetEta(0), isAccepted(kFALSE), deltaPhi(0), deltaEta(0), AJ(0), xJ(0), kTy(0) {}
    
    Int_t leadingHadronCutType;
    Int_t trigJetMinPtType;
    Int_t assJetMinPtType;
    
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
    
    void clear() {leadingHadronCutType=0; trigJetMinPtType=0; assJetMinPtType=0; trigJet=0; trigJetPt=0; trigJetPhi=0; trigJetEta=0;
                    assJet=0; assJetPt=0; assJetPhi=0; assJetEta=0; isAccepted=kFALSE; deltaPhi=0; deltaEta=0; AJ=0; xJ=0; kTy=0;}
  };

  AliAnalysisTaskEmcalDijetImbalance()                                          ;
  AliAnalysisTaskEmcalDijetImbalance(const char *name)                          ;
  virtual ~AliAnalysisTaskEmcalDijetImbalance()                                 ;

  void UserCreateOutputObjects()                                                ;
  
  // Setters
  void SetDeltaPhiCut(Double_t d)                           { fDeltaPhiMin = d; }
  void SetMaxPt(Double_t d)                                 { fMaxPt = d; }
  void SetPlotJetHistograms(Bool_t b)                       { fPlotJetHistograms = b; }
  void SetPlotDijetJetHistograms(Bool_t b)                  { fPlotDijetJetHistograms = b; }
  void SetPlotDijetImbalanceHistograms(Bool_t b)            { fPlotDijetImbalanceHistograms = b; }
  void SetDoMomentumBalance(Bool_t b)                       { fDoMomentumBalance = b; }
  void SetDoGeometricalMatching(Bool_t b, Double_t r, Double_t trackThresh, Double_t clusThresh)
    { fDoGeometricalMatching = b; fMatchingJetR = r; fTrackConstituentThreshold = trackThresh; fClusterConstituentThreshold = clusThresh;}
  void SetNDijetPtThresholds(Int_t n)                       { fNDijetPtThresholds = n; }
  void SetMinTrigJetPt(Double_t* arr)                       { fMinTrigJetPt = arr; }
  void SetMinAssJetPt(Double_t* arr)                        { fMinAssJetPt = arr; }
  void SetDijetLeadingHadronPt(Double_t pt)                 { fDijetLeadingHadronPt = pt; }
  void SetUseManualEvtCuts(Bool_t input)                    { fUseManualEventCuts = input;}

 protected:
  void                        ExecOnce()                                        ;
  Bool_t                      FillHistograms()                                  ;
  Bool_t                      Run()                                             ;
  Bool_t                      IsEventSelected()                                 ;

  // Analysis and plotting functions
  void                        GenerateHistoBins()                               ;
  void                        AllocateJetHistograms()                           ;
  void                        AllocateDijetJetHistograms()                      ;
  void                        AllocateDijetImbalanceHistograms()                ;
  void                        AllocateMomentumBalanceHistograms()               ;
  void                        AllocateGeometricalMatchingHistograms()           ;
  void                        FindDijet(AliJetContainer* jetCont, Int_t leadingHadronCutBin, Int_t trigJetMinPtBin, Int_t assJetMinPtBin);
  void                        DoMomentumBalance(TString histname)               ;
  void                        DoGeometricalMatching();
  void                        FindMatchingDijet(AliJetContainer* jetCont, Int_t assJetMinPtBin);
  void                        FillJetHistograms()                               ;
  void                        FillDijetJetHistograms(TString histname, Int_t isAccepted, Int_t IsAssJet, Double_t jetPt, Double_t jetPhi,
                                                     Double_t jetEta, Int_t nTracksJet, Double_t jetArea);
  void                        FillDijetImbalanceHistograms(TString histname)    ;
  void                        FillMomentumBalanceHistograms(TString histname, Double_t deltaPhi, Double_t trackPt, Double_t balancePt);
  void                        FillGeometricalMatchingHistograms();
  
  // Utility functions
  Double_t                    GetJetPt(AliJetContainer* jetCont, AliEmcalJet* jet);
  Double_t                    GetDeltaR(AliEmcalJet* jet1, AliEmcalJet* jet2);
  
  // Analysis parameters
  Double_t                    fDeltaPhiMin;                         ///< minimum delta phi between di-jets
  Int_t                       fNDijetPtThresholds;                  ///< number of pT thresholds on leading/subleading jets
  Double_t*                   fMinTrigJetPt;                        //[fNDijetPtThresholds] array of leading jet min pT's in a dijet pair
  Double_t*                   fMinAssJetPt;                         //[fNDijetPtThresholds] array of subleading jet min pT's in a dijet pair
  Double_t                    fDijetLeadingHadronPt;                ///< leading hadron pT threshold for leading jet in dijet
  Double_t                    fMatchingJetR;                        ///< jet R for matching study
  Double_t                    fTrackConstituentThreshold;           ///< constituent threshold for matching study
  Double_t                    fClusterConstituentThreshold;         ///< constituent threshold for matching study
  Dijet_t                     fDijet;                               //!<! dijet candidate (per event)
  Dijet_t                     fMatchingDijet;                       //!<! low-threshold matching dijet, for matching study

  // Analysis configuration and plotting options
  Bool_t                      fPlotJetHistograms;                   ///< Set whether to enable inclusive jet histograms
  Bool_t                      fPlotDijetJetHistograms;              ///< Set whether to enable dijet pair histograms
  Bool_t                      fPlotDijetImbalanceHistograms;        ///< Set whether to enable dijet imbalance histograms
  Bool_t                      fDoMomentumBalance;                   ///< Set whether to enable momentum balance study
  Bool_t                      fDoGeometricalMatching;               ///< Set whether to enable constituent study with geometrical matching

  // Plotting parameters
  Float_t                     fMaxPt;                               ///< Histogram pt limit
  Int_t                       fNCentHistBins;                       //!<! number of cent bins
  Double_t*                   fCentHistBins;                        //!<! cent bins
  
  // Event selection
  AliEventCuts                fEventCuts;                           ///< event selection utility
  TList                      *fEventCutList;                        //!<! Output list for event cut histograms
  Bool_t                      fUseManualEventCuts;                  ///< Flag to use manual event cuts
  
  // Hist manager
  THistManager                fHistManager;                         ///< Histogram manager

 private:
  AliAnalysisTaskEmcalDijetImbalance(const AliAnalysisTaskEmcalDijetImbalance&)           ; // not implemented
  AliAnalysisTaskEmcalDijetImbalance &operator=(const AliAnalysisTaskEmcalDijetImbalance&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalDijetImbalance, 3);
  /// \endcond
};
#endif
