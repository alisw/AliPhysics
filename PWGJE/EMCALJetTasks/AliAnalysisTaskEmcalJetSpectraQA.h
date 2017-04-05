/**
 * \file AliAnalysisTaskEmcalJetSpectraQA.h
 * \brief Declaration of class AliAnalysisTaskEmcalJetSpectraQA
 *
 * In this header file the class AliAnalysisTaskEmcalJetSpectraQA is declared.
 * This is a task used to do QA on jet spectra for EMCal jet analysis.
 *
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \date Feb 2, 2016
 */

#ifndef ALIANALYSISTASKEMCALJETSPECTRAQA_H
#define ALIANALYSISTASKEMCALJETSPECTRAQA_H

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "THistManager.h"
#include "AliTLorentzVector.h"
#include "AliAnalysisTaskEmcalJetLight.h"

/**
 * \class AliAnalysisTaskEmcalJetSpectraQA
 * \brief Implementation of a task to perform QA on jet spectra
 *
 * Implementation of a task that performs QA on jet spectra
 * for EMCal jet analysis.
 */
class AliAnalysisTaskEmcalJetSpectraQA : public AliAnalysisTaskEmcalJetLight {
 public:

  /**
   * \class AliEmcalJetInfo
   * \brief Class that encapsulates jets
   *
   * Implementation of a class that encapsulates jets.
   */
  class AliEmcalJetInfo : public AliTLorentzVector {
  public:
    AliEmcalJetInfo();
    AliEmcalJetInfo(const AliEmcalJet& jet);

    Double_t fArea;
    Double_t fMCPt;
    Double_t fNConstituents;
    Double_t fNEF;
    Double_t fCent;
    Double_t fEP;
    Double_t fCorrPt;
    Double_t fZ;
    Double_t fLeadingPt;
  };

  enum EHistoType_t {
    kTH2,
    kTHnSparse,
    kTTree
  };

  AliAnalysisTaskEmcalJetSpectraQA();
  AliAnalysisTaskEmcalJetSpectraQA(const char *name);
  virtual ~AliAnalysisTaskEmcalJetSpectraQA() {;}

  void                        UserCreateOutputObjects();

  void                        SetHistoType(EHistoType_t t)        { fHistoType             = t; }
  void                        SetJetEPaxis(Bool_t b)              { fJetEPaxis             = b; }
  void                        SetAreaAxis(Bool_t b)               { fAreaAxis              = b; }
  void                        SetPtBin(Float_t w, Float_t max)    { fPtBinWidth            = w; fMaxPt = max ; }
  void                        SetIsEmbedded(Bool_t i)             { fIsEmbedded            = i; }

  static AliAnalysisTaskEmcalJetSpectraQA* AddTaskEmcalJetSpectraQA(TString ntracks = "usedefault", TString nclusters = "usedefault", Double_t trackPtCut = 0.15, Double_t clusECut = 0.30, TString suffix = "");

 protected:
  virtual void                AllocateTHX(const AliJetContainer* jets);
  virtual void                AllocateTHnSparse(const AliJetContainer* jets);
  virtual void                AllocateTTree(const AliJetContainer* jets);

  virtual void                FillTHX(const AliEmcalJetInfo& jetInfo, const AliJetContainer* jets);
  virtual void                FillTHnSparse(const AliEmcalJetInfo& jetInfo, const AliJetContainer* jets);
  virtual void                FillTTree(const AliEmcalJetInfo& jetInfo, const AliJetContainer* jets);

  Bool_t                      FillHistograms();
  void                        FillJetHisto(const AliEmcalJetInfo& jetInfo, const AliJetContainer* jets);

  EHistoType_t                fHistoType;                   ///< histogram type
  Bool_t                      fJetEPaxis;                   ///< whether a EP-jet axis should be included in the THnSparse
  Bool_t                      fAreaAxis;                    ///< whether the area axis should be included
  Float_t                     fPtBinWidth;                  ///< Histogram pt bin width
  Float_t                     fMaxPt;                       ///< Histogram pt limit
  Bool_t                      fIsEmbedded;                  ///< Embedded data present
  THistManager                fHistManager;                 ///< Histogram manager

 private:
  AliAnalysisTaskEmcalJetSpectraQA(const AliAnalysisTaskEmcalJetSpectraQA&);            // not implemented
  AliAnalysisTaskEmcalJetSpectraQA &operator=(const AliAnalysisTaskEmcalJetSpectraQA&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalJetSpectraQA, 5);
  /// \endcond
};
#endif
