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
#include "AliAnalysisTaskEmcalJet.h"

/**
 * \class AliAnalysisTaskEmcalJetSpectraQA
 * \brief Implementation of a task to perform QA on jet spectra
 *
 * Implementation of a task that performs QA on jet spectra
 * for EMCal jet analysis.
 */
class AliAnalysisTaskEmcalJetSpectraQA : public AliAnalysisTaskEmcalJet {
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
  THistManager                fHistManager;                 ///< Histogram manager

 private:
  AliAnalysisTaskEmcalJetSpectraQA(const AliAnalysisTaskEmcalJetSpectraQA&);            // not implemented
  AliAnalysisTaskEmcalJetSpectraQA &operator=(const AliAnalysisTaskEmcalJetSpectraQA&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalJetSpectraQA, 4)
  /// \endcond
};
#endif
