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

  /**
   * \class AliEmcalJetInfoSummaryBase
   * \brief Class that encapsulates jets in a very compact structure
   *
   * Class that encapsulates jets in a very compact structure (46 bits)
   */
  class AliEmcalJetInfoSummaryBase {
  public:
    AliEmcalJetInfoSummaryBase() : fPt(0), fEta(0), fPhi(0), fNEF(0), fZLeading(0) {;}
    AliEmcalJetInfoSummaryBase(const AliEmcalJetInfo& source);

    virtual void Reset();
    virtual void Set(const AliEmcalJetInfo& source);

    /// Transverse momentum of the jet in GeV/c
    Double32_t  fPt        ; //[0,200,12]
    /// Eta of the jet
    Double32_t  fEta       ; //[-2,2,10]
    /// Phi of the jet
    Double32_t  fPhi       ; //[0,2*pi,10]
    // Fraction of neutral energy
    Double32_t  fNEF       ; //[0,1,6]
    // Momentum fraction of the leading particle
    Double32_t  fZLeading  ; //[0,1,8]

    /// \cond CLASSIMP
    ClassDef(AliEmcalJetInfoSummaryBase, 1);
    /// \endcond
  };

  /**
   * \class AliEmcalJetInfoSummaryPP
   * \brief Class that encapsulates jets in a very compact structure (pp analysis)
   *
   * Class that encapsulates jets in a very compact structure
   * for pp analysis (54 bits)
   */
  class AliEmcalJetInfoSummaryPP : public AliEmcalJetInfoSummaryBase {
  public:
    AliEmcalJetInfoSummaryPP() : AliEmcalJetInfoSummaryBase(), fNConstituents(0)  {;}
    AliEmcalJetInfoSummaryPP(const AliEmcalJetInfo& source);

    virtual void Reset();
    virtual void Set(const AliEmcalJetInfo& source);

    Char_t fNConstituents; ///< Number of constituents

    /// \cond CLASSIMP
    ClassDef(AliEmcalJetInfoSummaryPP, 1);
    /// \endcond
  };

  /**
   * \class AliEmcalJetInfoSummaryPbPb
   * \brief Class that encapsulates jets in a very compact structure (Pb-Pb analysis)
   *
   * Class that encapsulates jets in a very compact structure
   * for Pb-Pb analysis (96 bits)
   */
  class AliEmcalJetInfoSummaryPbPb : public AliEmcalJetInfoSummaryBase {
  public:
    AliEmcalJetInfoSummaryPbPb() : AliEmcalJetInfoSummaryBase(), fCent(0), fEP(0), fArea(0), fNConstituents(0), fCorrPt(0) {;}
    AliEmcalJetInfoSummaryPbPb(const AliEmcalJetInfo& source);

    virtual void Reset();
    virtual void Set(const AliEmcalJetInfo& source);

    Char_t      fCent          ; ///< Centrality
    /// Angle between the jet axis and the event plane
    Double32_t  fEP            ; //[0,pi,10]
    /// Jet area
    Double32_t  fArea          ; //[0,10,11]
    /// Number of constituents
    Double32_t  fNConstituents ; //[0,1000,9]
    /// Jet corrected pt
    Double32_t  fCorrPt        ; //[0,200,12]

    /// \cond CLASSIMP
    ClassDef(AliEmcalJetInfoSummaryPbPb, 1);
    /// \endcond
  };

  /**
   * \class AliEmcalJetInfoSummaryEmbedding
   * \brief Class that encapsulates jets in a very compact structure (embedding analysis)
   *
   * Class that encapsulates jets in a very compact structure
   * for embedding analysis (108 bits)
   */
  class AliEmcalJetInfoSummaryEmbedding : public AliEmcalJetInfoSummaryPbPb {
  public:
    AliEmcalJetInfoSummaryEmbedding() : AliEmcalJetInfoSummaryPbPb(), fMCPt(0) {;}
    AliEmcalJetInfoSummaryEmbedding(const AliEmcalJetInfo& source);

    virtual void Reset();
    virtual void Set(const AliEmcalJetInfo& source);

    Double32_t  fMCPt      ; //[0,200,12]

    /// \cond CLASSIMP
    ClassDef(AliEmcalJetInfoSummaryEmbedding, 1);
    /// \endcond
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
  void                        AllocateTHX(const AliJetContainer* jets);
  void                        AllocateTHnSparse(const AliJetContainer* jets);
  void                        AllocateTTree(const AliJetContainer* jets);

  void                        FillTHX(const AliEmcalJetInfo& jetInfo, const AliJetContainer* jets);
  void                        FillTHnSparse(const AliEmcalJetInfo& jetInfo, const AliJetContainer* jets);
  void                        FillTTree(const AliEmcalJetInfo& jetInfo, const AliJetContainer* jets);

  Bool_t                      FillHistograms();
  void                        FillJetHisto(const AliEmcalJetInfo& jetInfo, const AliJetContainer* jets);

  EHistoType_t                fHistoType;                   ///< histogram type
  Bool_t                      fJetEPaxis;                   ///< whether a EP-jet axis should be included in the THnSparse
  Bool_t                      fAreaAxis;                    ///< whether the area axis should be included
  THistManager                fHistManager;                 ///< Histogram manager

  AliEmcalJetInfoSummaryBase *fCurrentJetInfo;              //!<! Poitner to current jet info object attached to the tree branch
  THashList                  *fTreeList;                    //!<! Output tree list

 private:
  AliAnalysisTaskEmcalJetSpectraQA(const AliAnalysisTaskEmcalJetSpectraQA&);            // not implemented
  AliAnalysisTaskEmcalJetSpectraQA &operator=(const AliAnalysisTaskEmcalJetSpectraQA&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalJetSpectraQA, 3)
  /// \endcond
};
#endif
