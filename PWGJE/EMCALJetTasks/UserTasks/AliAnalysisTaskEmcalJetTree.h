/**
 * \file AliAnalysisTaskEmcalJetTree.h
 * \brief Declaration of class AliAnalysisTaskEmcalJetTree
 *
 * In this header file the class AliAnalysisTaskEmcalJetTree is declared.
 * This is a task used to create a tree with all found jets.
 *
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \date Apr 2, 2016
 */

#ifndef ALIANALYSISTASKEMCALJETTREE_H
#define ALIANALYSISTASKEMCALJETTREE_H

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <map>
#include <vector>
#include <string>
#include "AliAnalysisTaskEmcalJetSpectraQA.h"

/**
 * \class AliAnalysisTaskEmcalJetTreeBase
 * \brief Pure virtual base class for AliAnalysisTaskEmcalJetTree<T>
 *
 * This pure virtual class provides a basic interface for AliAnalysisTaskEmcalJetTree<T>
 */
class AliAnalysisTaskEmcalJetTreeBase : public AliAnalysisTaskEmcalJetSpectraQA {
 public:
  enum EAnalisysType_t {
    kJetPP,
    kJetPbPb,
    kJetEmbedding,
    kJetPPCharged,
    kJetPbPbCharged
  };

  /**
   * \class AliEmcalJetEventInfoSummaryPP
   * \brief Class that encapsulates event properties in a very compact structure (pp analysis)
   *
   * Class that encapsulates event properties in a very compact structure
   * for pp analysis. For the moment, a dummy class (0 bits)
   */
  class AliEmcalJetEventInfoSummaryPP {
  public:
    AliEmcalJetEventInfoSummaryPP() {;}
    AliEmcalJetEventInfoSummaryPP(Double_t, Double_t) {;}

    virtual ~AliEmcalJetEventInfoSummaryPP() {;}

    void Reset() {;}
    void Set(Double_t, Double_t) {;}

    /// \cond CLASSIMP
    ClassDef(AliEmcalJetEventInfoSummaryPP, 1);
    /// \endcond
  };

  /**
   * \class AliEmcalJetEventInfoSummaryPbPb
   * \brief Class that encapsulates event properties in a very compact structure (PbPb analysis)
   *
   * Class that encapsulates event properties in a very compact structure
   * for PbPb analysis (16 bits)
   */
  class AliEmcalJetEventInfoSummaryPbPb {
  public:
    AliEmcalJetEventInfoSummaryPbPb() : fCent(0), fEP(0) {;}
    AliEmcalJetEventInfoSummaryPbPb(Double_t cent, Double_t ep);

    virtual ~AliEmcalJetEventInfoSummaryPbPb() {;}

    void Reset();
    void Set(Double_t cent, Double_t ep);

    /// Eta of the jet
    Double32_t  fCent          ; //[-10,118,7]
    /// Phi of the jet
    Double32_t  fEP            ; //[0,2*pi,9]

    /// \cond CLASSIMP
    ClassDef(AliEmcalJetEventInfoSummaryPbPb, 2);
    /// \endcond
  };

  /**
   * \class AliEmcalJetInfoSummaryPP
   * \brief Class that encapsulates jets in a very compact structure (pp analysis)
   *
   * Class that encapsulates jets in a very compact structure
   * for pp analysis (52 bits)
   */
  class AliEmcalJetInfoSummaryPP {
  public:
    AliEmcalJetInfoSummaryPP() : fPt(0), fEta(0), fPhi(0), fNEF(0), fLeadingPt(0), fNConstituents(0)  {;}
    AliEmcalJetInfoSummaryPP(const AliEmcalJetInfo& source);

    virtual ~AliEmcalJetInfoSummaryPP() {;}

    void Reset();
    void Set(const AliEmcalJetInfo& source);

    /// Transverse momentum of the jet in GeV/c
    Double32_t  fPt            ; //[0,512,10]
    /// Eta of the jet
    Double32_t  fEta           ; //[-1.28,1.28,8]
    /// Phi of the jet
    Double32_t  fPhi           ; //[0,2*pi,9]
    // Fraction of neutral energy
    Double32_t  fNEF           ; //[0,1.28,7]
    // Momentum of leading particle
    Double32_t  fLeadingPt     ; //[0,512,10]
    /// Number of constituents
    Double32_t  fNConstituents ; //[0,256,8]

    /// \cond CLASSIMP
    ClassDef(AliEmcalJetInfoSummaryPP, 4);
    /// \endcond
  };

  /**
   * \class AliEmcalJetInfoSummaryPPCharged
   * \brief Class that encapsulates jets in a very compact structure (pp analysis)
   *
   * Class that encapsulates jets in a very compact structure
   * for pp analysis (45 bits)
   */
  class AliEmcalJetInfoSummaryPPCharged {
  public:
    AliEmcalJetInfoSummaryPPCharged() : fPt(0), fEta(0), fPhi(0), fLeadingPt(0), fNConstituents(0)  {;}
    AliEmcalJetInfoSummaryPPCharged(const AliEmcalJetInfo& source);

    virtual ~AliEmcalJetInfoSummaryPPCharged() {;}

    void Reset();
    void Set(const AliEmcalJetInfo& source);

    /// Transverse momentum of the jet in GeV/c
    Double32_t  fPt            ; //[0,512,10]
    /// Eta of the jet
    Double32_t  fEta           ; //[-1.28,1.28,8]
    /// Phi of the jet
    Double32_t  fPhi           ; //[0,2*pi,9]
    // Momentum of leading particle
    Double32_t  fLeadingPt     ; //[0,512,10]
    /// Number of constituents
    Double32_t  fNConstituents ; //[0,256,8]

    /// \cond CLASSIMP
    ClassDef(AliEmcalJetInfoSummaryPPCharged, 1);
    /// \endcond
  };

  /**
   * \class AliEmcalJetInfoSummaryPbPb
   * \brief Class that encapsulates jets in a very compact structure (Pb-Pb analysis)
   *
   * Class that encapsulates jets in a very compact structure
   * for Pb-Pb analysis (59 bits)
   */
  class AliEmcalJetInfoSummaryPbPb {
  public:
    AliEmcalJetInfoSummaryPbPb() : fPt(0), fEta(0), fPhi(0), fNEF(0), fLeadingPt(0), fArea(0), fNConstituents(0), fCorrPt(0) {;}
    AliEmcalJetInfoSummaryPbPb(const AliEmcalJetInfo& source);

    virtual ~AliEmcalJetInfoSummaryPbPb() {;}

    void Reset();
    void Set(const AliEmcalJetInfo& source);

    Double32_t  fPt            ; //[0,512,9]
    /// Eta of the jet
    Double32_t  fEta           ; //[-1.28,1.28,5]
    /// Phi of the jet
    Double32_t  fPhi           ; //[0,2*pi,7]
    // Fraction of neutral energy
    Double32_t  fNEF           ; //[0,1.28,5]
    // Momentum of leading particle
    Double32_t  fLeadingPt     ; //[0,512,9]
    /// Jet area
    Double32_t  fArea          ; //[0,1.28,6]
    /// Number of constituents
    Double32_t  fNConstituents ; //[0,1024,9]
    /// Jet corrected pt
    Double32_t  fCorrPt        ; //[-256,256,9]

    /// \cond CLASSIMP
    ClassDef(AliEmcalJetInfoSummaryPbPb, 5);
    /// \endcond
  };

  /**
   * \class AliEmcalJetInfoSummaryPbPbCharged
   * \brief Class that encapsulates jets in a very compact structure (Pb-Pb analysis)
   *
   * Class that encapsulates jets in a very compact structure
   * for Pb-Pb analysis (54 bits)
   */
  class AliEmcalJetInfoSummaryPbPbCharged {
  public:
    AliEmcalJetInfoSummaryPbPbCharged() : fPt(0), fEta(0), fPhi(0), fLeadingPt(0), fArea(0), fNConstituents(0), fCorrPt(0) {;}
    AliEmcalJetInfoSummaryPbPbCharged(const AliEmcalJetInfo& source);

    virtual ~AliEmcalJetInfoSummaryPbPbCharged() {;}

    void Reset();
    void Set(const AliEmcalJetInfo& source);

    Double32_t  fPt            ; //[0,512,9]
    /// Eta of the jet
    Double32_t  fEta           ; //[-1.28,1.28,5]
    /// Phi of the jet
    Double32_t  fPhi           ; //[0,2*pi,7]
    // Momentum of leading particle
    Double32_t  fLeadingPt     ; //[0,512,9]
    /// Jet area
    Double32_t  fArea          ; //[0,1.28,6]
    /// Number of constituents
    Double32_t  fNConstituents ; //[0,1024,9]
    /// Jet corrected pt
    Double32_t  fCorrPt        ; //[-256,256,9]

    /// \cond CLASSIMP
    ClassDef(AliEmcalJetInfoSummaryPbPbCharged, 2);
    /// \endcond
  };

  /**
   * \class AliEmcalJetInfoSummaryEmbedding
   * \brief Class that encapsulates jets in a very compact structure (embedding analysis)
   *
   * Class that encapsulates jets in a very compact structure
   * for embedding analysis (82 bits)
   */
  class AliEmcalJetInfoSummaryEmbedding {
  public:
    AliEmcalJetInfoSummaryEmbedding() : fPt(0), fEta(0), fPhi(0), fNEF(0), fLeadingPt(0), fArea(0), fNConstituents(0), fCorrPt(0), fMCPt(0) {;}
    AliEmcalJetInfoSummaryEmbedding(const AliEmcalJetInfo& source);

    virtual ~AliEmcalJetInfoSummaryEmbedding() {;}

    void Reset();
    void Set(const AliEmcalJetInfo& source);

    /// Transverse momentum of the jet in GeV/c
    Double32_t  fPt            ; //[0,512,10]
    /// Eta of the jet
    Double32_t  fEta           ; //[-1.28,1.28,8]
    /// Phi of the jet
    Double32_t  fPhi           ; //[0,2*pi,9]
    // Fraction of neutral energy
    Double32_t  fNEF           ; //[0,1.28,7]
    // Momentum of leading particle
    Double32_t  fLeadingPt     ; //[0,512,9]
    /// Jet area
    Double32_t  fArea          ; //[0,5.12,9]
    /// Number of constituents
    Double32_t  fNConstituents ; //[0,1024,9]
    /// Jet corrected pt
    Double32_t  fCorrPt        ; //[-512,512,11]
    /// Jet MC pt
    Double32_t  fMCPt          ; //[0,512,10]

    /// \cond CLASSIMP
    ClassDef(AliEmcalJetInfoSummaryEmbedding, 4);
    /// \endcond
  };

  AliAnalysisTaskEmcalJetTreeBase() : AliAnalysisTaskEmcalJetSpectraQA(), fTree(0) {;}
  AliAnalysisTaskEmcalJetTreeBase(const char *name);
  virtual ~AliAnalysisTaskEmcalJetTreeBase() {;}

  static AliAnalysisTaskEmcalJetTreeBase* CreateInstance(const char* name, EAnalisysType_t type = kJetPP);
  static AliAnalysisTaskEmcalJetTreeBase* AddTaskEmcalJetTree(TString ntracks = "usedefault", TString nclusters = "usedefault", Double_t trackPtCut = 0.15, Double_t clusECut = 0.30, EAnalisysType_t type = kJetPP, TString suffix = "");

 protected:
  void AllocateTTree(const AliJetContainer* jets) = 0;
  void FillTTree(const AliEmcalJetInfo& jetInfo, const AliJetContainer* jets) = 0;

  TTree   *fTree;    //!<! Output tree

private:
  AliAnalysisTaskEmcalJetTreeBase(const AliAnalysisTaskEmcalJetTreeBase&);            // not implemented
  AliAnalysisTaskEmcalJetTreeBase &operator=(const AliAnalysisTaskEmcalJetTreeBase&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalJetTreeBase, 1)
  /// \endcond
};

/**
 * \class AliAnalysisTaskEmcalJetTree
 * \brief Implementation of a task to generate a tree with all jets
 *
 * Implementation of a task that generates a tree with all jets
 * for EMCal jet analysis.
 */
template <class T, class U>
class AliAnalysisTaskEmcalJetTree : public AliAnalysisTaskEmcalJetTreeBase {
 public:

  AliAnalysisTaskEmcalJetTree();
  AliAnalysisTaskEmcalJetTree(const char *name);

  virtual ~AliAnalysisTaskEmcalJetTree() {;}

  void                        UserCreateOutputObjects();
  Bool_t                      FillHistograms();

 protected:
  void                        AllocateTTree(const AliJetContainer* jets);
  void                        FillTTree(const AliEmcalJetInfo& jetInfo, const AliJetContainer* jets);

  std::map<std::string, std::vector<T> > *fCurrentOutput;  //!<! This map contains the pointers of the tree branch objects
  U                                      *fCurrentEvent ;  //!<! This pointer holds the event object branch

 private:
  AliAnalysisTaskEmcalJetTree(const AliAnalysisTaskEmcalJetTree&);            // not implemented
  AliAnalysisTaskEmcalJetTree &operator=(const AliAnalysisTaskEmcalJetTree&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalJetTree, 1)
  /// \endcond
};

#endif
