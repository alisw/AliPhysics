/**
 * @file AliAnalysisTaskRhoDev.h
 * @brief Declaration of class AliAnalysisTaskRhoDev
 *
 * In this header file the class AliAnalysisTaskRhoDev is declared.
 *
 * @author Rosi Reed, Yale University
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * @date June 16, 2017
 */

/* Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#ifndef ALIANALYSISTASKRHODEV_H
#define ALIANALYSISTASKRHODEV_H

#include <utility>

#include "AliAnalysisTaskRhoBaseDev.h"

/** \class AliAnalysisTaskRhoDev
 * \brief Class for a task that calculates the UE
 *
 * Class for a task that calculates the average background
 * coming from the underlying event (UE) in jet analysis.
 * This task calculates the average background as the median
 * of the pt density of kt clusters. More details at: https://arxiv.org/pdf/0707.1378.pdf.
 * If scale function is given the scaled rho will be exported
 * with the name as "fOutRhoName".Apppend("_Scaled").
 * This is a development version. The stable version of this class
 * is AliAnalysisTaskRho.
 */
class AliAnalysisTaskRhoDev : public AliAnalysisTaskRhoBaseDev {

 public:
  AliAnalysisTaskRhoDev();
  AliAnalysisTaskRhoDev(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskRhoDev() {}

  void             UserCreateOutputObjects();

  void             SetExcludeLeadJets(UInt_t n)    { fNExclLeadJets = n    ; }
  void             SetRhoSparse(Bool_t b)          { fRhoSparse     = b    ; }
  void             SetExclJetOverlap(TString n)    { fExclJetOverlap= n    ; }

  static AliAnalysisTaskRhoDev* AddTaskRhoDev(
     TString        nTracks                        = "usedefault",
     Double_t       trackPtCut                     = 0.15,
     TString        nClusters                      = "usedefault",
     Double_t       clusECut                       = 0.30,
     TString        nRho                           = "Rho",
     Double_t       jetradius                      = 0.2,
     UInt_t         acceptance                     = AliEmcalJet::kTPCfid,
     AliJetContainer::EJetType_t jetType           = AliJetContainer::kChargedJet,
     AliJetContainer::ERecoScheme_t rscheme        = AliJetContainer::pt_scheme,
     Bool_t         histo                          = kTRUE,
     TString        suffix                         = ""
  );

 protected:
  void          CalculateRho();
  Bool_t        FillHistograms();
  Bool_t        VerifyContainers();

  std::pair<AliEmcalJet*, AliEmcalJet*>
                GetLeadingJets();

  UInt_t           fNExclLeadJets;                 ///< number of leading jets to be excluded from the median calculation
  Bool_t           fRhoSparse;                     ///< flag to run CMS method as described in https://arxiv.org/abs/1207.2392
  TString          fExclJetOverlap;                ///< name of the jet collection that should be used to reject jets that are considered "signal"

  Double_t         fOccupancyFactor;               //!<!occupancy correction factor for sparse events
  TH2F            *fHistOccCorrvsCent;             //!<!occupancy correction vs. centrality

  AliAnalysisTaskRhoDev(const AliAnalysisTaskRhoDev&);             // not implemented
  AliAnalysisTaskRhoDev& operator=(const AliAnalysisTaskRhoDev&);  // not implemented
  
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskRhoDev, 2);
  /// \endcond
};
#endif
