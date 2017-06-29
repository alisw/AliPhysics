/**
 * \file AliAnalysisTaskJetUEStudies.h
 * \brief Declaration of class AliAnalysisTaskJetUE
 *
 * In this header file the class AliAnalysisTaskJetUE is declared.
 * This is a task used to study the underlying event in jet analysis.
 *
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \date June 26, 2017
 */

#ifndef ALIANALYSISTASKJETUESTUDIES_H
#define ALIANALYSISTASKJETUESTUDIES_H

/* Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <map>
#include "THistManager.h"
#include "AliAnalysisTaskJetUE.h"

class TRandom;
class AliEmcalJet;

/**
 * \class AliAnalysisTaskJetUE
 * \brief Implementation of a task used to study the underlying event in jet analysis.
 *
 * Implementation of a task used to study the underlying event in jet analysis.
 */
class AliAnalysisTaskJetUEStudies : public AliAnalysisTaskJetUE {
 public:

  AliAnalysisTaskJetUEStudies();
  AliAnalysisTaskJetUEStudies(const char *name);
  virtual ~AliAnalysisTaskJetUEStudies() { delete fRandom; }

  void           UserCreateOutputObjects();

  void           AddAltRho(TString rhoName)    { fAlternativeRho.insert(std::make_pair(rhoName, nullptr));  }

  static AliAnalysisTaskJetUEStudies* AddTaskJetUEStudies(
      TString ntracks = "usedefault",
      TString nclusters = "usedefault",
      Double_t trackPtCut = 0.15,
      Double_t clusECut = 0.30,
      TString suffix = "");

 protected:
  void           ExecOnce();
  Bool_t         Run();
  Bool_t         FillHistograms();
  AliEmcalJet*   GetJetCone(Double_t radius, Double_t eta, Double_t phi);
  AliEmcalJet*   GetRandomCone(AliJetContainer* jetCont);
  AliEmcalJet*   GetRandomConePerp(AliJetContainer* jetCont, AliEmcalJet* leadJet);
  AliEmcalJet*   GetRandomConeExclLead(AliJetContainer* jetCont, AliEmcalJet* leadJet);

  template <class T, Int_t MAX_CONSTITUENTS>
  Int_t SumParticles(Double_t& pt, Double_t eta, Double_t phi, Double_t maxD2, std::map<std::string, T*>& CollArray, std::array<Int_t, MAX_CONSTITUENTS>& ConstList);

  std::map<TString, AliRhoParameter*>
                              fAlternativeRho;              ///< Alternative average background estimations
  THistManager                fHistManager;                 ///< Histogram manager

  std::set<TString>           fDefaultRhoNames;             //!<! Default rho names
  TRandom                    *fRandom;                      //!<! Random number generator

 private:
  AliAnalysisTaskJetUEStudies(const AliAnalysisTaskJetUEStudies&);            // not implemented
  AliAnalysisTaskJetUEStudies &operator=(const AliAnalysisTaskJetUEStudies&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskJetUEStudies, 1);
  /// \endcond
};
#endif
