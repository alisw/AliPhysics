#ifndef ALIANALYSISTASKRAWJETSPEC_H
#define ALIANALYSISTASKRAWJETSPEC_H
/**
 * \file AliAnalysisTaskRawJetSpec.h
 * \brief Declaration of class AliAnalysisTaskRawJetSpec
 *
 * In this header file the class AliAnalysisTaskRawJetSpec is declared.
 * This is a sample task that shows how to write a simple user analysis task
 * using the EMCal jet framework. It is also used to do automatic benchmark
 * tests of the software.
 *
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \date Apr 27, 2016
 */

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"

/**
 * \class AliAnalysisTaskRawJetSpec
 * \brief Implementation of a sample jet analysis task.
 *
 * This class in an implementation of a sample task for EMCal jet analysis.
 * It derives from AliAnalysisTaskEmcalJet.
 * It performs a simple analysis, producing track, cluster and jet spectra.
 * It also performs a QA of the cluster-track matching.
 * Note: if jets are not used this class can be simplified by deriving
 * from AliAnalysisTaskEmcal and removing the functions DoJetLoop()
 * and AllocateJetHistograms().
 */
class AliAnalysisTaskRawJetSpec : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskRawJetSpec();
  AliAnalysisTaskRawJetSpec(const char *name);
  virtual ~AliAnalysisTaskRawJetSpec();

  void                        UserCreateOutputObjects();
  void                        Terminate(Option_t *option);

  static AliAnalysisTaskRawJetSpec* AddTaskRawJetSpec(
      const char *ntracks            = "usedefault",
      const char *nclusters          = "usedefault",
      const char* ncells             = "usedefault",
      const char *suffix             = "");


 protected:
  void                        ExecOnce() ;
  Bool_t                      FillHistograms()                                  ;
  Bool_t                      Run()                                             ;

  void                        AllocateJetHistograms()                           ;
  void                        DoJetLoop()                                       ;

  THistManager                fHistManager                                      ;///< Histogram manager

 private:
  AliAnalysisTaskRawJetSpec(const AliAnalysisTaskRawJetSpec&)           ; // not implemented
  AliAnalysisTaskRawJetSpec &operator=(const AliAnalysisTaskRawJetSpec&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskRawJetSpec, 7);
  /// \endcond
};
#endif
