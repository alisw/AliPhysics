#ifndef ALIANALYSISTASKJETVN_H
#define ALIANALYSISTASKJETVN_H
/**
 * \file AliAnalysisTaskEmcalJetSample.h
 * \brief Declaration of class AliAnalysisTaskEmcalJetSample
 *
 * In this header file the class AliAnalysisTaskEmcalJetSample is declared.
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
 * \class AliAnalysisTaskEmcalJetSample
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
class AliAnalysisTaskJetVn : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskJetVn()                                               ;
  AliAnalysisTaskJetVn(const char *name)                               ;
  virtual ~AliAnalysisTaskJetVn()                                      ;

  void                        UserCreateOutputObjects()                         ;
  void                        Terminate(Option_t *option)                       ;

  static AliAnalysisTaskJetVn* AddTaskJetVn(
      const char *ntracks            = "usedefault",
      const char *nclusters          = "usedefault",
      const char* ncells             = "usedefault",
      const char *suffix             = "");

 protected:
  void                        ExecOnce()                                        ;
  Bool_t                      FillHistograms()                                  ;
  Bool_t                      Run()                                             ;

  void                        AllocateJetHistograms()                           ;
  void                        AllocateTrackHistograms()                         ;

  void                        DoJetLoop()                                       ;
  void                        DoTrackLoop()                                     ;

  THistManager                fHistManager                                      ;///< Histogram manager

 private:
  AliAnalysisTaskJetVn(const AliAnalysisTaskJetVn&)           ; // not implemented
  AliAnalysisTaskJetVn &operator=(const AliAnalysisTaskJetVn&); // not implemented

  // Will's functions
  Double_t                    CalculateEventPlaneVZERO(Int_t n);
  Double_t                    CalculateEventPlaneVZEROA(Int_t n);
  Double_t                    CalculateEventPlaneVZEROC(Int_t n);
  //Double_t                    CalculateEventPlaneTPC(Int_t n);
  //Double_t                    CalculateEventPlaneTPCPositiveEta(Int_t n);
  //Double_t                    CalculateEventPlaneTPCNegativeEta(Int_t n);
  //Double_t                    CalculateEventPlaneResolution(Double_t psiVZERO, Double_t psiTPCPositiveEta, Double_t psiTPCNegativeEta);
  Double_t                    CalculateSubEventPlaneResolution(Int_t n, Double_t psiI, Double_t psiJ);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskJetVn, 1);
  /// \endcond
};
#endif
