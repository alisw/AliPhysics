#ifndef ALIANALYSISTASKEMCALSAMPLE_H
#define ALIANALYSISTASKEMCALSAMPLE_H
/**
 * \file AliAnalysisTaskEmcalSample.h
 * \brief Declaration of class AliAnalysisTaskEmcalSample
 *
 * In this header file the class AliAnalysisTaskEmcalSample is declared.
 * This is a sample task that shows how to write a simple user analysis task
 * using the EMCal framework.
 *
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \date Apr 28, 2016
 */

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcal.h"
#include "THistManager.h"

/**
 * \class AliAnalysisTaskEmcalSample
 * \brief Implementation of a sample jet analysis task.
 *
 * This class in an implementation of a sample task for EMCal jet analysis.
 * It derives from AliAnalysisTaskEmcalJet.
 * It performs a simple analysis, producing track and cluster spectra.
 * It also performs a QA of the cluster-track matching.
 */
class AliAnalysisTaskEmcalSample : public AliAnalysisTaskEmcal {
 public:

  AliAnalysisTaskEmcalSample()                                               ;
  AliAnalysisTaskEmcalSample(const char *name)                               ;
  virtual ~AliAnalysisTaskEmcalSample()                                      ;

  void                        UserCreateOutputObjects()                         ;
  void                        Terminate(Option_t *option)                       ;

 protected:
  void                        ExecOnce()                                        ;
  Bool_t                      FillHistograms()                                  ;
  Bool_t                      Run()                                             ;

  void                        AllocateTrackHistograms()                         ;
  void                        AllocateClusterHistograms()                       ;
  void                        AllocateCellHistograms()                          ;

  void                        DoTrackLoop()                                     ;
  void                        DoClusterLoop()                                   ;
  void                        DoCellLoop()                                      ;

  THistManager                fHistManager                                      ;///< Histogram manager

 private:
  AliAnalysisTaskEmcalSample(const AliAnalysisTaskEmcalSample&)           ; // not implemented
  AliAnalysisTaskEmcalSample &operator=(const AliAnalysisTaskEmcalSample&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalSample, 2);
  /// \endcond
};
#endif
