#ifndef AliAnalysisTaskEmcalPIDinJet_H
#define AliAnalysisTaskEmcalPIDinJet_H
/**
 * \file AliAnalysisTaskEmcalPIDinJet.h
 * \brief Declaration of class AliAnalysisTaskEmcalPIDinJet
 *
 * In this header file the class AliAnalysisTaskEmcalPIDinJet is declared.
 * This is a sample task that shows how to write a simple user analysis task
 * using the EMCal jet framework. It is also used to do automatic benchmark
 * tests of the software.
 * AliAnalysisTaskEmcalPIDinJet is based on AliAnalysisTaskEmcalJetSample.
 *
 * \author Yasaki Keisuke <keisuke.yasaki@cern.ch>, Tsukuba University
 * \date Mar 20, 2023

 * Copyright(c) 1998-2023, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                              
**/

#include "AliAnalysisTaskSE.h"
#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"
#include "AliInputEventHandler.h"
#include "AliPIDResponse.h"

class AliPIDResponse;

/**
 * \class AliAnalysisTaskEmcalPIDinJet
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
class AliAnalysisTaskEmcalPIDinJet : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskEmcalPIDinJet()                                               ;
  AliAnalysisTaskEmcalPIDinJet(const char *name)                               ;
  virtual ~AliAnalysisTaskEmcalPIDinJet()                                      ;

  void                        UserCreateOutputObjects()                         ;
  void                        Terminate(Option_t *option)                       ;

  static AliAnalysisTaskEmcalPIDinJet* AddTaskEmcalPIDinJet(
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
  AliInputEventHandler        inputHandler;

 private:
  AliAnalysisTaskEmcalPIDinJet(const AliAnalysisTaskEmcalPIDinJet&)           ; // not implemented
  AliAnalysisTaskEmcalPIDinJet &operator=(const AliAnalysisTaskEmcalPIDinJet&); // not implemented
  AliPIDResponse*         fPIDResponse;
  
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalPIDinJet, 7);
  /// \endcond
};
#endif
