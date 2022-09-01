#ifndef ALIANALYSISTASKEPCalibForJet_H
#define ALIANALYSISTASKEPCalibForJet_H
/**
 * \file AliAnalysisTaskEPCalibForJet.h
 * \brief Declaration of class AliAnalysisTaskEPCalibForJet
 *
 * \author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * \date Apr 27, 2016
 */

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"

/**
 * \class AliAnalysisTaskEPCalibForJet
 * \brief Implementation of a 1le jet analysis task.
 */
class AliAnalysisTaskEPCalibForJet : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskEPCalibForJet()                                               ;
  AliAnalysisTaskEPCalibForJet(const char *name)                               ;
  virtual ~AliAnalysisTaskEPCalibForJet()                                      ;

  void                        UserCreateOutputObjects()                         ;
  void                        Terminate(Option_t *option)                       ;

  static AliAnalysisTaskEPCalibForJet* AddTaskEPCalibForJet(
      const char *ntracks            = "usedefault",
      const char *nclusters          = "usedefault",
      const char* ncells             = "usedefault",
      const char *suffix             = "");

 protected:
  void                        ExecOnce()                                        ;
  Bool_t                      FillHistograms()                                  ;
  Bool_t                      Run()                                             ;

  void                        AllocateJetHistograms()                           ;
  void                        DoJetLoop()                                       ;

  THistManager                fHistManager                                      ;///< Histogram manager

 private:
  AliAnalysisTaskEPCalibForJet(const AliAnalysisTaskEPCalibForJet&)           ; // not implemented
  AliAnalysisTaskEPCalibForJet &operator=(const AliAnalysisTaskEPCalibForJet&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEPCalibForJet, 4);
  /// \endcond
};
#endif