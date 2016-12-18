#ifndef ALIANALYSISTASKJETSUBSTRUCTURE_H
#define ALIANALYSISTASKJETSUBSTRUCTURE_H
/**
 * \file AliAnalysisTaskJetSubstructure.h
 * \brief Declaration of class AliAnalysisTaskJetSubstructure
 *
 * Very quick analysis task to get a feeling for
 * jet shapes in pp collisions
 *
 * \author Redmer Alexander Bertens <rbertens@cern.ch>
 * \date Sep 26, 2016
 */

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"

class AliFJWrapper;

class AliAnalysisTaskJetSubstructure : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskJetSubstructure()                                               ;
  AliAnalysisTaskJetSubstructure(const char *name)                               ;
  virtual ~AliAnalysisTaskJetSubstructure()                                      ;

  void                        UserCreateOutputObjects()                         ;
  void                        Terminate(Option_t *option)                       ;

 protected:
  void                        ExecOnce()                                        ;
  Bool_t                      FillHistograms()                                  ;
  Bool_t                      Run()                                             ;

  void                        AllocateJetHistograms()                           ;
  void                        AllocateTrackHistograms()                         ;
  void                        AllocateClusterHistograms()                       ;
  void                        AllocateCellHistograms()                          ;
  void                        AllocateJetSubstructureHistograms()               ;

  void                        DoJetLoop()                                       ;
  void                        DoJetSubstructureLoop()                           ;
  void                        DoTrackLoop()                                     ;
  void                        DoClusterLoop()                                   ;
  void                        DoCellLoop()                                      ;

  void                        AnalyzeJets()                                     ;

  THistManager                fHistManager                                      ;///< Histogram manager
  AliFJWrapper                *fAliFJWrapper                                    ;// fastjet wrapper

 private:
  AliAnalysisTaskJetSubstructure(const AliAnalysisTaskJetSubstructure&)           ; // not implemented
  AliAnalysisTaskJetSubstructure &operator=(const AliAnalysisTaskJetSubstructure&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskJetSubstructure, 1);
  /// \endcond
};
#endif
