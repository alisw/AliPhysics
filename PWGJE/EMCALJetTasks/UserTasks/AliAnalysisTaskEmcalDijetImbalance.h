#ifndef AliAnalysisTaskEmcalDijetImbalance_H
#define AliAnalysisTaskEmcalDijetImbalance_H
/**
 * \file AliAnalysisTaskEmcalDijetImbalance.h
 * \brief Declaration of class AliAnalysisTaskEmcalDijetImbalance
 *
 * In this header file the class AliAnalysisTaskEmcalDijetImbalance is declared.
 * Based on AliAnalysisTaskEmcalJetSample. 
 *
 * \author James Mulligan <james.mulligan@yale.edu>, Yale University
 * \date Jun 27, 2016
 */

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"

/**
 * \class AliAnalysisTaskEmcalDijetImbalance
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
class AliAnalysisTaskEmcalDijetImbalance : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskEmcalDijetImbalance()                                               ;
  AliAnalysisTaskEmcalDijetImbalance(const char *name)                               ;
  virtual ~AliAnalysisTaskEmcalDijetImbalance()                                      ;

  void                        UserCreateOutputObjects()                         ;
  void                        Terminate(Option_t *option)                       ;
  
  void SetDeltaPhiCut(Double_t d)         { fDeltaPhiMin = d; }
  void SetTrigJetMinPt(Double_t d)        { fTrigJetMinPt = d; }
  void SetAssJetMinPt(Double_t d)         { fAssJetMinPt = d; }

 protected:
  void                        ExecOnce()                                        ;
  Bool_t                      FillHistograms()                                  ;
  Bool_t                      Run()                                             ;

  void                        AllocateJetHistograms()                           ;
  void                        AllocateTrackHistograms()                         ;
  void                        AllocateClusterHistograms()                       ;
  void                        AllocateCellHistograms()                          ;
  void                        AllocateDijetHistograms()                         ;

  void                        DoJetLoop()                                       ;
  void                        DoTrackLoop()                                     ;
  void                        DoClusterLoop()                                   ;
  void                        DoCellLoop()                                      ;
  
  Double_t                    fDeltaPhiMin;                   //! minimum delta phi between di-jets
  Double_t                    fTrigJetMinPt;                  //! Pt threshold for trigger (full) jet
  Double_t                    fAssJetMinPt;                   //! Pt threshold for associated (charged) jet (note: unscaled)

  THistManager                fHistManager                                      ;///< Histogram manager

 private:
  AliAnalysisTaskEmcalDijetImbalance(const AliAnalysisTaskEmcalDijetImbalance&)           ; // not implemented
  AliAnalysisTaskEmcalDijetImbalance &operator=(const AliAnalysisTaskEmcalDijetImbalance&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskEmcalDijetImbalance, 1);
  /// \endcond
};
#endif
