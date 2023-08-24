#ifndef AliAnalysisTaskJetOutlierRemoval_H
#define AliAnalysisTaskJetOutlierRemoval_H
/**
 * \file AliAnalysisTaskJetOutlierRemoval.h
 *
 * \author Nicolas Schmidt <nicolas.schmidt@cern.ch>, ORNL
 * \date Nov 20, 2019
 */

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"

/**
 * \class AliAnalysisTaskJetOutlierRemoval
 * \brief Implementation of a sample jet analysis task.
 *
 * This class in an implementation of a sample task for EMCal jet analysis.
 * It derives from AliAnalysisTaskEmcalJet.
 * It performs a simple analysis, returning the maximum jet pT from MC.
 */
class AliAnalysisTaskJetOutlierRemoval : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskJetOutlierRemoval()                                               ;
  AliAnalysisTaskJetOutlierRemoval(const char *name)                               ;
  virtual ~AliAnalysisTaskJetOutlierRemoval()                                      ;

  void                        UserCreateOutputObjects()                         ;
  void                        Terminate(Option_t *option)                       ;

  static AliAnalysisTaskJetOutlierRemoval* AddTask_GammaOutlierRemoval();

  Double_t GetMaxJetPt() {return fMaxJetPt;}


  Double_t Get_Jet_Radius(){
      AliJetContainer* jetCont = 0;
      TIter next(&fJetCollArray);
      Double_t radius = -1;
      while ((jetCont = static_cast<AliJetContainer*>(next()))) {
         radius = jetCont->GetJetRadius();
      }
      return radius;
 }

 protected:
  void                        ExecOnce();
  Bool_t                      FillHistograms();
  Bool_t                      Run();

  void                        DoJetLoop();

  Double_t                    fMaxJetPt;                          // Number of reconstructed jets

 private:
  AliAnalysisTaskJetOutlierRemoval(const AliAnalysisTaskJetOutlierRemoval&)           ;
  AliAnalysisTaskJetOutlierRemoval &operator=(const AliAnalysisTaskJetOutlierRemoval&);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskJetOutlierRemoval, 1);
  /// \endcond
};
#endif
