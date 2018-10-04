#ifndef AliAnalysisTaskConvJet_H
#define AliAnalysisTaskConvJet_H
/**
 * \file AliAnalysisTaskConvJet.h
 * \brief Declaration of class AliAnalysisTaskConvJet
 *
 * In this header file the class AliAnalysisTaskConvJet is declared.
 * This is a sample task that shows how to write a simple user analysis task
 * using the EMCal jet framework. It is also used to do automatic benchmark
 * tests of the software.
 *
 * \author Lizette Lamers <lizette.jacqueline.lamers@cern.ch>, Utrecht University
 * \author Mike Sas <mike.sas@cern.ch>, Utrecht University
 * \date Okt 4, 2018
 */

/* Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include "AliAnalysisTaskEmcalJet.h"
#include "THistManager.h"

/**
 * \class AliAnalysisTaskConvJet
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
class AliAnalysisTaskConvJet : public AliAnalysisTaskEmcalJet {
 public:

  AliAnalysisTaskConvJet()                                               ;
  AliAnalysisTaskConvJet(const char *name)                               ;
  virtual ~AliAnalysisTaskConvJet()                                      ;

  void                        UserCreateOutputObjects()                         ;
  void                        Terminate(Option_t *option)                       ;

  static AliAnalysisTaskConvJet* AddTask_GammaConvJet(
      const char *ntracks            = "usedefault",
      const char *nclusters          = "usedefault",
      const char* ncells             = "usedefault",
      const char *suffix             = "");

  Double_t GetNJets() {return fNJets;}
  vector<Double_t> GetVectorJetPt()  {return fVectorJetPt;}
  vector<Double_t> GetVectorJetEta() {return fVectorJetEta;}
  vector<Double_t> GetVectorJetPhi() {return fVectorJetPhi;}
  vector<Double_t> GetVectorJetR()   {return fVectorJetR;}

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
  void                        ExecOnce()                                        ;
  Bool_t                      FillHistograms()                                  ;
  Bool_t                      Run()                                             ;

  void                        DoJetLoop()                                       ;

  Double_t                    fNJets                                            ; //
  vector<Double_t>            fVectorJetPt                                      ;
  vector<Double_t>            fVectorJetEta                                     ;
  vector<Double_t>            fVectorJetPhi                                     ;
  vector<Double_t>            fVectorJetR                                       ;

 private:
  AliAnalysisTaskConvJet(const AliAnalysisTaskConvJet&)           ; // not implemented
  AliAnalysisTaskConvJet &operator=(const AliAnalysisTaskConvJet&); // not implemented

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskConvJet, 7);
  /// \endcond
};
#endif
