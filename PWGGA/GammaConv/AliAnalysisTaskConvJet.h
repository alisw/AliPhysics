#ifndef AliAnalysisTaskConvJet_H
#define AliAnalysisTaskConvJet_H
/**
 * \file AliAnalysisTaskConvJet.h
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
 * It performs a simple analysis, producing jet spectra.
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

  Double_t                    fNJets                                            ;
  vector<Double_t>            fVectorJetPt                                      ;
  vector<Double_t>            fVectorJetEta                                     ;
  vector<Double_t>            fVectorJetPhi                                     ;
  vector<Double_t>            fVectorJetR                                       ;

 private:
  AliAnalysisTaskConvJet(const AliAnalysisTaskConvJet&)           ;
  AliAnalysisTaskConvJet &operator=(const AliAnalysisTaskConvJet&);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskConvJet, 1);
  /// \endcond
};
#endif
