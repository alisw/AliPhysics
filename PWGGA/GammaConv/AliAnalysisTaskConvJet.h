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
  std::vector<Double_t> GetVectorJetPt()  {return fVectorJetPt;}
  std::vector<Double_t> GetVectorJetPx()  {return fVectorJetPx;}
  std::vector<Double_t> GetVectorJetPy()  {return fVectorJetPy;}
  std::vector<Double_t> GetVectorJetPz()  {return fVectorJetPz;}
  std::vector<Double_t> GetVectorJetEta() {return fVectorJetEta;}
  std::vector<Double_t> GetVectorJetPhi() {return fVectorJetPhi;}
  std::vector<Double_t> GetVectorJetArea() {return fVectorJetR;}

  Double_t GetTrueNJets() {return fTrueNJets;}
  std::vector<Double_t> GetTrueVectorJetPt()  {return fTrueVectorJetPt;}
  std::vector<Double_t> GetTrueVectorJetPx()  {return fTrueVectorJetPx;}
  std::vector<Double_t> GetTrueVectorJetPy()  {return fTrueVectorJetPy;}
  std::vector<Double_t> GetTrueVectorJetPz()  {return fTrueVectorJetPz;}
  std::vector<Double_t> GetTrueVectorJetEta() {return fTrueVectorJetEta;}
  std::vector<Double_t> GetTrueVectorJetPhi() {return fTrueVectorJetPhi;}
  std::vector<Double_t> GetTrueVectorJetArea()   {return fTrueVectorJetR;}

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

  Double_t                    fNJets;                          // Number of reconstructed jets
  std::vector<Double_t>       fVectorJetPt;                    // Vector for the pt of the reconstructed jets
  std::vector<Double_t>       fVectorJetPx;                    // Vector for the px of the reconstructed jets
  std::vector<Double_t>       fVectorJetPy;                    // Vector for the py of the reconstructed jets
  std::vector<Double_t>       fVectorJetPz;                    // Vector for the pz of the reconstructed jets
  std::vector<Double_t>       fVectorJetEta;                   // Vector for the eta of the reconstructed jets
  std::vector<Double_t>       fVectorJetPhi;                   // Vector for the phi of the reconstructed jets
  std::vector<Double_t>       fVectorJetR;                     // Vector for the radius of the reconstructed jets

  Double_t                    fTrueNJets;                      // Number of true jets
  std::vector<Double_t>       fTrueVectorJetPt;                // Vector for the pt of the true jets
  std::vector<Double_t>       fTrueVectorJetPx;                // Vector for the px of the true jets
  std::vector<Double_t>       fTrueVectorJetPy;                // Vector for the py of the true jets
  std::vector<Double_t>       fTrueVectorJetPz;                // Vector for the pz of the true jets
  std::vector<Double_t>       fTrueVectorJetEta;               // Vector for the eta of the true jets
  std::vector<Double_t>       fTrueVectorJetPhi;               // Vector for the phi of the true jets
  std::vector<Double_t>       fTrueVectorJetR;                 // Vector for the radius of the true jets

 private:
  AliAnalysisTaskConvJet(const AliAnalysisTaskConvJet&)           ;
  AliAnalysisTaskConvJet &operator=(const AliAnalysisTaskConvJet&);

  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskConvJet, 11);
  /// \endcond
};
#endif
