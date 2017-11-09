/**
 * @file AliAnalysisTaskJetUE.h
 * @brief Declaration of class AliAnalysisTaskJetUE
 *
 * In this header file the class AliAnalysisTaskJetUE is declared.
 *
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * @date June 26, 2017
 */

/* Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#ifndef ALIANALYSISTASKJETUE_H
#define ALIANALYSISTASKJETUE_H

class TString;
class TF1;
class TH1F;
class TH2F;
class TH3F;
class AliRhoParameter;

#include <map>
#include <list>
#include <string>

#include "AliAnalysisTaskEmcalJetLight.h"


/** \class AliAnalysisTaskJetUE
 * \brief Base class for a task that studies the UE
 *
 * Base class for tasks used to study the underlying event (UE) in jet analysis.
 */
class AliAnalysisTaskJetUE : public AliAnalysisTaskEmcalJetLight {
 public:
  AliAnalysisTaskJetUE();
  AliAnalysisTaskJetUE(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskJetUE() {}

  void             SetBackToBackJetPtFraction(Double_t f)    { fBackToBackJetPtFraction = f    ; }
  void             SetMaxMomentumThridJet(Double_t pt)       { fMaxMomentumThridJet     = pt   ; }
  void             SetPtBin(Float_t w, Float_t max)          { fPtBinWidth              = w    ; fMaxPt = max ; }

 protected:
  virtual void                        CalculateEventProperties();
  void                                SortJets();
  Bool_t                              IsB2BEvent(std::string jetCollName = "Signal");
  Bool_t                              AreJetsOverlapping(AliEmcalJet* jet1, AliEmcalJet* jet2);

  Double_t                            fBackToBackJetPtFraction;       ///< Minimum pt fraction of the back-to-back jet
  Double_t                            fMaxMomentumThridJet;           ///< Maximum pt of any additional jet in the event (other than the back-to-back fraction
  Float_t                             fPtBinWidth;                    ///<  Histogram pt bin width
  Float_t                             fMaxPt;                         ///<  Histogram pt limit

  // Event properties
  Int_t                               fNtracks;                       //!<!number of tracks
  Int_t                               fNclusters;                     //!<!number of clusters
  std::map<std::string, Int_t>        fNjets;                         //!<!number of jets

  std::map<std::string, Double_t>     fTotJetArea;                    //!<!total area covered by jets

  AliVParticle                       *fLeadingParticle;               //!<!leading particle
  AliVCluster                        *fLeadingCluster;                //!<!leading cluster
  std::map<std::string, AliEmcalJet*> fLeadingJet;                    //!<!leading jet

  std::map<std::string, std::list<AliEmcalJet*> >
                                      fSortedJets;                    //!<!jets sorted by momentum

  AliAnalysisTaskJetUE(const AliAnalysisTaskJetUE&);             // not implemented
  AliAnalysisTaskJetUE& operator=(const AliAnalysisTaskJetUE&);  // not implemented
  
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskJetUE, 1);
  /// \endcond
};
#endif
