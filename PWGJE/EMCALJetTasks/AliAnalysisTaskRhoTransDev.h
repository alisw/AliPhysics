/**
 * @file AliAnalysisTaskRhoTransDev.h
 * @brief Declaration of class AliAnalysisTaskRhoTransDev
 *
 * In this header file the class AliAnalysisTaskRhoTransDev is declared.
 *
 * @author Rosi Reed, Yale University
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * @date June 25, 2017
 */

/* Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#ifndef ALIANALYSISTASKRHOTRANSDEV_H
#define ALIANALYSISTASKRHOTRANSDEV_H

#include <map>

#include "AliAnalysisTaskRhoBaseDev.h"

/** \class AliAnalysisTaskRhoTransDev
 * \brief Class for a task that calculates the UE
 *
 * Class for a task that calculates the average background
 * coming from the underlying event (UE) in jet analysis.
 * In each event, tha task finds the leading jet and sum
 * the particle pT in a cone of R = 0.4 perpendicular
 * to the leading jet axis. In addition another histogram is filled
 * with the UE estimated in events that fulfill the back-to-back
 * condition: a back-to-back jet, satisfying |∆φ| > 5/6π and
 * carrying at least a given fraction (default: 60%)
 * of the transverse momentum of the leading jet, is required
 * to be present in the event. Additionally, all the other jets in the event
 * should have pT smaller than a certain threshold (default: 12 GeV/c).
 * This implements the techniques described in the analysis note
 * of the ALICE inclusive jet spectrum in pp collisions at 2.76 TeV,
 * section 7.3: https://aliceinfo.cern.ch/Notes/sites/aliceinfo.cern.ch.Notes/files/notes/analysis/rma/2013-Mar-29-analysis_note-ppJet_note.pdf
 * If scale function is given the scaled rho will be exported
 * with the name as "fOutRhoName".Apppend("_Scaled").
 */
class AliAnalysisTaskRhoTransDev : public AliAnalysisTaskRhoBaseDev {

 public:
  AliAnalysisTaskRhoTransDev();
  AliAnalysisTaskRhoTransDev(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskRhoTransDev() {}

  void             UserCreateOutputObjects();

  static AliAnalysisTaskRhoTransDev* AddTaskRhoTransDev(
     TString        nTracks                        = "usedefault",
     Double_t       trackPtCut                     = 0.15,
     TString        nClusters                      = "usedefault",
     Double_t       clusECut                       = 0.30,
     TString        nRho                           = "Rho",
     Double_t       jetradius                      = 0.2,
     UInt_t         acceptance                     = AliEmcalJet::kTPCfid,
     AliJetContainer::EJetType_t jetType           = AliJetContainer::kChargedJet,
     AliJetContainer::ERecoScheme_t rscheme        = AliJetContainer::pt_scheme,
     Bool_t         histo                          = kTRUE,
     TString        suffix                         = ""
  );

 protected:
  void          CalculateRho();
  Bool_t        FillHistograms();
  Bool_t        VerifyContainers();

  Double_t      GetPerpPtDensity(AliEmcalContainer* cont, AliVParticle* leadingJet);

  TH2                                *fHistB2BRhoVsCent;                 //!<!rho vs. centrality

  std::map<std::string, TH2*>         fHistB2BRhoVsLeadJetPt;            //!<!rho vs. leading jet pt
  TH2                                *fHistB2BRhoVsLeadTrackPt;          //!<!rho vs. leading track pt
  TH2                                *fHistB2BRhoVsNtrack;               //!<!rho vs. no. of tracks
  TH2                                *fHistB2BRhoVsLeadClusterE;         //!<!rho vs. leading cluster energy
  TH2                                *fHistB2BRhoVsNcluster;             //!<!rho vs. no. of clusters
  TH2                                *fHistB2BRhoScaledVsCent;           //!<!rhoscaled vs. centrality
  TH2                                *fHistB2BRhoScaledVsNtrack;         //!<!rhoscaled vs. no. of tracks
  TH2                                *fHistB2BRhoScaledVsNcluster;       //!<!rhoscaled vs. no. of clusters

  AliAnalysisTaskRhoTransDev(const AliAnalysisTaskRhoTransDev&);             // not implemented
  AliAnalysisTaskRhoTransDev& operator=(const AliAnalysisTaskRhoTransDev&);  // not implemented
  
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskRhoTransDev, 1);
  /// \endcond
};
#endif
