/**
 * @file AliAnalysisTaskRhoBaseDev.h
 * @brief Declaration of class AliAnalysisTaskRhoBaseDev
 *
 * In this header file the class AliAnalysisTaskRhoBaseDev is declared.
 *
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * @date June 16, 2017
 */

/* Copyright(c) 1998-2017, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#ifndef ALIANALYSISTASKRHOBASEDEV_H
#define ALIANALYSISTASKRHOBASEDEV_H

class TString;
class TF1;
class TH1F;
class TH2F;
class TH3F;
class AliRhoParameter;

#include <map>
#include <string>

#include "AliAnalysisTaskJetUE.h"


/** \class AliAnalysisTaskRhoBaseDev
 * \brief Base class for a task that calculates the UE
 *
 * Base class for a task that calculates the average background
 * coming from the underlying event (UE) in jet analysis.
 * This is a development version. The stable version of this class
 * is AliAnalysisTaskRhoBase.
 */
class AliAnalysisTaskRhoBaseDev : public AliAnalysisTaskJetUE {
 public:
  AliAnalysisTaskRhoBaseDev();
  AliAnalysisTaskRhoBaseDev(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskRhoBaseDev() {}

  void                   UserCreateOutputObjects();

  void                   SetOutRhoName(const char *name)                       { fOutRhoName           = name    ;
                                                                                 fOutRhoScaledName     = Form("%s_Scaled",name);     }
  void                   SetScaleFunction(TF1* sf)                             { fScaleFunction        = sf      ;                   }
  void                   SetRhoFunction(TF1* rf)                               { fRhoFunction          = rf      ;                   }
  TF1*                   LoadRhoFunction(const char* path, const char* name);
  void                   SetAttachToEvent(Bool_t a)                            { fAttachToEvent        = a       ;                   }

  const char*            GetOutRhoName() const                                 { return fOutRhoName.Data()       ;                   }
  const char*            GetOutRhoScaledName() const                           { return fOutRhoScaledName.Data() ;                   }

  static AliAnalysisTaskRhoBaseDev* AddTaskRhoBaseDev(
     TString        nTracks                        = "usedefault",
     TString        nClusters                      = "usedefault",
     TString        nRho                           = "Rho",
     Double_t       jetradius                      = 0.2,
     UInt_t         acceptance                     = AliEmcalJet::kTPCfid,
     AliJetContainer::EJetType_t jetType           = AliJetContainer::kChargedJet,
     AliJetContainer::ERecoScheme_t rscheme        = AliJetContainer::pt_scheme,
     Bool_t         histo                          = kTRUE,
     TString        suffix                         = ""
  );

 protected:
  void                                ExecOnce();
  Bool_t                              Run();
  Bool_t                              FillHistograms();

  virtual Bool_t                      VerifyContainers() { return kTRUE; }

  virtual void                        CalculateRho();

  virtual Double_t                    GetRhoFactor(Double_t cent);
  virtual Double_t                    GetScaleFactor(Double_t cent);

  TString                             fOutRhoName;                    ///< name of output rho object
  TString                             fOutRhoScaledName;              ///< name of output scaled rho object
  TF1                                *fRhoFunction;                   ///< pre-computed rho as a function of centrality
  TF1                                *fScaleFunction;                 ///< pre-computed scale factor as a function of centrality
  Bool_t                              fAttachToEvent;                 ///< whether or not attach rho to the event objects list

  Bool_t                              fTaskConfigured;                //!<!kTRUE if the task is properly configured

  // Exported background density
  AliRhoParameter                    *fOutRho;                        //!<!output rho object
  AliRhoParameter                    *fOutRhoScaled;                  //!<!output scaled rho object

  // Histograms
  TH2                                *fHistRhoVsCent;                 //!<!rho vs. centrality

  std::map<std::string, TH2*>         fHistRhoVsLeadJetPt;            //!<!rho vs. leading jet pt
  std::map<std::string, TH2*>         fHistLeadJetPtVsCent;           //!<!leading jet pt vs. centrality
  std::map<std::string, TH2*>         fHistLeadJetPtDensityVsCent;    //!<!leading jet area vs. centrality
  std::map<std::string, TH2*>         fHistTotJetAreaVsCent;          //!<!total area covered by jets vs. centrality
  std::map<std::string, TH2*>         fHistLeadJetNconstVsCent;       //!<!leading jet constituents vs. cent
  std::map<std::string, TH2**>        fHistLeadJetNconstVsPt;         //!<!leading jet constituents vs. pt
  std::map<std::string, TH2*>         fHistNjetVsCent;                //!<!no. of jets vs. centrality
  std::map<std::string, TH2*>         fHistNjetVsNtrack;              //!<!no. of jets vs. no. of tracks

  TH2                                *fHistRhoVsLeadTrackPt;          //!<!rho vs. leading track pt
  TH2                                *fHistRhoVsNtrack;               //!<!rho vs. no. of tracks
  TH2                                *fHistLeadTrackPtVsCent;         //!<!leading track pt vs. centrality
  TH2                                *fHistNtrackVsCent;              //!<!no. of tracks vs. centrality

  TH2                                *fHistRhoVsLeadClusterE;         //!<!rho vs. leading cluster energy
  TH2                                *fHistRhoVsNcluster;             //!<!rho vs. no. of clusters
  TH2                                *fHistLeadClusterEVsCent;        //!<!leading cluster energy vs. centrality
  TH2                                *fHistNclusterVsCent;            //!<!no. of cluster vs. centrality

  TH2                                *fHistRhoScaledVsCent;           //!<!rhoscaled vs. centrality
  TH2                                *fHistRhoScaledVsNtrack;         //!<!rhoscaled vs. no. of tracks
  TH2                                *fHistRhoScaledVsNcluster;       //!<!rhoscaled vs. no. of clusters

  AliAnalysisTaskRhoBaseDev(const AliAnalysisTaskRhoBaseDev&);             // not implemented
  AliAnalysisTaskRhoBaseDev& operator=(const AliAnalysisTaskRhoBaseDev&);  // not implemented
  
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskRhoBaseDev, 2);
  /// \endcond
};
#endif
