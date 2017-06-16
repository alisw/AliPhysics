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

#include "AliAnalysisTaskEmcalJetLight.h"


/** \class AliAnalysisTaskRhoBaseDev
 * \brief Base class for a task that calculates the UE
 *
 * Base class for a task that calculates the average background
 * coming from the underlying event (UE) in jet analysis.
 * This is a development version. The stable version of this class
 * is AliAnalysisTaskRhoBase.
 */
class AliAnalysisTaskRhoBaseDev : public AliAnalysisTaskEmcalJetLight {
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
  void                   SetHistoBins(Int_t nbins, Double_t min, Double_t max) { fNbins = nbins; fMinBinPt = min; fMaxBinPt = max    ; }

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
  void                   ExecOnce();
  Bool_t                 Run();
  Bool_t                 FillHistograms();

  virtual Double_t       GetRhoFactor(Double_t cent);
  virtual Double_t       GetScaleFactor(Double_t cent);

  TString                fOutRhoName;                    ///< name of output rho object
  TString                fOutRhoScaledName;              ///< name of output scaled rho object
  TF1                   *fRhoFunction;                   ///< pre-computed rho as a function of centrality
  TF1                   *fScaleFunction;                 ///< pre-computed scale factor as a function of centrality
  Bool_t                 fAttachToEvent;                 ///< whether or not attach rho to the event objects list
  Int_t                  fNbins;                         ///< no. of pt bins
  Double_t               fMinBinPt;                      ///< min pt in histograms
  Double_t               fMaxBinPt;                      ///< max pt in histograms

  AliRhoParameter       *fOutRho;                        //!<!output rho object
  AliRhoParameter       *fOutRhoScaled;                  //!<!output scaled rho object

  TH2                   *fHistJetPtvsCent;               //!<!jet pt vs. centrality
  TH2                   *fHistJetAreavsCent;             //!<!jet area vs. centrality
  TH2                   *fHistJetRhovsCent;              //!<!jet pt/area vs. centrality
  TH2                   *fHistNjetvsCent;                //!<!no. of jets vs. centrality
  TH2                   *fHistJetPtvsNtrack;             //!<!jet pt vs. no. of tracks
  TH2                   *fHistJetAreavsNtrack;           //!<!jet area vs. no. of tracks
  TH2                   *fHistNjetvsNtrack;              //!<!no. of jets vs. no. of tracks
  TH2                   *fHistJetNconstVsPt[4];          //!<!jet no. of constituents vs. pt
  TH2                   *fHistJetRhovsEta[4];            //!<!rho vs. eta
  TH2                   *fHistRhovsCent;                 //!<!rho vs. centrality
  TH2                   *fHistRhoScaledvsCent;           //!<!rhoscaled vs. centrality
  TH2                   *fHistRhovsNtrack;               //!<!rho vs. no. of tracks
  TH2                   *fHistRhoScaledvsNtrack;         //!<!rhoscaled vs. no. of tracks
  TH2                   *fHistRhovsNcluster;             //!<!rho vs. no. of clusters
  TH2                   *fHistRhoScaledvsNcluster;       //!<!rhoscaled vs. no. of clusters

  AliAnalysisTaskRhoBaseDev(const AliAnalysisTaskRhoBaseDev&);             // not implemented
  AliAnalysisTaskRhoBaseDev& operator=(const AliAnalysisTaskRhoBaseDev&);  // not implemented
  
  /// \cond CLASSIMP
  ClassDef(AliAnalysisTaskRhoBaseDev, 1);
  /// \endcond
};
#endif
