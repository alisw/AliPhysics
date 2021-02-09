/************************************************************************************
 * Copyright (C) 2017, Copyright Holders of the ALICE Collaboration                 *
 * All rights reserved.                                                             *
 *                                                                                  *
 * Redistribution and use in source and binary forms, with or without               *
 * modification, are permitted provided that the following conditions are met:      *
 *     * Redistributions of source code must retain the above copyright             *
 *       notice, this list of conditions and the following disclaimer.              *
 *     * Redistributions in binary form must reproduce the above copyright          *
 *       notice, this list of conditions and the following disclaimer in the        *
 *       documentation and/or other materials provided with the distribution.       *
 *     * Neither the name of the <organization> nor the                             *
 *       names of its contributors may be used to endorse or promote products       *
 *       derived from this software without specific prior written permission.      *
 *                                                                                  *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND  *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED    *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE           *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY              *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES       *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;     *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND      *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT       *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS    *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                     *
 ************************************************************************************/
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


/** 
 * @class AliAnalysisTaskRhoBaseDev
 * @brief Base class for a task that calculates the UE
 * @author Salvatore Aiola <salvatore.aiola@cern.ch>, Yale University
 * @date June 16, 2017
 *
 * Base class for a task that calculates the average background
 * coming from the underlying event (UE) in jet analysis.
 * This is a development version. The stable version of this class
 * is AliAnalysisTaskRhoBase.
 */
class AliAnalysisTaskRhoBaseDev : public AliAnalysisTaskJetUE {
 public:
  /**
   * @brief Default constructor. Needed by ROOT I/O
   */
  AliAnalysisTaskRhoBaseDev();

  /**
   * @brief Standard constructor. Should be used by the user.
   *
   * @param[in] name  Name of the task
   * @param[in] histo If kTRUE, the task will also produce QA histograms
   */
  AliAnalysisTaskRhoBaseDev(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskRhoBaseDev() {}

   /**
    * @brief Creating user output
    * 
    * Performing run-independent initialization.
    * Here the histograms should be instantiated.
    */
  void                   UserCreateOutputObjects();

  void                   SetOutRhoName(const char *name)                       { fOutRhoName           = name    ;
                                                                                 fOutRhoScaledName     = Form("%s_Scaled",name);     }
  void                   SetScaleFunction(TF1* sf)                             { fScaleFunction        = sf      ;                   }
  void                   SetRhoFunction(TF1* rf)                               { fRhoFunction          = rf      ;                   }

  /**
   * @brief Load the scale function from a file.
   * @param path Path to the file
   * @param name Name of the object inside the file
   * @return On success, the pointer to the TF1 object
   */
  TF1*                   LoadRhoFunction(const char* path, const char* name);
  void                   SetAttachToEvent(Bool_t a)                            { fAttachToEvent        = a       ;                   }

  const char*            GetOutRhoName() const                                 { return fOutRhoName.Data()       ;                   }
  const char*            GetOutRhoScaledName() const                           { return fOutRhoScaledName.Data() ;                   }

   /**
    * @brief Create an instance of this class and add it to the analysis manager
    * @param trackName name of the track collection
    * @param clusName name of the calorimeter cluster collection
    * @param nRho name of the output rho object
    * @param jetradius Radius of the kt jets used to calculate the background
    * @param acceptance Fiducial acceptance of the kt jets
    * @param jetType Jet type (full/charged)
    * @param rscheme Recombination scheme
    * @param histo If kTRUE the task will also produce QA histograms
    * @param suffix additional suffix that can be added at the end of the task name
    * @return pointer to the new AliAnalysisTaskRhoDev task
    */
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
  /**
   * Init the analysis.
   */
  void                                ExecOnce();

  /**
   * @brief Run the analysis.
   */
  Bool_t                              Run();

  /**
   * @brief Fill histograms.
   */
  Bool_t                              FillHistograms();

  virtual Bool_t                      VerifyContainers() { return kTRUE; }

  /**
    * @brief Calculates the average background using a given parametrization
    * as a function of centrality.
    */
  virtual void                        CalculateRho();

  /**
   * @brief Access rho per centrality.
   * @param cent Centrality percentile.
   * @return The average background from the parametrization loaded in memory.
   */
  virtual Double_t                    GetRhoFactor(Double_t cent);

  /**
   * @brief Get scale factor.
   * @param cent Centrality percentile.
   * @return The scale factor (charged -> full) from the parametrization loaded in memory.
   */
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
  
  ClassDef(AliAnalysisTaskRhoBaseDev, 2);
};
#endif
