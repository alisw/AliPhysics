/************************************************************************************
 * Copyright (C) 2012, Copyright Holders of the ALICE Collaboration                 *
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
#ifndef ALIANALYSISTASKRHOBASE_H
#define ALIANALYSISTASKRHOBASE_H

class TString;
class TF1;
class TH1F;
class TH2F;
class TH3F;
class AliRhoParameter;

#include "AliAnalysisTaskEmcalJet.h"

/**
 * @class AliAnalysisTaskRhoBase
 * @brief Base class for rho calculation.
 * @ingroup PWGJEBASE
 * @author Salvatore Aiola
 * @since May 17, 2012
 * 
 * Calculates parameterized rho for given centrality independent of input.
 */
class AliAnalysisTaskRhoBase : public AliAnalysisTaskEmcalJet {
 public:
  /**
   * @brief Dummy constructor, for ROOT I/O
   */
  AliAnalysisTaskRhoBase();
  /**
   * @brief Construct a new Ali Analysis Task Rho Base object
   * 
   * @param name Name of the rho task
   * @param histo If true QA/debug histograms are created
   */
  AliAnalysisTaskRhoBase(const char *name, Bool_t histo=kFALSE);

  /**
   * @brief Destructor
   */
  virtual ~AliAnalysisTaskRhoBase() {}

  /**
   * @brief User create output objects, called at the beginning of the analysis.
   */
  void                   UserCreateOutputObjects();

  void                   SetOutRhoName(const char *name)                       { fOutRhoName           = name    ;
                                                                                 fOutRhoScaledName     = Form("%s_Scaled",name);     }
  void                   SetCompareRhoName(const char *name)                   { fCompareRhoName       = name    ;                   }
  void                   SetCompareRhoScaledName(const char *name)             { fCompareRhoScaledName = name    ;                   }
  void                   SetScaleFunction(TF1* sf)                             { fScaleFunction        = sf      ;                   }
  void                   SetRhoFunction(TF1* rf)                               { fRhoFunction          = rf      ;                   }
  /**
   * @brief Load the scale function from a file.
   * @param path Path of the file from which to load the scale function
   * @param name Name of the scale function
   * @return Scale function 
   */
  TF1*                   LoadRhoFunction(const char* path, const char* name);
  void                   SetInEventSigmaRho(Double_t s)                        { fInEventSigmaRho      = s       ;                   }
  void                   SetAttachToEvent(Bool_t a)                            { fAttachToEvent        = a       ;                   }
  void                   SetSmallSystem(Bool_t setter = kTRUE)                 { fIsPbPb               = !setter ;                   }

  const char*            GetOutRhoName() const                                 { return fOutRhoName.Data()       ;                   }
  const char*            GetOutRhoScaledName() const                           { return fOutRhoScaledName.Data() ;                   }

 protected:
  /**
   * @brief Init the analysis.
   */
  void                   ExecOnce();

  /**
   * @brief Run the analysis.
   * @return Always true
   */
  Bool_t                 Run();
  /**
   * @brief Fill histograms.
   * @return Always true
   */
  Bool_t                 FillHistograms();

  /**
   * @brief Access rho valie per centrality.
   * @param cent Centrality class
   * @return Rho value for the given centrality class
   */
  virtual Double_t       GetRhoFactor(Double_t cent);

  /**
   * @brief Get scale factor per centrality.
   * @param cent Centrality class 
   * @return Scale factor for a given centrality 
   */
  virtual Double_t       GetScaleFactor(Double_t cent);

  TString                fOutRhoName;                    ///< name of output rho object
  TString                fOutRhoScaledName;              ///< name of output scaled rho object
  TString                fCompareRhoName;                ///< name of rho object to compare
  TString                fCompareRhoScaledName;          ///< name of scaled rho object to compare
  TF1                   *fRhoFunction;                   ///< pre-computed rho as a function of centrality
  TF1                   *fScaleFunction;                 ///< pre-computed scale factor as a function of centrality
  Double_t               fInEventSigmaRho;               ///< in-event sigma rho
  Bool_t                 fAttachToEvent;                 ///< whether or not attach rho to the event objects list
  Bool_t                 fIsPbPb;                        ///< different histogram ranges for pp/pPb and PbPb
  
  AliRhoParameter       *fOutRho;                        //!<! output rho object
  AliRhoParameter       *fOutRhoScaled;                  //!<! output scaled rho object
  AliRhoParameter       *fCompareRho;                    //!<! rho object to compare
  AliRhoParameter       *fCompareRhoScaled;              //!<! scaled rho object to compare

  TH2F                  *fHistJetPtvsCent;               //!<! jet pt vs. centrality
  TH2F                  *fHistJetAreavsCent;             //!<! jet area vs. centrality
  TH2F                  *fHistJetRhovsCent;              //!<! jet pt/area vs. centrality
  TH2F                  *fHistNjetvsCent;                //!<! no. of jets vs. centrality
  TH2F                  *fHistJetPtvsNtrack;             //!<! jet pt vs. no. of tracks
  TH2F                  *fHistJetAreavsNtrack;           //!<! jet area vs. no. of tracks
  TH2F                  *fHistNjetvsNtrack;              //!<! no. of jets vs. no. of tracks
  TH2F                  *fHistNjUEoverNjVsNj[12];        //!<! ratio no. of jets below rho*A+sigma_rho over. no. of jets vs. no. of jets
  TH2F                  *fHistJetNconstVsPt[4];          //!<! jet no. of constituents vs. pt
  TH2F                  *fHistJetRhovsEta[4];            //!<! rho vs. eta
  TH2F                  *fHistRhovsCent;                 //!<! rho vs. centrality
  TH2F                  *fHistRhoScaledvsCent;           //!<! rhoscaled vs. centrality
  TH2F                  *fHistDeltaRhovsCent;            //!<! delta rho vs. centrality
  TH2F                  *fHistDeltaRhoScalevsCent;       //!<! delta rhoscaled vs. centrality

  TH3F                  *fHistRhovsNtrackvsV0Mult;       //!<! rho vs. no. of tracks vs V0mult
  TH3F                  *fHistRhoScaledvsNtrackvsV0Mult; //!rhoscaled vs. no. of tracks vs V0mult
  TH2F                  *fHistDeltaRhovsNtrack;          //!delta rho vs. no. of tracks
  TH2F                  *fHistDeltaRhoScalevsNtrack;     //!delta rho scaled vs. no. of tracks
 
  TH2F                  *fHistRhovsNcluster;             //!rho vs. no. of clusters
  TH2F                  *fHistRhoScaledvsNcluster;       //!rhoscaled vs. no. of clusters

  AliAnalysisTaskRhoBase(const AliAnalysisTaskRhoBase&);             // not implemented
  AliAnalysisTaskRhoBase& operator=(const AliAnalysisTaskRhoBase&);  // not implemented
  
  ClassDef(AliAnalysisTaskRhoBase, 11); // Rho base task
};
#endif
