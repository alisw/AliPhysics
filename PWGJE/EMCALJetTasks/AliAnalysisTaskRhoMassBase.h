/************************************************************************************
 * Copyright (C) 2014, Copyright Holders of the ALICE Collaboration                 *
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
#ifndef ALIANALYSISTASKRHOMASSBASE_H
#define ALIANALYSISTASKRHOMASSBASE_H

class TString;
class TF1;
class TH1F;
class TH2F;
class AliRhoParameter;

#include "AliAnalysisTaskEmcalJet.h"

/**
 * @class AliAnalysisTaskRhoMassBase
 * @brief  Base class for rho mass calculation.
 * @ingroup PWGJEBASE
 * @author Martha Verweij
 * @since Dec 19, 2014
 *
 * Calculates parameterized rho mass for given centrality independent of input.
 * Similar to AliAnalysisTaskRhoBase
 */
class AliAnalysisTaskRhoMassBase : public AliAnalysisTaskEmcalJet {
 public:
  /**
   * @brief Constructor.
   */
  AliAnalysisTaskRhoMassBase();

  /**
   * @brief Constructor.
   * @param name of the rho task
   * @param histo If true QA/debug histograms are created
   */
  AliAnalysisTaskRhoMassBase(const char *name, Bool_t histo=kFALSE);

  /**
   * @brief Destructor
   * 
   */
  virtual ~AliAnalysisTaskRhoMassBase() {}

  /**
   * @brief User create output objects, called at the beginning of the analysis.
   */
  void                   UserCreateOutputObjects();

  void                   SetOutRhoMassName(const char *name)                   { fOutRhoMassName           = name ; 
                                                                                 fOutRhoMassScaledName     = Form("%s_Scaled",name) ; }
  void                   SetCompareRhoMassName(const char *name)               { fCompareRhoMassName       = name ;                   }
  void                   SetCompareRhoMassScaledName(const char *name)         { fCompareRhoMassScaledName = name ;                   }
  void                   SetScaleFunction(TF1* sf)                             { fScaleFunction        = sf   ;                       }
  void                   SetRhoMassFunction(TF1* rf)                           { fRhoMassFunction      = rf   ;                       }
  void                   SetAttachToEvent(Bool_t a)                            { fAttachToEvent        = a    ;                       }
  void                   SetSmallSystem(Bool_t setter = kTRUE)                 {fIsPbPb = !setter; }


  const TString&         GetOutRhoMassName() const                             { return fOutRhoMassName;                              }
  const TString&         GetOutRhoMassScaledName() const                       { return fOutRhoMassScaledName;                        } 

 protected:

  /**
   * @brief Init the analysis.
   */
  void                   ExecOnce();

  /**
   * @brief  Run the analysis.
   * @return Always true
   */
  Bool_t                 Run();

  /**
   * @brief Fill histograms.
   * @return Always true
   */
  Bool_t                 FillHistograms();

  /**
   * @brief Return rho per centrality.
   * @param cent Centrality class
   * @return Rho value for the centrality class
   */
  virtual Double_t       GetRhoMassFactor(Double_t cent);

  /**
   * @brief Get scale factor.
   * @param cent Centrality class
   * @return Double_t Scale factor for centrality class
   */
  virtual Double_t       GetScaleFactor(Double_t cent);

  TString                fOutRhoMassName;                ///< name of output rho mass object
  TString                fOutRhoMassScaledName;          ///< name of output scaled rho mass object
  TString                fCompareRhoMassName;            ///< name of rho mass object to compare
  TString                fCompareRhoMassScaledName;      ///< name of scaled rho mass object to compare
  TF1                   *fRhoMassFunction;               ///< pre-computed rho mass as a function of centrality
  TF1                   *fScaleFunction;                 ///< pre-computed scale factor as a function of centrality
  Bool_t                 fAttachToEvent;                 ///< whether or not attach rho mass to the event objects list
  Bool_t                 fIsPbPb;                        ///< different histogram ranges for pp/pPb and PbPb
  
  AliRhoParameter       *fOutRhoMass;                    //!<! output rho object
  AliRhoParameter       *fOutRhoMassScaled;              //!<! output scaled rho object
  AliRhoParameter       *fCompareRhoMass;                //!<! rho object to compare
  AliRhoParameter       *fCompareRhoMassScaled;          //!<! scaled rho object to compare

  TH2F                  *fHistJetMassvsCent;             //!<! jet mass vs. centrality
  TH2F                  *fHistRhoMassvsCent;             //!<! rho mass vs. centrality
  TH2F                  *fHistRhoMassScaledvsCent;       //!<! rho mass scaled vs. centrality
  TH2F                  *fHistDeltaRhoMassvsCent;        //!<! delta rho mass vs. centrality
  TH2F                  *fHistDeltaRhoMassScalevsCent;   //!<! delta rho mass scaled vs. centrality

  TH2F                  *fHistRhoMassvsNtrack;           //!<! rho mass vs. no. of tracks
  TH2F                  *fHistRhoMassScaledvsNtrack;     //!<! rho mass scaled vs. no. of tracks
  TH2F                  *fHistDeltaRhoMassvsNtrack;      //!<! delta rho mass vs. no. of tracks
  TH2F                  *fHistDeltaRhoMassScalevsNtrack; //!<! delta rho mass scaled vs. no. of tracks
 
  TH2F                  *fHistRhoMassvsNcluster;         //!<! rho mass vs. no. of clusters
  TH2F                  *fHistRhoMassScaledvsNcluster;   //!<! rho mass scaled vs. no. of clusters

  TH2F                  *fHistGammaVsNtrack;             //!<! Gamma(<E>/<M>) vs Ntrack

  AliAnalysisTaskRhoMassBase(const AliAnalysisTaskRhoMassBase&);             // not implemented
  AliAnalysisTaskRhoMassBase& operator=(const AliAnalysisTaskRhoMassBase&);  // not implemented
  
  ClassDef(AliAnalysisTaskRhoMassBase, 2); // Rho mass base task
};
#endif
