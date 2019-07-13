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
#ifndef ALIANALYSISTASKRHOMASSSPARSE_H
#define ALIANALYSISTASKRHOMASSSPARSE_H

#include "AliAnalysisTaskRhoMassBase.h"

/**
 * @class AliAnalysisTaskRhoMassSparse
 * @brief Calculation of rho mass from a collection of jets.
 * @ingroup PWGJEBASE
 * @author Martha Verweij
 * @since Dec 19, 2014
 * 
 * If scale function is given the scaled rho will be exported
 * with the name as "fOutRhoMassName".Apppend("_Scaled").
 */
class AliAnalysisTaskRhoMassSparse : public AliAnalysisTaskRhoMassBase {

 public:

  /**
   * @brief Constructor.
   */
  AliAnalysisTaskRhoMassSparse();

  /**
   * @brief Constructor.
   * @param name Name of the rho task
   * @param histo If true QA/debug histograms are created
   */
  AliAnalysisTaskRhoMassSparse(const char *name, Bool_t histo=kFALSE);

  /**
   * @brief Destructor
   */
  virtual ~AliAnalysisTaskRhoMassSparse() {}

  enum JetRhoMassType {
    kMd     = 0,            ///< rho_m from arXiv:1211.2811
    kMdP    = 1,            ///< rho_m using P instead of pT
    kMd4    = 2             ///< rho_m using addition of 4-vectors
  };

  /**
   * @brief User create output objects, called at the beginning of the analysis.
   */
  void             UserCreateOutputObjects();

  void             SetExcludeLeadJets(UInt_t n)     { fNExclLeadJets  = n   ; }
  void             SetRhoMassType(JetRhoMassType t) { fJetRhoMassType = t   ; }
  void             SetPionMassForClusters(Bool_t b) { fPionMassClusters = b ; }
  void             SetRhoCMS(Bool_t cms)           { fRhoCMS = cms ; }
  Bool_t           IsJetOverlapping(AliEmcalJet* jet1, AliEmcalJet* jet2);
  Bool_t           IsJetSignal(AliEmcalJet* jet1);


 protected:
  /**
   * @brief Run the analysis.
   * @return Always true
   */
  Bool_t           Run();

  Double_t         GetSumMConstituents(AliEmcalJet *jet);
  Double_t         GetSumPtConstituents(AliEmcalJet *jet);
  
  /**
   * @brief Get md as defined in http://arxiv.org/pdf/1211.2811.pdf
   * @param jet Jet for which md is calculated
   * @return the md value 
   */
  Double_t         GetMd(AliEmcalJet *jet);

  UInt_t           fNExclLeadJets;                 ///< number of leading jets to be excluded from the median calculation
  Bool_t           fRhoCMS;                        ///< flag to run CMS method
  JetRhoMassType   fJetRhoMassType;                ///< method for rho_m calculation
  Bool_t           fPionMassClusters;              ///< assume pion mass for clusters

  TH2F            *fHistMdAreavsCent;              //!<! Md/Area vs cent for all kt clusters
  TH2F            *fHistOccCorrvsCent;             //!<! occupancy correction vs. centrality

  AliAnalysisTaskRhoMassSparse(const AliAnalysisTaskRhoMassSparse&);             // not implemented
  AliAnalysisTaskRhoMassSparse& operator=(const AliAnalysisTaskRhoMassSparse&);  // not implemented
  
  ClassDef(AliAnalysisTaskRhoMassSparse, 1); // Rho_m task
};
#endif
