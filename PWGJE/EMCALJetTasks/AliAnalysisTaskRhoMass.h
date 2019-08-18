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
#ifndef ALIANALYSISTASKRHOMASS_H
#define ALIANALYSISTASKRHOMASS_H

#include "AliAnalysisTaskRhoMassBase.h"

/**
 * @class AliAnalysisTaskRhoMass
 * @brief Calculation of rho mass from a collection of jets.
 * @ingroup PWGJEBASE
 * @author Martha Verweij
 * @since Dec. 19, 2014
 * 
 * If scale function is given the scaled rho will be exported
 * with the name as "fOutRhoMassName".Apppend("_Scaled").
 */
class AliAnalysisTaskRhoMass : public AliAnalysisTaskRhoMassBase {

 public:
  /**
   * @brief Constructor.
   */
  AliAnalysisTaskRhoMass();

  /**
   * @brief Constructor.
   * @param name Name of the rho task
   * @param histo If true QA/debug histograms are created
   */
  AliAnalysisTaskRhoMass(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskRhoMass() {}

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
   * @return Double_t md value
   */
  Double_t         GetMd(AliEmcalJet *jet);

  UInt_t           fNExclLeadJets;                 ///< number of leading jets to be excluded from the median calculation
  JetRhoMassType   fJetRhoMassType;                ///< method for rho_m calculation
  Bool_t           fPionMassClusters;              ///< assume pion mass for clusters

  TH2F            *fHistMdAreavsCent;              //!<! Md/Area vs cent for all kt clusters

  AliAnalysisTaskRhoMass(const AliAnalysisTaskRhoMass&);             // not implemented
  AliAnalysisTaskRhoMass& operator=(const AliAnalysisTaskRhoMass&);  // not implemented
  
  ClassDef(AliAnalysisTaskRhoMass, 2); // Rho_m task
};
#endif
