/**********************************************************************************
 * Copyright (C) 2016, Copyright Holders of the ALICE Collaboration                *
 * All rights reserved.                                                            *
 *                                                                                 *
 * Redistribution and use in source and binary forms, with or without              *
 * modification, are permitted provided that the following conditions are met:     *
 *   * Redistributions of source code must retain the above copyright              *
 *     notice, this list of conditions and the following disclaimer.               *
 *   * Redistributions in binary form must reproduce the above copyright           *
 *     notice, this list of conditions and the following disclaimer in the         *
 *     documentation and/or other materials provided with the distribution.        *
 *   * Neither the name of the <organization> nor the                              *
 *     names of its contributors may be used to endorse or promote products        *
 *     derived from this software without specific prior written permission.       *
 *                                                                                 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED   *
 * WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE          *
 * DISCLAIMED. IN NO EVENT SHALL ALICE COLLABORATION BE LIABLE FOR ANY             *
 * DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES      *
 * (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;    *
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND     *
 * ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT      *
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS   *
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.                    *
 * *********************************************************************************/
#ifndef ALIANALYSISTASKRHOSPARSE_H
#define ALIANALYSISTASKRHOSPARSE_H

#include "AliAnalysisTaskRhoBase.h"

/**
 * @class AliAnalysisTaskRhoSparse
 * @brief Calculation of background rho for sparse events (pp and pPb).
 * @ingroup PWGJEBASE
 * @authors R.Reed, S.Aiola, M.Connors, E. Epple, Yale University
 * @date Oct 11, 2018
 */
class AliAnalysisTaskRhoSparse : public AliAnalysisTaskRhoBase {

 public:
  /**
   * @brief Dummy constructor, for ROOT I/O
   */
  AliAnalysisTaskRhoSparse();

  /**
   * @brief Main constructor
   * @param name Name of the task
   * @param histo If true then QA/debug histograms are filled
   */
  AliAnalysisTaskRhoSparse(const char *name, Bool_t histo=kFALSE);

  /**
   * @brief Destructor
   */
  virtual ~AliAnalysisTaskRhoSparse() {}

  /**
   * Rho Sparse AddTask.
   */
  static AliAnalysisTaskRhoSparse* AddTaskRhoSparse(
  		const char    *nTracks     = "usedefault",
  		const char    *nClusters   = "usedefault",
  		const char    *nRho        = "Rho",
  		Double_t       jetradius   = 0.2,
  		UInt_t         acceptance  = AliEmcalJet::kTPC,
  		AliJetContainer::EJetType_t jetType    = AliJetContainer::kChargedJet,
  		AliJetContainer::ERecoScheme_t rscheme = AliJetContainer::pt_scheme,
  		const Bool_t   histo       = kFALSE,
  		const char    *nJetsSig    = "",
  		const char    *cutType     = "TPC",
  		Double_t       jetptcut    = 0.0,
  		Double_t       jetareacut  = 0.01,
  		Double_t       emcareacut  = 0,
  		const char    *suffix      = ""
  		);

  /**
   * @brief Create output histograms
   */
  void             UserCreateOutputObjects();
  void             SetExcludeLeadJets(UInt_t n)           { fNExclLeadJets = n       ; }
  void             SetExcludeOverlapJets(Bool_t input)    { fExcludeOverlaps = input ; }
  void             SetRhoCMS(Bool_t cms)                  { fRhoCMS = cms            ; }
  void             SetAreaCalculationDetails(Bool_t inputTPCArea, Bool_t inputExcludeJetArea){ fUseTPCArea = inputTPCArea ; fExcludeAreaExcludedJets = inputExcludeJetArea; }

  /**
   * @brief Check whether two jets are overlapping
   * 
   * Overlapping jets must have at least one track in common
   * 
   * @param jet1 First jet
   * @param jet2 Second jet
   * @return Bool_t True if the jets are overlapping, false otherwise
   */
  Bool_t           IsJetOverlapping(AliEmcalJet* jet1, AliEmcalJet* jet2);

  /**
   * @brief Select jet as signal jet
   * 
   * Signal jets must have a pt > 5 GeV/c
   * 
   * @param jet1 Jet to be tested
   * @return Bool_t True if jet is classified as signal jet, false otherwise
   */
  Bool_t           IsJetSignal(AliEmcalJet* jet1);

 protected:
  /**
   * @brief Run the analysis.
   * @return always true
   */
  Bool_t           Run();

  UInt_t           fNExclLeadJets;                                    ///< number of leading jets to be excluded from the median calculation
  Bool_t           fExcludeOverlaps;                                  ///< exclude background jets that overlap (share at least one track) with anti-KT signal jets
  Bool_t           fRhoCMS;                                           ///< flag to run CMS method
  Bool_t           fUseTPCArea;                                       ///< use the full TPC area for the denominator of the occupancy calculation
  Bool_t           fExcludeAreaExcludedJets;                          ///<
  TH2F            *fHistOccCorrvsCent;            				            //!<! occupancy correction vs. centrality

  AliAnalysisTaskRhoSparse(const AliAnalysisTaskRhoSparse&);           ///< not implemented
  AliAnalysisTaskRhoSparse& operator=(const AliAnalysisTaskRhoSparse&);///< not implemented
  
  ClassDef(AliAnalysisTaskRhoSparse, 2);                               ///< Rho task
};
#endif
