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
#ifndef ALIANALYSISTASKRHO_H
#define ALIANALYSISTASKRHO_H

#include "AliAnalysisTaskRhoBase.h"

/**
 * @class
 * @brief Calculation of rho from a collection of jets.
 * @ingroup
 * @author R.Reed
 * @author S.Aiola
 * @since May 10, 2012
 * 
 * If scale function is given the scaled rho will be exported
 * with the name as "fOutRhoName".Apppend("_Scaled").
 */
class AliAnalysisTaskRho : public AliAnalysisTaskRhoBase {

 public:
  AliAnalysisTaskRho();
  AliAnalysisTaskRho(const char *name, Bool_t histo=kFALSE);
  virtual ~AliAnalysisTaskRho() {}

  void             SetExcludeLeadJets(UInt_t n)    { fNExclLeadJets = n    ; }

  static AliAnalysisTaskRho* AddTaskRhoNew (
    const char    *nTracks                        = "usedefault",
    const char    *nClusters                      = "usedefault",
    const char    *nRho                           = "Rho",
    Double_t       jetradius                      = 0.2,
    UInt_t         acceptance                     = AliEmcalJet::kTPCfid,
    AliJetContainer::EJetType_t jetType           = AliJetContainer::kChargedJet,
    const Bool_t   histo                          = kFALSE,
    AliJetContainer::ERecoScheme_t rscheme        = AliJetContainer::pt_scheme,
    const char    *suffix                         = ""
);

 protected:
  Bool_t           Run();

  UInt_t           fNExclLeadJets;                 ///< number of leading jets to be excluded from the median calculation

  AliAnalysisTaskRho(const AliAnalysisTaskRho&);             // not implemented
  AliAnalysisTaskRho& operator=(const AliAnalysisTaskRho&);  // not implemented
  
  ClassDef(AliAnalysisTaskRho, 10); // Rho task
};
#endif
