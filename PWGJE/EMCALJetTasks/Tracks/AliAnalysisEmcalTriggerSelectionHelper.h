/************************************************************************************
 * Copyright (C) 2019, Copyright Holders of the ALICE Collaboration                 *
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
#ifndef ALIANALYSISEMCALTRIGGERSELECTIONHELPER_H
#define ALIANALYSISEMCALTRIGGERSELECTIONHELPER_H

#include <string>
#include <vector>
#include <TObject.h>
#include "AliEmcalStringView.h"

namespace PWGJE {

namespace EMCALJetTasks {

class AliAnalysisEmcalTriggerSelectionHelperImpl {
public:
  enum TriggerCluster_t {
    kTrgClusterANY,
    kTrgClusterCENT,
    kTrgClusterCENTNOTRD,
    kTrgClusterCALO,
    kTrgClusterCALOFAST,
    kTrgClusterCENTBOTH,
    kTrgClusterOnlyCENT,
    kTrgClusterOnlyCENTNOTRD,
    kTrgClusterCALOBOTH,
    kTrgClusterOnlyCALO,
    kTrgClusterOnlyCALOFAST,
    kTrgClusterN
  };

  AliAnalysisEmcalTriggerSelectionHelperImpl() {}
  virtual ~AliAnalysisEmcalTriggerSelectionHelperImpl() {}

#if (defined(__CINT_) && !defined(__CLING__)) || (defined(__MAKECINT__) && !defined(__ROOTCLING__))
  // ROOT5 function headers
  std::vector<PWGJE::EMCALJetTasks::AliAnalysisEmcalTriggerSelectionHelperImpl::TriggerCluster_t> GetTriggerClusterIndices(EMCAL_STRINGVIEW triggerstring) const;
  std::vector<PWGJE::EMCALJetTasks::AliAnalysisEmcalTriggerSelectionHelperImpl::TriggerCluster_t> GetTriggerClustersANY() const { return {kTrgClusterANY}; }
#else
  // ROOT6 function headers
  std::vector<TriggerCluster_t> GetTriggerClusterIndices(EMCAL_STRINGVIEW triggerstring) const;
  std::vector<TriggerCluster_t> GetTriggerClustersANY() const { return {kTrgClusterANY}; }
#endif
  bool IsSelectEmcalTriggers(EMCAL_STRINGVIEW triggerstring) const;
  std::string MatchTrigger(EMCAL_STRINGVIEW striggerstring, EMCAL_STRINGVIEW triggerselectionstring, bool useMuonCalo = false) const;
  std::string GetNameTriggerCluster(TriggerCluster_t clust) const;

  ClassDef(AliAnalysisEmcalTriggerSelectionHelperImpl, 1);
};

class AliAnalysisEmcalTriggerSelectionHelper : public TObject, public AliAnalysisEmcalTriggerSelectionHelperImpl {
public:
  AliAnalysisEmcalTriggerSelectionHelper() : TObject(), AliAnalysisEmcalTriggerSelectionHelperImpl() {}
  virtual ~AliAnalysisEmcalTriggerSelectionHelper() {}

  ClassDef(AliAnalysisEmcalTriggerSelectionHelper, 1);
}; 

}

}
#endif
