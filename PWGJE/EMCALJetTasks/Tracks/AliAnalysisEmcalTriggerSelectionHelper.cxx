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
#include <array>
#include <algorithm>
#include <iostream>

#include "AliAnalysisEmcalTriggerSelectionHelper.h"
#include "AliEmcalTriggerStringDecoder.h"
#include "AliLog.h"

ClassImp(PWGJE::EMCALJetTasks::AliAnalysisEmcalTriggerSelectionHelperImpl)
ClassImp(PWGJE::EMCALJetTasks::AliAnalysisEmcalTriggerSelectionHelper)

using namespace PWGJE::EMCALJetTasks;

bool AliAnalysisEmcalTriggerSelectionHelperImpl::IsSelectEmcalTriggers(EMCAL_STRINGVIEW triggerstring) const {
  const std::array<std::string, 9> kEMCALTriggers = {
    "EJE", "EJ1", "EJ2", "DJ1", "DJ2", "EG1", "EG2", "DG1", "DG2"
  };
  bool isEMCAL = false;
  for(auto emcaltrg : kEMCALTriggers) {
    if(triggerstring.find(emcaltrg) != std::string::npos) {
      isEMCAL = true;
      break;
    }
  }
  return isEMCAL;
}

std::string AliAnalysisEmcalTriggerSelectionHelperImpl::MatchTrigger(EMCAL_STRINGVIEW triggerstring, EMCAL_STRINGVIEW triggerselectionstring, bool useMuonCalo) const{
  auto triggerclasses = PWG::EMCAL::Triggerinfo::DecodeTriggerString(triggerstring.data());
  std::string result;
  for(const auto &t : triggerclasses) {
    // Use CENT cluster for downscaling
    if(t.BunchCrossing() != "B") continue;
    if(useMuonCalo){
      if(t.Triggercluster() != "CALO") continue;
    } else {
      if(t.Triggercluster() != "CENT") continue;
    }
    if(t.Triggerclass().find(triggerselectionstring.data()) == std::string::npos) continue; 
    result = t.ExpandClassName();
    break;
  }
  return result;
}

std::vector<AliAnalysisEmcalTriggerSelectionHelper::TriggerCluster_t> AliAnalysisEmcalTriggerSelectionHelperImpl::GetTriggerClusterIndices(EMCAL_STRINGVIEW triggerstring) const {
  // decode trigger string in order to determine the trigger clusters
  AliDebugGeneralStream("AliAnalysisEmcalTriggerSelectionHelperImpl::GetTriggerClusterIndices", 4) << "Triggerstring: " << triggerstring.data() << std::endl;
  std::vector<TriggerCluster_t> result;
  result.emplace_back(kTrgClusterANY);      // cluster ANY always included 

  // Data - separate trigger clusters
  std::vector<std::string> clusternames;
  auto triggerinfos = PWG::EMCAL::Triggerinfo::DecodeTriggerString(triggerstring.data());
  for(auto t : triggerinfos) {
    if(std::find(clusternames.begin(), clusternames.end(), t.Triggercluster()) == clusternames.end()) clusternames.emplace_back(t.Triggercluster());
  }
  bool isCENT = (std::find(clusternames.begin(), clusternames.end(), "CENT") != clusternames.end()),
       isCENTNOTRD = (std::find(clusternames.begin(), clusternames.end(), "CENTNOTRD") != clusternames.end()),
       isCALO = (std::find(clusternames.begin(), clusternames.end(), "CALO") != clusternames.end()),
       isCALOFAST = (std::find(clusternames.begin(), clusternames.end(), "CALOFAST") != clusternames.end()),
       isCENTNOPMD = (std::find(clusternames.begin(), clusternames.end(), "CENTNOPMD") != clusternames.end()),
       isALL = (std::find(clusternames.begin(), clusternames.end(), "ALL") != clusternames.end()),
       isALLNOTRD = (std::find(clusternames.begin(), clusternames.end(), "ALLNOTRD") != clusternames.end());
  AliDebugGeneralStream("AliAnalysisEmcalTriggerSelectionHelperImpl::GetTriggerClusterIndices", 4) << "Selected trigger clusters: CENT: " << (isCENT ? "yes" : "no") << ", CENTNOTRD: " << (isCENTNOTRD ? "yes" : "no") << ", CALO: " << (isCALO ? "yes" : "no") << ", CALOFAST: " << (isCALOFAST ? "yes" :  "no") << std::endl;
  if(isCENT || isCENTNOTRD) {
    if(isCENT) {
      result.emplace_back(kTrgClusterCENT);
      if(isCENTNOTRD) {
        result.emplace_back(kTrgClusterCENTNOTRD);
        result.emplace_back(kTrgClusterCENTBOTH);
      } else result.emplace_back(kTrgClusterOnlyCENT);
    } else {
      result.emplace_back(kTrgClusterCENTNOTRD);
      result.emplace_back(kTrgClusterOnlyCENTNOTRD);
    }
  }
  if(isCALO || isCALOFAST) {
    if(isCALO) {
      result.emplace_back(kTrgClusterCALO);
      if(isCALOFAST) {
        result.emplace_back(kTrgClusterCALOFAST);
        result.emplace_back(kTrgClusterCALOBOTH);
      } else result.emplace_back(kTrgClusterOnlyCALO);
    } else {
      result.emplace_back(kTrgClusterCALOFAST);
      result.emplace_back(kTrgClusterOnlyCALOFAST);
    }
  }
  if(isALL || isALLNOTRD) {
    if(isALL) {
      result.emplace_back(kTrgClusterALL);
      if(isALLNOTRD) {
        result.emplace_back(kTrgClusterALLNOTRD);
        result.emplace_back(kTrgClusterALLBOTH);
      } else result.emplace_back(kTrgClusterOnlyALL);
    } else {
      result.emplace_back(kTrgClusterALLNOTRD);
      result.emplace_back(kTrgClusterOnlyALLNOTRD);
    }
  }
  if(isCENTNOPMD) result.emplace_back(kTrgClusterCENTNOPMD);
  return result;
}

std::string AliAnalysisEmcalTriggerSelectionHelperImpl::GetNameTriggerCluster(TriggerCluster_t clust) const{
  const std::array<std::string, kTrgClusterN> kNamesTriggerCluster = {{"ANY", "CENT", "CENTNOTRD", "CALO", "CALOFAST", 
                                                                       "CENTBOTH", "OnlyCENT", "OnlyCENTNOTRD", "CALOBOTH", 
                                                                       "OnluCALO", "OnlyCALOFAST", "CENTNOPMD", "ALL", "ALLNOTRD", 
                                                                       "ALLBOTH", "OnlyALL", "OnlyALLNOTRD"}};
  return kNamesTriggerCluster[clust];
}