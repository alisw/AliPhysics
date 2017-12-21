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
#include <functional>
#include <iostream>
#include <string>

#include "TFormula.h"

#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEmcalESDHybridTrackCuts.h"
#include "AliLog.h"

using namespace PWG::EMCAL;

AliEmcalESDHybridTrackCuts::AliEmcalESDHybridTrackCuts():
  AliVCuts(),
  fLocalInitialized(kFALSE),
  fHybridTrackDefinition(kDef2010),
  fHybridTrackCutsGlobal(nullptr),
  fHybridTrackCutsConstrained(nullptr),
  fHybridTrackCutsNoItsRefit(nullptr)
{

}

AliEmcalESDHybridTrackCuts::AliEmcalESDHybridTrackCuts(const char *name, HybridDefinition_t hybriddef):
  AliVCuts(name, ""),
  fLocalInitialized(kFALSE),
  fHybridTrackDefinition(hybriddef),
  fHybridTrackCutsGlobal(nullptr),
  fHybridTrackCutsConstrained(nullptr),
  fHybridTrackCutsNoItsRefit(nullptr)
{

}

AliEmcalESDHybridTrackCuts::~AliEmcalESDHybridTrackCuts(){
  if(fHybridTrackCutsGlobal) delete fHybridTrackCutsGlobal;
  if(fHybridTrackCutsConstrained) delete fHybridTrackCutsConstrained;
  if(fHybridTrackCutsNoItsRefit) delete fHybridTrackCutsNoItsRefit;
}

bool AliEmcalESDHybridTrackCuts::IsSelected(TObject *o){
  AliDebugStream(1) << "AliEmcalESDHybridTrackCuts::IsSelected(): Called" << std::endl;
  if(!fLocalInitialized) Init();
  if(auto esdtrack = dynamic_cast<AliESDtrack *>(o)) {
    bool selected[3] = {false, false, false};
    if(fHybridTrackCutsGlobal && fHybridTrackCutsGlobal->AcceptTrack(esdtrack)) selected[0] = true;
    if(fHybridTrackCutsConstrained && fHybridTrackCutsConstrained->AcceptTrack(esdtrack)) selected[1] = true;
    if(fHybridTrackCutsNoItsRefit && fHybridTrackCutsNoItsRefit->AcceptTrack(esdtrack)) selected[2]= true;
    return selected[0] || selected[1] || selected[2];
  }
  AliErrorStream() << "No ESD track" << std::endl;
  return false;
}

void AliEmcalESDHybridTrackCuts::Init(){
  switch(fHybridTrackDefinition){
  case kDef2010: InitHybridTracks2010(); break;
  case kDef2011: InitHybridTracks2011(); break;
  default:
    AliErrorStream() << "No matching initialization found for requested hybrid track definition" <<std::endl;
  };
  fLocalInitialized = true;
}

void AliEmcalESDHybridTrackCuts::InitHybridTracks2010() {
  std::cout << "Initializing hybrid track cuts based on the 2010 definition" << std::endl;
  auto baseCutsFactory = [](std::string cutsname) -> AliESDtrackCuts * {
    auto trackcuts = new AliESDtrackCuts(cutsname.data());    
    auto f1NClustersTPCLinearPtDep = new TFormula("f1NClustersTPCLinearPtDep","70.+30./20.*x");
    trackcuts->SetMinNClustersTPCPtDep(f1NClustersTPCLinearPtDep,20.);
    trackcuts->SetMinNClustersTPC(70);
    trackcuts->SetMaxChi2PerClusterTPC(4);
    trackcuts->SetRequireTPCStandAlone(kTRUE); //cut on NClustersTPC and chi2TPC Iter1
    trackcuts->SetAcceptKinkDaughters(kFALSE);
    trackcuts->SetRequireTPCRefit(kTRUE);
    trackcuts->SetMaxFractionSharedTPCClusters(0.4);
    // ITS
    trackcuts->SetRequireITSRefit(kTRUE);
    //accept secondaries
    trackcuts->SetMaxDCAToVertexXY(2.4);
    trackcuts->SetMaxDCAToVertexZ(3.2);
    trackcuts->SetDCAToVertex2D(kTRUE);
    //reject fakes
    trackcuts->SetMaxChi2PerClusterITS(36);
    trackcuts->SetMaxChi2TPCConstrainedGlobal(36);

    trackcuts->SetRequireSigmaToVertex(kFALSE);

    trackcuts->SetEtaRange(-0.9,0.9);
    trackcuts->SetPtRange(0.15, 1E+15);

    return trackcuts;
  };

  fHybridTrackCutsGlobal = baseCutsFactory("JetCuts10001006");
  fHybridTrackCutsGlobal->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kAny);

  // the complement to the one with SPD requirement: tracks with ITS refit but no SPD hit
  //  AliESDtrackCuts* esdTrackCutsHG1 = CreateTrackCutsPWGJE(10011006);
  fHybridTrackCutsConstrained = baseCutsFactory("JetCuts10011006");
  fHybridTrackCutsConstrained->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);

  if(this->fSelectNonRefitTracks) {
    // all complementary hybrid tracks: no SPD requirement, no ITS refit requirement
    //  AliESDtrackCuts* esdTrackCutsGCOnly = CreateTrackCutsPWGJE(10041006);
    AliInfoStream() << "Create selection for non-refit tracks" << std::endl;
    fHybridTrackCutsNoItsRefit = baseCutsFactory("JetCuts10041006");
    fHybridTrackCutsConstrained->SetClusterRequirementITS(AliESDtrackCuts::kSPD, AliESDtrackCuts::kNone);
    fHybridTrackCutsNoItsRefit->SetRequireITSRefit(kFALSE);
  }
}

void AliEmcalESDHybridTrackCuts::InitHybridTracks2011() {
  std::cout << "Initialiing hybrid tracks based on the 2011 definition" <<std::endl;
  // first the global tracks we want to take
  fHybridTrackCutsGlobal = AliESDtrackCuts::GetStandardITSTPCTrackCuts2011(kFALSE); 
  fHybridTrackCutsGlobal->SetName("Global Hybrid tracks, loose DCA");
  fHybridTrackCutsGlobal->SetMaxDCAToVertexXY(2.4);
  fHybridTrackCutsGlobal->SetMaxDCAToVertexZ(3.2);
  fHybridTrackCutsGlobal->SetDCAToVertex2D(kTRUE);
  fHybridTrackCutsGlobal->SetMaxChi2TPCConstrainedGlobal(36);
  fHybridTrackCutsGlobal->SetMaxFractionSharedTPCClusters(0.4);

  // then the global constrainted tracks
  fHybridTrackCutsConstrained = new AliESDtrackCuts(*fHybridTrackCutsGlobal);
  fHybridTrackCutsConstrained->SetName("Global Constraint Hybrid tracks, loose DCA no it requirement");
  fHybridTrackCutsConstrained->SetClusterRequirementITS(AliESDtrackCuts::kSPD,AliESDtrackCuts::kOff);
  fHybridTrackCutsConstrained->SetRequireITSRefit(kTRUE);

  if(fSelectNonRefitTracks){
    AliInfoStream() << "Create selection for non-refit tracks" << std::endl;
    fHybridTrackCutsNoItsRefit = new AliESDtrackCuts(*fHybridTrackCutsConstrained);
    fHybridTrackCutsNoItsRefit->SetRequireITSRefit(kFALSE);
  }
}
