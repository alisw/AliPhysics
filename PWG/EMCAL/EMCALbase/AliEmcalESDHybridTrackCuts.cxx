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
#include <algorithm>
#include <functional>
#include <iostream>
#include <string>

#include "TFormula.h"

#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliEmcalESDHybridTrackCuts.h"
#include "AliEmcalTrackSelResultHybrid.h"
#include "AliLog.h"

ClassImp(PWG::EMCAL::AliEmcalESDHybridTrackCuts)

using namespace PWG::EMCAL;

AliEmcalESDHybridTrackCuts::AliEmcalESDHybridTrackCuts():
  AliEmcalCutBase(),
  fLocalInitialized(kFALSE),
  fHybridTrackDefinition(kDef2010),
  fHybridTrackCutsGlobal(nullptr),
  fHybridTrackCutsConstrained(nullptr),
  fHybridTrackCutsNoItsRefit(nullptr),
  fRequireTPCTRDClusters(false),
  fMinClustersTPCTRD(0),
  fPtDepParamClusterCut(0.)
{

}

AliEmcalESDHybridTrackCuts::AliEmcalESDHybridTrackCuts(const char *name, HybridDefinition_t hybriddef):
  AliEmcalCutBase(name, ""),
  fLocalInitialized(kFALSE),
  fHybridTrackDefinition(hybriddef),
  fHybridTrackCutsGlobal(nullptr),
  fHybridTrackCutsConstrained(nullptr),
  fHybridTrackCutsNoItsRefit(nullptr),
  fRequireTPCTRDClusters(false),
  fMinClustersTPCTRD(0),
  fPtDepParamClusterCut(0.)
{

}

AliEmcalESDHybridTrackCuts::~AliEmcalESDHybridTrackCuts(){
  if(fHybridTrackCutsGlobal) delete fHybridTrackCutsGlobal;
  if(fHybridTrackCutsConstrained) delete fHybridTrackCutsConstrained;
  if(fHybridTrackCutsNoItsRefit) delete fHybridTrackCutsNoItsRefit;
}

AliEmcalTrackSelResultPtr AliEmcalESDHybridTrackCuts::IsSelected(TObject *o){
  AliDebugStream(1) << "AliEmcalESDHybridTrackCuts::IsSelected(): Called" << std::endl;
  if(!fLocalInitialized) Init();
  AliDebugStream(1) << "Global cuts:      " << (fHybridTrackCutsGlobal ? "yes" : "no") << std::endl;
  AliDebugStream(1) << "Constrained cuts: " << (fHybridTrackCutsConstrained ? "yes" : "no") << std::endl;
  AliDebugStream(1) << "Non-refit cuts:   " << (fHybridTrackCutsNoItsRefit ? "yes" : "no") << std::endl;
  if(auto esdtrack = dynamic_cast<AliESDtrack *>(o)) {
    AliEmcalTrackSelResultHybrid::HybridType_t tracktype = AliEmcalTrackSelResultHybrid::kUndefined;
    if(fHybridTrackCutsGlobal && fHybridTrackCutsGlobal->AcceptTrack(esdtrack)){
      // Temporary hack for TPC+TRD number of cluster cut which is not yet available in AliESDtrackCuts
      bool isSelected = true;
      auto tracklength = GetTPCTRDNumberOfClusters(esdtrack);
      auto tracklengthcut = GetPtDepCutTPCTRDNumberOfClusters(esdtrack);
      AliDebugStream(3) << "Global: Combined track length: " << tracklength << "(" << esdtrack->GetTPCCrossedRows() << "/" << static_cast<Int_t>(esdtrack->GetTRDncls()) << "), cut " << tracklengthcut << std::endl;
      if(fRequireTPCTRDClusters && (static_cast<Double_t>(tracklength) < tracklengthcut)) isSelected = false;
      if(isSelected){
        AliDebugStream(2) << "Track selected as global hybrid track" << std::endl;
        tracktype = AliEmcalTrackSelResultHybrid::kHybridGlobal;
      }
    } else {
      if(fHybridTrackCutsConstrained && fHybridTrackCutsConstrained->AcceptTrack(esdtrack)){
        // Temporary hack for TPC+TRD number of cluster cut which is not yet available in AliESDtrackCuts
        bool isSelected = true;
        auto tracklength = GetTPCTRDNumberOfClusters(esdtrack);
        auto tracklengthcut = GetPtDepCutTPCTRDNumberOfClusters(esdtrack);
        AliDebugStream(3) << "Constrained: Combined track length: " << tracklength << "(" << esdtrack->GetTPCCrossedRows() << "/" << static_cast<Int_t>(esdtrack->GetTRDncls()) << "), cut " << tracklengthcut << std::endl;
        if(fRequireTPCTRDClusters && (tracklength < tracklengthcut)) isSelected = false;
        if(isSelected){
          AliDebugStream(2) << "Track selected as constrained hybrid track" << std::endl;
          if(IsActiveITSModule(esdtrack, 0) || IsActiveITSModule(esdtrack, 1)) tracktype = AliEmcalTrackSelResultHybrid::kHybridConstrainedFake;
          else tracktype = AliEmcalTrackSelResultHybrid::kHybridConstrainedTrue;
        }
      } else if(fHybridTrackCutsNoItsRefit && fHybridTrackCutsNoItsRefit->AcceptTrack(esdtrack)) {
        // Temporary hack for TPC+TRD number of cluster cut which is not yet available in AliESDtrackCuts
        bool isSelected = true;
        if(fRequireTPCTRDClusters && (GetTPCTRDNumberOfClusters(esdtrack) < GetPtDepCutTPCTRDNumberOfClusters(esdtrack))) isSelected = false;
        if(isSelected){
          AliDebugStream(2) << "Track selected as non-refit hybrid track" << std::endl;
          tracktype = AliEmcalTrackSelResultHybrid::kHybridConstrainedNoITSrefit;
        }
      } else {
        AliDebugStream(2) << "Track not selected as hybrid track" << std::endl;
      }
    }
    AliEmcalTrackSelResultPtr result(esdtrack, tracktype != AliEmcalTrackSelResultHybrid::kUndefined);
    if(result) result.SetUserInfo(new AliEmcalTrackSelResultHybrid(tracktype));
    
    return result;
  }
  AliErrorStream() << "No ESD track" << std::endl;
  return AliEmcalTrackSelResultPtr(nullptr, kFALSE);
}

void AliEmcalESDHybridTrackCuts::Init(){
  switch(fHybridTrackDefinition){
  case kDef2010: InitHybridTracks2010(); break;
  case kDef2011: InitHybridTracks2011(); break;
  case kDef2018TRD: InitHybridTracks2018TRD(); break;
  default:
    AliErrorStream() << "No matching initialization found for requested hybrid track definition" <<std::endl;
  };
  fLocalInitialized = true;
}

Int_t AliEmcalESDHybridTrackCuts::GetTPCTRDNumberOfClusters(const AliVTrack *const trk) const {
  auto esdtrack = static_cast<const AliESDtrack *>(trk);
  return static_cast<Int_t>(trk->GetTPCCrossedRows()) + static_cast<Int_t>(esdtrack->GetTRDntracklets() * 20);
}

Double_t AliEmcalESDHybridTrackCuts::GetPtDepCutTPCTRDNumberOfClusters(const AliVTrack *trk) const {
  return std::min(static_cast<Double_t>(fMinClustersTPCTRD), static_cast<Double_t>(fMinClustersTPCTRD) - fPtDepParamClusterCut / trk->Pt());
}

Bool_t AliEmcalESDHybridTrackCuts::IsActiveITSModule(const AliESDtrack *const trk, int layer) const {
  int det, status;
  Float_t xloc, zloc;
  trk->GetITSModuleIndexInfo(layer, det, status, xloc, zloc);
  // dead status:
  // - dead (2)
  // - skipped (3)
  // - outinz (4)
  // - holeinz (7)
  return !(status == 2 || status == 3 || status == 4 || status == 7);
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

void AliEmcalESDHybridTrackCuts::InitHybridTracks2018TRD(){
  std::cout << "Initialiing hybrid tracks based on the 2018 definition" <<std::endl;
  InitHybridTracks2011();

  // Deactivate TPC crossed rows cut
  if(fHybridTrackCutsGlobal) fHybridTrackCutsGlobal->SetMinNCrossedRowsTPC(30);
  if(fHybridTrackCutsConstrained) fHybridTrackCutsConstrained->SetMinNCrossedRowsTPC(30);
  if(fHybridTrackCutsNoItsRefit) fHybridTrackCutsNoItsRefit->SetMinNCrossedRowsTPC(30);

  // Set min. number of TPC+TRD clusters cut
  fRequireTPCTRDClusters = true;
  fMinClustersTPCTRD = 120;
  fPtDepParamClusterCut = 10.;
}