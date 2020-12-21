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
#include <iostream>

#include <TClonesArray.h>
#include <TBits.h>
#include <TObjArray.h>

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliEmcalAODFilterBitCuts.h"
#include "AliEmcalAODHybridTrackCuts.h"
#include "AliEmcalAODTPCOnlyTrackCuts.h"
#include "AliEmcalCutBase.h"
#include "AliEmcalTrackSelResultCombined.h"
#include "AliEmcalTrackSelResultHybrid.h"
#include "AliEmcalTrackSelectionAOD.h"
#include "AliEmcalVCutsWrapper.h"
#include "AliESDtrackCuts.h"
#include "AliLog.h"
#include "AliPicoTrack.h"

/// \cond CLASSIMP
ClassImp(AliEmcalTrackSelectionAOD)
ClassImp(PWG::EMCAL::TestAliEmcalTrackSelectionAOD)
/// \endcond

using namespace PWG::EMCAL;

AliEmcalTrackSelectionAOD::AliEmcalTrackSelectionAOD() :
	AliEmcalTrackSelection()
{
}

AliEmcalTrackSelectionAOD::AliEmcalTrackSelectionAOD(AliVCuts* cuts, UInt_t filterbits):
	AliEmcalTrackSelection()
{
  if(cuts) AddTrackCuts(cuts);
  if(filterbits) {
    auto filterbitcuts = new PWG::EMCAL::AliEmcalAODFilterBitCuts("filterbitcuts", "AOD filter bit cuts");
    filterbitcuts->SetFilterBits(filterbits, true);
    AddTrackCuts(filterbitcuts);
  }
}

AliEmcalTrackSelectionAOD::AliEmcalTrackSelectionAOD(ETrackFilterType_t type, const char* period):
  AliEmcalTrackSelection()
{
  GenerateTrackCuts(type, period);
}

void AliEmcalTrackSelectionAOD::GenerateTrackCuts(ETrackFilterType_t type, const char* period)
{
  switch (type) {
  case kHybridTracks:
    {
      auto hybridcuts = new PWG::EMCAL::AliEmcalAODHybridTrackCuts("hybridcuts");
      Char_t hybridbits[2];
      GetHybridFilterBits(hybridbits, period);
      hybridcuts->SetHybridFilterBits(hybridbits[0], hybridbits[1]);
      AddTrackCuts(hybridcuts);
      break;
    }

  case kTPCOnlyTracks:
    {
      AddTrackCuts(new PWG::EMCAL::AliEmcalAODTPCOnlyTrackCuts("tpconlycuts", "hybrid track cuts for TPC only tracks"));
      break;
    }

  case AliEmcalTrackSelection::kITSPureTracks:
    {
      if (fListOfCuts) fListOfCuts->Clear();
      AddTrackCuts(AliESDtrackCuts::GetStandardITSSATrackCuts2010());
      break;
    }

  case kHybridTracks2010wNoRefit:
    {
      auto trackcuts = new PWG::EMCAL::AliEmcalAODHybridTrackCuts("hybridcuts2010_wNoRefit");
      Char_t hybridbits[2];
      GetHybridFilterBits(hybridbits, period);
      trackcuts->SetHybridFilterBits(hybridbits[0], hybridbits[1]);
      AddTrackCuts(trackcuts);
      break;
    }

  case kHybridTracks2010woNoRefit:
    {
      auto trackcuts = new PWG::EMCAL::AliEmcalAODHybridTrackCuts("hybridcuts2010_woNoRefit");
      Char_t hybridbits[2];
      GetHybridFilterBits(hybridbits, period);
      trackcuts->SetHybridFilterBits(hybridbits[0], hybridbits[1]);
      trackcuts->SetSelectNonITSrefitTracks(kFALSE);
      AddTrackCuts(trackcuts);
      break;
    }

  case kHybridTracks2011wNoRefit:
    {
      auto trackcuts = new PWG::EMCAL::AliEmcalAODHybridTrackCuts("hybridcuts2011_wNoRefit");
      Char_t hybridbits[2];
      GetHybridFilterBits(hybridbits, period);
      trackcuts->SetHybridFilterBits(hybridbits[0], hybridbits[1]);
      AddTrackCuts(trackcuts);
      break;
    }

  case kHybridTracks2011woNoRefit:
    {
      auto trackcuts = new PWG::EMCAL::AliEmcalAODHybridTrackCuts("hybridcuts2011_woNoRefit");
      Char_t hybridbits[2];
      GetHybridFilterBits(hybridbits, period);
      trackcuts->SetHybridFilterBits(hybridbits[0], hybridbits[1]);
      trackcuts->SetSelectNonITSrefitTracks(kFALSE);
      AddTrackCuts(trackcuts);
      break;
    }

  default:
    break;
  }
}

PWG::EMCAL::AliEmcalTrackSelResultPtr AliEmcalTrackSelectionAOD::IsTrackAccepted(AliVTrack * const trk)
{
  AliAODTrack *aodt = dynamic_cast<AliAODTrack*>(trk);
  if (!aodt){
    AliPicoTrack *picotrack = dynamic_cast<AliPicoTrack*>(trk);
    if(picotrack) {
      aodt = dynamic_cast<AliAODTrack *>(picotrack->GetTrack());
    }
    else {
      AliError("Track neither AOD track nor pico track");
      return PWG::EMCAL::AliEmcalTrackSelResultPtr(nullptr, kFALSE);
    }
  }
  if(!aodt){
    AliError("Failed getting AOD track");
    return PWG::EMCAL::AliEmcalTrackSelResultPtr(nullptr, kFALSE);
  }

  TBits trackbitmap(64);
  trackbitmap.ResetAllBits();
  UInt_t cutcounter(0);
  std::vector<PWG::EMCAL::AliEmcalTrackSelResultPtr> selectionStatus(fListOfCuts->GetEntries());
  if (fListOfCuts) {
    for (auto cutIter : *fListOfCuts){ // @suppress("Symbol is not resolved")
      PWG::EMCAL::AliEmcalCutBase *trackCuts = static_cast<PWG::EMCAL::AliEmcalCutBase*>(static_cast<AliEmcalManagedObject *>(cutIter)->GetObject());
      PWG::EMCAL::AliEmcalTrackSelResultPtr cutresults = trackCuts->IsSelected(aodt);
      if (cutresults) trackbitmap.SetBitNumber(cutcounter);
      selectionStatus.push_back(cutresults);
      cutcounter++;
    }
  }

  PWG::EMCAL::AliEmcalTrackSelResultPtr result(aodt, kFALSE, new PWG::EMCAL::AliEmcalTrackSelResultCombined(selectionStatus));
  if (fSelectionModeAny){
    // In case of ANY one of the cuts need to be fulfilled (equivalent to one but set)
    result.SetSelectionResult(trackbitmap.CountBits() > 0 || cutcounter == 0);
  }
  else {
    // In case of ALL all of the cuts need to be fulfilled (equivalent to all bits set)
    result.SetSelectionResult(trackbitmap.CountBits() == cutcounter);
  }
  return result;
}

void AliEmcalTrackSelectionAOD::AddFilterBit(UInt_t filterbits){
  PWG::EMCAL::AliEmcalAODFilterBitCuts *filtercuts = nullptr;
  // Find existing filter bit cuts
  for(auto c : *fListOfCuts) { // @suppress("Symbol is not resolved")
    if(auto vcutswrapper = dynamic_cast<PWG::EMCAL::AliEmcalVCutsWrapper*>(c)) {
      if(auto aodcuts = dynamic_cast<PWG::EMCAL::AliEmcalAODFilterBitCuts *>(vcutswrapper->GetCutObject())){
        filtercuts = aodcuts;
        break;
      }
    }
  }
  if(filtercuts) filtercuts->SetFilterBits(filterbits, false);
  else {
    filtercuts = new PWG::EMCAL::AliEmcalAODFilterBitCuts("filterbitcuts", "AOD filter bit cuts");
    filtercuts->SetFilterBits(filterbits, kTRUE);
    AddTrackCuts(filtercuts);
  }
}

Bool_t AliEmcalTrackSelectionAOD::GetHybridFilterBits(Char_t bits[], TString period)
{
  period.ToLower();
  if (period == "lhc10b" || period == "lhc10c" || period == "lhc10d" ||
      period == "lhc10e" || period == "lhc10h" ||
      period == "lhc11h" || period == "lhc12a" || period == "lhc12b" ||
      period == "lhc12c" || period == "lhc12d" || period == "lhc12e" ||
      period == "lhc12f" || period == "lhc12g" || period == "lhc12h" ||
      period == "lhc12i" || period == "lhc13b" || period == "lhc13c" ||
      period == "lhc13d" || period == "lhc13e" || period == "lhc13f" ||
      period == "lhc13g" ||

      (period.Length() == 6 && (period.BeginsWith("lhc15") || period.BeginsWith("lhc16") || period.BeginsWith("lhc17") || period.BeginsWith("lhc18"))) // all Run-2 data, excluding MC productions
  ) {
    bits[0] = 8;
    bits[1] = 9;
  }

  else if (period == "lhc10f7a" || period == "lhc12a15e" || period.BeginsWith("lhc12a17") ||
      period == "lhc14j4b" || period == "lhc14j4c" || period == "lhc14j4d" ||
      period == "lhc14j4e" ||
      period == "lhc15i2b" || period == "lhc15i2c" || period == "lhc15i2d" ||
      period == "lhc15i2e" ||
      period == "lhc15g6b" || period == "lhc15g6c" || period == "lhc15g6d" ||
      period == "lhc15g6e" ||
      period == "lhc13b4" || period == "lhc13b4_fix" || period == "lhc13b4_plus" || period == "lhc14k1a" || period == "lhc14k1b" || period == "lhc13e4" ||
      period.BeginsWith("lhc14a1") || period.BeginsWith("lhc13b2_efix") ||
      period.BeginsWith("lhc15g6") || period.BeginsWith("lhc16c2") || period.BeginsWith("lhc16e1") || period.BeginsWith("lhc17f8") || 
      period.BeginsWith("lhc18b8") || period.BeginsWith("lhc18f5") || period.BeginsWith("lhc18g2") || period.BeginsWith("lhc17g8a") || 
      period.BeginsWith("lhc19a1") || period.BeginsWith("lhc19f4")) {
    bits[0] = 8;
    bits[1] = 9;
  }

  else if (period == "lhc11a" || period == "lhc10hold" || period == "lhc11c" || period == "lhc11d") {
    bits[0] = 8;
    bits[1] = 4;
  }

  else if (period.Contains("lhc12a15a") || period == "lhc12a15f" ||
      period == "lhc12a15g" || period.BeginsWith("lhc11a1")) {
    bits[0] = 8;
    bits[1] = 4;
  }

  else {
    ::Error("AliEmcalTrackSelectionAOD::GetHybridFilterBits", "Could not find period %s! Hybrid tracks will be selected, but will not be able to distinguish between global and constrained.", period.Data());
    bits[0] = -1;
    bits[1] = -1;
    return kFALSE;
  }

  return kTRUE;
}

namespace PWG {

namespace EMCAL {

TestAliEmcalTrackSelectionAOD::TestAliEmcalTrackSelectionAOD() :
  TObject(),
  fTrackSelHybrid2010wRefit(nullptr),
  fTrackSelHybrid2010woRefit(nullptr),
  fTrackSelHybrid2011(nullptr),
  fTrackSelTPConly(nullptr)
{

}

TestAliEmcalTrackSelectionAOD::~TestAliEmcalTrackSelectionAOD(){
  if(fTrackSelHybrid2010woRefit) delete fTrackSelHybrid2010woRefit;
  if(fTrackSelHybrid2010wRefit) delete fTrackSelHybrid2010wRefit;
  if(fTrackSelHybrid2011) delete fTrackSelHybrid2011;
  if(fTrackSelTPConly) delete fTrackSelTPConly;
}

void TestAliEmcalTrackSelectionAOD::Init() {
  fTrackSelHybrid2010woRefit = new AliEmcalTrackSelectionAOD(AliEmcalTrackSelection::kHybridTracks2010woNoRefit, "lhc10hold");
  fTrackSelHybrid2010wRefit = new AliEmcalTrackSelectionAOD(AliEmcalTrackSelection::kHybridTracks2010wNoRefit, "lhc10hold");
  fTrackSelHybrid2011 = new AliEmcalTrackSelectionAOD(AliEmcalTrackSelection::kHybridTracks2010woNoRefit, "lhc11h");
  fTrackSelTPConly = new AliEmcalTrackSelectionAOD(AliEmcalTrackSelection::kTPCOnlyTracks);
}

bool TestAliEmcalTrackSelectionAOD::RunAllTests() const {
 return TestHybridDef2010wRefit() && TestHybridDef2010woRefit() && TestHybridDef2011() && TestTPConly();
}

bool TestAliEmcalTrackSelectionAOD::TestHybridDef2010wRefit() const {
  AliInfoStream() << "Running test for 2010 Definition with non-refit tracks" << std::endl;
  AliAODTrack testCat1WithRefit, testCat2WithRefit, testCat2WithoutRefit, testNoHybrid;
  testCat1WithRefit.SetIsHybridGlobalConstrainedGlobal();
  testCat2WithRefit.SetIsHybridGlobalConstrainedGlobal();
  testCat2WithoutRefit.SetIsHybridGlobalConstrainedGlobal();
  testCat1WithRefit.SetStatus(AliVTrack::kITSrefit);
  testCat2WithRefit.SetStatus(AliVTrack::kITSrefit);
  testCat1WithRefit.SetFilterMap(BIT(8));
  testCat2WithRefit.SetFilterMap(BIT(4));
  testCat2WithoutRefit.SetFilterMap(BIT(4));

  int nfailure = 0;
  auto result_cat1_wrefit = fTrackSelHybrid2010wRefit->IsTrackAccepted(&testCat1WithRefit);
  if(!result_cat1_wrefit){
    AliErrorStream() << "Hybrid track CAT1 not selected" << std::endl;
    nfailure++;
  } else {
    auto hybridcat = FindHybridSelectionResult(result_cat1_wrefit);
    if(!hybridcat){
      AliErrorStream() << "No hybrid selection result found for CAT1 hybrid track" << std::endl;
      nfailure++;
    } else {
      if(hybridcat->GetHybridTrackType() != AliEmcalTrackSelResultHybrid::kHybridGlobal) {
        AliErrorStream() << "Incorrect hybrid track type for CAT1 hybrid track: " << hybridcat->GetHybridTrackType() << std::endl;
        nfailure++;
      }
    }
  }

  auto result_cat2_wrefit = fTrackSelHybrid2010wRefit->IsTrackAccepted(&testCat2WithRefit);
  if(!result_cat2_wrefit){
    AliErrorStream() << "Hybrid track CAT2 not selected" << std::endl;
    nfailure++;
  } else {
    auto hybridcat = FindHybridSelectionResult(result_cat2_wrefit);
    if(!hybridcat){
      AliErrorStream() << "No hybrid selection result found for CAT2 hybrid track" << std::endl;
      nfailure++;
    } else {
      if(hybridcat->IsHybridTrackConstrained()) {
        AliErrorStream() << "Incorrect hybrid track type for CAT2 hybrid track: " << hybridcat->GetHybridTrackType() << std::endl;
        nfailure++;
      }
    }
  }

  auto result_cat2_worefit = fTrackSelHybrid2010wRefit->IsTrackAccepted(&testCat2WithoutRefit);
  if(!result_cat2_worefit){
    AliErrorStream() << "Hybrid track CAT3 not selected" << std::endl;
    nfailure++;
  } else {
    auto hybridcat = FindHybridSelectionResult(result_cat2_worefit);
    if(!hybridcat){
      AliErrorStream() << "No hybrid selection result found for CAT3 hybrid track" << std::endl;
      nfailure++;
    } else {
      if(hybridcat->GetHybridTrackType() != AliEmcalTrackSelResultHybrid::kHybridConstrainedNoITSrefit) {
        AliErrorStream() << "Incorrect hybrid track type for CAT3 hybrid track: " << hybridcat->GetHybridTrackType() << std::endl;
        nfailure++;
      }
    }
  }

  auto result_nohybrid = fTrackSelHybrid2010wRefit->IsTrackAccepted(&testNoHybrid);
  if(result_nohybrid){
    AliErrorStream() << "Non-hybrid track selected as hybrid track" << std::endl;
    nfailure++;
  }

  return nfailure == 0;
}

bool TestAliEmcalTrackSelectionAOD::TestHybridDef2010woRefit() const {
  AliInfoStream() << "Running test for 2010 Definition without non-refit tracks" << std::endl;
  AliAODTrack testCat1WithRefit, testCat2WithRefit, testCat2WithoutRefit, testNoHybrid;
  testCat1WithRefit.SetIsHybridGlobalConstrainedGlobal();
  testCat2WithRefit.SetIsHybridGlobalConstrainedGlobal();
  testCat2WithoutRefit.SetIsHybridGlobalConstrainedGlobal();
  testCat1WithRefit.SetStatus(AliVTrack::kITSrefit);
  testCat2WithRefit.SetStatus(AliVTrack::kITSrefit);
  testCat1WithRefit.SetFilterMap(BIT(8));
  testCat2WithRefit.SetFilterMap(BIT(4));
  testCat2WithoutRefit.SetFilterMap(BIT(4));

  int nfailure = 0;
  auto result_cat1_wrefit = fTrackSelHybrid2010woRefit->IsTrackAccepted(&testCat1WithRefit);
  if(!result_cat1_wrefit){
    AliErrorStream() << "Hybrid track CAT1 not selected" << std::endl;
    nfailure++;
  } else {
    auto hybridcat = FindHybridSelectionResult(result_cat1_wrefit);
    if(!hybridcat){
      AliErrorStream() << "No hybrid selection result found for CAT1 hybrid track" << std::endl;
      nfailure++;
    } else {
      if(hybridcat->GetHybridTrackType() != AliEmcalTrackSelResultHybrid::kHybridGlobal) {
        AliErrorStream() << "Incorrect hybrid track type for CAT1 hybrid track: " << hybridcat->GetHybridTrackType() << std::endl;
        nfailure++;
      }
    }
  }

  auto result_cat2_wrefit = fTrackSelHybrid2010woRefit->IsTrackAccepted(&testCat2WithRefit);
  if(!result_cat2_wrefit){
    AliErrorStream() << "Hybrid track CAT2 not selected" << std::endl;
    nfailure++;
  } else {
    auto hybridcat = FindHybridSelectionResult(result_cat2_wrefit);
    if(!hybridcat){
      AliErrorStream() << "No hybrid selection result found for CAT2 hybrid track" << std::endl;
      nfailure++;
    } else {
      if(hybridcat->IsHybridTrackConstrained()) {
        AliErrorStream() << "Incorrect hybrid track type for CAT2 hybrid track: " << hybridcat->GetHybridTrackType() << std::endl;
        nfailure++;
      }
    }
  }
  
  auto result_cat2_worefit = fTrackSelHybrid2010woRefit->IsTrackAccepted(&testCat2WithoutRefit);
  if(result_cat2_worefit){
    AliErrorStream() << "CAT2 track without refit selected as hybrid track in track selection excluding non-refit tracks" << std::endl;
    nfailure++;
  }
 
  auto result_nohybrid = fTrackSelHybrid2010woRefit->IsTrackAccepted(&testNoHybrid);
  if(result_nohybrid){
    AliErrorStream() << "Non-hybrid track selected as hybrid track" << std::endl;
    nfailure++;
  }

  return nfailure == 0;
}

bool TestAliEmcalTrackSelectionAOD::TestHybridDef2011() const {
  AliInfoStream() << "Running test for 2011 Definition" << std::endl;
  AliAODTrack testCat1, testCat2, testNoHybrid;
  testCat1.SetIsHybridGlobalConstrainedGlobal();
  testCat2.SetIsHybridGlobalConstrainedGlobal();
  testCat1.SetStatus(AliVTrack::kITSrefit);
  testCat2.SetStatus(AliVTrack::kITSrefit);
  testCat1.SetFilterMap(BIT(8));
  testCat2.SetFilterMap(BIT(9));

  int nfailure = 0;
  
  auto result_cat1 = fTrackSelHybrid2011->IsTrackAccepted(&testCat1);
  if(!result_cat1){
    AliErrorStream() << "Hybrid track CAT1 not selected" << std::endl;
    nfailure++;
  } else {
    auto hybridcat = FindHybridSelectionResult(result_cat1);
    if(!hybridcat){
      AliErrorStream() << "No hybrid selection result found for CAT1 hybrid track" << std::endl;
      nfailure++;
    } else {
      if(hybridcat->GetHybridTrackType() != AliEmcalTrackSelResultHybrid::kHybridGlobal) {
        AliErrorStream() << "Incorrect hybrid track type for CAT1 hybrid track: " << hybridcat->GetHybridTrackType() << std::endl;
        nfailure++;
      }
    }
  }

  auto result_cat2 = fTrackSelHybrid2011->IsTrackAccepted(&testCat2);
  if(!result_cat2){
    AliErrorStream() << "Hybrid track CAT2 not selected" << std::endl;
    nfailure++;
  } else {
    auto hybridcat = FindHybridSelectionResult(result_cat2);
    if(!hybridcat){
      AliErrorStream() << "No hybrid selection result found for CAT2 hybrid track" << std::endl;
      nfailure++;
    } else {
      if(hybridcat->IsHybridTrackConstrained()) {
        AliErrorStream() << "Incorrect hybrid track type for CAT2 hybrid track: " << hybridcat->GetHybridTrackType() << std::endl;
        nfailure++;
      }
    }
  }

  auto result_nohybrid = fTrackSelHybrid2011->IsTrackAccepted(&testNoHybrid);
  if(result_nohybrid){
    AliErrorStream() << "Non-hybrid track selected as hybrid track" << std::endl;
    nfailure++;
  }
  return nfailure == 0;
}

bool TestAliEmcalTrackSelectionAOD::TestTPConly() const {
  AliInfoStream() << "Running test for TPC-only tracks" << std::endl;
  AliAODTrack testtrackTrue, testtrackFalse;
  testtrackTrue.SetIsHybridTPCConstrainedGlobal(true);

  int nfailure = 0;
  auto result_true = fTrackSelTPConly->IsTrackAccepted(&testtrackTrue);
  if(!result_true) {
    AliErrorStream() << "TPC-only constrained track rejected" << std::endl;
    nfailure++;
  }

  auto result_false = fTrackSelTPConly->IsTrackAccepted(&testtrackFalse);
  if(result_false) {
    AliErrorStream() << "Non-TPC-only track selected as constrained track" << std::endl;
    nfailure++;
  }

  return nfailure == 0;
}

const AliEmcalTrackSelResultHybrid *TestAliEmcalTrackSelectionAOD::FindHybridSelectionResult(const AliEmcalTrackSelResultPtr &data) const {
  if(!data.GetUserInfo()) return nullptr;
  if(auto hybridinfo = dynamic_cast<const AliEmcalTrackSelResultHybrid *>(data.GetUserInfo())) return hybridinfo;
  if(auto combinedinfo = dynamic_cast<const AliEmcalTrackSelResultCombined *>(data.GetUserInfo())) {
    for(int i = 0; i < combinedinfo->GetNumberOfSelectionResults(); i++) {
      try{
        auto res = FindHybridSelectionResult((*combinedinfo)[i]);
        if(res) return res;
      } catch (AliEmcalTrackSelResultCombined::IndexException &e) {
        // just go on
      }
    }
  }  
  return nullptr;
}

}

}
