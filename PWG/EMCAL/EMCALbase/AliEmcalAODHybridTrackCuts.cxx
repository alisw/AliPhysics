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
#include "AliAODTrack.h"
#include "AliEmcalAODHybridTrackCuts.h"
#include "AliEmcalTrackSelResultHybrid.h"
#include "AliLog.h"
#include <iostream>

/// \cond CLASSIMP
ClassImp(PWG::EMCAL::AliEmcalAODHybridTrackCuts)
/// \endcond

using namespace PWG::EMCAL;

AliEmcalAODHybridTrackCuts::AliEmcalAODHybridTrackCuts():
  AliEmcalCutBase(),
  fSelectNonITSrefitTracks(kTRUE)
{
  fHybridFilterBits[0] = -1;
  fHybridFilterBits[1] = -1;
}

AliEmcalAODHybridTrackCuts::AliEmcalAODHybridTrackCuts(const char *name):
  AliEmcalCutBase(name,""),
  fSelectNonITSrefitTracks(kTRUE)
{
  
}

AliEmcalTrackSelResultPtr AliEmcalAODHybridTrackCuts::IsSelected(TObject *o){
  AliAODTrack *aodtrack = dynamic_cast<AliAODTrack *>(o);
  if(!aodtrack) return AliEmcalTrackSelResultPtr(nullptr, kFALSE);
  bool selectionresult = aodtrack->IsHybridGlobalConstrainedGlobal();
  // Reject non-ITSrefit tracks if requested
  if((fSelectNonITSrefitTracks == false) && (!(aodtrack->GetStatus() & AliVTrack::kITSrefit))) selectionresult = false;
  AliEmcalTrackSelResultPtr result(aodtrack, selectionresult);
  // Create user object defining the hybrid track type (only in case the object is selected as hybrid track)
  if(selectionresult){
    AliEmcalTrackSelResultHybrid::HybridType_t tracktype = AliEmcalTrackSelResultHybrid::kHybridGlobal;
    if(fHybridFilterBits[0] > -1 && fHybridFilterBits[1] > -1) {
      if(aodtrack->TestFilterBit(BIT(fHybridFilterBits[0]))) tracktype = AliEmcalTrackSelResultHybrid::kHybridGlobal;
      else if(aodtrack->TestFilterBit(BIT(fHybridFilterBits[1]))){
        if(aodtrack->GetStatus() & AliVTrack::kITSrefit) tracktype = AliEmcalTrackSelResultHybrid::kHybridConstrainedTrue; // no module map -> set all complementary hybrid tracks to true complementary hybrid tracks
        else tracktype = AliEmcalTrackSelResultHybrid::kHybridConstrainedNoITSrefit;
      }
    }
    result.SetUserInfo(new AliEmcalTrackSelResultHybrid(tracktype));
  }
  return result;
}

TestAliEmcalAODHybridTrackCuts::TestAliEmcalAODHybridTrackCuts():
  TObject(),
  fDef2010wRefit(nullptr),
  fDef2010woRefit(nullptr),
  fDef2011(nullptr)
{

}

TestAliEmcalAODHybridTrackCuts::~TestAliEmcalAODHybridTrackCuts(){
  if(fDef2010wRefit) delete fDef2010wRefit;
  if(fDef2010woRefit) delete fDef2010woRefit;
  if(fDef2011) delete fDef2011;
}

void TestAliEmcalAODHybridTrackCuts::Init(){
  fDef2010wRefit = new AliEmcalAODHybridTrackCuts("def2010wNoRefit");
  fDef2010wRefit->SetSelectNonITSrefitTracks(kTRUE);
  fDef2010wRefit->SetHybridFilterBits(8,4);
  fDef2010woRefit = new AliEmcalAODHybridTrackCuts("def2010woNoRefit");
  fDef2010woRefit->SetSelectNonITSrefitTracks(kFALSE);
  fDef2010woRefit->SetHybridFilterBits(8,4);
  fDef2011 = new AliEmcalAODHybridTrackCuts("def2011");
  fDef2011->SetHybridFilterBits(8,9);
}

bool TestAliEmcalAODHybridTrackCuts::RunAllTests() const {
  return TestDef2010wRefit() && TestDef2010woRefit() && TestDef2011();
}

bool TestAliEmcalAODHybridTrackCuts::TestDef2010wRefit() const {
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
  auto result_cat1_wrefit = fDef2010wRefit->IsSelected(&testCat1WithRefit);
  if(!result_cat1_wrefit){
    AliErrorStream() << "Track CAT1 not selected as hybrid track" << std::endl;
    nfailure++;   // track not selected - failure
  } else {
    auto tracktype = dynamic_cast<const AliEmcalTrackSelResultHybrid *>(result_cat1_wrefit.GetUserInfo());
    if(!tracktype){ 
      AliErrorStream() << "Hybrid track information not found for track CAT1" << std::endl;
      nfailure++;          // no hybrid track type found - failure
    } else {
      if(tracktype->GetHybridTrackType() != AliEmcalTrackSelResultHybrid::kHybridGlobal){
        AliErrorStream() << "Track not selected as hybrid track CAT1: " << int(tracktype->GetHybridTrackType()) << std::endl;
        nfailure++; // wrong hybrid track type
      }
    }
  }

  auto result_cat2_wrefit = fDef2010wRefit->IsSelected(&testCat2WithRefit);
  if(!result_cat2_wrefit){
    AliErrorStream() << "Track CAT2 not selected as hybrid track" << std::endl;
    nfailure++;  // track not selected - failure
  } else {
    auto tracktype = dynamic_cast<const AliEmcalTrackSelResultHybrid *>(result_cat2_wrefit.GetUserInfo());
    if(!tracktype){
      AliErrorStream() << "Hybrid track information not found for track CAT2" << std::endl;
      nfailure++;          // no hybrid track type found - failure
    } else {
      if(tracktype->IsHybridTrackConstrained()){
        AliErrorStream() << "Track not selected as hybrid track CAT2: " << int(tracktype->GetHybridTrackType()) << std::endl;
        nfailure++; // wrong hybrid track type
      }
    }
  }

  auto result_cat2_worefit = fDef2010wRefit->IsSelected(&testCat2WithoutRefit);
  if(!result_cat2_worefit){
    AliErrorStream() << "Track CAT3 not selected as hybrid track" << std::endl;
    nfailure++;  // track not selected - failure
  } else {
    auto tracktype = dynamic_cast<const AliEmcalTrackSelResultHybrid *>(result_cat2_worefit.GetUserInfo());
    if(!tracktype){
      AliErrorStream() << "Hybrid track information not found for track CAT3" << std::endl;
      nfailure++;          // no hybrid track type found - failure
    } else {
      if(tracktype->GetHybridTrackType() != AliEmcalTrackSelResultHybrid::kHybridConstrainedNoITSrefit){
        AliErrorStream() << "Track not selected as hybrid track CAT3: " << int(tracktype->GetHybridTrackType()) << std::endl;
        nfailure++; // wrong hybrid track type
      }
    }
  }

  auto result_nohybrid = fDef2010wRefit->IsSelected(&testNoHybrid);
  if(result_nohybrid || result_nohybrid.GetUserInfo()){
    AliErrorStream() << "Non-hybrid track selected or user object attached " << std::endl;
    nfailure++;  // Non-hybrid track selected as hybrid track
  }

  return nfailure == 0;
}

bool TestAliEmcalAODHybridTrackCuts::TestDef2010woRefit() const {
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
  auto result_cat1_wrefit = fDef2010woRefit->IsSelected(&testCat1WithRefit);
  if(!result_cat1_wrefit) {
    AliErrorStream() << "Track CAT1 not selected as hybrid track" << std::endl;
    nfailure++;   // track not selected - failure
  } else {
    auto tracktype = dynamic_cast<const AliEmcalTrackSelResultHybrid *>(result_cat1_wrefit.GetUserInfo());
    if(!tracktype) {
      AliErrorStream() << "Hybrid track information not found for track CAT1" << std::endl;
      nfailure++;          // no hybrid track type found - failure
    } else {
      if(tracktype->GetHybridTrackType() != AliEmcalTrackSelResultHybrid::kHybridGlobal){
        AliErrorStream() << "Track not selected as hybrid track CAT1: " << int(tracktype->GetHybridTrackType()) << std::endl;
        nfailure++; // wrong hybrid track type
      } 
    }
  }

  auto result_cat2_wrefit = fDef2010woRefit->IsSelected(&testCat2WithRefit);
  if(!result_cat2_wrefit){
    AliErrorStream() << "Track CAT2 not selected as hybrid track" << std::endl;
    nfailure++;  // track not selected - failure
  } else {
    auto tracktype = dynamic_cast<const AliEmcalTrackSelResultHybrid *>(result_cat2_wrefit.GetUserInfo());
    if(!tracktype){
      AliErrorStream() << "Hybrid track information not found for track CAT2" << std::endl;
      nfailure++;          // no hybrid track type found - failure
    } else {
      if(tracktype->IsHybridTrackConstrained()){
        AliErrorStream() << "Track not selected as hybrid track CAT2: " << int(tracktype->GetHybridTrackType()) << std::endl;
        nfailure++; // wrong hybrid track type
      }
    }
  }

  auto result_cat2_worefit = fDef2010woRefit->IsSelected(&testCat2WithoutRefit);
  if(result_cat2_worefit || result_cat2_worefit.GetUserInfo()) {
    AliErrorStream() << "Non-refit track selected or user object attached " << std::endl;
    nfailure++;  // hybrid tracks without refit rejected in this test
  }

  auto result_nohybrid = fDef2010woRefit->IsSelected(&testNoHybrid);
  if(result_nohybrid || result_nohybrid.GetUserInfo()){
    AliErrorStream() << "Non-hybrid track selected or user object attached " << std::endl;
    nfailure++;  // Non-hybrid track selected as hybrid track
  }

  return nfailure == 0;
}

bool TestAliEmcalAODHybridTrackCuts::TestDef2011() const {
  AliInfoStream() << "Running test for 2011 Definition" << std::endl;
  AliAODTrack testCat1, testCat2, testNoHybrid;
  testCat1.SetIsHybridGlobalConstrainedGlobal(kTRUE);
  testCat2.SetIsHybridGlobalConstrainedGlobal(kTRUE);
  testCat1.SetStatus(AliVTrack::kITSrefit);
  testCat2.SetStatus(AliVTrack::kITSrefit);
  testCat1.SetFilterMap(BIT(8));
  testCat2.SetFilterMap(BIT(9));

  int nfailure = 0;

  auto result_cat1 = fDef2011->IsSelected(&testCat1);
  if(!result_cat1){
    AliErrorStream() << "Track CAT1 not selected as hybrid track" << std::endl;
    nfailure++;      // track not selected - failure
  } else {
    auto tracktype = dynamic_cast<const AliEmcalTrackSelResultHybrid *>(result_cat1.GetUserInfo());
    if(!tracktype){
      AliErrorStream() << "Hybrid track information not found for track CAT1" << std::endl;
      nfailure++;
    } else {
      if(tracktype->GetHybridTrackType() != AliEmcalTrackSelResultHybrid::kHybridGlobal){ 
        AliErrorStream() << "Track not selected as hybrid track CAT1: " << int(tracktype->GetHybridTrackType()) << std::endl;
        nfailure++;
      }
    }
  }

  auto result_cat2 = fDef2011->IsSelected(&testCat2);
  if(!result_cat2){
    AliErrorStream() << "Track CAT2 not selected as hybrid track" << std::endl;
    nfailure++;      // track not selected - failure
  } else {
    auto tracktype = dynamic_cast<const AliEmcalTrackSelResultHybrid *>(result_cat2.GetUserInfo());
    if(!tracktype){
      AliErrorStream() << "Hybrid track information not found for track CAT2" << std::endl;
      nfailure++;
    } else {
      if(tracktype->IsHybridTrackConstrained()){
        AliErrorStream() << "Track not selected as hybrid track CAT2: " << int(tracktype->GetHybridTrackType()) << std::endl;
        nfailure++;
      }
    }
  }

  auto result_nohybrid = fDef2011->IsSelected(&testNoHybrid);
  if(result_nohybrid || result_nohybrid.GetUserInfo()){
    AliErrorStream() << "Non-hybrid track selected or user object attached " << std::endl;
    nfailure++;
  }

  return nfailure == 0;
}
