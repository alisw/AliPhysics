// clang-format off
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
// clang-format on
#include <iostream>

#include <THistManager.h>
#include <TObjArray.h>

#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliAnalysisManager.h"
#include "AliESDtrackCuts.h"
#include "AliEmcalTrackSelectionAOD.h"
#include "AliLog.h"

#include "TestAliEmcalTrackSelection.h"

/// \cond CLASSIMP
ClassImp(PWG::EMCAL::TestAliEmcalTrackSelection);
ClassImp(PWG::EMCAL::TestImplAliEmcalTrackSelection);
/// \endcond

namespace PWG {
namespace EMCAL {
TestAliEmcalTrackSelection::TestAliEmcalTrackSelection()
  : AliAnalysisTaskEmcalLight(), fPeriod(), fTestSuite(nullptr), fTestResults(nullptr) {}

TestAliEmcalTrackSelection::TestAliEmcalTrackSelection(const char* name)
  : AliAnalysisTaskEmcalLight(name, kTRUE), fPeriod(), fTestSuite(nullptr), fTestResults(nullptr) {}

TestAliEmcalTrackSelection::~TestAliEmcalTrackSelection() {
  if (fTestSuite)
    delete fTestSuite;
}

void TestAliEmcalTrackSelection::UserCreateOutputObjects() {
  AliAnalysisTaskEmcalLight::UserCreateOutputObjects();

  GenerateTestSuite();

  fTestResults = new THistManager("testresults");
  for (auto test : *fTestSuite)
    fTestResults->CreateTH1(Form("TestStatus%s", test->GetName()), Form("Test Status %s", test->GetName()), 2, -0.5,
                            1.5);
  for (auto hist : *fTestResults->GetListOfHistograms())
    fOutput->Add(hist);

  PostData(1, fOutput);
}

bool TestAliEmcalTrackSelection::Run() {
  for (auto test : *fTestSuite) {
    fTestResults->FillTH1(Form("TestStatus%s", test->GetName()),
                          EvaluateTest(static_cast<TestImplAliEmcalTrackSelection*>(test)) ? 1 : 0);
  }
  return kTRUE;
}

void TestAliEmcalTrackSelection::AddTestImpl(TestImplAliEmcalTrackSelection* test) {
  if (!fTestSuite) {
    fTestSuite = new TObjArray;
    fTestSuite->SetOwner(true);
  }
  AliInfoStream() << "Adding test " << test->GetName() << std::endl;
  fTestSuite->Add(test);
}

void TestAliEmcalTrackSelection::GenerateTestSuite() {
  AddTestImpl(new TestImplAliEmcalTrackSelectionHybrid("hybrid", fPeriod));
  AddTestImpl(new TestImplAliEmcalTrackSelectionTPConly("tpconly", fPeriod));
  AddTestImpl(new TestImplAliEmcalTrackSelectionITSpure("itspure", fPeriod));
}

bool TestAliEmcalTrackSelection::EvaluateTest(TestImplAliEmcalTrackSelection* test) {
  AliInfoStream() << "Evaluating test " << test->GetName() << std::endl;
  return test->RunTest(static_cast<AliAODEvent*>(fInputEvent));
}

TestAliEmcalTrackSelection* TestAliEmcalTrackSelection::AddTestAliEmcalTrackSelection(const char* name) {
  AliAnalysisManager* mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    std::cerr << "TestAliEmcalTrackSelection::AddTestEmcalTrackSelection: No "
                 "analysis manager found. Not adding test!\n";
    return nullptr;
  }

  TestAliEmcalTrackSelection* test = new TestAliEmcalTrackSelection(name);
  mgr->AddTask(test);

  TString outputdir(mgr->GetCommonFileName());
  outputdir += ":TestResults" + TString(name);

  mgr->ConnectInput(test, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(
    test, 1,
    mgr->CreateContainer(Form("TestResults%s", name), TList::Class(), AliAnalysisManager::kOutputContainer, outputdir));

  return test;
}

TestImplAliEmcalTrackSelection::TestImplAliEmcalTrackSelection() : TNamed(), fTrackSelection(nullptr) {}

TestImplAliEmcalTrackSelection::TestImplAliEmcalTrackSelection(const char* name)
  : TNamed(name, ""), fTrackSelection(nullptr) {}

TestImplAliEmcalTrackSelection::~TestImplAliEmcalTrackSelection() {
  if (fTrackSelection)
    delete fTrackSelection;
}

bool TestImplAliEmcalTrackSelection::RunTest(const AliAODEvent* const ev) {
  int nfailure(0);
  for (int itrk = 0; itrk < ev->GetNumberOfTracks(); itrk++) {
    auto track = static_cast<AliAODTrack*>(ev->GetTrack(itrk));
    auto sel = fTrackSelection->IsTrackAccepted(track);
    auto truth = IsTrueTrack(track);
    if (sel != truth)
      nfailure++;
  }

  return nfailure == 0;
}

TestImplAliEmcalTrackSelectionITSpure::TestImplAliEmcalTrackSelectionITSpure()
  : TestImplAliEmcalTrackSelection(), fRefCuts(nullptr) {}

TestImplAliEmcalTrackSelectionITSpure::TestImplAliEmcalTrackSelectionITSpure(const char* name, const char* period)
  : TestImplAliEmcalTrackSelection(name), fRefCuts(nullptr) {
  fTrackSelection = new AliEmcalTrackSelectionAOD;
  fTrackSelection->GenerateTrackCuts(AliEmcalTrackSelection::kITSPureTracks, period);
  fRefCuts = AliESDtrackCuts::GetStandardITSPureSATrackCuts2010(true);
}

TestImplAliEmcalTrackSelectionITSpure::~TestImplAliEmcalTrackSelectionITSpure() {
  if (fRefCuts)
    delete fRefCuts;
}

bool TestImplAliEmcalTrackSelectionITSpure::IsTrueTrack(const AliAODTrack* const trk) const {
  if (!(trk->GetStatus() & AliVTrack::kITSpureSA))
    return false; // addditional test - make sure track is a true stand alone
  // track
  if (!(trk->GetStatus() & AliVTrack::kITSrefit))
    return false;

  if (!fRefCuts->AcceptVTrack(trk))
    return false;
  return true;
}

TestImplAliEmcalTrackSelectionHybrid::TestImplAliEmcalTrackSelectionHybrid(const char* name, const char* period)
  : TestImplAliEmcalTrackSelection(name) {
  fTrackSelection = new AliEmcalTrackSelectionAOD;
  fTrackSelection->GenerateTrackCuts(AliEmcalTrackSelection::kHybridTracks, period);
}

bool TestImplAliEmcalTrackSelectionHybrid::IsTrueTrack(const AliAODTrack* const track) const {
  return track->IsHybridGlobalConstrainedGlobal();
}

TestImplAliEmcalTrackSelectionTPConly::TestImplAliEmcalTrackSelectionTPConly(const char* name, const char* period)
  : TestImplAliEmcalTrackSelection(name) {
  fTrackSelection = new AliEmcalTrackSelectionAOD;
  fTrackSelection->GenerateTrackCuts(AliEmcalTrackSelection::kTPCOnlyTracks, period);
}

bool TestImplAliEmcalTrackSelectionTPConly::IsTrueTrack(const AliAODTrack* const track) const {
  return track->IsHybridTPCConstrainedGlobal();
}
}
}
