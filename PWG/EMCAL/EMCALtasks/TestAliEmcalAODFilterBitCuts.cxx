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
#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliAODTrack.h"
#include "AliEmcalAODFilterBitCuts.h"
#include "AliLog.h"
#include <THistManager.h>
#include <TList.h>
#include <TObjArray.h>
#include <TString.h>
#include <iostream>

#include "TestAliEmcalAODFilterBitCuts.h"

/// \cond CLASSIMP
ClassImp(PWG::EMCAL::TestAliEmcalAODFilterBitCuts)
ClassImp(PWG::EMCAL::TestImplAliEmcalAODFilterBitCuts)
ClassImp(PWG::EMCAL::TestImplAliEmcalAODFilterBitCutsTPCconstrained)
ClassImp(PWG::EMCAL::TestImplAliEmcalAODFilterBitCutsHybrid)
/// \endcond

namespace PWG {
namespace EMCAL {

TestAliEmcalAODFilterBitCuts::TestAliEmcalAODFilterBitCuts() :
    AliAnalysisTaskEmcalLight(),
    fTestSuite(nullptr),
    fTestStatus(nullptr)
{

}

TestAliEmcalAODFilterBitCuts::TestAliEmcalAODFilterBitCuts(const char *name) :
    AliAnalysisTaskEmcalLight(name, kTRUE),
    fTestSuite(nullptr),
    fTestStatus(nullptr)
{
  GenerateTestSuite();
}

TestAliEmcalAODFilterBitCuts::~TestAliEmcalAODFilterBitCuts() {
  if(fTestSuite) delete fTestSuite;
}

void TestAliEmcalAODFilterBitCuts::AddTestImpl(TestImplAliEmcalAODFilterBitCuts *test) {
  if(!fTestSuite) fTestSuite = new TObjArray;
  fTestSuite->SetOwner(kTRUE);
  AliInfoStream() << "Adding test: " << test->GetName() << std::endl;
  fTestSuite->Add(test);
}

void TestAliEmcalAODFilterBitCuts::UserCreateOutputObjects(){
  AliAnalysisTaskEmcalLight::UserCreateOutputObjects();

  fTestStatus = new THistManager("testStatus");
  for(auto t : *fTestSuite) fTestStatus->CreateTH1(Form("TestResult%s", t->GetName()), Form("Result Test %s", t->GetName()), 2, -0.5, 1.5);
  for(auto h : *fTestStatus->GetListOfHistograms()) fOutput->Add(h);

  PostData(1, fOutput);
}

bool TestAliEmcalAODFilterBitCuts::Run(){
  AliAODEvent *ev = dynamic_cast<AliAODEvent *>(fInputEvent);
  if(!ev) return false;
  for(auto test : *fTestSuite) EvaluateTest(static_cast<TestImplAliEmcalAODFilterBitCuts *>(test), ev);
  return true;
}

void TestAliEmcalAODFilterBitCuts::EvaluateTest(TestImplAliEmcalAODFilterBitCuts *test, const AliAODEvent *const ev) {
  AliInfoStream() << "Running test: " << test->GetName() << std::endl;
  fTestStatus->FillTH1(Form("TestResult%s", test->GetName()), test->RunTest(ev) ? 1 : 0);
}

void TestAliEmcalAODFilterBitCuts::GenerateTestSuite(){
  AddTestImpl(new TestImplAliEmcalAODFilterBitCutsTPCconstrained("TPCconstrained"));
  AddTestImpl(new TestImplAliEmcalAODFilterBitCutsHybrid("Hybridall"));
}

TestAliEmcalAODFilterBitCuts *TestAliEmcalAODFilterBitCuts::AddTestAliEmcalAODFilterBitCuts(const char *testname){
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    std::cerr << "TestAliEmcalAODFilterBitCuts: No Analysis manager\n";
    return nullptr;
  }

  TestAliEmcalAODFilterBitCuts *test = new TestAliEmcalAODFilterBitCuts(testname);
  mgr->AddTask(test);

  TString outname(mgr->GetCommonFileName());
  outname += TString::Format(":%s", testname);

  mgr->ConnectInput(test, 0, mgr->GetCommonInputContainer());
  mgr->ConnectOutput(test, 1, mgr->CreateContainer(Form("Histos%s", testname), TList::Class(), AliAnalysisManager::kOutputContainer, outname));
  return test;
}

TestImplAliEmcalAODFilterBitCuts::TestImplAliEmcalAODFilterBitCuts():
    TNamed(),
    fTestObject(nullptr)
{
}

TestImplAliEmcalAODFilterBitCuts::TestImplAliEmcalAODFilterBitCuts(const char *name):
    TNamed(name, ""),
    fTestObject(nullptr)
{
}

TestImplAliEmcalAODFilterBitCuts::~TestImplAliEmcalAODFilterBitCuts(){
  if(fTestObject) delete fTestObject;
}

bool TestImplAliEmcalAODFilterBitCuts::RunTest(const AliAODEvent *const ev) {
  int nerror(0);
  for(int itrk = 0; itrk < ev->GetNumberOfTracks(); itrk++){
    AliAODTrack *trk = static_cast<AliAODTrack *>(ev->GetTrack(itrk));
    bool truth = this->IsTrueTrack(trk), test = fTestObject->IsSelected(trk);
    if(truth != test) nerror++;
  }
  return nerror != 0;
}

TestImplAliEmcalAODFilterBitCutsTPCconstrained::TestImplAliEmcalAODFilterBitCutsTPCconstrained(const char *name):
    TestImplAliEmcalAODFilterBitCuts(name)
{
  fTestObject = new AliEmcalAODFilterBitCuts(Form("CutsTest%s", name), "Test cuts");
  fTestObject->SetStatusBits(static_cast<ULong_t>(AliAODTrack::kIsTPCConstrained));
}

bool TestImplAliEmcalAODFilterBitCutsTPCconstrained::IsTrueTrack(const AliAODTrack *track) const {
  return track->IsTPCConstrained();
}

TestImplAliEmcalAODFilterBitCutsHybrid::TestImplAliEmcalAODFilterBitCutsHybrid(const char *name):
    TestImplAliEmcalAODFilterBitCuts(name)
{
  fTestObject = new AliEmcalAODFilterBitCuts(Form("CutsTest%s", name), "Test cuts");
  fTestObject->SetStatusBits(static_cast<ULong_t>(AliAODTrack::kIsHybridGCG | AliAODTrack::kIsHybridTPCCG));
  fTestObject->SetSelectionMode(AliEmcalAODFilterBitCuts::kSelAny);
}

bool TestImplAliEmcalAODFilterBitCutsHybrid::IsTrueTrack(const AliAODTrack *track) const {
  return track->IsHybridGlobalConstrainedGlobal() || track->IsHybridTPCConstrainedGlobal();
}

} /* namespace EMCAL */
} /* namespace PWG */
