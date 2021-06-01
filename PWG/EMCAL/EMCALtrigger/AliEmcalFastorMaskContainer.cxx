/************************************************************************************
 * Copyright (C) 2021, Copyright Holders of the ALICE Collaboration                 *
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
#include <iostream>
#include <iterator>
#include <memory>
#include <random>
#include <string>

#include <TBrowser.h>
#include <TFile.h>
#include <TH2.h>
#include <TSystem.h>    // For I/O unit test

#include "AliEMCALTriggerMappingV1.h"
#include "AliEMCALTriggerMappingV2.h"
#include "AliEmcalFastorMaskContainer.h"

ClassImp(PWG::EMCAL::AliEmcalFastorMaskContainer)
ClassImp(PWG::EMCAL::AliEmcalFastorMaskContainer::MaskedFastor)
ClassImp(PWG::EMCAL::TestAliEmcalFastorMaskContainer)

using namespace PWG::EMCAL;

AliEmcalFastorMaskContainer::AliEmcalFastorMaskContainer():
  TNamed(),
  fMaskedFastors(),
  fRunType(RunType_t::kRuntypeUndefined)
{
  fMaskedFastors.SetOwner(true);
}

AliEmcalFastorMaskContainer::AliEmcalFastorMaskContainer(const char *name):
  TNamed(name, "Offline masked FastORs"),
  fMaskedFastors(),
  fRunType(RunType_t::kRuntypeUndefined)
{
  fMaskedFastors.SetOwner(true);
}

void AliEmcalFastorMaskContainer::Browse(TBrowser *b) {
  if(fRunType == RunType_t::kRuntypeUndefined) {
    Error("AliEmcalFastorMaskContainer::Browse", "Run type not defined, cannot draw bad channel mask");
    return;
  }
  std::unique_ptr<AliEMCALTriggerMapping> triggermapping;
  switch(fRunType) {
    case RunType_t::kRun1: triggermapping = std::unique_ptr<AliEMCALTriggerMapping>(new AliEMCALTriggerMappingV1); break;
    case RunType_t::kRun2: triggermapping = std::unique_ptr<AliEMCALTriggerMapping>(new AliEMCALTriggerMappingV2); break;
    case RunType_t::kRuntypeUndefined:
    default:
      Error("AliEmcalFastorMaskContainer::Browse", "Trigger mapping unknown");
  };
  if(!triggermapping) {
    Error("AliEmcalFastorMaskContainer::Browse", "Trigger mapping not set, cannot draw bad channel mask");
    return;
  }

  int ncol = 48, nrow = fRunType == RunType_t::kRun1 ? 64 : 104;
  auto plot = new TH2C("TriggerMask", "Offline trigger mask (bad and dead)", ncol, -0.5, ncol - 0.5, nrow, -0.5, nrow - 0.5);
  plot->SetStats(false);
  plot->SetXTitle("col (#eta)");
  plot->SetYTitle("row (#phi)");
  int col, row;
  for(auto dead : GetMaskDead()) {
    triggermapping->GetPositionInEMCALFromAbsFastORIndex(dead, col, row);
    int binx = plot->GetXaxis()->FindBin(col),
        biny = plot->GetYaxis()->FindBin(row);
    plot->SetBinContent(binx, biny, 1);
  }
  for(auto bad : GetMaskBad()) {
    triggermapping->GetPositionInEMCALFromAbsFastORIndex(bad, col, row);
    int binx = plot->GetXaxis()->FindBin(col),
        biny = plot->GetYaxis()->FindBin(row);
    plot->SetBinContent(binx, biny, 2);
  }
  b->Add(plot);
  plot->Draw("col");
}

std::vector<int> AliEmcalFastorMaskContainer::GetMaskAll() const {
  std::vector<int> result;
  for(auto fastor : fMaskedFastors) {
    auto fastorInfo = dynamic_cast<MaskedFastor *>(fastor);
    if(!fastorInfo) continue;
    if(fastorInfo->IsDead() || fastorInfo->IsBad()) result.emplace_back(fastorInfo->GetAbsID());
  }
  std::sort(result.begin(), result.end(), std::less<int>());
  return result;
}

std::vector<int> AliEmcalFastorMaskContainer::GetMaskDead() const {
  std::vector<int> result;
  for(auto fastor : fMaskedFastors) {
    auto fastorInfo = dynamic_cast<MaskedFastor *>(fastor);
    if(!fastorInfo) continue;
    if(fastorInfo->IsDead()) result.emplace_back(fastorInfo->GetAbsID());
  }
  std::sort(result.begin(), result.end(), std::less<int>());
  return result;
}

std::vector<int> AliEmcalFastorMaskContainer::GetMaskBad() const {
  std::vector<int> result;
  for(auto fastor : fMaskedFastors) {
    auto fastorInfo = dynamic_cast<MaskedFastor *>(fastor);
    if(!fastorInfo) continue;
    if(fastorInfo->IsBad()) result.emplace_back(fastorInfo->GetAbsID());
  }
  std::sort(result.begin(), result.end(), std::less<int>());
  return result;
}

AliEmcalFastorMaskContainer::MaskStatus_t AliEmcalFastorMaskContainer::GetFastorMaskStatus(int absID) const {
  auto fastorInfo = FindInContainer(absID);
  if(fastorInfo) return fastorInfo->GetMaskStatus();
  return MaskStatus_t::kStatusUndefined;
}

void AliEmcalFastorMaskContainer::AddDeadFastors(std::vector<int> absIDs) {
  for(auto absID : absIDs) AddDeadFastor(absID);
}

void AliEmcalFastorMaskContainer::AddBadFastors(std::vector<int> absIDs) {
  for(auto absID : absIDs) AddBadFastor(absID);
}

void AliEmcalFastorMaskContainer::SetFastorStatus(int absID, MaskStatus_t status) {
  auto fastorInfo = FindInContainer(absID);
  if(fastorInfo) {
    if(status == MaskStatus_t::kStatusUndefined) fMaskedFastors.Remove(fastorInfo);
    else fastorInfo->SetMaskStatus(status);
  } 
  else {
    if(status != MaskStatus_t::kStatusUndefined) fMaskedFastors.Add(new MaskedFastor(absID, status));
  }
}

AliEmcalFastorMaskContainer::MaskedFastor *AliEmcalFastorMaskContainer::FindInContainer(int absID) const {
  MaskedFastor test(absID, MaskStatus_t::kStatusUndefined);
  auto found = fMaskedFastors.FindObject(&test);
  if(found) {
    return dynamic_cast<MaskedFastor *>(found);
  }
  return nullptr;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/// Content of TestAliEmcalFastorMaskContainer
////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool TestAliEmcalFastorMaskContainer::RunAllTests() const {
  bool status = true;
  std::cout << "Running test: In-memory" <<std::endl;
  if(!TestInMemory()) {
    std::cout << "Failed test: In-memory" << std::endl;
    status = false;
  }
  std::cout << "Running test: I/O" <<std::endl;
  if(!TestIO()) {
    std::cout << "Failed test: I/O" << std::endl;
    status = false;
  } 
  std::cout << "Running test: Change status" <<std::endl;
  if(!TestChangeStatus()) {
    std::cout << "Failed test: Change status" << std::endl;
    status = false;
  } 
  std::cout << "Running test: Clean" <<std::endl;
  if(!TestClean()) {
    std::cout << "Failed test: Clean" << std::endl;
    status = false;
  }
  return status;
}

bool TestAliEmcalFastorMaskContainer::TestInMemory() const {
  std::vector<int> badchannels = GetTestBadChannels(),
                   deadchannels = GetTestDeadChannels(),
                   allchannels;
  std::merge(deadchannels.begin(), deadchannels.end(), badchannels.begin(), badchannels.end(), std::back_inserter(allchannels));
  AliEmcalFastorMaskContainer testcontainer("testcontainer");
  testcontainer.AddBadFastors(badchannels);
  testcontainer.AddDeadFastors(deadchannels);

  auto readbadchannels = testcontainer.GetMaskBad(),
       readdeadchannels = testcontainer.GetMaskDead(),
       readallchannels = testcontainer.GetMaskAll();
  bool status=true;
  if(!CheckEqual(badchannels, readbadchannels)) {
    std::cout << "Test In-memory: Fail equalness bad" << std::endl;
    printVector(badchannels, "Bad channels, true");
    printVector(readbadchannels, "Bad channels, test");
    status = false;
  }
  if(!CheckEqual(deadchannels, readdeadchannels)) {
    std::cout << "Test In-memory: Fail equalness dead" << std::endl;
    printVector(deadchannels, "Dead channels, true");
    printVector(readdeadchannels, "Dead channels, test");
    status = false;
  }
  if(!CheckEqual(allchannels, readallchannels)){
    std::cout << "Test In-memory: Fail equalness all" << std::endl;
    printVector(allchannels, "All channels, true");
    printVector(readallchannels, "All channels, test");
    status = false;
  }
  return status;
}

bool TestAliEmcalFastorMaskContainer::TestIO() const {
  std::vector<int> badchannels = GetTestBadChannels(),
                   deadchannels = GetTestDeadChannels(),
                   allchannels;
  std::merge(deadchannels.begin(), deadchannels.end(), badchannels.begin(), badchannels.end(), std::back_inserter(allchannels));

  std::string filename(Form("%s/TestFastorMaskContainer.root", gSystem->WorkingDirectory())), containername("testcontainer");
  {
    // First step: prepare file and write
    AliEmcalFastorMaskContainer testcontainer(containername.data());
    testcontainer.AddBadFastors(badchannels);
    testcontainer.AddDeadFastors(deadchannels);
    std::unique_ptr<TFile> testwriter(TFile::Open(filename.data(), "RECREATE"));
    testwriter->cd();
    testcontainer.Write(testcontainer.GetName(), TObject::kSingleKey);
  }

  std::vector<int> readbadchannels, readdeadchannels, readallchannels;
  { 
    std::unique_ptr<TFile> testreader(TFile::Open(filename.data(), "READ"));
    if(!testreader || testreader->IsZombie()){
      std::cout << "File entry not found" << std::endl;
      return false;
    } 
    auto testcontainer = dynamic_cast<AliEmcalFastorMaskContainer *>(testreader->Get(containername.data()));
    if(!testcontainer) {
      std::cout << "No container " << testcontainer << " found in file" << std::endl;
      return false;
    }
    readbadchannels = testcontainer->GetMaskBad();
    readdeadchannels = testcontainer->GetMaskDead();
    readallchannels = testcontainer->GetMaskAll();
  }
  // cleanup after test
  if(!gSystem->AccessPathName(filename.data()), kFileExists) gSystem->Exec(Form("rm %s", filename.data()));
  bool status=true;
  if(!CheckEqual(badchannels, readbadchannels)) {
    std::cout << "Test I/O: Fail equalness bad" << std::endl;
    printVector(badchannels, "Bad channels, true");
    printVector(readbadchannels, "Bad channels, test");
    status = false;
  }
  if(!CheckEqual(deadchannels, readdeadchannels)) {
    std::cout << "Test I/O: Fail equalness dead" << std::endl;
    printVector(deadchannels, "Dead channels, true");
    printVector(readdeadchannels, "Dead channels, test");
    status = false;
  }
  if(!CheckEqual(allchannels, readallchannels)){
    std::cout << "Test I/O: Fail equalness all" << std::endl;
    printVector(allchannels, "All channels, true");
    printVector(readallchannels, "All channels, test");
    status = false;
  }
  return status;
}

bool TestAliEmcalFastorMaskContainer::TestChangeStatus() const {
  std::vector<int> badchannels = GetTestBadChannels(),
                   deadchannels = GetTestDeadChannels(),
                   swapToDead = GetRandomSubsample(badchannels, 5),
                   swapToBad = GetRandomSubsample(deadchannels, 5),
                   expectDead,
                   expectBad;
  std::set_difference(badchannels.begin(), badchannels.end(), swapToDead.begin(), swapToDead.end(), std::back_inserter(expectBad));
  for(auto en : swapToBad) expectBad.emplace_back(en);
  std::set_difference(deadchannels.begin(), deadchannels.end(), swapToBad.begin(), swapToBad.end(), std::back_inserter(expectDead));
  for(auto en : swapToDead) expectDead.emplace_back(en);

  AliEmcalFastorMaskContainer testcontainer("testcontainer");
  testcontainer.AddDeadFastors(deadchannels);
  testcontainer.AddBadFastors(badchannels);
  for(auto en : swapToDead) testcontainer.SetFastorStatus(en, AliEmcalFastorMaskContainer::MaskStatus_t::kDeadFastOR);
  for(auto en : swapToBad) testcontainer.SetFastorStatus(en, AliEmcalFastorMaskContainer::MaskStatus_t::kBadFastOR);

  auto readbadchannels = testcontainer.GetMaskBad(),
       readdeadchannels = testcontainer.GetMaskDead();
  bool status=true;
  if(!CheckEqual(expectBad, readbadchannels)) {
    std::cout << "Test I/O: Fail equalness bad" << std::endl;
    printVector(expectBad, "Dead channels, true");
    printVector(readbadchannels, "Dead channels, test");
    status = false;
  }
  if(!CheckEqual(expectDead, readdeadchannels)) {
    std::cout << "Test I/O: Fail equalness dead" << std::endl;
    printVector(expectDead, "Dead channels, true");
    printVector(readdeadchannels, "Dead channels, test");
    status = false;
  }       
  return status;
}

bool TestAliEmcalFastorMaskContainer::TestClean() const {
  std::vector<int> badchannels = GetTestBadChannels(),
                   deadchannels = GetTestDeadChannels();
  AliEmcalFastorMaskContainer testcontainer("testcontainer");
  testcontainer.AddBadFastors(badchannels);
  testcontainer.AddDeadFastors(deadchannels);
  testcontainer.Clear();
  auto allread = testcontainer.GetMaskAll();
  if(allread.size() != 0) {
    std::cout << "Test Clean: size all not 0 (" << allread.size() << ")" << std::endl;
  }
  return true;
}

std::vector<int> TestAliEmcalFastorMaskContainer::GetTestBadChannels() const {
  return {87, 119, 165, 204, 215, 238, 276, 279, 303, 421, 516, 534, 615, 713, 769, 
          776, 825, 910, 935, 971, 974, 975, 1124, 1156, 1159, 1215, 1234, 1267, 1329, 
          1340, 1414, 1499, 1533, 1534, 1538, 1564, 1612, 1735, 1737, 1745, 1747, 1772, 
          1803, 1812, 1828, 1858, 1876, 1937, 1969, 1999, 2053, 2200, 2235, 2239, 2374, 
          2414, 2441, 2538, 2708, 2730, 2765, 2767, 2809, 2872, 2877, 2899, 2937, 2969, 
          3022, 3079, 3081, 3247, 3249, 3282, 3326, 3411, 3691, 3709, 3795, 3853, 3873, 
          3986, 4108, 4135, 4217, 4243, 4251, 4270, 4345, 4411, 4646, 4655, 4676, 4780, 
          4880, 4888, 4898, 4928, 4935, 4969};
}

std::vector<int> TestAliEmcalFastorMaskContainer::GetTestDeadChannels() const {
  return {35, 146, 211, 232, 535, 607, 626, 640, 660, 667, 704, 775, 788, 948, 1041, 1055, 
          1109, 1313, 1331, 1570, 1723, 1878, 1907, 1942, 2044, 2265, 2308, 2328, 2369, 2381, 
          2544, 2612, 2735, 2858, 3014, 3020, 3308, 3402, 3702, 3904, 3916, 4407, 4413, 4485, 
          4626, 4719, 4745, 4851, 4884, 4979};
}

bool TestAliEmcalFastorMaskContainer::CheckEqual(const std::vector<int> &list1, const std::vector<int> &list2) const {
  if(list1.size() != list2.size()) return false;
  return std::is_permutation(list1.begin(), list1.end(), list2.begin());
}

std::vector<int> TestAliEmcalFastorMaskContainer::GetRandomSubsample(std::vector<int> inputlist, int length) const {
  std::vector<int> result(length);
  std::random_device randomhandler;
  std::mt19937 reshuffler(randomhandler());
  std::shuffle(inputlist.begin(), inputlist.end(), reshuffler);
  std::copy(inputlist.begin(), inputlist.begin()+5, result.begin());
  std::sort(result.begin(), result.end(), std::less<int>());
  return result;
}

void TestAliEmcalFastorMaskContainer::printVector(std::vector<int> data, const char *header) const {
  std::sort(data.begin(), data.end(), std::less<int>());
  std::cout << header << " (size " << data.size() << "):" << std::endl;
  std::cout << "__________________________________" << std::endl;
  std::copy(data.begin(), data.end(), std::ostream_iterator<int>(std::cout, " "));
  std::cout << "\n";
  std::cout << "\n";
}