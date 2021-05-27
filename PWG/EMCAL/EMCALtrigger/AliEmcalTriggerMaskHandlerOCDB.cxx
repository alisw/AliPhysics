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
#include <bitset>
#include <functional>
#include <iostream>

#include <TH2.h>

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliEMCALGeometry.h"
#include "AliEMCALTriggerDCSConfig.h"
#include "AliEMCALTriggerTRUDCSConfig.h"
#include "AliEMCALTriggerSTUDCSConfig.h"
#include "AliEmcalTriggerMaskHandlerOCDB.h"

ClassImp(PWG::EMCAL::AliEmcalTriggerMaskHandlerOCDB)

using namespace PWG::EMCAL;

AliEmcalTriggerMaskHandlerOCDB *AliEmcalTriggerMaskHandlerOCDB::fgInstance = nullptr;

AliEmcalTriggerMaskHandlerOCDB *AliEmcalTriggerMaskHandlerOCDB::Instance() {
  if(!fgInstance) fgInstance = new AliEmcalTriggerMaskHandlerOCDB;
  return fgInstance;
}

std::vector<int> AliEmcalTriggerMaskHandlerOCDB::GetMaskedFastorIndicesL0(int runnumber) {
  std::vector<int> maskedfastors;
  UpdateCache(runnumber);
  auto geo = GetGeometry(runnumber);

  std::function<int (int, int)> ChannelMaskHandler = geo->GetTriggerMappingVersion() == 2 ? GetChannelForMaskRun2 : GetChannelForMaskRun1;
  int itru = 0;
  for(auto truconfigobj : *fCurrentConfig->GetTRUArr()) {
    auto truconfig = dynamic_cast<AliEMCALTriggerTRUDCSConfig *>(truconfigobj);
    int localtru = itru % 32, detector = itru >= 32 ? 1 : 0,
        globaltru = geo->GetTriggerMapping()->GetTRUIndexFromSTUIndex(localtru, detector);
    for(int ipos = 0; ipos < 6; ipos++) {
      auto regmask = truconfig->GetMaskReg(ipos);
      std::bitset<16> bitsregmask(regmask);
      for(int ibit = 0; ibit < 16; ibit++) {
        if(bitsregmask.test(ibit)) {
          auto channel = ChannelMaskHandler(ipos, ibit);
          int absfastor;
          geo->GetTriggerMapping()->GetAbsFastORIndexFromTRU(globaltru, channel, absfastor);
          maskedfastors.push_back(absfastor);
        }
      }
    }
    itru++;
  }
  std::sort(maskedfastors.begin(), maskedfastors.end(), std::less<int>());
  return maskedfastors;
}

std::vector<AliEmcalTriggerMaskHandlerOCDB::FastORPosition> AliEmcalTriggerMaskHandlerOCDB::GetMaskedFastorPositionsL0(int runnumber) {
  std::vector<FastORPosition> maskedfastors;
  auto geo = GetGeometry(runnumber);
  for(auto absfastor : GetMaskedFastorIndicesL0(runnumber)) {
    int row, col;
    geo->GetTriggerMapping()->GetPositionInEMCALFromAbsFastORIndex(absfastor, col, row); 
    maskedfastors.push_back({col, row});
  }
  std::sort(maskedfastors.begin(), maskedfastors.end(), std::less<FastORPosition>());
  return maskedfastors;
}

std::vector<int> AliEmcalTriggerMaskHandlerOCDB::GetMaskedFastorIndicesL1(int runnumber) {
  std::vector<int> maskedfastors;
  // Find masked FastORs at L0
  auto maskedL0 = GetMaskedFastorIndicesL0(runnumber);
  std::copy(maskedL0.begin(), maskedL0.end(), maskedfastors.begin());

  // Look for TRUs which are masked in the L1 region
  // For each STU which is masked at L1 mask also the FastORs 
  // if not yet masked at L0
  AliEMCALGeometry *geo = nullptr;
  for(auto truID : GetGlobalMaskedTRUIndices()) {
    if(!geo) geo = GetGeometry(runnumber);
    for(int ichannel = 0; ichannel < 96; ichannel++) {
      int absfastor;
      geo->GetTriggerMapping()->GetAbsFastORIndexFromTRU(truID, ichannel, absfastor);
      if(std::find(maskedfastors.begin(), maskedfastors.end(), absfastor)  != maskedfastors.end()) maskedfastors.push_back(absfastor);
    }
  }

  std::sort(maskedfastors.begin(), maskedfastors.end(), std::less<int>());
  return maskedfastors;
}

std::vector<AliEmcalTriggerMaskHandlerOCDB::FastORPosition> AliEmcalTriggerMaskHandlerOCDB::GetMaskedFastorPositionsL1(int runnumber) {
  std::vector<FastORPosition> maskedfastors;
  auto geo = GetGeometry(runnumber);
  for(auto absfastor : GetMaskedFastorIndicesL1(runnumber)) {
    int row, col;
    geo->GetTriggerMapping()->GetPositionInEMCALFromAbsFastORIndex(absfastor, col, row); 
    maskedfastors.push_back({col, row});
  }
  std::sort(maskedfastors.begin(), maskedfastors.end(), std::less<FastORPosition>());

  return maskedfastors;
}

std::vector<int> AliEmcalTriggerMaskHandlerOCDB::GetGlobalMaskedTRUIndices(int runnumber) {
  std::vector<int> truindices;
  UpdateCache(runnumber);
  auto geo = GetGeometry(runnumber);

  auto emcalstu = fCurrentConfig->GetSTUDCSConfig(false),
       dcalstu = fCurrentConfig->GetSTUDCSConfig(false);

  if(emcalstu) {
    std::bitset<32> mask(emcalstu->GetRegion());
    for(int itru = 0; itru < 32; itru++) {
      if(!mask.test(itru)) {
        // TRU excluded from region, add to masked TRUs
        truindices.push_back(geo->GetTRUIndexFromSTUIndex(itru, 0));
      } 
    }
  }

  if(dcalstu) {
    std::bitset<32> mask(dcalstu->GetRegion());
    for(int itru = 0; itru < 14; itru++) {
      if(!mask.test(itru)) {
        // TRU excluded from region, add to masked TRUs
        truindices.push_back(geo->GetTRUIndexFromSTUIndex(itru, 1));
      } 
    }
  }

  std::sort(truindices.begin(), truindices.end(), std::less<int>());
  return truindices;
}

TH2 *AliEmcalTriggerMaskHandlerOCDB::MonitorMaskedFastORsL0(int runnumber) {
  std::string runstring = runnumber < 0 ? "Current" : Form("%d", runnumber);
  std::string histname = Form("maskedFastorsL0%s", runstring.data()), 
              histtitle = Form("Masked fastors at L0 for run %s", runstring.data());
  auto outputhist = PrepareHistogram(histname.data(), histtitle.data());
  auto maskedfastors = GetMaskedFastorPositionsL0(runnumber);
  FillMaskedFastors(outputhist, maskedfastors);
  return outputhist;
}

TH2 *AliEmcalTriggerMaskHandlerOCDB::MonitorMaskedFastORsL1(int runnumber) {
  std::string runstring = runnumber < 0 ? "Current" : Form("%d", runnumber);
  std::string histname = Form("maskedFastorsL1%s", runstring.data()), 
              histtitle = Form("Masked fastors at L1 for run %s", runstring.data());
  auto outputhist = PrepareHistogram(histname.data(), histtitle.data());
  auto maskedfastors = GetMaskedFastorPositionsL1(runnumber);
  FillMaskedFastors(outputhist, maskedfastors);
  return outputhist;
}

int AliEmcalTriggerMaskHandlerOCDB::GetChannelForMaskRun1(int mask, int bitnumber) {
  return mask * 16 + bitnumber;
}

int AliEmcalTriggerMaskHandlerOCDB::GetChannelForMaskRun2(int mask, int bitnumber) {
  const int kChannelMap[6][16] = {{ 8, 9,10,11,20,21,22,23,32,33,34,35,44,45,46,47},   // Channels in mask0
                                  {56,57,58,59,68,69,70,71,80,81,82,83,92,93,94,95},   // Channels in mask1
                                  { 4, 5, 6, 7,16,17,18,19,28,29,30,31,40,41,42,43},   // Channels in mask2
                                  {52,53,54,55,64,65,66,67,76,77,78,79,88,89,90,91},   // Channels in mask3
                                  { 0, 1, 2, 3,12,13,14,15,24,25,26,27,36,37,38,39},   // Channels in mask4
                                  {48,49,50,51,60,61,62,63,72,73,74,75,84,85,86,87}};  // Channels in mask5
  return kChannelMap[mask][bitnumber];
}

AliEMCALGeometry *AliEmcalTriggerMaskHandlerOCDB::GetGeometry(int runnumber) {
  if(runnumber >= 0) return AliEMCALGeometry::GetInstanceFromRunNumber(runnumber);
  auto geo = AliEMCALGeometry::GetInstance();
  if(!geo) throw GeometryNotSetException();
  return geo;
}

void AliEmcalTriggerMaskHandlerOCDB::ClearCache() {
  if(fCurrentConfig) delete fCurrentConfig;
  fCurrentConfig = nullptr;
  fCurrentRunnumber = -1;
}

void AliEmcalTriggerMaskHandlerOCDB::UpdateCache(int runnumber) {
  if((runnumber == fCurrentRunnumber) && fCurrentConfig) return;
  ClearCache();
  auto cdb = InitCDB(runnumber);
  auto trgobject = cdb->Get("EMCAL/CalibTrigger")->GetObject();
  if(!trgobject) throw OCDBNotInitializedException();
  fCurrentConfig = dynamic_cast<AliEMCALTriggerDCSConfig *>(trgobject);
  if(!fCurrentConfig) throw OCDBNotInitializedException();
  if(runnumber != fCurrentRunnumber) fCurrentRunnumber = runnumber;
}


AliCDBManager *AliEmcalTriggerMaskHandlerOCDB::InitCDB(int runnumber) {
  AliCDBManager *cdb = AliCDBManager::Instance();
  if(!cdb) throw OCDBNotInitializedException();
  if(!cdb->IsDefaultStorageSet()) throw OCDBNotInitializedException();
  if(runnumber >= 0) {
    // Run number handled by handler
    auto currentrun = cdb->GetRun();
    if(currentrun != runnumber) {
      // Changing run numbre
      cdb->SetRun(runnumber);
    }
  } else {
    // Run number not specified - run number handled outside
    if(cdb->GetRun() <= 0) throw OCDBNotInitializedException();
  }
  return cdb;
}

TH2 *AliEmcalTriggerMaskHandlerOCDB::PrepareHistogram(const char *name, const char *title) {
  TH2 *maskhist = new TH2I(name, title, 48, -0.5, 47.5, 104, -0.5, 103.5);
  maskhist->SetDirectory(nullptr);
  maskhist->GetXaxis()->SetTitle("#eta (col)");
  maskhist->GetYaxis()->SetTitle("#phi (row)");
  return maskhist;
}

void AliEmcalTriggerMaskHandlerOCDB::FillMaskedFastors(TH2 *outputhist, const std::vector<FastORPosition> &fastors) const {
  for(auto fastor : fastors) {
    outputhist->Fill(fastor.column, fastor.row);
  }
}