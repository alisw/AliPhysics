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

#include <TObjArray.h>

#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliLog.h"
#include "AliTriggerClass.h"
#include "AliTriggerConfiguration.h"

#include "AliEmcalDownscaleFactorsOCDB.h"

/// \cond CLASSIMP
ClassImp(PWG::EMCAL::AliEmcalDownscaleFactorsOCDB)
/// \endcond

using namespace PWG::EMCAL;

AliEmcalDownscaleFactorsOCDB *AliEmcalDownscaleFactorsOCDB::fgDownscaleFactors = nullptr;

AliEmcalDownscaleFactorsOCDB::AliEmcalDownscaleFactorsOCDB() :
  fCurrentRun(0),
  fDownscaleFactors()
{
}

AliEmcalDownscaleFactorsOCDB *AliEmcalDownscaleFactorsOCDB::Instance(){
  if(!fgDownscaleFactors) {
    fgDownscaleFactors = new AliEmcalDownscaleFactorsOCDB;
  }
  return fgDownscaleFactors;
}

void AliEmcalDownscaleFactorsOCDB::SetRun(Int_t runnumber){
  if(runnumber == fCurrentRun) return;

  fCurrentRun = runnumber;
  fDownscaleFactors.clear();
  AliInfoStream() << "Loading downscale factors for run " << fCurrentRun << std::endl;

  AliCDBManager *mgr = AliCDBManager::Instance();
  AliCDBEntry *trgcdb = mgr->Get("GRP/CTP/Config");
  AliTriggerConfiguration *trgconf = static_cast<AliTriggerConfiguration *>(trgcdb->GetObject());
  for(auto e : trgconf->GetClasses()){
    AliTriggerClass *trgcls = static_cast<AliTriggerClass *>(e);
    Double_t downscalefactor;
    trgcls->GetDownscaleFactor(downscalefactor);
    fDownscaleFactors.insert(std::pair<TString, Double_t>(trgcls->GetName(), downscalefactor));
  }
}

double AliEmcalDownscaleFactorsOCDB::GetDownscaleFactorForTriggerClass(const TString &trigger) const {
  double downscale = 1.;
  std::map<TString, Double_t>::const_iterator found = fDownscaleFactors.find(trigger);
  if(found != fDownscaleFactors.end()) downscale = found->second;
  return downscale;
}

std::vector<TString> AliEmcalDownscaleFactorsOCDB::GetTriggerClasses() const {
  std::vector<TString> triggerclasses;
  for(const auto &en : fDownscaleFactors) triggerclasses.push_back(en.first);
  return triggerclasses;
}
