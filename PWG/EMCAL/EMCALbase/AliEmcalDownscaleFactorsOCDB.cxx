#include <iostream>

#include <TObjArray.h>

/**************************************************************************
 * Copyright(c) 1998-2016, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliLog.h"
#include "AliTriggerClass.h"
#include "AliTriggerConfiguration.h"

#include "AliEmcalDownscaleFactorsOCDB.h"

/// \cond CLASSIMP
ClassImp(AliEmcalDownscaleFactorsOCDB)
/// \endcond

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
