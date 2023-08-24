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
#include <algorithm>
#include <cfloat>
#include <iostream>
#include <vector>

#include <TClonesArray.h>
#include <TH2.h>
#include <TLorentzVector.h>
#include <TRandom.h>

#include "AliEMCALTriggerConstants.h"
#include "AliEMCALTriggerPatchInfo.h"
#include "AliEmcalTriggerOfflineSelection.h"
#include "AliLog.h"
#include "AliVEvent.h"

ClassImp(PWGJE::EMCALJetTasks::AliEmcalTriggerOfflineSelection)

using namespace PWGJE::EMCALJetTasks;

const TString AliEmcalTriggerOfflineSelection::fgkTriggerNames[AliEmcalTriggerOfflineSelection::kTrgn] = {
    "EMC7", "EG1", "EG2", "EJ1", "EJ2", "DMC7", "DG1", "DG2", "DJ1", "DJ2"
};

AliEmcalTriggerOfflineSelection::AliEmcalTriggerOfflineSelection():
    TObject(),
    fEnergyDefinition(kFEEEnergy),
    fNameClusterContainer(),
    fResolution(0),
    fUseSmearedEnergy(false)
{
  for(int itrg = 0; itrg < kTrgn; itrg++) fOfflineEnergyThreshold[itrg] = 100000.;  // unimplemented triggers get very high threshold assinged, so that the result is automatically false
  memset(fAcceptanceMaps, 0, sizeof(TH2 *) * kTrgn);
}

AliEmcalTriggerOfflineSelection::~AliEmcalTriggerOfflineSelection(){
  for(int itrg = 0; itrg < kTrgn; itrg++){
    if(fAcceptanceMaps[itrg]) delete fAcceptanceMaps[itrg];
  }
}

Bool_t AliEmcalTriggerOfflineSelection::IsOfflineSelected(EmcalTriggerClass trgcls, const AliVEvent * const data) const {
  if(fOfflineEnergyThreshold[trgcls] < 0) return true;
  AliDebugStream(1) << "Applying offline threshold " << fOfflineEnergyThreshold[trgcls] << " for trigger class " << GetTriggerName(trgcls) << std::endl;
  if(UseClusters()){
    return ApplyClusterTrigger(trgcls, data);
  }
  return ApplyPatchTrigger(trgcls, static_cast<TClonesArray *>(data->FindListObject("EmcalTriggers")));
}

bool AliEmcalTriggerOfflineSelection::ApplyPatchTrigger(EmcalTriggerClass trgcls, const TClonesArray * const triggerpatches) const {
  AliDebugStream(1) << "Using patch trigger with energy definition " << fEnergyDefinition << std::endl;
  bool isSingleShower = IsSingleShower(trgcls), isDCAL = IsDCAL(trgcls);
  int nfound = 0;
  AliEMCALTriggerPatchInfo *patch = NULL;
  std::vector<double> patchefficiencies;
  for(auto patchIter : *triggerpatches){
    patch = static_cast<AliEMCALTriggerPatchInfo *>(patchIter);
    if(!patch->IsOfflineSimple()) continue;
    if((isDCAL && !patch->IsDCalPHOS()) || (!isDCAL && patch->IsDCalPHOS())) continue;      // reject patches in opposite detector
    if(isSingleShower){
     if(!patch->IsGammaLowSimple()) continue;
    } else {
      if(!patch->IsJetLowSimple()) continue;
    }
    AliDebugStream(1) << "Patch energy: " << patch->GetPatchE() << ", smeared " << patch->GetSmearedEnergy() << std::endl;
    double energy(0);
    // No switch as only cases for patches are handled in this class
    if(fEnergyDefinition == kFEEEnergy){
      if(fUseSmearedEnergy){
        AliDebugStream(1) << "Using smeared energy" << std::endl;
        energy = patch->GetSmearedEnergy();
      } else {
        AliDebugStream (1) << "Using default energy" << std::endl;
        energy = patch->GetPatchE();
      }
    }
    else if(fEnergyDefinition == kFEETransverseEnergy){
      if(fUseSmearedEnergy){
        TLorentzVector cm = patch->GetLorentzVectorCM();
        cm.SetE(patch->GetSmearedEnergy());
        energy = cm.Et();
      } else {
        energy = patch->GetPatchET();
      }
    }
    else if(fEnergyDefinition == kFEEADC) energy = patch->GetADCOfflineAmp();
    else if(fEnergyDefinition == kFEETransverseADC) energy = patch->GetPatchET() / EMCALTrigger::kEMCL1ADCtoGeV;
    double threshold(fOfflineEnergyThreshold[trgcls]);
    if(TMath::Abs(fResolution) > DBL_EPSILON){
      // Simulation effect of the finite online energy resolution:
      // The threshold is smeared using a gaussinan using the
      // energy resolution as parameter for the smearing. A patch is
      // selected in case the energy is above the smeared threshold.
      // The approach is equivalent to smearing the patch energy itself.
      AliDebugStream(1) << "Smearing energy with resolution " << fResolution << std::endl;
      threshold = gRandom->Gaus(threshold, fResolution);
      AliDebugStream(1) << "Thresholds: original(" << fOfflineEnergyThreshold[trgcls] << "), smeared(" << threshold << ")" << std::endl;
    }
    if(energy > threshold){
      AliDebugStream(2) << GetTriggerName(trgcls) << " patch above threshold (" << patch->GetPatchE() << " | " << fOfflineEnergyThreshold[trgcls] << ")" <<  std::endl;
      if(fAcceptanceMaps[trgcls]){
        // Handle azimuthal inefficiencies of the trigger observed online:
        // For each patch provide an efficiency. Once all patches in the event
        // is determined, only the patch with the maximum efficiency is chosen.
        // The event is selected in this case if the sample value is below the
        // efficiency for the chosen patch.
        double peff = fAcceptanceMaps[trgcls]->GetBinContent(patch->GetColStart(), patch->GetRowStart());
        patchefficiencies.push_back(peff);
        AliDebugStream(2) << "Spatial Efficiency " << peff
            << " for trigger patch at position (" << patch->GetColStart()
            << "," << patch->GetRowStart() << ")" << std::endl;
      } else{
        patchefficiencies.push_back(1.);
      }
    }
  }
  if(patchefficiencies.size()){
    std::sort(patchefficiencies.begin(), patchefficiencies.end(), std::greater<double>());
    double sample = gRandom->Uniform(0., 1.);
    if(sample < patchefficiencies[0]){
      AliDebugStream(1) << "Event selected for trigger class " << GetTriggerName(trgcls) << ", " << nfound << " good patch(es) found" << std::endl;
      return true;
    }
  }
  return false;
}

bool AliEmcalTriggerOfflineSelection::ApplyClusterTrigger(EmcalTriggerClass trgcls, const AliVEvent * const ev) const {
  AliDebugStream(1) << "Applying cluster trigger with energy definition " << fEnergyDefinition << std::endl;
  int ntrigger = 0;
  TClonesArray *clusters = static_cast<TClonesArray *>(ev->FindListObject(fNameClusterContainer));
  double vertex[3]; ev->GetPrimaryVertex()->GetXYZ(vertex);
  for(auto o : *clusters){
    AliVCluster *c = dynamic_cast<AliVCluster *>(o);
    if(!c->IsEMCAL()) continue;
    // don't apply a time cut in case of simulations
    double energy = c->GetNonLinCorrEnergy();
    if(fEnergyDefinition == kClusterTransverseEnergy) {
      TLorentzVector vec;
      c->GetMomentum(vec, vertex);
      vec.SetE(c->GetNonLinCorrEnergy());
      energy = vec.Et();
    }
    double threshold(fOfflineEnergyThreshold[trgcls]);
    if(TMath::Abs(fResolution) > DBL_EPSILON){
      // Simulation effect of the finite online energy resolution:
      // The threshold is smeared using a gaussinan using the
      // energy resolution as parameter for the smearing. A patch is
      // selected in case the energy is above the smeared threshold.
      // The approach is equivalent to smearing the patch energy itself.
      threshold = gRandom->Gaus(threshold, fResolution);
    }
    if(energy > threshold) ntrigger++;
  }
  return ntrigger > 0;
}

Bool_t AliEmcalTriggerOfflineSelection::IsSingleShower(EmcalTriggerClass cls){
 return ((cls == kTrgEG1) || (cls == kTrgEG2) || (cls == kTrgEL0) || (cls == kTrgDG1) || (cls == kTrgDG2) || (cls == kTrgDL0));

}

Bool_t AliEmcalTriggerOfflineSelection::IsDCAL(EmcalTriggerClass cls){
  return ((cls == kTrgDL0) || (cls == kTrgDG1) || (cls == kTrgDG2) || (cls == kTrgDJ1) || (cls == kTrgDJ2));
}

Bool_t AliEmcalTriggerOfflineSelection::UseClusters() const {
  return fEnergyDefinition == kClusterEnergy || fEnergyDefinition == kClusterTransverseEnergy;
}

Bool_t AliEmcalTriggerOfflineSelection::UsePatches() const {
  return fEnergyDefinition == kFEEEnergy || fEnergyDefinition == kFEETransverseEnergy || fEnergyDefinition == kFEEADC || fEnergyDefinition == kFEETransverseADC;
}
