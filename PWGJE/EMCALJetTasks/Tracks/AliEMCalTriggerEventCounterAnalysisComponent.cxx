/**************************************************************************
 * Copyright(c) 1998-2014, ALICE Experiment at CERN, All rights reserved. *
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
/*
 * Analysis component counting events for different trigger classes. Task needs
 * to be grouped with a global event selection
 *
 *   Author: Markus Fasel
 */
#include <map>
#include <string>

#include "AliVEvent.h"
#include "AliVVertex.h"

#include "AliEMCalTriggerAnaTriggerDecision.h"
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerEventData.h"
#include "AliEMCalHistoContainer.h"
#include "AliEMCalTriggerEventCounterAnalysisComponent.h"

ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerEventCounterAnalysisComponent)

namespace EMCalTriggerPtAnalysis {

//______________________________________________________________________________
AliEMCalTriggerEventCounterAnalysisComponent::AliEMCalTriggerEventCounterAnalysisComponent():
  AliEMCalTriggerTracksAnalysisComponent(),
  fUsePatches(kFALSE)
{
  /*
   * Default (I/O) constructor, not to be used
   */
}

//______________________________________________________________________________
AliEMCalTriggerEventCounterAnalysisComponent::AliEMCalTriggerEventCounterAnalysisComponent(const char *name):
  AliEMCalTriggerTracksAnalysisComponent(name),
  fUsePatches(kFALSE)
{
  /*
   * Main constructor
   */
}

//______________________________________________________________________________
void AliEMCalTriggerEventCounterAnalysisComponent::CreateHistos() {
  /*
   * Create event counter histograms
   */
  AliEMCalTriggerTracksAnalysisComponent::CreateHistos();

  // Create trigger definitions
  std::map<std::string, std::string> triggerCombinations;
  const char *triggernames[11] = {"MinBias", "EMCJHigh", "EMCJLow", "EMCGHigh",
      "EMCGLow", "EMCHighBoth", "EMCHighGammaOnly", "EMCHighJetOnly",
      "EMCLowBoth", "EMCLowGammaOnly", "EMCLowJetOnly"};
  // Define names and titles for different triggers in the histogram container
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[0], "min. bias events"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[1], "jet-triggered events (high threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[2], "jet-triggered events (low threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[3], "gamma-triggered events (high threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[4], "gamma-triggered events (low threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[5], "jet and gamma triggered events (high threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[6], "exclusively gamma-triggered events (high threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[7], "exclusively jet-triggered events (high threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[8], "jet and gamma triggered events (low threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[9], "exclusively gamma-triggered events (low threshold)"));
  triggerCombinations.insert(std::pair<std::string,std::string>(triggernames[10], "exclusively-triggered events (low threshold)"));

  AliEMCalTriggerBinningDimension *vertexbinning = fBinning->GetBinning("zvertex");

  for(std::map<std::string,std::string>::iterator it = triggerCombinations.begin(); it != triggerCombinations.end(); ++it){
    const std::string name = it->first, &title = it->second;
    // Create event-based histogram
    fHistos->CreateTH1(Form("hEventHist%s", name.c_str()), Form("Event-based data for %s events; pileup rejection; z_{V} (cm)", title.c_str()), vertexbinning->GetNumberOfBins(), vertexbinning->GetBinLimits());
  }

  // Make correlation histogram for different trigger classes
  const TAxis *triggeraxis[5]; memset(triggeraxis, 0, sizeof(const TAxis *) * 5);
  const char *binlabels[2] = {"OFF", "ON"};
  TAxis mytrgaxis[5];
  for(int itrg = 0; itrg < 5; ++itrg){
    DefineAxis(mytrgaxis[itrg], triggernames[itrg], triggernames[itrg], 2, -0.5, 1.5, binlabels);
    triggeraxis[itrg] = mytrgaxis+itrg;
  }
  fHistos->CreateTHnSparse("hEventTriggers", "Trigger type per event", 5, triggeraxis);
}

//______________________________________________________________________________
void AliEMCalTriggerEventCounterAnalysisComponent::Process(const AliEMCalTriggerEventData* const data) {
  /*
   * Do event counting
   */
  if(!fTriggerDecision) return;

  double vz = data->GetRecEvent()->GetPrimaryVertex()->GetZ();
  double triggerCorrelation[5]; memset(triggerCorrelation, 0, sizeof(double) * 5);

  if(fComponentDebugLevel > 2)
    fTriggerDecision->Print();

  if(fTriggerDecision->IsMinBias()){
    triggerCorrelation[0] = 1.;
    fHistos->FillTH1("hEventHistMinBias", vz);
  }
  if(fTriggerDecision->IsTriggered(AliEMCalTriggerAnaTriggerDecision::kTAEMCJHigh, fUsePatches)){
    triggerCorrelation[2] = 1.;
    fHistos->FillTH1("hEventHistEMCJHigh", vz);
    // Check whether also the gamma high-threshold trigger fired
    if(fTriggerDecision->IsTriggered(AliEMCalTriggerAnaTriggerDecision::kTAEMCGHigh, fUsePatches)){
      fHistos->FillTH1("hEventHistEMCHighBoth", vz);
    } else {
      fHistos->FillTH1("hEventHistEMCHighJetOnly", vz);
    }
  }
  if(fTriggerDecision->IsTriggered(AliEMCalTriggerAnaTriggerDecision::kTAEMCJLow, fUsePatches)){
    triggerCorrelation[1] = 1.;
    fHistos->FillTH1("hEventHistEMCJLow", vz);
    // Check whether also the gamma high-threshold trigger fired
    if(fTriggerDecision->IsTriggered(AliEMCalTriggerAnaTriggerDecision::kTAEMCGLow, fUsePatches)){
      fHistos->FillTH1("hEventHistEMCLowBoth", vz);
    } else {
      fHistos->FillTH1("hEventHistEMCLowJetOnly", vz);
    }
  }
  if(fTriggerDecision->IsTriggered(AliEMCalTriggerAnaTriggerDecision::kTAEMCGHigh, fUsePatches)){
    triggerCorrelation[3] = 1.;
    fHistos->FillTH1("hEventHistEMCGHigh", vz);
    if(!fTriggerDecision->IsTriggered(AliEMCalTriggerAnaTriggerDecision::kTAEMCJHigh, fUsePatches))
      fHistos->FillTH1("hEventHistEMCHighGammaOnly", vz);
  }
  if(fTriggerDecision->IsTriggered(AliEMCalTriggerAnaTriggerDecision::kTAEMCGLow, fUsePatches)){
    triggerCorrelation[4] = 1.;
    fHistos->FillTH1("hEventHistEMCGLow", vz);
    if(!fTriggerDecision->IsTriggered(AliEMCalTriggerAnaTriggerDecision::kTAEMCJLow, fUsePatches))
      fHistos->FillTH1("hEventHistEMCLowGammaOnly", vz);
  }

  fHistos->FillTHnSparse("hEventTriggers", triggerCorrelation);
}

//______________________________________________________________________________
void AliEMCalTriggerEventCounterAnalysisComponent::DefineAxis(TAxis& axis, const char* name,
    const char* title, int nbins, double min, double max,
    const char** labels) const {
  /*
   * Define an axis with number of bins from min to max
   *
   * @param axis: Axis to be defined
   * @param name: Name of the axis
   * @param title: Title of the axis
   * @param nbins: Number of bins
   * @param min: lower limit of the axis
   * @param max: upper limit of the axis
   * @param labels (@optional): array of bin labels
   */
  axis.Set(nbins, min, max);
  axis.SetName(name);
  axis.SetTitle(title);
  if(labels){
    for(int ib = 1; ib <= axis.GetNbins(); ++ib)
      axis.SetBinLabel(ib, labels[ib-1]);
  }
}

} /* namespace EMCalTriggerPtAnalysis */
