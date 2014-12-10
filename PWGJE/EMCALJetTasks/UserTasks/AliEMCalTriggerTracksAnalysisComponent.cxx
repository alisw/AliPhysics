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
 * Base class for anaysis components. Inheriting classes have to implement the
 * functions CreateHistos and Process.
 *
 *   Author: Markus Fasel
 */
#include <TAxis.h>

#include "AliEMCalTriggerAnaTriggerDecision.h"
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerTracksAnalysisComponent.h"

ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerTracksAnalysisComponent)

namespace EMCalTriggerPtAnalysis {

//______________________________________________________________________________
AliEMCalTriggerTracksAnalysisComponent::AliEMCalTriggerTracksAnalysisComponent() :
  TNamed(),
  fHistos(NULL),
  fBinning(NULL),
  fKineCuts(NULL),
  fTriggerDecision(NULL)
{
  /*
   * Dummy (I/O) constructor
   */
}

//______________________________________________________________________________
AliEMCalTriggerTracksAnalysisComponent::~AliEMCalTriggerTracksAnalysisComponent() {
  /*
   * Release histogram container
   */
  if(fHistos) delete fHistos;
}

//______________________________________________________________________________
AliEMCalTriggerTracksAnalysisComponent::AliEMCalTriggerTracksAnalysisComponent(const char* name) :
  TNamed(name,""),
  fHistos(NULL),
  fBinning(NULL),
  fKineCuts(NULL),
  fTriggerDecision(NULL)
{
  /*
   * Main constructor, to be called by the user
   *
   * @param name: component name
   */
}

//______________________________________________________________________________
void AliEMCalTriggerTracksAnalysisComponent::CreateHistos() {
  /*
   * Create Container for histograms. Inheriting classes overwrite this method, in which they call
   * this and add the histograms of their choise.
   */
  fHistos = new AliEMCalHistoContainer(Form("Histos%s", GetName()));
  fHistos->ReleaseOwner();
}

//______________________________________________________________________________
TAxis* AliEMCalTriggerTracksAnalysisComponent::DefineAxis(const char* name, const AliEMCalTriggerBinningDimension* binning) {
  /*
   * Create and define axis
   *
   * @param name: Name of the axis
   * @param binning: binning information
   * @return: the new axis
   */
  TAxis *result = new TAxis(binning->GetNumberOfBins(), binning->GetBinLimits());
  result->SetName(name);
  return result;
}

//______________________________________________________________________________
TAxis* AliEMCalTriggerTracksAnalysisComponent::DefineAxis(const char* name, int nbins, double min, double max) {
  /*
   * Create and define axis
   *
   * @param name: Name of the axis
   * @param nbins: number of bins
   * @param min: min. range
   * @param max: max. range
   * @return: the new axis
   */
  TAxis *result = new TAxis(nbins, min, max);
  result->SetName(name);
  return result;
}

//______________________________________________________________________________
void AliEMCalTriggerTracksAnalysisComponent::GetMachingTriggerNames(std::vector<std::string>& triggernames, Bool_t usePatches) {
  /*
   * Get a set of names of trigger strings that is matching with the trigger decision.
   *
   * @param triggernames: output container for selected trigger names
   * @param usePatches: determines whether we use the trigger decision from patches
   */
  triggernames.clear();
  if(!fTriggerDecision) return;
  if(fTriggerDecision->IsMinBias()) triggernames.push_back("MinBias");
  if(fTriggerDecision->IsTriggered(AliEMCalTriggerAnaTriggerDecision::kTAEMCJHigh, usePatches)){
    triggernames.push_back("EMCJHigh");
    if(fTriggerDecision->IsTriggered(AliEMCalTriggerAnaTriggerDecision::kTAEMCGHigh, usePatches))
      triggernames.push_back("EMCHighBoth");
    else
      triggernames.push_back("EMCHighJetOnly");
  }
  if(fTriggerDecision->IsTriggered(AliEMCalTriggerAnaTriggerDecision::kTAEMCJLow, usePatches)){
    triggernames.push_back("EMCJLow");
    if(fTriggerDecision->IsTriggered(AliEMCalTriggerAnaTriggerDecision::kTAEMCGLow, usePatches))
      triggernames.push_back("EMCLowBoth");
    else
      triggernames.push_back("EMCLowJetOnly");
  }
  if(fTriggerDecision->IsTriggered(AliEMCalTriggerAnaTriggerDecision::kTAEMCGHigh, usePatches)){
    triggernames.push_back("EMCGHigh");
    if(!fTriggerDecision->IsTriggered(AliEMCalTriggerAnaTriggerDecision::kTAEMCJHigh, usePatches))
      triggernames.push_back("EMCHighGammaOnly");
  }
  if(fTriggerDecision->IsTriggered(AliEMCalTriggerAnaTriggerDecision::kTAEMCGLow, usePatches)){
    triggernames.push_back("EMCGLow");
    if(!fTriggerDecision->IsTriggered(AliEMCalTriggerAnaTriggerDecision::kTAEMCJLow, usePatches))
      triggernames.push_back("EMCLowGammaOnly");
  }
}

} /* namespace EMCalTriggerPtAnalysis */

