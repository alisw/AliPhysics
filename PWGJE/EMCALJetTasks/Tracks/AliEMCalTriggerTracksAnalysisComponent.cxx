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
#include <iostream>

#include <TAxis.h>
#include "AliEMCalTriggerTracksAnalysisComponent.h"
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerAnaTriggerDecision.h"

/// \cond CLASSIMP
ClassImp(EMCalTriggerPtAnalysis::AliEMCalTriggerTracksAnalysisComponent)
/// \endcond

namespace EMCalTriggerPtAnalysis {

/**
 * Dummy (I/O) constructor, not to be used
 */
AliEMCalTriggerTracksAnalysisComponent::AliEMCalTriggerTracksAnalysisComponent() :
  TNamed(),
  fHistos(NULL),
  fBinning(NULL),
  fKineCuts(NULL),
  fTriggerDecision(NULL),
  fWeightHandler(NULL),
  fComponentDebugLevel(0)
{
}

/**
 * Destructor, release histogram container
 */
AliEMCalTriggerTracksAnalysisComponent::~AliEMCalTriggerTracksAnalysisComponent() {
  if(fHistos) delete fHistos;
}

/**
 * Main constructor, to be called by the user. Initializes all fields with NULL.
 * \param name: component name
 */
AliEMCalTriggerTracksAnalysisComponent::AliEMCalTriggerTracksAnalysisComponent(const char* name) :
  TNamed(name,""),
  fHistos(NULL),
  fBinning(NULL),
  fKineCuts(NULL),
  fTriggerDecision(NULL),
  fWeightHandler(NULL),
  fComponentDebugLevel(0)
{
}

/**
 * Create Container for histograms. Inheriting classes overwrite this method, in which they call
 * this and add the histograms of their choise.
 */
void AliEMCalTriggerTracksAnalysisComponent::CreateHistos() {
  fHistos = new AliEMCalHistoContainer(Form("Histos%s", GetName()));
  fHistos->ReleaseOwner();
}

/**
 * Create and define axis
 *
 * \param name Name of the axis
 * \param binning binning information
 * \return the new axis
 */
TAxis* AliEMCalTriggerTracksAnalysisComponent::DefineAxis(const char* name, const AliEMCalTriggerBinningDimension* binning) {
  TAxis *result = new TAxis(binning->GetNumberOfBins(), binning->GetBinLimits());
  result->SetName(name);
  return result;
}

/**
 * Create and define axis
 * \param name Name of the axis
 * \param nbins number of bins
 * \param min min. range
 * \param max max. range
 * \return the new axis
 */
TAxis* AliEMCalTriggerTracksAnalysisComponent::DefineAxis(const char* name, int nbins, double min, double max) {
  TAxis *result = new TAxis(nbins, min, max);
  result->SetName(name);
  return result;
}

/**
 * Get a set of names of trigger strings that is matching with the trigger decision.
 *
 * \param triggernames: output container for selected trigger names
 * \param usePatches: determines whether we use the trigger decision from patches
 */
void AliEMCalTriggerTracksAnalysisComponent::GetMachingTriggerNames(std::vector<std::string>& triggernames, ETriggerMethod_t method) {
  triggernames.clear();
  if(!fTriggerDecision) return;
  if(fTriggerDecision->IsMinBias()) triggernames.push_back("MinBias");
  if(fTriggerDecision->IsTriggered(kTAEMCJHigh, method)){
    triggernames.push_back("EMCJHigh");
    if(fTriggerDecision->IsTriggered(kTAEMCGHigh, method))
      triggernames.push_back("EMCHighBoth");
    else
      triggernames.push_back("EMCHighJetOnly");
  }
  if(fTriggerDecision->IsTriggered(kTAEMCJLow, method)){
    triggernames.push_back("EMCJLow");
    if(fTriggerDecision->IsTriggered(kTAEMCGLow, method))
      triggernames.push_back("EMCLowBoth");
    else
      triggernames.push_back("EMCLowJetOnly");
  }
  if(fTriggerDecision->IsTriggered(kTAEMCGHigh, method)){
    triggernames.push_back("EMCGHigh");
    if(!fTriggerDecision->IsTriggered(kTAEMCJHigh, method))
      triggernames.push_back("EMCHighGammaOnly");
  }
  if(fTriggerDecision->IsTriggered(kTAEMCGLow, method)){
    triggernames.push_back("EMCGLow");
    if(!fTriggerDecision->IsTriggered(kTAEMCJLow, method))
      triggernames.push_back("EMCLowGammaOnly");
  }
}

/**
 * Helper function to print the names of the selected trigger classes. For debugging purposes.
 * \param triggernames Selected trigger names
 * \param componentName Name of the component responsible for the printout
 */
void EMCalTriggerPtAnalysis::AliEMCalTriggerTracksAnalysisComponent::PrintTriggerNames(
    const std::vector<std::string>& triggernames, const std::string& componentName) const {

    std::cout << componentName << ", Triggers found: " << std::endl;
    std::cout << "==========================================" << std::endl;
    if(!triggernames.size()) return;
    std::vector<std::string>::const_iterator trgiter = triggernames.begin();
    std::cout << trgiter->c_str();
    trgiter++;
    while(trgiter != triggernames.end()){
      std::cout << ", " << trgiter->c_str();
      trgiter++;
    }
    std::cout << std::endl;
}

} /* namespace EMCalTriggerPtAnalysis */

