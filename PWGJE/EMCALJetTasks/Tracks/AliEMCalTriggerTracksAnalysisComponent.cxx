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
#include <TClass.h>
#include <TBinning.h>
#include <THistManager.h>

#include "AliEMCalTriggerAnaClassManager.h"
#include "AliEMCalTriggerAnaTriggerClass.h"
#include "AliEMCalTriggerTracksAnalysisComponent.h"
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerAnaTriggerDecision.h"

ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerTracksAnalysisComponent)

using namespace PWGJE::EMCALJetTasks;

/**
 * Dummy (I/O) constructor, not to be used
 */
AliEMCalTriggerTracksAnalysisComponent::AliEMCalTriggerTracksAnalysisComponent() :
  TNamed(),
  fHistos(NULL),
  fTriggerClassManager(NULL),
  fBinning(NULL),
  fKineCuts(NULL),
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
  fTriggerClassManager(NULL),
  fBinning(NULL),
  fKineCuts(NULL),
  fWeightHandler(NULL),
  fComponentDebugLevel(0)
{
}

/**
 * Create Container for histograms. Inheriting classes overwrite this method, in which they call
 * this and add the histograms of their choise.
 */
void AliEMCalTriggerTracksAnalysisComponent::CreateHistos() {
  fHistos = new THistManager(Form("Histos%s", GetName()));
  fHistos->ReleaseOwner();
}

/**
 * Create and define axis
 *
 * \param name Name of the axis
 * \param binning binning information
 * \return the new axis
 */
TAxis* AliEMCalTriggerTracksAnalysisComponent::DefineAxis(const char* name, const TBinning &binning) {
  TArrayD binedges;
  binning.CreateBinEdges(binedges);
  TAxis *result = new TAxis(binedges.GetSize() - 1, binedges.GetArray());
  result->SetName(name);
  return result;
}

/**
 * Get a set of names of trigger strings that is matching with the trigger decision.
 * \param triggernames: output container for selected trigger names
 * \throw TriggerHandlerNotFoundException in case no trigger handler is available
 */
void AliEMCalTriggerTracksAnalysisComponent::GetMachingTriggerNames(std::vector<std::string>& triggernames) const {
  triggernames.clear();
  if(!fTriggerClassManager) throw TriggerManagerNotFoundException(this->IsA()->GetName());
  for(TIter trgiter = TIter(fTriggerClassManager->GetSelectedTriggerClasses()).Begin(); trgiter != TIter::End(); ++ trgiter){
    triggernames.push_back((static_cast<AliEMCalTriggerAnaTriggerClass *>(*trgiter))->GetName());
  }
}

/**
 * Get trigger names and titles for the event
 * \param triggers Map with trigger names and titles
 * \throw TriggerHandlerNotFoundException in case no trigger handler is available
 */
void AliEMCalTriggerTracksAnalysisComponent::GetAllTriggerNamesAndTitles(std::map<std::string, std::string> &triggers) const {
  triggers.clear();
  if(!fTriggerClassManager) {
    throw TriggerManagerNotFoundException(this->IsA()->GetName());
  }
  for(TIter trgiter = TIter(fTriggerClassManager->GetAllTriggerClasses()).Begin(); trgiter != TIter::End(); ++ trgiter){
    AliEMCalTriggerAnaTriggerClass *trgcls = static_cast<AliEMCalTriggerAnaTriggerClass *>(*trgiter);
    triggers.insert(std::pair<std::string, std::string>(trgcls->GetName(), trgcls->GetTitle()));
  }
}

/**
 * Helper function to print the names of the selected trigger classes. For debugging purposes.
 * \param triggernames Selected trigger names
 * \param componentName Name of the component responsible for the printout
 */
void AliEMCalTriggerTracksAnalysisComponent::PrintTriggerNames(
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
