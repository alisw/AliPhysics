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
#include <map>
#include <string>

#include <TArrayD.h>
#include <TAxis.h>
#include <TBinning.h>
#include <THnSparse.h>
#include <THistManager.h>

#include "AliVEvent.h"
#include "AliVVertex.h"

#include "AliEMCalTriggerAnaClassManager.h"
#include "AliEMCalTriggerBinningComponent.h"
#include "AliEMCalTriggerEventData.h"
#include "AliEMCalTriggerEventCounterAnalysisComponent.h"

ClassImp(PWGJE::EMCALJetTasks::AliEMCalTriggerEventCounterAnalysisComponent)

using namespace PWGJE::EMCALJetTasks;

/**
 * Default (I/O) constructor, not to be used
 */
AliEMCalTriggerEventCounterAnalysisComponent::AliEMCalTriggerEventCounterAnalysisComponent():
  AliEMCalTriggerTracksAnalysisComponent()
{
}

/**
 * Main constructor
 */
AliEMCalTriggerEventCounterAnalysisComponent::AliEMCalTriggerEventCounterAnalysisComponent(const char *name):
  AliEMCalTriggerTracksAnalysisComponent(name)
{
}

/**
 * Create event counter histograms
 */
void AliEMCalTriggerEventCounterAnalysisComponent::CreateHistos() {
  AliEMCalTriggerTracksAnalysisComponent::CreateHistos();

  // Create trigger definitions
  std::map<std::string, std::string> triggerCombinations;
  GetAllTriggerNamesAndTitles(triggerCombinations);

  if(fComponentDebugLevel > 0){
    std::cout << "Event counter component - Found the following triggers:" << std::endl;
    for(std::map<std::string, std::string>::iterator it = triggerCombinations.begin(); it != triggerCombinations.end(); it++){
      std::cout << it->first << ", " << it->second << std::endl;
    }
  }

  const TBinning *vertexbinning = fBinning->GetBinning("zvertex");

  for(std::map<std::string,std::string>::iterator it = triggerCombinations.begin(); it != triggerCombinations.end(); ++it){
    const std::string name = it->first, &title = it->second;
    // Create event-based histogram
    fHistos->CreateTH1(Form("hEventHist%s", name.c_str()), Form("Event-based data for %s events; pileup rejection; z_{V} (cm)", title.c_str()), *vertexbinning);
  }

  // Make correlation histogram for different trigger classes
  const TAxis **triggeraxis = new const TAxis *[triggerCombinations.size()];
  memset(triggeraxis, 0, sizeof(const TAxis *) * triggerCombinations.size());
  const char *binlabels[2] = {"OFF", "ON"};
  TAxis *mytrgaxis = new TAxis[triggerCombinations.size()];
  std::map<std::string, std::string>::iterator trgiter = triggerCombinations.begin();
  for(int itrg = 0; itrg < triggerCombinations.size(); ++itrg){
    DefineAxis(mytrgaxis[itrg], trgiter->first.c_str(), trgiter->first.c_str(), 2, -0.5, 1.5, binlabels);
    triggeraxis[itrg] = mytrgaxis+itrg;
    trgiter++;
  }
  fHistos->CreateTHnSparse("hEventTriggers", "Trigger type per event", triggerCombinations.size(), triggeraxis);
  delete[] mytrgaxis;
  delete[] triggeraxis;
}

/**
 * Do event counting
 *  -# Select trigger class and fill vertex distribution
 *  -# Fill also correlation histogram
 */
void AliEMCalTriggerEventCounterAnalysisComponent::Process(const AliEMCalTriggerEventData* const data) {
  if(!fTriggerClassManager) return;

  double vz = data->GetRecEvent()->GetPrimaryVertex()->GetZ();
  TArrayD triggerCorrelation(fTriggerClassManager->GetAllTriggerClasses()->GetEntries());
  memset(triggerCorrelation.GetArray(), 0, sizeof(double) * triggerCorrelation.GetSize());

  std::vector<std::string> triggernames;
  this->GetMachingTriggerNames(triggernames);

  THnSparse *correlationhist = dynamic_cast<THnSparse *>(fHistos->FindObject("hEventTriggers"));
  for(std::vector<std::string>::iterator it = triggernames.begin(); it != triggernames.end(); it++){
    fHistos->FillTH1(Form("hEventHist%s", it->c_str()) ,vz);
    if(correlationhist){
      int idim = FindAxis(correlationhist, it->c_str());
      if(idim >= 0) triggerCorrelation[idim] = 1.;
    }
  }
  if(correlationhist) correlationhist->Fill(triggerCorrelation.GetArray());
}

/**
 * Define an axis with number of bins from min to max
 *
 * \param axis Axis to be defined
 * \param name Name of the axis
 * \param title Title of the axis
 * \param nbins Number of bins
 * \param min lower limit of the axis
 * \param max upper limit of the axis
 * \param labels array of bin labels (optional)
 */
void AliEMCalTriggerEventCounterAnalysisComponent::DefineAxis(TAxis& axis, const char* name,
    const char* title, int nbins, double min, double max,
    const char** labels) const {
  axis.Set(nbins, min, max);
  axis.SetName(name);
  axis.SetTitle(title);
  if(labels){
    for(int ib = 1; ib <= axis.GetNbins(); ++ib)
      axis.SetBinLabel(ib, labels[ib-1]);
  }
}

Int_t AliEMCalTriggerEventCounterAnalysisComponent::FindAxis(THnSparse* hist, const char* title) const {
  Int_t naxis = hist->GetNdimensions();
  Int_t result = -1;
  for(int idim = 0; idim < naxis; idim++){
    if(!TString(hist->GetAxis(idim)->GetName()).CompareTo(title)){
      result = idim;
      break;
    }
  }
  return result;
}
