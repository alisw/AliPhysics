#include <vector>

#include "AliAnalysisTaskEmcalEmbeddingHelper.h"
#include "AliAnalysisTaskEmcalEmbeddingHelperData.h"
#include "AliAODEvent.h"
#include "AliESDEvent.h"
#include "AliVAODHeader.h"
#include "AliEmcalList.h"
#include "AliAnalysisManager.h"

#include "TClonesArray.h"
#include "AliEmcalTrackSelectionAOD.h"
#include "AliTrackContainer.h"
#include "TObjArray.h"


/// \cond CLASSIMP
ClassImp(AliAnalysisTaskEmcalEmbeddingHelperData);
/// \endcond

AliAnalysisTaskEmcalEmbeddingHelperData::AliAnalysisTaskEmcalEmbeddingHelperData():
  AliAnalysisTaskEmcalEmbeddingHelper(),
  fPtSelection(),
  fTrackContainer(0x0)
{
}


AliAnalysisTaskEmcalEmbeddingHelperData::AliAnalysisTaskEmcalEmbeddingHelperData(const char *name):
  AliAnalysisTaskEmcalEmbeddingHelper(name),
  fPtSelection(),
  fTrackContainer(0x0)
{
}

AliAnalysisTaskEmcalEmbeddingHelperData::~AliAnalysisTaskEmcalEmbeddingHelperData()
{
}

void AliAnalysisTaskEmcalEmbeddingHelperData::ExecOnce()
{
//

}

void AliAnalysisTaskEmcalEmbeddingHelperData::UserCreateOutputObjects() 
{

  // run function from embedding helper base class
  AliAnalysisTaskEmcalEmbeddingHelper::UserCreateOutputObjects();

  TString histName = "fHistPtSelection";
  TString histTitle = "pT selection;Result;p_{T}^{lead}";
  auto histPtSelection= fHistManager.CreateTH2(histName, histTitle, 2, 0, 2, 200, 0, 100);
  histPtSelection->GetXaxis()->SetBinLabel(1,"Accepted");
  histPtSelection->GetXaxis()->SetBinLabel(2,"Rejected");

  fOutput->Add(histPtSelection);

  PostData(1, fOutput);

  // set track container array
  if(fTrackContainer) fTrackContainer->SetArray(InputEvent());
  else {
    AliFatal("must set track container to use this task");
  }

}

/**
 * Override inherited function and add check on track with max pT.
 * Rejects event if track within set pT range is not found
 */

Bool_t AliAnalysisTaskEmcalEmbeddingHelperData::CheckIsEmbeddedEventSelected()
{
  fTrackContainer->NextEvent(InputEvent());

  if(!AliAnalysisTaskEmcalEmbeddingHelper::CheckIsEmbeddedEventSelected()) return kFALSE;

  Int_t nTracksInRange = 0;
  Double_t ptMax = 0;
  for(auto track : fTrackContainer->accepted()) {
    Double_t pt = track->Pt();
    if(pt>ptMax) ptMax = pt;
    //now check set windows
    Int_t npt = fPtSelection.size();
    for(Int_t j=0;j<npt;j++) {
      if(pt>fPtSelection[j].first && pt<fPtSelection[j].second) { 
        // found track within pt
        nTracksInRange++;
      }
    }
  }

  TString histName = "fHistPtSelection";
  if(nTracksInRange==0) {
    fHistManager.FillTH2(histName, 1.5, ptMax );
    return kFALSE;
  }
  else {
    fHistManager.FillTH2(histName, 0.5, ptMax );
  }


  return kTRUE;
}


void AliAnalysisTaskEmcalEmbeddingHelperData::RetrieveTaskPropertiesFromYAMLConfig() {

  AliAnalysisTaskEmcalEmbeddingHelper::RetrieveTaskPropertiesFromYAMLConfig();

  // get pt selection from YAML file
  // defined in 2 arrays, pt

  std::vector <double> ptSelectionMin;
  std::vector <double> ptSelectionMax;
  bool resMin = fYAMLConfig.GetProperty( "ptSelectionMin", ptSelectionMin, false);
  bool resMax = fYAMLConfig.GetProperty( "ptSelectionMax", ptSelectionMax, false);
  if (resMin && resMax) {
    if (ptSelectionMin.size() != ptSelectionMax.size()) {
      AliErrorStream() << "Passed pt selection min with size" << ptSelectionMin.size() << " entries and pt selection max with size " << ptSelectionMax.size() << " - need same size arrays\n";
    }
    else {
      for(std::size_t i = 0; i < ptSelectionMin.size(); i++) {
        AliDebugStream(1) << "Setting pt selection [" << ptSelectionMin.at(i) << ", " << ptSelectionMax.at(i) << "]\n";
        fPtSelection.push_back(std::make_pair(ptSelectionMin.at(i), ptSelectionMax.at(i)));
      }
    }
  }
  else if( (resMin && !resMax) || (!resMin && resMax)) {
    AliErrorStream() << "Only set one pt selection array - should set both\n";
  }
}

/**
 *
 * The main implementation of the AddTask. 
 * Same structure as AddTaskEmbeddingHelper
 *
 */

AliAnalysisTaskEmcalEmbeddingHelperData *AliAnalysisTaskEmcalEmbeddingHelperData::AddTaskEmcalEmbeddingHelperData()
{  

  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalEmbeddingHelperData", "No analysis manager to connect to.");
    return 0;
  }

  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  AliVEventHandler* handler = mgr->GetInputEventHandler();
  if (!handler)
  {
    ::Error("AddTaskEmcalEmbeddingHelperData", "This task requires an input event handler");
    return 0;
  }

  TString name = "AliAnalysisTaskEmcalEmbeddingHelperData";

  AliAnalysisTaskEmcalEmbeddingHelperData * mgrTask = static_cast<AliAnalysisTaskEmcalEmbeddingHelperData *>(mgr->GetTask(name.Data()));
  if (mgrTask) return mgrTask;

  // Create the task that manages
  AliAnalysisTaskEmcalEmbeddingHelperData * embeddingHelper = new AliAnalysisTaskEmcalEmbeddingHelperData(name.Data());

  // set embedded track container
  AliTrackContainer *contEmb = new AliTrackContainer("tracks");
  contEmb->SetTrackCutsPeriod("LHC17p");
  embeddingHelper->SetEmbeddedTrackContainer(contEmb);

  //-------------------------------------------------------
  // Final settings, pass to manager and set the containers
  //-------------------------------------------------------

  mgr->AddTask(embeddingHelper);

  // Create containers for input/output
  AliAnalysisDataContainer* cInput = mgr->GetCommonInputContainer();

  TString outputContainerName(name);
  outputContainerName += "_histos";

  AliAnalysisDataContainer * cOutput = mgr->CreateContainer(outputContainerName.Data(),
      TList::Class(),
      AliAnalysisManager::kOutputContainer,
      Form("%s", AliAnalysisManager::GetCommonFileName()));

  mgr->ConnectInput(embeddingHelper, 0, cInput);
  mgr->ConnectOutput(embeddingHelper, 1, cOutput);

  return embeddingHelper;
}
