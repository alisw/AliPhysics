#include <TList.h>
#include <TMath.h>
#include <TTree.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TNtuple.h>
#include <TROOT.h>
#include <TInterpreter.h>

#include "AliAnalysisTaskCorrelationsDev.h"
#include "AliAnalyseLeadingTrackUE.h"
#include "AliPIDResponse.h"

#include "AliInputEventHandler.h"
#include "AliAnalysisCuts.h"

#include "AliLog.h"
#include "AliVParticle.h"
#include "AliVTrack.h"

#include "AliCentrality.h"
#include "AliMultSelection.h"

ClassImp(AliAnalysisTaskCorrelationsDev)

//____________________________________________________________________
AliAnalysisTaskCorrelationsDev:: AliAnalysisTaskCorrelationsDev(const char* name):
AliAnalysisTaskSE(name),
// general configuration
fAnalyseUE(0x0),
fDebug(0),
fListOfHistos(0x0), 
fEventCuts(0),
fSelectBit(0),
fZVertex(-1),
fCentralityMethod("V0M"),
fUseUncheckedCentrality(kFALSE),
fUseNewCentralityFramework(kFALSE),
fTrackEtaCut(10),
fTrackEtaCutMin(-1.),
fPtMin(0),
fFilterBit(0),
fTrackStatus(0),
fOutputTree(0),
fOutputContainer(0),
fTreeEventConfig(0),
fTreeMaxTracks(0),
fTreeTrackConfig(0)
{
  // Default constructor

  DefineOutput(1, TList::Class());
  DefineOutput(2, TTree::Class());
}

AliAnalysisTaskCorrelationsDev::~AliAnalysisTaskCorrelationsDev() 
{ 
  // destructor
  
  if (fListOfHistos) 
    delete fListOfHistos;
}

//____________________________________________________________________
void  AliAnalysisTaskCorrelationsDev::UserCreateOutputObjects()
{
  // Create the output container
  
  AliInfo(Form("Particle selection: filter bit = %d | eta_min = %f | eta_max = %f | pt_min = %f | track_status = %d", fFilterBit, fTrackEtaCutMin, fTrackEtaCut, fPtMin, fTrackStatus));
  
  // Initialize class with main algorithms, event and track selection. 
  fAnalyseUE = new AliAnalyseLeadingTrackUE();
  fAnalyseUE->SetParticleSelectionCriteria(fFilterBit, kFALSE, fTrackEtaCut, fTrackEtaCutMin, fPtMin);
  fAnalyseUE->SetTrackStatus(fTrackStatus);
  fAnalyseUE->SetDebug(fDebug); 
  fAnalyseUE->SetEventSelection(fSelectBit);
  
  // create output tree
  if (fTreeEventConfig->GetEntriesFast() > 0)
  {
    TString varlist = "trigger";
    for (int i=0; i<fTreeEventConfig->GetEntriesFast(); i++) {
      varlist += ":";
      varlist += fTreeEventConfig->At(i)->GetName();
    }
    
    varlist += ":ntracks";
    
    for (int j=0; j<fTreeMaxTracks; j++)
      for (int i=0; i<fTreeTrackConfig->GetEntriesFast(); i++)
        varlist += Form(":%s_%d", fTreeTrackConfig->At(i)->GetName(), j);
    
    varlist.ReplaceAll("(", "");
    varlist.ReplaceAll(")", "");
    varlist.ReplaceAll("[", "");
    varlist.ReplaceAll("]", "");
    varlist.ReplaceAll(".", "_");
    varlist.ReplaceAll("->", "_");
      
    fOutputTree = new TNtuple("fOutputTree", "fOutputTree", varlist);
    fOutputContainer = new Float_t[fOutputTree->GetNvar()];

    Printf("Created tree with %d branches: %s", fOutputTree->GetNvar(), varlist.Data());

    PostData(2, fOutputTree);
  }

  // Initialize output list of containers
  if (fListOfHistos != NULL){
    delete fListOfHistos;
    fListOfHistos = NULL;
  }
  if (!fListOfHistos){
    fListOfHistos = new TList();
    fListOfHistos->SetOwner(kTRUE); 
  }
  
  fListOfHistos->Add(new TH2F("tofSignalvsMult", ";multiplicity;tofSignal", 10, 0.5, 10.5, 1000, 0, 40000));
  fListOfHistos->Add(new TH2F("tofSignalvsTOFMult", ";multiplicity;tofSignal", 10, 0.5, 10.5, 1000, 0, 40000));
  
  PostData(1, fListOfHistos);
  
  AddSettingsTree();
}

//____________________________________________________________________
void  AliAnalysisTaskCorrelationsDev::AddSettingsTree()
{
  // Write settings to output list
  TTree *settingsTree   = new TTree("AnalysisSettings","Analysis Settings");
  settingsTree->Branch("fZVertex", &fZVertex,"ZVertex/D");
  //settingsTree->Branch("fCentralityMethod", fCentralityMethod.Data(),"CentralityMethod/C");
  settingsTree->Branch("fTrackEtaCut", &fTrackEtaCut, "TrackEtaCut/D");
  settingsTree->Branch("fTrackEtaCutMin", &fTrackEtaCutMin, "TrackEtaCutMin/D");
  settingsTree->Branch("fPtMin", &fPtMin, "PtMin/D");
  settingsTree->Branch("fFilterBit", &fFilterBit,"FilterBit/I");
  settingsTree->Branch("fTrackStatus", &fTrackStatus,"TrackStatus/I");
  settingsTree->Branch("fSelectBit", &fSelectBit,"EventSelectionBit/I");
  settingsTree->Branch("fUseUncheckedCentrality", &fUseUncheckedCentrality,"UseUncheckedCentrality/O");
  settingsTree->Branch("fUseNewCentralityFramework", &fUseNewCentralityFramework,"fUseNewCentralityFramework/O");
  
  settingsTree->Fill();
  fListOfHistos->Add(settingsTree);
}  

//____________________________________________________________________
void  AliAnalysisTaskCorrelationsDev::UserExec(Option_t */*option*/)
{
  // exec (per event)
  fAnalyseUE->NextEvent();

  // skip not selected events here (the AOD is not updated for those)
  if (fSelectBit != 0 && !(fInputHandler->IsEventSelected() & fSelectBit))
    return;
  
  // Trigger selection ************************************************
  if (fSelectBit != 0 && !fAnalyseUE->TriggerSelection(fInputHandler)) 
    return;

  // Vertex selection *************************************************
  if (fZVertex > 0 && !fAnalyseUE->VertexSelection(fInputEvent, 1, fZVertex)) 
    return;
  
  // general event selection
  if (fEventCuts && !fEventCuts->IsSelected(fInputEvent))
    return;
  
  Double_t centrality = GetCentrality(fInputEvent, 0);

  TObjArray* tracks = fAnalyseUE->GetAcceptedParticles(fInputEvent, 0, kTRUE, -1, kTRUE);
  
  Int_t totalTracks = tracks->GetEntriesFast();
  Int_t tofTracks = 0;
  for (int i=0; i<totalTracks; i++)
  {
    AliVTrack* track = dynamic_cast<AliVTrack*> (tracks->At(i));
    if (track->GetTOFsignal() < 90000)
      tofTracks++;
  }
  
  for (int i=0; i<totalTracks; i++)
  {
    AliVTrack* track = dynamic_cast<AliVTrack*> (tracks->At(i));
    
    Float_t tofSignal = track->GetTOFsignal();
    
    if (tofSignal < 90000) {
      ((TH2F*) fListOfHistos->FindObject("tofSignalvsMult"))->Fill(totalTracks, track->GetTOFsignal());
      ((TH2F*) fListOfHistos->FindObject("tofSignalvsTOFMult"))->Fill(tofTracks, track->GetTOFsignal());
    }
  }

  if (fDebug > 1)
    Printf("Tracks: %d %d", totalTracks, tofTracks);
  
  if (fOutputTree)
  {
    for (int i=0; i<fOutputTree->GetNvar(); i++)
      fOutputContainer[i] = 0;
    
    Int_t containerCounter = 0;
    fOutputContainer[containerCounter++] = fInputHandler->IsEventSelected();
    
    // event properties
    for (int i=0; i<fTreeEventConfig->GetEntriesFast(); i++)
    {
      TString command(fTreeEventConfig->At(i)->GetName());
      
      TObject* input = fInputEvent;
      if (command.BeginsWith("_pid->"))
      {
        input = fInputHandler->GetPIDResponse();
        command.Remove(0, strlen("_pid->"));
      }
      
      Float_t value = GetValueInterpreted(input, command);
      fOutputContainer[containerCounter++] = value;
      
      if (fDebug > 5)
        Printf("E %d %s %f", i, command.Data(), value);
    }
    
    // tracks
    fOutputContainer[containerCounter++] = tracks->GetEntriesFast();
    Int_t nTracksToWrite = TMath::Min(fTreeMaxTracks, tracks->GetEntriesFast());
    
    for (int itrack=0; itrack<nTracksToWrite; itrack++)
    {
      for (int i=0; i<fTreeTrackConfig->GetEntriesFast(); i++)
      {
        Float_t value = GetValueInterpreted(dynamic_cast<AliVTrack*> (tracks->At(itrack)), fTreeTrackConfig->At(i)->GetName());
        fOutputContainer[containerCounter++] = value;
      
        if (fDebug > 5)
          Printf("T %d %d %s %f", itrack, i, fTreeTrackConfig->At(i)->GetName(), value);
      }
    }
    
    fOutputTree->Fill(fOutputContainer);
  }

  delete tracks;
}

Float_t AliAnalysisTaskCorrelationsDev::GetValueInterpreted(TObject* source, const char* command)
{
  // get value from the interpreter to allow configuration

  // with some float <-> int casting magic...
  Int_t error = 0;
  Int_t tmpInt = (Int_t) gROOT->ProcessLine(Form("Float_t tmpFloat = ((%s*) %p)->%s; Int_t tmpInt; memcpy(&tmpInt, &tmpFloat, sizeof(Float_t)); tmpInt;", source->ClassName(), (void*) source, command), &error);
  
  Float_t value;
  if (error == TInterpreter::kNoError)
    memcpy(&value, &tmpInt, sizeof(Float_t));
  else
    value = -999;
  
  return value;
}

Double_t AliAnalysisTaskCorrelationsDev::GetCentrality(AliVEvent* inputEvent, TObject* mc)
{
  // return centrality
  
  if (fCentralityMethod.Length() == 0)
    return 0;
  
  Double_t centrality = 0;
  
  if (fUseNewCentralityFramework) 
  {
    AliMultSelection *multSelection = (AliMultSelection*) inputEvent->FindListObject("MultSelection");
    if (!multSelection)
      AliFatal("MultSelection not found in input event");
      
    if (fUseUncheckedCentrality)
      centrality = multSelection->GetMultiplicityPercentile(fCentralityMethod, kFALSE);
    else
      centrality = multSelection->GetMultiplicityPercentile(fCentralityMethod, kTRUE);
    
    // error handling
    if (centrality > 100)
      centrality = -1;
  }
  else 
  {
    AliCentrality *centralityObj = 0;

    if (fCentralityMethod == "TRACKS_MANUAL")
    {
      // for pp
      TObjArray* tracks = fAnalyseUE->GetAcceptedParticles(inputEvent, 0, kTRUE, -1, kTRUE);
      centrality = tracks->GetEntriesFast();
      if (centrality > 40)
        centrality = 41;
//       Printf("%d %f", tracks->GetEntriesFast(), centrality);

      delete tracks;
    }
    else
    {
      centralityObj = inputEvent->GetCentrality();
      
      if (centralityObj)
      {
        if (fUseUncheckedCentrality)
          centrality = centralityObj->GetCentralityPercentileUnchecked(fCentralityMethod);
        else
          centrality = centralityObj->GetCentralityPercentile(fCentralityMethod);
      }
      else
        centrality = -1;
    }
  }
  AliInfo(Form("Centrality is %f", centrality));
  
  return centrality;
}
