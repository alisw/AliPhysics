///
/// \file AliFemtoMultCorrAnalysis.cxx
/// \author Jeremi Niedziela

#include "AliFemtoMultCorrAnalysis.h"
#include "AliFemtoPicoEvent.h"
#include "AliFemtoTrackCut.h"
#include "AliFemtoV0Cut.h"

#include <iostream>
#include <iterator>

#ifdef __ROOT__
/// \cond CLASSIMP
ClassImp(AliFemtoMultCorrAnalysis);
/// \endcond
#endif

AliFemtoMultCorrAnalysis::AliFemtoMultCorrAnalysis():
  AliFemtoSimpleAnalysis(),
fEventCut(nullptr),
fNeventsPassed(0)
{
  fMotherCut[0] = nullptr;
  fMotherCut[1] = nullptr;
  fDaughterCut[0] = nullptr;
  fDaughterCut[1] = nullptr;
  fDaughterCut[2] = nullptr;
  
  fMultCorrFctn             = new TH2D("multCorrFctn","multCorrFctn",500,0,500,1000,0,1.0);
  fMultCorrTimesMultFctn    = new TH2D("multCorrTimesMultFctn","multCorrTimesMultFctn",500,0,500,1000,0,50);
  
  fNmotherNdaughterCorrFctn = new TH2D("NmotherNdaughterCorrFctn","NmotherNdaughterCorrFctn",500,0,500,1000,0,0.01);
  fNmotherNdaughterRootsCorrFctn = new TH2D("NmotherNdaughterRootsCorrFctn","NmotherNdaughterRootsCorrFctn",500,0,500,1000,0,1.0);
}

AliFemtoMultCorrAnalysis::~AliFemtoMultCorrAnalysis()
{
  if(fEventCut)        delete fEventCut;
  if(fMotherCut[0])    delete fMotherCut[0];
  if(fMotherCut[1])    delete fMotherCut[1];
  if(fDaughterCut[0])    delete fDaughterCut[0];
  if(fDaughterCut[1])    delete fDaughterCut[1];
  if(fDaughterCut[2])    delete fDaughterCut[2];
  
  if(fMultCorrFctn) delete fMultCorrFctn;
  if(fMultCorrTimesMultFctn) delete fMultCorrTimesMultFctn;
  if(fNmotherNdaughterCorrFctn) delete fNmotherNdaughterCorrFctn;
  if(fNmotherNdaughterRootsCorrFctn) delete fNmotherNdaughterRootsCorrFctn;
}

void AliFemtoMultCorrAnalysis::ProcessEvent(const AliFemtoEvent* currentEvent)
{
  EventBegin(currentEvent);
  
  bool tmpPassEvent = fEventCut->Pass(currentEvent);
  fEventCut->FillCutMonitor(currentEvent, tmpPassEvent);
  if (!tmpPassEvent) {
    EventEnd(currentEvent);  // cleanup for EbyE
    return;
  }
  
  AliFemtoParticleCollection *motherCollection[2];
  AliFemtoParticleCollection *daughterCollection[3];
  
  motherCollection[0] = new AliFemtoParticleCollection();
  motherCollection[1] = new AliFemtoParticleCollection();
  daughterCollection[0] = new AliFemtoParticleCollection();
  daughterCollection[1] = new AliFemtoParticleCollection();
  daughterCollection[2] = new AliFemtoParticleCollection();
  
  FillHbtParticleCollection(fMotherCut[0],   currentEvent,motherCollection[0]);
  FillHbtParticleCollection(fMotherCut[1],   currentEvent,motherCollection[1]);
  FillHbtParticleCollection(fDaughterCut[0], currentEvent,daughterCollection[0]);
  FillHbtParticleCollection(fDaughterCut[1], currentEvent,daughterCollection[1]);
  FillHbtParticleCollection(fDaughterCut[2], currentEvent,daughterCollection[2]);
  
  int mothersMult = motherCollection[0]->size() + motherCollection[1]->size();
  int mothersMultProduct = motherCollection[0]->size() * motherCollection[1]->size();
  int daughtersMult = daughterCollection[0]->size() + daughterCollection[1]->size() + daughterCollection[2]->size();
  int daughtersMultProduct = daughterCollection[0]->size() * daughterCollection[1]->size() * daughterCollection[2]->size();
  int totalMult = currentEvent->NumberOfTracks();
  
  if(daughtersMult > 0){
    fMultCorrFctn->Fill(totalMult, mothersMult/(double)(daughtersMult));
    fMultCorrTimesMultFctn->Fill(totalMult, mothersMult*totalMult/(double)(daughtersMult));
  }
  
  if(daughtersMultProduct > 0){
    fNmotherNdaughterCorrFctn->Fill(totalMult, mothersMultProduct/(double)(daughtersMultProduct));
    fNmotherNdaughterRootsCorrFctn->Fill(totalMult, sqrt(mothersMultProduct)/(double)pow(daughtersMultProduct,1/3.));
    cout<<mothersMultProduct/(double)(daughtersMultProduct)<<"\t"<<sqrt(mothersMultProduct)/(double)pow(daughtersMultProduct,1/3.)<<endl;
  }
  
  delete motherCollection[0];
  delete motherCollection[1];
  delete daughterCollection[0];
  delete daughterCollection[1];
  delete daughterCollection[2];
  
  EventEnd(currentEvent);
  fNeventsPassed++;
}

void AliFemtoMultCorrAnalysis::FillHbtParticleCollection(AliFemtoParticleCut *cut,
                                                         const AliFemtoEvent *hbtEvent,
                                                         AliFemtoParticleCollection *outputCollection)
{
  if(cut->Type() == hbtTrack) {
    for (const auto &track : *(hbtEvent->TrackCollection())){
      const bool track_passes = ((AliFemtoTrackCut*)cut)->Pass(track);
      cut->FillCutMonitor(track, track_passes);
      if (track_passes) {
        outputCollection->push_back(new AliFemtoParticle(track, cut->Mass()));
      }
    }
  }
  else if(cut->Type() == hbtV0){
    for (const auto &track : *(hbtEvent->V0Collection())){
      const bool track_passes = ((AliFemtoV0Cut*)cut)->Pass(track);
      cut->FillCutMonitor(track, track_passes);
      if (track_passes) {
        outputCollection->push_back(new AliFemtoParticle(track, cut->Mass()));
      }
    }
  }
}


TList* AliFemtoMultCorrAnalysis::GetOutputList()
{
  TList *outputList = new TList();
  outputList->Add(fMultCorrFctn);
  outputList->Add(fMultCorrTimesMultFctn);
  outputList->Add(fNmotherNdaughterCorrFctn);
  outputList->Add(fNmotherNdaughterRootsCorrFctn);
  
  return outputList;
}
