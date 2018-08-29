///
/// \file AliFemtoTrioAnalysis.cxx
/// \author Jeremi Niedziela

#include "AliFemtoTrioAnalysis.h"
#include "AliFemtoPicoEvent.h"

#include <iostream>
#include <iterator>

#ifdef __ROOT__
/// \cond CLASSIMP
ClassImp(AliFemtoTrioAnalysis);
/// \endcond
#endif

AliFemtoEventCut*    copyTheCut(AliFemtoEventCut*);
AliFemtoParticleCut* copyTheCut(AliFemtoParticleCut*);
AliFemtoCorrFctn*    copyTheCorrFctn(AliFemtoCorrFctn*);

extern void FillHbtParticleCollection(AliFemtoParticleCut* partCut,
                                      const AliFemtoEvent* currentEvent,
                                      AliFemtoParticleCollection* partCollection,
                                      bool performSharedDaughterCut=kFALSE);

AliFemtoTrioAnalysis::AliFemtoTrioAnalysis():
AliFemtoSimpleAnalysis(),
fTrioFctnCollection(NULL),
fEventCut(NULL),
fFirstParticleCut(NULL),
fSecondParticleCut(NULL),
fThirdParticleCut(NULL),
fNeventsProcessed(0),
fNeventsPassed(0),
fPicoEvent(nullptr),
fPerformSharedDaughterCut(kFALSE),
fDoEventMixing(false),
fCollection1type(AliFemtoTrio::kUnknown),
fCollection2type(AliFemtoTrio::kUnknown),
fCollection3type(AliFemtoTrio::kUnknown)
{
  // Default constructor
  fTrioFctnCollection = new AliFemtoTrioFctnCollection();
  for(int i=0;i<3;i++){
    fMixingBuffer[i] = nullptr;
  }
}

AliFemtoTrioAnalysis::~AliFemtoTrioAnalysis()
{
  cout << " AliFemtoTrioAnalysis::~AliFemtoTrioAnalysis()" << endl;
  
  if(fEventCut)           delete fEventCut;
  if(fFirstParticleCut)   delete fFirstParticleCut;
  if(fSecondParticleCut)  delete fSecondParticleCut;
  if(fThirdParticleCut)   delete fThirdParticleCut;

  if (fTrioFctnCollection) {
    for (AliFemtoTrioFctnIterator iter = fTrioFctnCollection->begin(); iter != fTrioFctnCollection->end(); iter++) {
      if(*iter) delete *iter;
    }
    if(fTrioFctnCollection) delete fTrioFctnCollection;
  }
  
  for(int i=0;i<3;i++){
    if (fMixingBuffer[i]) {delete fMixingBuffer[i];}
  }
}

void AliFemtoTrioAnalysis::ProcessEvent(const AliFemtoEvent* currentEvent)
{
  fPicoEvent = nullptr;
  EventBegin(currentEvent);
  bool tmpPassEvent = fEventCut->Pass(currentEvent);
  
  if (!tmpPassEvent) {
    fEventCut->FillCutMonitor(currentEvent, tmpPassEvent);
    EventEnd(currentEvent);  // cleanup for EbyE
    return;
  }
  
  fPicoEvent = new AliFemtoPicoEvent();
  AliFemtoParticleCollection *collection1 = fPicoEvent->FirstParticleCollection();
  AliFemtoParticleCollection *collection2 = fPicoEvent->SecondParticleCollection();
  AliFemtoParticleCollection *collection3 = fPicoEvent->ThirdParticleCollection();
  
  if (!collection1 || !collection2 || !collection3){
    cout << "E-AliFemtoTrioAnalysis::ProcessEvent: new PicoEvent is missing particle collections!\n";
    EventEnd(currentEvent);  // cleanup for EbyE
    if(fPicoEvent) delete fPicoEvent;
    return;
  }
  
  FillHbtParticleCollection(fFirstParticleCut, currentEvent,collection1,fPerformSharedDaughterCut);
  FillHbtParticleCollection(fSecondParticleCut,currentEvent,collection2,fPerformSharedDaughterCut);
  FillHbtParticleCollection(fThirdParticleCut, currentEvent,collection3,fPerformSharedDaughterCut);

  fEventCut->FillCutMonitor(currentEvent, tmpPassEvent);

  if (!tmpPassEvent){
    EventEnd(currentEvent);
    if(fPicoEvent) delete fPicoEvent;
    return;
  }
  
  // add real events
  bool mixing = false;
  AddParticles(collection1,collection2,collection3,mixing);
  // add mixed events (if enough entries in the mixing buffer)
  if(fDoEventMixing){
    if(fNeventsPassed>2){
      mixing = true;
      collection1 = fMixingBuffer[0]->FirstParticleCollection();
      collection2 = fMixingBuffer[1]->SecondParticleCollection();
      collection3 = fMixingBuffer[2]->ThirdParticleCollection();
      AddParticles(collection1,collection2,collection3,mixing);
      
      collection1 = fMixingBuffer[0]->FirstParticleCollection();
      collection2 = fMixingBuffer[2]->SecondParticleCollection();
      collection3 = fMixingBuffer[1]->ThirdParticleCollection();
      AddParticles(collection1,collection2,collection3,mixing);
      
      collection1 = fMixingBuffer[1]->FirstParticleCollection();
      collection2 = fMixingBuffer[2]->SecondParticleCollection();
      collection3 = fMixingBuffer[0]->ThirdParticleCollection();
      AddParticles(collection1,collection2,collection3,mixing);
      
      collection1 = fMixingBuffer[1]->FirstParticleCollection();
      collection2 = fMixingBuffer[0]->SecondParticleCollection();
      collection3 = fMixingBuffer[2]->ThirdParticleCollection();
      AddParticles(collection1,collection2,collection3,mixing);
      
      collection1 = fMixingBuffer[2]->FirstParticleCollection();
      collection2 = fMixingBuffer[0]->SecondParticleCollection();
      collection3 = fMixingBuffer[1]->ThirdParticleCollection();
      AddParticles(collection1,collection2,collection3,mixing);
      
      collection1 = fMixingBuffer[2]->FirstParticleCollection();
      collection2 = fMixingBuffer[1]->SecondParticleCollection();
      collection3 = fMixingBuffer[0]->ThirdParticleCollection();
      AddParticles(collection1,collection2,collection3,mixing);
      
    }
    // delete the oldest event, shift others and save the current one in the buffer
    if(fMixingBuffer[2]) delete fMixingBuffer[2];
    
    fMixingBuffer[2] = fMixingBuffer[1];
    fMixingBuffer[1] = fMixingBuffer[0];
    fMixingBuffer[0] = fPicoEvent;
  }
  else{
    if(fPicoEvent) delete fPicoEvent;
  }
  EventEnd(currentEvent);
  fNeventsPassed++;
}

//_________________________
void AliFemtoTrioAnalysis::AddParticles(AliFemtoParticleCollection *collection1,
                                            AliFemtoParticleCollection *collection2,
                                            AliFemtoParticleCollection *collection3,
                                            bool mixing)
{
  AliFemtoTrio *trio = new AliFemtoTrio();
  
  for (AliFemtoParticleConstIterator iPart1 = collection1->begin();iPart1 != collection1->end();++iPart1){
    trio->SetTrack1((AliFemtoParticle *)(*iPart1),fCollection1type);
    
    for (AliFemtoParticleConstIterator  iPart2 = collection2->begin();iPart2 != collection2->end();++iPart2){
      trio->SetTrack2((AliFemtoParticle *)(*iPart2),fCollection2type);
      
      for (AliFemtoParticleConstIterator  iPart3 = collection3->begin();iPart3 != collection3->end();++iPart3){
        trio->SetTrack3((AliFemtoParticle *)(*iPart3),fCollection3type);
        
        for (AliFemtoTrioFctnIterator iFun = fTrioFctnCollection->begin();iFun != fTrioFctnCollection->end();++iFun){
          AliFemtoTrioFctn *triofctn = *iFun;
    
          if(mixing)  triofctn->AddMixedTrio(trio);
          else        triofctn->AddRealTrio(trio);
        }
      }
    }
  }
  if(trio) delete trio;
}


AliFemtoEventCut*      AliFemtoTrioAnalysis::EventCut(){return fEventCut;}
AliFemtoParticleCut*   AliFemtoTrioAnalysis::FirstParticleCut(){return fFirstParticleCut;}
AliFemtoParticleCut*   AliFemtoTrioAnalysis::SecondParticleCut(){return fSecondParticleCut;}
AliFemtoParticleCut*   AliFemtoTrioAnalysis::ThirdParticleCut(){return fThirdParticleCut;}

void AliFemtoTrioAnalysis::AddTrioFctn(AliFemtoTrioFctn* function){fTrioFctnCollection->push_back(function);}

void AliFemtoTrioAnalysis::SetEventCut(AliFemtoEventCut* cut){fEventCut = cut;cut->SetAnalysis(this);}
void AliFemtoTrioAnalysis::SetFirstParticleCut(AliFemtoParticleCut* cut){fFirstParticleCut = cut;cut->SetAnalysis(this);}
void AliFemtoTrioAnalysis::SetSecondParticleCut(AliFemtoParticleCut* cut){fSecondParticleCut = cut;cut->SetAnalysis(this);}
void AliFemtoTrioAnalysis::SetThirdParticleCut(AliFemtoParticleCut* cut){fThirdParticleCut = cut;cut->SetAnalysis(this);}

void AliFemtoTrioAnalysis::SetCollection1type(AliFemtoTrio::EPart type){fCollection1type=type;}
void AliFemtoTrioAnalysis::SetCollection2type(AliFemtoTrio::EPart type){fCollection2type=type;}
void AliFemtoTrioAnalysis::SetCollection3type(AliFemtoTrio::EPart type){fCollection3type=type;}

void AliFemtoTrioAnalysis::SetDoEventMixing(bool mix){fDoEventMixing = mix;}


void AliFemtoTrioAnalysis::SetV0SharedDaughterCut(bool perform){fPerformSharedDaughterCut = perform;}
bool AliFemtoTrioAnalysis::V0SharedDaughterCut(){return fPerformSharedDaughterCut;}


AliFemtoString AliFemtoTrioAnalysis::Report(){return "Report";}
int AliFemtoTrioAnalysis::GetNeventsProcessed(){return fNeventsProcessed;}
TList* AliFemtoTrioAnalysis::ListSettings(){return nullptr;}

void AliFemtoTrioAnalysis::EventBegin(const AliFemtoEvent* TheEventToBegin){};
void AliFemtoTrioAnalysis::EventEnd(const AliFemtoEvent* TheEventToWrapUp){};

TList* AliFemtoTrioAnalysis::GetOutputList()
{
  TList *outputList = new TList();
  
  TList *p1Cut = fFirstParticleCut->GetOutputList();
  TListIter nextp1(p1Cut);
  while (TObject *obj = nextp1.Next()) {outputList->Add(obj);}
  if(p1Cut) delete p1Cut;
  
  if (fSecondParticleCut != fFirstParticleCut) {
    TList *p2Cut = fSecondParticleCut->GetOutputList();
    TIter nextp2(p2Cut);
    while (TObject *obj = nextp2()) {outputList->Add(obj);}
    if(p2Cut) delete p2Cut;
  }
  
  if (fThirdParticleCut != fFirstParticleCut) {
    TList *p3Cut = fThirdParticleCut->GetOutputList();
    TIter nextp3(p3Cut);
    while (TObject *obj = nextp3()) {outputList->Add(obj);}
    if(p3Cut) delete p3Cut;
  }
  
  TList *eventCut = fEventCut->GetOutputList();
  TIter nextEvent(eventCut);
  while (TObject *obj = nextEvent()) {outputList->Add(obj);}
  if(eventCut) delete eventCut;
  
  for (AliFemtoTrioFctnIterator iter = fTrioFctnCollection->begin();iter != fTrioFctnCollection->end();++iter){
    TList *trioFunctionList = (*iter)->GetOutputList();
    TIter nextTrioFctn(trioFunctionList);
    while (TObject *obj = nextTrioFctn()) {outputList->Add(obj);}
    if(trioFunctionList) delete trioFunctionList;
  }
  
//  TH1D *nEvents = new TH1D("nEvents","nEvents",999999999,0,999999999);
//  nEvents->Fill(fNeventsPassed);
//  outputList->Add(nEvents);
  
  return outputList;
}
