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
fNeventsPassed(0),
fPicoEvent(nullptr),
fPerformSharedDaughterCut(kFALSE),
fMinSizePart1Collection(1),
fMinSizePart2Collection(1),
fMinSizePart3Collection(1),
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

  // increment number of events processed
  AddEventProcessed();
  
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

  

  fEventCut->FillCutMonitor(currentEvent, tmpPassEvent);

  bool coll_1_size_passes = kFALSE;
  bool coll_2_size_passes = kFALSE;
  bool coll_3_size_passes = kFALSE;
  
  if(fFirstParticleCut && fSecondParticleCut && fThirdParticleCut)
    {
      FillHbtParticleCollection(fFirstParticleCut, currentEvent,collection1,fPerformSharedDaughterCut);
      FillHbtParticleCollection(fSecondParticleCut, currentEvent,collection2,fPerformSharedDaughterCut);
      FillHbtParticleCollection(fThirdParticleCut, currentEvent,collection3,fPerformSharedDaughterCut);

      coll_1_size_passes = (collection1->size() >= fMinSizePart1Collection);
      coll_2_size_passes = (collection2->size() >= fMinSizePart2Collection);
      coll_3_size_passes = (collection3->size() >= fMinSizePart3Collection);

      fEventCut->FillCutMonitor(collection1, collection2);
    }
  if(fFirstParticleCut && fSecondParticleCut && fThirdParticleCut==NULL)
    {
      FillHbtParticleCollection(fFirstParticleCut, currentEvent,collection1,fPerformSharedDaughterCut);
      FillHbtParticleCollection(fSecondParticleCut, currentEvent,collection2,fPerformSharedDaughterCut);
      FillHbtParticleCollection(fSecondParticleCut, currentEvent,collection3,fPerformSharedDaughterCut);

      coll_1_size_passes = (collection1->size() >= fMinSizePart1Collection);
      coll_2_size_passes = (collection2->size() >= fMinSizePart2Collection);
      coll_3_size_passes = kTRUE;

      fEventCut->FillCutMonitor(collection1, collection2);
    }
  if(fFirstParticleCut && fSecondParticleCut==NULL && fThirdParticleCut==NULL)
    {
      FillHbtParticleCollection(fFirstParticleCut, currentEvent,collection1,fPerformSharedDaughterCut);
      FillHbtParticleCollection(fFirstParticleCut, currentEvent,collection2,fPerformSharedDaughterCut);
      FillHbtParticleCollection(fFirstParticleCut, currentEvent,collection3,fPerformSharedDaughterCut);

      coll_1_size_passes = (collection1->size() >= fMinSizePart1Collection);
      coll_2_size_passes = kTRUE;
      coll_3_size_passes = kTRUE;

      fEventCut->FillCutMonitor(collection1, collection1);
    }
  
  // passes only if the sizes are above min
  tmpPassEvent = tmpPassEvent && coll_1_size_passes && coll_2_size_passes && coll_3_size_passes;

  if (!tmpPassEvent){
    EventEnd(currentEvent);
    if(fPicoEvent) delete fPicoEvent;
    return;
  }

  
  // add real events
  if(fFirstParticleCut && fSecondParticleCut && fThirdParticleCut)
    {
      AddParticles(collection1,collection2,collection3);
    }
  if(fFirstParticleCut && fSecondParticleCut && fThirdParticleCut==NULL)
    {
      AddParticles(collection1,collection2);
    }
  if(fFirstParticleCut && fSecondParticleCut==NULL && fThirdParticleCut==NULL)
    {
      AddParticles(collection1);
    }
  
  bool mixing = false;
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

  //we have to swap particles in trio in order to avoid any implicit ordering of particles in collections
  int swpart3= fNeventsProcessed % 3;
  bool swpart2= fNeventsProcessed % 2;
  
  // All three collections (particles) are the same
  if(collection1 && collection2==NULL && collection3==NULL)
    {
      
      // We have three loops - outer, middle, and inner
      AliFemtoParticleConstIterator tStartOuterLoop = collection1->begin(),
	tEndOuterLoop = collection1->end(),
	tStartMiddleLoop,
	tEndMiddleLoop  = collection1->end(),
	tStartInnerLoop,
	tEndInnerLoop = collection1->end();
      
      tEndOuterLoop--; tEndOuterLoop--;   //   Outer loop goes to next-to-next-to-last particle (2nd from the last)
      tEndMiddleLoop--;                            //   Middle loop goes to next-to-last particle (1st from the last) 
      tEndInnerLoop = collection1->end() ;     //   Inner loop goes to last particle

      
      // Begin the outer loop
      for (AliFemtoParticleConstIterator tPartIter1 = tStartOuterLoop;  tPartIter1 != tEndOuterLoop; ++tPartIter1)
	{
	  // start middle loop at the particle
	  // after the current outer loop position, (loops until end)
	  tStartMiddleLoop = tPartIter1;
	  tStartMiddleLoop++;


	  // Begin the middle loop
	  for (AliFemtoParticleConstIterator tPartIter2 = tStartMiddleLoop;  tPartIter2 != tEndMiddleLoop; ++tPartIter2)
	    {
	      // start middle loop at the particle
	      // after the current outer loop position, (loops until end)
	      tStartInnerLoop = tPartIter2;
	      tStartInnerLoop++;
	      
	      // Begin the inner loop
	      for (AliFemtoParticleConstIterator tPartIter3 = tStartInnerLoop;  tPartIter3 != tEndInnerLoop; ++tPartIter3)
		{	  
		  // Swap between all three particles to avoid biased ordering		  
		  if(swpart3 == 0)
		    {
		      trio->SetTrack1(*tPartIter1,fCollection1type);
		      trio->SetTrack2(*tPartIter2,fCollection1type);
		      trio->SetTrack3(*tPartIter3,fCollection1type);
		      swpart3++;
		    }
		  else if(swpart3 == 1)
		    {
		      trio->SetTrack1(*tPartIter3,fCollection1type);
		      trio->SetTrack2(*tPartIter1,fCollection1type);
		      trio->SetTrack3(*tPartIter2,fCollection1type);
		      swpart3++;
		    }
		  else
		    {
		      trio->SetTrack1(*tPartIter2,fCollection1type);
		      trio->SetTrack2(*tPartIter3,fCollection1type);
		      trio->SetTrack3(*tPartIter1,fCollection1type);
		      swpart3 = 0;
		    }

		  
		  for (AliFemtoTrioFctnIterator iFun = fTrioFctnCollection->begin();iFun != fTrioFctnCollection->end();++iFun)
		    {
		      AliFemtoTrioFctn *triofctn = *iFun;
		      
		      if(mixing)  triofctn->AddMixedTrio(trio);
		      else        triofctn->AddRealTrio(trio);
		    }

		} // end of inner loop

	    } // end of middle loop

	} //end of outer loop

    }




  
  // Two collections (particles) are the same, the third one is different
  // collection 1 is different than collection 2 and collection 3
  // collections 2 and 3 are the same
  if(collection1 && collection2 && collection3==NULL)
    {
      for (AliFemtoParticleConstIterator iPart1 = collection1->begin();iPart1 != collection1->end();++iPart1)
	{
	  trio->SetTrack1(*iPart1,fCollection1type);
	  // Setup iterator ranges
	  //
	  // The outer loop alway starts at beginning of particle collection 1.
	  // * If we are iterating over both particle collections, then the loop simply
	  // runs through both from beginning to end.
	  // * If we are only iterating over one particle collection, the inner loop
	  // loops over all particles between the outer iterator and the end of the
	  // collection. The outer loop must skip the last entry of the list.
	  AliFemtoParticleConstIterator tStartOuterLoop = collection2->begin(),
	    tEndOuterLoop = collection2->end(),
	    tStartInnerLoop,
	    tEndInnerLoop;
	  
	  
	  
	  tEndOuterLoop--;                             //   Outer loop goes to next-to-last particle
	  tEndInnerLoop = collection2->end() ;     //   Inner loop goes to last particle
	  

	  // Begin the outer loop
	  for (AliFemtoParticleConstIterator tPartIter2 = tStartOuterLoop;  tPartIter2 != tEndOuterLoop; ++tPartIter2)
	    {

	      // start inner loop at the particle
	      // after the current outer loop position, (loops until end)
	      tStartInnerLoop = tPartIter2;
	      tStartInnerLoop++;
	    


	      // Begin the inner loop
	      for (AliFemtoParticleConstIterator tPartIter3 = tStartInnerLoop; tPartIter3 != tEndInnerLoop; ++tPartIter3)
		{


		  // Swap between first and second particles to avoid biased ordering
		  trio->SetTrack2(swpart2 ? *tPartIter2 : *tPartIter3,fCollection2type);
		  trio->SetTrack3(swpart2 ? *tPartIter3 : *tPartIter2,fCollection2type);
		  swpart2 = !swpart2;

		  for (AliFemtoTrioFctnIterator iFun = fTrioFctnCollection->begin();iFun != fTrioFctnCollection->end();++iFun)
		    {
		      AliFemtoTrioFctn *triofctn = *iFun;
		    
		      if(mixing)  triofctn->AddMixedTrio(trio);
		      else        triofctn->AddRealTrio(trio);
		    }

		  
		}    // loop over second particle
	    }      // loop over first particle
	}
    }

  // The three collections are different (we take three different particles)
  if(collection1 && collection2 && collection3)
    {
      //cout<<"Collections not  the same"<<endl;
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

void AliFemtoTrioAnalysis::SetMinSizePart1Collection(unsigned int minSize)
{ fMinSizePart1Collection = minSize; }
void AliFemtoTrioAnalysis::SetMinSizePart2Collection(unsigned int minSize)
{ fMinSizePart2Collection = minSize; }
void AliFemtoTrioAnalysis::SetMinSizePart3Collection(unsigned int minSize)
{ fMinSizePart3Collection = minSize; }


void AliFemtoTrioAnalysis::SetDoEventMixing(bool mix){fDoEventMixing = mix;}


void AliFemtoTrioAnalysis::SetV0SharedDaughterCut(bool perform){fPerformSharedDaughterCut = perform;}
bool AliFemtoTrioAnalysis::V0SharedDaughterCut(){return fPerformSharedDaughterCut;}


AliFemtoString AliFemtoTrioAnalysis::Report(){return "Report";}
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
  
  if (fSecondParticleCut && fSecondParticleCut != fFirstParticleCut) {
    TList *p2Cut = fSecondParticleCut->GetOutputList();
    TIter nextp2(p2Cut);
    while (TObject *obj = nextp2()) {outputList->Add(obj);}
    if(p2Cut) delete p2Cut;
  }
  
  if (fThirdParticleCut && fThirdParticleCut != fFirstParticleCut) {
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
