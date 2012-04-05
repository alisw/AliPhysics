///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoCutMonitorHandler: a handler for cut monitors                 //
// You add cut monitors to the collection which are stored in two        //
// separate collections - one which stores characteristics of the        //
// entities (tracks, particles, pairs, events) that pass the respective  //
// cuts and the other for the ones that fail the cut.                    //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TList.h>
#include "AliFemtoCutMonitorHandler.h"
#include "AliFemtoTypes.h"

#ifdef __ROOT__
ClassImp(AliFemtoCutMonitorHandler)
#endif
// ---------------------------------------------------------------------------
AliFemtoCutMonitorHandler::AliFemtoCutMonitorHandler():
  fCollectionsEmpty(0), fPassColl(0), fFailColl(0)
{
  // Default constructor
  cout << " *** AliFemtoCutMonitorHandler::AliFemtoCutMonitorHandler() " << endl;
  fCollectionsEmpty = 0;
  fPassColl = new AliFemtoCutMonitorCollection();
  fFailColl = new AliFemtoCutMonitorCollection();
}
// ---------------------------------------------------------------------------
AliFemtoCutMonitorHandler::AliFemtoCutMonitorHandler(const AliFemtoCutMonitorHandler& aHan):
  fCollectionsEmpty(0), fPassColl(0), fFailColl(0)
{
  // Copy constructor
  fCollectionsEmpty = aHan.fCollectionsEmpty;
  fPassColl = new AliFemtoCutMonitorCollection();
  AliFemtoCutMonitorIterator iter;
  for (iter=aHan.fPassColl->begin(); iter!=aHan.fPassColl->end(); iter++){
    fPassColl->push_back(*iter);
  }
  fFailColl = new AliFemtoCutMonitorCollection();
  for (iter=aHan.fFailColl->begin(); iter!=aHan.fFailColl->end(); iter++){
    fFailColl->push_back(*iter);
  }
}

// ---------------------------------------------------------------------------
AliFemtoCutMonitorHandler::~AliFemtoCutMonitorHandler() { 
  // Default destructor
  delete fPassColl;
  delete fFailColl;
}   
//__________________________
AliFemtoCutMonitorHandler& AliFemtoCutMonitorHandler::operator=(const AliFemtoCutMonitorHandler& aHan)
{
  // assignment operator
  if (this == &aHan)
    return *this;

  AliFemtoCutMonitorIterator iter;
  if (fPassColl) {
    fPassColl->clear();
    delete fPassColl;
  }
  if (fFailColl) {
    fFailColl->clear();
    delete fFailColl;
  }
  fPassColl = new AliFemtoCutMonitorCollection();
  for (iter=aHan.fPassColl->begin(); iter!=aHan.fPassColl->end(); iter++){
    fPassColl->push_back(*iter);
  }
  fFailColl = new AliFemtoCutMonitorCollection();
  for (iter=aHan.fFailColl->begin(); iter!=aHan.fFailColl->end(); iter++){
    fFailColl->push_back(*iter);
  }
  return *this;
}

// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoEvent* event, bool pass) { 
  // fill event cut monitors
  if (fCollectionsEmpty) return;
  AliFemtoCutMonitorIterator iter;
  AliFemtoCutMonitor* tCM;
  if ( pass) {
    for (iter=fPassColl->begin(); iter!=fPassColl->end(); iter++){
      tCM = *iter;
      tCM->Fill(event);
    }
  } else {
    for (iter=fFailColl->begin(); iter!=fFailColl->end(); iter++){
      tCM = *iter;
      tCM->Fill(event);
    }
  }
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoTrack* track, bool pass) { 
  // Fill track cut monitors
  if (fCollectionsEmpty) return;
  AliFemtoCutMonitorIterator iter;
  AliFemtoCutMonitor* tCM;
  if ( pass) {
    for (iter=fPassColl->begin(); iter!=fPassColl->end(); iter++){
      tCM = *iter;
      tCM->Fill(track);
    }
  } else {
    for (iter=fFailColl->begin(); iter!=fFailColl->end(); iter++){
      tCM = *iter;
      tCM->Fill(track);
    }
  }
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoV0* v0, bool pass) { 
  // fill V0 cut monitors
  if (fCollectionsEmpty) return;
  AliFemtoCutMonitorIterator iter;
  AliFemtoCutMonitor* tCM;
  if ( pass) {
    for (iter=fPassColl->begin(); iter!=fPassColl->end(); iter++){
      tCM = *iter;
      tCM->Fill(v0);
    }
  } else {
    for (iter=fFailColl->begin(); iter!=fFailColl->end(); iter++){
      tCM = *iter;
      tCM->Fill(v0);
    }
  }
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoKink* kink, bool pass) { 
  // fill kink cut monitors
  if (fCollectionsEmpty) return;
  AliFemtoCutMonitorIterator iter;
  AliFemtoCutMonitor* tCM;
  if ( pass) {
    for (iter=fPassColl->begin(); iter!=fPassColl->end(); iter++){
      tCM = *iter;
      tCM->Fill(kink);
    }
  } else {
    for (iter=fFailColl->begin(); iter!=fFailColl->end(); iter++){
      tCM = *iter;
      tCM->Fill(kink);
    }
  }
}
// ---------------------------------Gael/12/04/02-----------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoPair* pair, bool pass) { 
  // fill pair cut monitors
  if (fCollectionsEmpty) return;
  AliFemtoCutMonitorIterator iter;
  AliFemtoCutMonitor* tCM;
  if ( pass) {
    for (iter=fPassColl->begin(); iter!=fPassColl->end(); iter++){
      tCM = *iter;
      tCM->Fill(pair);
    }
  } else {
    for (iter=fFailColl->begin(); iter!=fFailColl->end(); iter++){
      tCM = *iter;
      tCM->Fill(pair);
    }
  }
}
// ---------------------------------Gael/19/06/02-----------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoParticleCollection* partColl) {
  // fill particle collection cut monitor 
  if (fCollectionsEmpty) return;
  AliFemtoCutMonitorIterator iter;
  AliFemtoCutMonitor* tCM;
  
  for (iter=fPassColl->begin(); iter!=fPassColl->end(); iter++){
    tCM = *iter;
    tCM->Fill(partColl);
  }
}
// ------------------------------------Gael/19/06/02-------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoEvent* event,const AliFemtoParticleCollection* partColl) {
  // Fill event particle collection
  //cout<<"In AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoEvent* event, AliFemtoPicoEvent* picoEvent)"<<endl;
  if (fCollectionsEmpty) return;
  AliFemtoCutMonitorIterator iter;
  AliFemtoCutMonitor* tCM;
  
  for (iter=fPassColl->begin(); iter!=fPassColl->end(); iter++){
    tCM = *iter;
    tCM->Fill(event,partColl);
  }
}

void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoParticleCollection* partColl1, const AliFemtoParticleCollection* partColl2) {
  // Fill event particle collection
  //cout<<"***In AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoEvent* event, AliFemtoPicoEvent* picoEvent)"<<endl;
  if (fCollectionsEmpty) return;
  AliFemtoCutMonitorIterator iter;
  AliFemtoCutMonitor* tCM;
  
  for (iter=fPassColl->begin(); iter!=fPassColl->end(); iter++){
    tCM = *iter;
    tCM->Fill(partColl1,partColl2);
  }
}

// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::Finish() { 
  // Perform finish operations on cut monitors
  AliFemtoCutMonitorIterator iter;
  for (iter=fPassColl->begin(); iter!=fPassColl->end(); iter++){
    (*iter)->Finish();
  }
  for (iter=fFailColl->begin(); iter!=fFailColl->end(); iter++){
    (*iter)->Finish();
  }
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::AddCutMonitor(AliFemtoCutMonitor* cutMoni1, AliFemtoCutMonitor* cutMoni2) { 
  // Add cut monitors to collections
  fPassColl->push_back(cutMoni1);
  fFailColl->push_back(cutMoni2);
  fCollectionsEmpty=false;
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::AddCutMonitor(AliFemtoCutMonitor* cutMoni) { 
  // make a copy of the cut monitor
  cout << " make a copy of the cutmonitor and push both into the collections " << endl;
  cout << " not yet implemented" << endl;
  fPassColl->push_back(cutMoni);
  cout << " only pass collection pushed" << endl;
  fCollectionsEmpty=false;
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::AddCutMonitorPass(AliFemtoCutMonitor* cutMoni) { 
  // add monitors to pass
  fPassColl->push_back(cutMoni);
  fCollectionsEmpty=false;
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::AddCutMonitorFail(AliFemtoCutMonitor* cutMoni) { 
  // add monitors to fail
  fFailColl->push_back(cutMoni);
  fCollectionsEmpty=false;
}
// ---------------------------------------------------------------------------
AliFemtoCutMonitor* AliFemtoCutMonitorHandler::PassMonitor(int n) { 
  // return pass monitor number n
  AliFemtoCutMonitorIterator iter = fPassColl->begin();
  if ( (int)fPassColl->size() <= n ) return NULL;
  for ( int i=0; i<n; i++)
    iter++;
  return *iter;
}
// ---------------------------------------------------------------------------
AliFemtoCutMonitor* AliFemtoCutMonitorHandler::FailMonitor(int n) { 
  // return fail monitor number n
  AliFemtoCutMonitorIterator iter = fFailColl->begin();
  if ( (int)fFailColl->size() <= n ) return NULL;
  for ( int i=0; i<n; i++)
    iter++;
  return *iter;
}
//_____________________________________________________________________________
TList *AliFemtoCutMonitorHandler::GetOutputList()
{
  TList *tOutputList = new TList();

  for (unsigned int ipass=0; ipass<fPassColl->size(); ipass++) {
    TList *tLp = PassMonitor(ipass)->GetOutputList();

    TIter nextLp(tLp);
    while (TObject *obj = nextLp()) {
      tOutputList->Add(obj);
    }
    
    delete tLp;
  }

  for (unsigned int ipass=0; ipass<fFailColl->size(); ipass++) {
    TList *tLf = FailMonitor(ipass)->GetOutputList();

    TIter nextLf(tLf);
    while (TObject *obj = nextLf()) {
      tOutputList->Add(obj);
    }
    
    delete tLf;
  }

  return tOutputList;
}
//_____________________________________________________________________________
void AliFemtoCutMonitorHandler::EventBegin(const AliFemtoEvent* aEvent) 
{ 
  if (fCollectionsEmpty) return;

  AliFemtoCutMonitorIterator iter;
  AliFemtoCutMonitor* tCM;

  for (iter=fPassColl->begin(); iter!=fPassColl->end(); iter++){
    tCM = *iter;
    tCM->EventBegin(aEvent);
  }
  
  for (iter=fFailColl->begin(); iter!=fFailColl->end(); iter++){
    tCM = *iter;
    tCM->EventBegin(aEvent);
  }
}
//_____________________________________________________________________________
void AliFemtoCutMonitorHandler::EventEnd(const AliFemtoEvent* aEvent) 
{ 
  if (fCollectionsEmpty) return;

  AliFemtoCutMonitorIterator iter;
  AliFemtoCutMonitor* tCM;

  for (iter=fPassColl->begin(); iter!=fPassColl->end(); iter++){
    tCM = *iter;
    tCM->EventEnd(aEvent);
  }
  
  for (iter=fFailColl->begin(); iter!=fFailColl->end(); iter++){
    tCM = *iter;
    tCM->EventEnd(aEvent);
  }
}


 
 
