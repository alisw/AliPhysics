
#include "Infrastructure/AliFemtoCutMonitorHandler.h"
#include "Infrastructure/AliFemtoTypes.h"

#ifdef __ROOT__
ClassImp(AliFemtoCutMonitorHandler)
#endif
// ---------------------------------------------------------------------------
AliFemtoCutMonitorHandler::AliFemtoCutMonitorHandler() {
  cout << " *** AliFemtoCutMonitorHandler::AliFemtoCutMonitorHandler() " << endl;
  fCollectionsEmpty = 0;
  fPassColl = new AliFemtoCutMonitorCollection();
  fFailColl = new AliFemtoCutMonitorCollection();
}
// ---------------------------------------------------------------------------
AliFemtoCutMonitorHandler::~AliFemtoCutMonitorHandler() { 
  delete fPassColl;
  delete fFailColl;
}   
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoEvent* event, bool pass) { 
  if (fCollectionsEmpty) return;
  AliFemtoCutMonitorIterator iter;
  AliFemtoCutMonitor* CM;
  if ( pass) {
    for (iter=fPassColl->begin(); iter!=fPassColl->end(); iter++){
      CM = *iter;
      CM->Fill(event);
    }
  } else {
    for (iter=fFailColl->begin(); iter!=fFailColl->end(); iter++){
      CM = *iter;
      CM->Fill(event);
    }
  }
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoTrack* track, bool pass) { 
  if (fCollectionsEmpty) return;
  AliFemtoCutMonitorIterator iter;
  AliFemtoCutMonitor* CM;
  if ( pass) {
    for (iter=fPassColl->begin(); iter!=fPassColl->end(); iter++){
      CM = *iter;
      CM->Fill(track);
    }
  } else {
    for (iter=fFailColl->begin(); iter!=fFailColl->end(); iter++){
      CM = *iter;
      CM->Fill(track);
    }
  }
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoV0* v0, bool pass) { 
  if (fCollectionsEmpty) return;
  AliFemtoCutMonitorIterator iter;
  AliFemtoCutMonitor* CM;
  if ( pass) {
    for (iter=fPassColl->begin(); iter!=fPassColl->end(); iter++){
      CM = *iter;
      CM->Fill(v0);
    }
  } else {
    for (iter=fFailColl->begin(); iter!=fFailColl->end(); iter++){
      CM = *iter;
      CM->Fill(v0);
    }
  }
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoKink* kink, bool pass) { 
  if (fCollectionsEmpty) return;
  AliFemtoCutMonitorIterator iter;
  AliFemtoCutMonitor* CM;
  if ( pass) {
    for (iter=fPassColl->begin(); iter!=fPassColl->end(); iter++){
      CM = *iter;
      CM->Fill(kink);
    }
  } else {
    for (iter=fFailColl->begin(); iter!=fFailColl->end(); iter++){
      CM = *iter;
      CM->Fill(kink);
    }
  }
}
// ---------------------------------Gael/12/04/02-----------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoPair* pair, bool pass) { 
  if (fCollectionsEmpty) return;
  AliFemtoCutMonitorIterator iter;
  AliFemtoCutMonitor* CM;
  if ( pass) {
    for (iter=fPassColl->begin(); iter!=fPassColl->end(); iter++){
      CM = *iter;
      CM->Fill(pair);
    }
  } else {
    for (iter=fFailColl->begin(); iter!=fFailColl->end(); iter++){
      CM = *iter;
      CM->Fill(pair);
    }
  }
}
// ---------------------------------Gael/19/06/02-----------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoParticleCollection* partColl) {
  if (fCollectionsEmpty) return;
  AliFemtoCutMonitorIterator iter;
  AliFemtoCutMonitor* CM;
  
  for (iter=fPassColl->begin(); iter!=fPassColl->end(); iter++){
    CM = *iter;
    CM->Fill(partColl);
  }
}
// ------------------------------------Gael/19/06/02-------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoEvent* event,const AliFemtoParticleCollection* partColl) {
  
  cout<<"In AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoEvent* event, AliFemtoPicoEvent* picoEvent)"<<endl;
  if (fCollectionsEmpty) return;
  AliFemtoCutMonitorIterator iter;
  AliFemtoCutMonitor* CM;
  
  for (iter=fPassColl->begin(); iter!=fPassColl->end(); iter++){
    CM = *iter;
    CM->Fill(event,partColl);
  }
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::Finish() { 
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
  fPassColl->push_back(cutMoni1);
  fFailColl->push_back(cutMoni2);
  fCollectionsEmpty=false;
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::AddCutMonitor(AliFemtoCutMonitor* cutMoni) { 
  cout << " make a copy of the cutmonitor and push both into the collections " << endl;
  cout << " not yet implemented" << endl;
  fPassColl->push_back(cutMoni);
  cout << " only pass collection pushed" << endl;
  fCollectionsEmpty=false;
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::AddCutMonitorPass(AliFemtoCutMonitor* cutMoni) { 
  fPassColl->push_back(cutMoni);
  fCollectionsEmpty=false;
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::AddCutMonitorFail(AliFemtoCutMonitor* cutMoni) { 
  fFailColl->push_back(cutMoni);
  fCollectionsEmpty=false;
}
// ---------------------------------------------------------------------------
AliFemtoCutMonitor* AliFemtoCutMonitorHandler::PassMonitor(int n) { 
  AliFemtoCutMonitorIterator iter = fPassColl->begin();
  if ( (int)fPassColl->size() <= n ) return NULL;
  for ( int i=0; i<n; i++)
    iter++;
  return *iter;
}
// ---------------------------------------------------------------------------
AliFemtoCutMonitor* AliFemtoCutMonitorHandler::FailMonitor(int n) { 
  AliFemtoCutMonitorIterator iter = fFailColl->begin();
  if ( (int)fFailColl->size() <= n ) return NULL;
  for ( int i=0; i<n; i++)
    iter++;
  return *iter;
}



 
 
