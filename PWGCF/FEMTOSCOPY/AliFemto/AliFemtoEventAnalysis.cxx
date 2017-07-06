///
/// \file AliFemtoEventAnalysis.cxx
///

#include "AliFemtoEventAnalysis.h"
#include "AliFemtoTrackCut.h"
#include "AliFemtoV0Cut.h"
#include "AliFemtoPicoEvent.h"

#include <string>
#include <iostream>
#include <iterator>

#ifdef __ROOT__
/// \cond CLASSIMP
ClassImp(AliFemtoEventAnalysis);
/// \endcond
#endif

AliFemtoEventCut*    copyTheCut(AliFemtoEventCut*);
AliFemtoParticleCut* copyTheCut(AliFemtoParticleCut*);
AliFemtoCorrFctn*    copyTheCorrFctn(AliFemtoCorrFctn*);

extern void FillHbtParticleCollection(AliFemtoParticleCut* partCut,
                                      AliFemtoEvent* hbtEvent,
                                      AliFemtoParticleCollection* partCollection,
                                      bool performSharedDaughterCut=kFALSE);

template <class TrackCollectionType, class TrackCutType>
void DoFillParticleCollection(TrackCutType *cut,
                              TrackCollectionType *track_collection,
                              AliFemtoParticleCollection *output)
{
  // lets's just name the iterator type
  typedef typename TrackCollectionType::iterator TrackCollectionIterType;
  
  for (TrackCollectionIterType pIter = track_collection->begin();
       pIter != track_collection->end();
       pIter++) {
    const Bool_t track_passes = cut->Pass(*pIter);
    cut->FillCutMonitor(*pIter, track_passes);
    if (track_passes) {
      output->push_back(new AliFemtoParticle(*pIter, cut->Mass()));
    }
  }
}

AliFemtoEventAnalysis::AliFemtoEventAnalysis():
fCorrFctnCollection(NULL),
fEventCut(NULL),
fFirstParticleCut(NULL),
fSecondParticleCut(NULL),
fNeventsProcessed(0),
fPerformSharedDaughterCut(kFALSE)
{
  // Default constructor
  fCorrFctnCollection = new AliFemtoCorrFctnCollection;
}
//____________________________
AliFemtoEventAnalysis::AliFemtoEventAnalysis(const AliFemtoEventAnalysis& a):
AliFemtoAnalysis(),
fCorrFctnCollection(NULL),
fEventCut(NULL),
fFirstParticleCut(NULL),
fSecondParticleCut(NULL),
fNeventsProcessed(0),
fPerformSharedDaughterCut(a.fPerformSharedDaughterCut)
{
  /// Copy constructor
  
  const char msg_template[] = " AliFemtoEventAnalysis::AliFemtoEventAnalysis(const AliFemtoEventAnalysis& a) - %s",
  warn_template[] = " WARNING [AliFemtoEventAnalysis::AliFemtoEventAnalysis(const AliFemtoEventAnalysis& a)] %s";
  
  fCorrFctnCollection = new AliFemtoCorrFctnCollection;
  
  // Clone the event cut
  AliFemtoEventCut *ev_cut = a.fEventCut->Clone();
  if (ev_cut) {
    SetEventCut(ev_cut);
    cout << " event cut set " << endl;
    
  } else {
    cerr << TString::Format(warn_template, "Could not clone event cut") << endl;
    // TODO: handle uncloned event cut
  }
  
  fFirstParticleCut = a.fFirstParticleCut->Clone();
  if (fFirstParticleCut) {
    SetFirstParticleCut(fFirstParticleCut);
    cout << TString::Format(msg_template, "first particle cut set") << endl;
  }
  else {
    cerr << TString::Format(warn_template, "Could not clone first particle cut") << endl;
    // TODO: handle uncloned track cut
  }
  
  
  fSecondParticleCut = (a.fFirstParticleCut == a.fSecondParticleCut)
  ? fFirstParticleCut
  : a.fSecondParticleCut->Clone();
  if (fSecondParticleCut) {
    SetSecondParticleCut(fSecondParticleCut);
    cout << TString::Format(msg_template, "second particle cut set") << endl;
    
  }
  else{
    cerr << TString::Format(warn_template, "Could not clone second particle cut") << endl;
    // TODO: handle uncloned track cut
  }
  
  AliFemtoCorrFctnIterator iter;
  
  cout << TString::Format(msg_template, "looking for correlation functions") << endl;
  
  for (iter = a.fCorrFctnCollection->begin(); iter != a.fCorrFctnCollection->end(); ++iter) {
    AliFemtoCorrFctn* fctn = (*iter)->Clone();
    if (fctn) {
      AddCorrFctn(fctn);
    } else {
      cout << TString::Format(msg_template, "correlation function not found") << endl;
    }
  }
  cout << TString::Format(msg_template, "analysis copied") << endl;
}
//____________________________
AliFemtoEventAnalysis::~AliFemtoEventAnalysis()
{
  cout << " AliFemtoEventAnalysis::~AliFemtoEventAnalysis()" << endl;
  
  if (fFirstParticleCut == fSecondParticleCut) fSecondParticleCut = NULL;
  
  delete fEventCut;
  delete fFirstParticleCut;
  delete fSecondParticleCut;

  if (fCorrFctnCollection) {
    for (AliFemtoCorrFctnIterator iter = fCorrFctnCollection->begin(); iter != fCorrFctnCollection->end(); iter++) {
      delete *iter;
    }
    delete fCorrFctnCollection;
  }
}
//______________________
AliFemtoEventAnalysis& AliFemtoEventAnalysis::operator=(const AliFemtoEventAnalysis& aAna)
{
  if (this == &aAna) return *this;
  if (fFirstParticleCut == fSecondParticleCut) fSecondParticleCut = NULL;
  
  delete fEventCut;
  delete fFirstParticleCut;
  delete fSecondParticleCut;
  
  if (fCorrFctnCollection) {
    for (AliFemtoCorrFctnIterator iter = fCorrFctnCollection->begin(); iter != fCorrFctnCollection->end(); ++iter) {
      delete *iter;
    }
    fCorrFctnCollection->clear();
  } else {
    cerr << " WARNING [AliFemtoEventAnalysis::operator=()] fCorrFctnCollection was NULL, this should not happen." << endl;
    fCorrFctnCollection = new AliFemtoCorrFctnCollection;
  }
  
  fEventCut = aAna.fEventCut->Clone();
  fFirstParticleCut = aAna.fFirstParticleCut->Clone();
  fSecondParticleCut = (aAna.fFirstParticleCut == aAna.fSecondParticleCut)
  ? fFirstParticleCut
  : aAna.fSecondParticleCut->Clone();
  
  if (fEventCut) {
    SetEventCut(fEventCut);
  } else {
    cerr << " WARNING [AliFemtoEventAnalysis::operator=()] Could not clone event cut." << endl;
  }
  
  if (fFirstParticleCut) {
    SetFirstParticleCut(fFirstParticleCut);
  } else {
    cerr << " WARNING [AliFemtoEventAnalysis::operator=()] Could not clone first particle cut." << endl;
  }
  
  if (fSecondParticleCut) {
    SetSecondParticleCut(fSecondParticleCut);
  } else {
    cerr << " WARNING [AliFemtoEventAnalysis::operator=()] Could not clone second particle cut." << endl;
  }
  
  for (AliFemtoCorrFctnIterator iter = aAna.fCorrFctnCollection->begin(); iter != aAna.fCorrFctnCollection->end(); ++iter) {
    AliFemtoCorrFctn* fctn = (*iter)->Clone();
    if (fctn) AddCorrFctn(fctn);
  }
  
  fPerformSharedDaughterCut = aAna.fPerformSharedDaughterCut;
  
  return *this;
}

AliFemtoCorrFctn* AliFemtoEventAnalysis::CorrFctn(int n)
{
  if (n < 0 || n > (int)fCorrFctnCollection->size()) return NULL;
  
  AliFemtoCorrFctnIterator iter = fCorrFctnCollection->begin();
  std::advance(iter, n);
  return *iter;
}

AliFemtoString AliFemtoEventAnalysis::Report()
{
  /// Create a simple report from the analysis execution
  
  cout << "AliFemtoEventAnalysis - constructing Report..."<<endl;
  string temp = "-----------\nHbt Analysis Report:\n";
  temp += "\nEvent Cuts:\n";
  temp += fEventCut->Report();
  temp += "\nParticle Cuts - First Particle:\n";
  temp += fFirstParticleCut->Report();
  temp += "\nParticle Cuts - Second Particle:\n";
  temp += fSecondParticleCut->Report();
  temp += "\nCorrelation Functions:\n";
  
  if (fCorrFctnCollection->empty()) {
    cout << "AliFemtoEventAnalysis-Warning : no correlations functions in this analysis " << endl;
  }
  for (AliFemtoCorrFctnIterator iter = fCorrFctnCollection->begin(); iter != fCorrFctnCollection->end(); ++iter) {
    temp += (*iter)->Report();
    temp += "\n";
  }
  temp += "-------------\n";
  AliFemtoString returnThis=temp;
  return returnThis;
}

void AliFemtoEventAnalysis::ProcessEvent(const AliFemtoEvent* hbtEvent)
{
  fPicoEvent = NULL;
  fNeventsProcessed++;

  EventBegin(hbtEvent);
  bool tmpPassEvent = fEventCut->Pass(hbtEvent);
  
  if (!tmpPassEvent) {
    fEventCut->FillCutMonitor(hbtEvent, tmpPassEvent);
    EventEnd(hbtEvent);  // cleanup for EbyE
    return;
  }
  fPicoEvent = new AliFemtoPicoEvent;
  
  AliFemtoParticleCollection  *collection1 = fPicoEvent->FirstParticleCollection(),
                              *collection2 = fPicoEvent->SecondParticleCollection();
  
  if (collection1 == NULL || collection2 == NULL) {
    cout << "E-AliFemtoEventAnalysis::ProcessEvent: new PicoEvent is missing particle collections!\n";
    EventEnd(hbtEvent);  // cleanup for EbyE
    delete fPicoEvent;
    return;
  }
  
  FillHbtParticleCollection(fFirstParticleCut,
                            (AliFemtoEvent*)hbtEvent,
                            fPicoEvent->FirstParticleCollection(),
                            fPerformSharedDaughterCut);
  
  FillHbtParticleCollection(fSecondParticleCut,
                            (AliFemtoEvent*)hbtEvent,
                            fPicoEvent->SecondParticleCollection(),
                            fPerformSharedDaughterCut);
  
  fEventCut->FillCutMonitor(collection1, collection2); //MJ!
  fEventCut->FillCutMonitor(hbtEvent, tmpPassEvent);
  
  if (!tmpPassEvent) {
    EventEnd(hbtEvent);
    delete fPicoEvent;
    return;
  }
  
  AddParticles("first", collection1);
  AddParticles("second", collection2);
  
  
  cout << "AliFemtoEventAnalysis::ProcessEvent() - reals done ";
  
  EventEnd(hbtEvent);
}

//_________________________
void AliFemtoEventAnalysis::AddParticles(const char* typeIn, AliFemtoParticleCollection *partCollection)
{
  const string type = typeIn;
  
  for (AliFemtoCorrFctnIterator tCorrFctnIter = fCorrFctnCollection->begin();
       tCorrFctnIter != fCorrFctnCollection->end();
       ++tCorrFctnIter)
  {
    AliFemtoCorrFctn *tCorrFctn = (*tCorrFctnIter);
    
    for (AliFemtoParticleConstIterator  tPartIter = partCollection->begin();
         tPartIter != partCollection->end();
         ++tPartIter)
    {
      AliFemtoParticle *particle = *tPartIter;
      if (type == "first")      tCorrFctn->AddFirstParticle(particle);
      else if(type == "second") tCorrFctn->AddSecondParticle(particle);
    }
    
    tCorrFctn->CalculateAnglesForEvent();
  }
}

void AliFemtoEventAnalysis::EventBegin(const AliFemtoEvent* ev)
{
  fFirstParticleCut->EventBegin(ev);
  fSecondParticleCut->EventBegin(ev);
  
  for (AliFemtoCorrFctnIterator iter = fCorrFctnCollection->begin();
       iter != fCorrFctnCollection->end();
       ++iter){
    (*iter)->EventBegin(ev);
  }
}

void AliFemtoEventAnalysis::EventEnd(const AliFemtoEvent* ev)
{
  fFirstParticleCut->EventEnd(ev);
  fSecondParticleCut->EventEnd(ev);
  
  for (AliFemtoCorrFctnIterator iter = fCorrFctnCollection->begin();
       iter != fCorrFctnCollection->end();
       ++iter){
    (*iter)->EventEnd(ev);
  }
}

void AliFemtoEventAnalysis::Finish()
{
  for (AliFemtoCorrFctnIterator iter = fCorrFctnCollection->begin();
       iter != fCorrFctnCollection->end();
       ++iter){
    (*iter)->Finish();
  }
}

TList* AliFemtoEventAnalysis::GetOutputList()
{
  TList *tOutputList = new TList();
  
  TList *p1Cut = fFirstParticleCut->GetOutputList();
  
  TListIter nextp1(p1Cut);
  while (TObject *obj = nextp1.Next()) {
    tOutputList->Add(obj);
  }
  delete p1Cut;
  
  if (fSecondParticleCut != fFirstParticleCut) {
    TList *p2Cut = fSecondParticleCut->GetOutputList();
    
    TIter nextp2(p2Cut);
    while (TObject *obj = nextp2()) {
      tOutputList->Add(obj);
    }
    delete p2Cut;
  }
  
  TList *eventCut = fEventCut->GetOutputList();
  
  TIter nextevent(eventCut);
  while (TObject *obj = nextevent()) {
    tOutputList->Add(obj);
  }
  delete eventCut;
  
  for (AliFemtoCorrFctnIterator iter = fCorrFctnCollection->begin();
       iter != fCorrFctnCollection->end();
       ++iter) {
    
    TList *tListCf = (*iter)->GetOutputList();
    
    TIter nextListCf(tListCf);
    while (TObject *obj = nextListCf()) {
      tOutputList->Add(obj);
    }
    delete tListCf;
  }
  
  return tOutputList;
}
