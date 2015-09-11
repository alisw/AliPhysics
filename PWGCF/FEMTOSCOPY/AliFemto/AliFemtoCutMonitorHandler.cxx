///
/// \file AliFemtoCutMonitorHandler.cxx
///

#include <TList.h>
#include "AliFemtoCutMonitorHandler.h"
#include "AliFemtoTypes.h"

#include <iterator>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoCutMonitorHandler)
  /// \endcond
#endif

// ---------------------------------------------------------------------------
AliFemtoCutMonitorHandler::AliFemtoCutMonitorHandler():
  fCollectionsEmpty(kTRUE),
  fPassColl(NULL),
  fFailColl(NULL)
{
  // Default constructor
  fPassColl = new AliFemtoCutMonitorCollection();
  fFailColl = new AliFemtoCutMonitorCollection();
}
// ---------------------------------------------------------------------------
AliFemtoCutMonitorHandler::AliFemtoCutMonitorHandler(const AliFemtoCutMonitorHandler& aHan):
  fCollectionsEmpty(aHan.fCollectionsEmpty),
  fPassColl(NULL),
  fFailColl(NULL)
{
  // Copy constructor
  fPassColl = new AliFemtoCutMonitorCollection(aHan.fPassColl->begin(), aHan.fPassColl->end());
  fFailColl = new AliFemtoCutMonitorCollection(aHan.fFailColl->begin(), aHan.fFailColl->end());
}

// ---------------------------------------------------------------------------
AliFemtoCutMonitorHandler::~AliFemtoCutMonitorHandler()
{
  // Default destructor
  delete fPassColl;
  delete fFailColl;
}
//__________________________
AliFemtoCutMonitorHandler& AliFemtoCutMonitorHandler::operator=(const AliFemtoCutMonitorHandler& aHan)
{
  // assignment operator
  if (this == &aHan) {
    return *this;
  }

  if (fPassColl) {
    fPassColl->clear();
    fPassColl->insert(fPassColl->begin(), aHan.fPassColl->begin(), aHan.fPassColl->end());
  } else {
    fPassColl = new AliFemtoCutMonitorCollection(aHan.fPassColl->begin(), aHan.fPassColl->end());
  }

  if (fFailColl) {
    fFailColl->clear();
    fFailColl->insert(fFailColl->begin(), aHan.fFailColl->begin(), aHan.fFailColl->end());
  } else {
    fFailColl = new AliFemtoCutMonitorCollection(aHan.fFailColl->begin(), aHan.fFailColl->end());
  }

  return *this;
}

// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoEvent* event, bool pass)
{
  // fill event cut monitors
  if (fCollectionsEmpty) return;

  AliFemtoCutMonitorCollection *output = pass ? fPassColl : fFailColl;

  AliFemtoCutMonitorIterator iter;
  for (iter = output->begin(); iter != output->end(); ++iter) {
    (*iter)->Fill(event);
  }
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoTrack* track, bool pass)
{
  // Fill track cut monitors
  if (fCollectionsEmpty) return;

  AliFemtoCutMonitorCollection *output = pass ? fPassColl : fFailColl;

  AliFemtoCutMonitorIterator iter;
  for (iter = output->begin(); iter != output->end(); ++iter) {
    (*iter)->Fill(track);
  }
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoV0* v0, bool pass)
{
  // fill V0 cut monitors
  if (fCollectionsEmpty) return;

  AliFemtoCutMonitorCollection *output = pass ? fPassColl : fFailColl;

  AliFemtoCutMonitorIterator iter;
  for (iter = output->begin(); iter != output->end(); ++iter) {
    (*iter)->Fill(v0);
  }
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoXi* xi, bool pass)
{
  // fill V0 cut monitors
  if (fCollectionsEmpty) return;

  AliFemtoCutMonitorCollection *output = pass ? fPassColl : fFailColl;

  AliFemtoCutMonitorIterator iter;
  for (iter = output->begin(); iter != output->end(); ++iter) {
    (*iter)->Fill(xi);
  }
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoKink* kink, bool pass)
{
  // fill kink cut monitors
  if (fCollectionsEmpty) return;

  AliFemtoCutMonitorCollection *output = pass ? fPassColl : fFailColl;

  AliFemtoCutMonitorIterator iter;
  for (iter = output->begin(); iter != output->end(); ++iter) {
    (*iter)->Fill(kink);
  }
}
// ---------------------------------Gael/12/04/02-----------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoPair* pair, bool pass)
{
  // fill pair cut monitors
  if (fCollectionsEmpty) return;

  AliFemtoCutMonitorCollection *output = pass ? fPassColl : fFailColl;

  AliFemtoCutMonitorIterator iter;
  for (iter = output->begin(); iter != output->end(); ++iter) {
    (*iter)->Fill(pair);
  }
}
// ---------------------------------Gael/19/06/02-----------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoParticleCollection* partColl)
{
  // fill particle collection cut monitor
  if (fCollectionsEmpty) return;

  AliFemtoCutMonitorIterator iter;
  for (iter = fPassColl->begin(); iter != fPassColl->end(); ++iter) {
    (*iter)->Fill(partColl);
  }
}
// ------------------------------------Gael/19/06/02-------------------------
void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoEvent* event,
                                               const AliFemtoParticleCollection* partColl)
{
  // Fill event particle collection
  if (fCollectionsEmpty) return;

  AliFemtoCutMonitorIterator iter;
  for (iter = fPassColl->begin(); iter != fPassColl->end(); ++iter){
    (*iter)->Fill(event, partColl);
  }
}

void AliFemtoCutMonitorHandler::FillCutMonitor(const AliFemtoParticleCollection* partColl1,
                                               const AliFemtoParticleCollection* partColl2)
{
  // Fill event particle collection
  if (fCollectionsEmpty) return;

  AliFemtoCutMonitorIterator iter;
  for (iter = fPassColl->begin(); iter != fPassColl->end(); ++iter) {
    (*iter)->Fill(partColl1, partColl2);
  }
}

// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::Finish()
{
  // Perform finish operations on cut monitors
  AliFemtoCutMonitorIterator iter;
  for (iter = fPassColl->begin(); iter != fPassColl->end(); ++iter) {
    (*iter)->Finish();
  }
  for (iter = fFailColl->begin(); iter != fFailColl->end(); ++iter) {
    (*iter)->Finish();
  }
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::AddCutMonitor(AliFemtoCutMonitor* cutMoni1,
                                              AliFemtoCutMonitor* cutMoni2)
{
  // Add cut monitors to collections
  fPassColl->push_back(cutMoni1);
  fFailColl->push_back(cutMoni2);
  fCollectionsEmpty = false;
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::AddCutMonitor(AliFemtoCutMonitor* cutMoni)
{
  // make a copy of the cut monitor
  cout << " make a copy of the cutmonitor and push both into the collections\n"
          " not yet implemented" << endl;
  fPassColl->push_back(cutMoni);
  cout << " only pass collection pushed" << endl;
  fCollectionsEmpty = false;
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::AddCutMonitorPass(AliFemtoCutMonitor* cutMoni)
{
  // add monitors to pass
  fPassColl->push_back(cutMoni);
  fCollectionsEmpty = false;
}
// ---------------------------------------------------------------------------
void AliFemtoCutMonitorHandler::AddCutMonitorFail(AliFemtoCutMonitor* cutMoni)
{
  // add monitors to fail
  fFailColl->push_back(cutMoni);
  fCollectionsEmpty = false;
}
// ---------------------------------------------------------------------------
AliFemtoCutMonitor* AliFemtoCutMonitorHandler::PassMonitor(int n)
{
  // return pass monitor number n
  if (static_cast<int>(fPassColl->size()) <= n) {
    return NULL;
  }
  AliFemtoCutMonitorIterator iter = fPassColl->begin();
  std::advance(iter, n);
  return *iter;
}
// ---------------------------------------------------------------------------
AliFemtoCutMonitor* AliFemtoCutMonitorHandler::FailMonitor(int n)
{
  // return fail monitor number n
  if (static_cast<int>(fFailColl->size()) <= n) {
    return NULL;
  }
  AliFemtoCutMonitorIterator iter = fFailColl->begin();
  std::advance(iter, n);
  return *iter;
}
//_____________________________________________________________________________
TList *AliFemtoCutMonitorHandler::GetOutputList()
{
  TList *tOutputList = new TList();

  for (unsigned int ipass = 0; ipass < fPassColl->size(); ipass++) {
    TList *tLp = PassMonitor(ipass)->GetOutputList();

    TIter nextLp(tLp);
    while (TObject *obj = nextLp()) {
      tOutputList->Add(obj);
    }

    delete tLp;
  }

  for (unsigned int ipass = 0; ipass < fFailColl->size(); ipass++) {
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

  for (iter = fPassColl->begin(); iter != fPassColl->end(); ++iter) {
    (*iter)->EventBegin(aEvent);
  }

  for (iter = fFailColl->begin(); iter != fFailColl->end(); ++iter) {
    (*iter)->EventBegin(aEvent);
  }
}
//_____________________________________________________________________________
void AliFemtoCutMonitorHandler::EventEnd(const AliFemtoEvent* aEvent)
{
  if (fCollectionsEmpty) return;

  AliFemtoCutMonitorIterator iter;

  for (iter = fPassColl->begin(); iter != fPassColl->end(); ++iter) {
    (*iter)->EventEnd(aEvent);
  }

  for (iter = fFailColl->begin(); iter != fFailColl->end(); ++iter) {
    (*iter)->EventEnd(aEvent);
  }
}
