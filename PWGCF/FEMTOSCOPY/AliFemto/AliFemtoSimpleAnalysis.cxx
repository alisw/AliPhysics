///
/// \file AliFemtoSimpleAnalysis.cxx
///

#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoTrackCut.h"
#include "AliFemtoV0Cut.h"
#include "AliFemtoKinkCut.h"
#include "AliFemtoXiCut.h"
#include "AliFemtoXiTrackCut.h"
#include "AliFemtoPicoEvent.h"

#include <string>
#include <iostream>
#include <iterator>

#ifdef __ROOT__
  /// \cond CLASSIMP
  ClassImp(AliFemtoSimpleAnalysis);
  /// \endcond
#endif

AliFemtoEventCut*    copyTheCut(AliFemtoEventCut*);
AliFemtoParticleCut* copyTheCut(AliFemtoParticleCut*);
AliFemtoPairCut*     copyTheCut(AliFemtoPairCut*);
AliFemtoCorrFctn*    copyTheCorrFctn(AliFemtoCorrFctn*);


/// Generalized particle collection filler function - called by
/// FillParticleCollection()
///
/// This function loops over the track_collection, calling the cut's Pass
/// method on each track. If it passes, a new AliFemtoParticle is constructed
/// from the track's pointer (or more generally, the item returned by
/// dereferencing the container's iterator) and the cut's expected mass.
///
/// This templated function accepts a track cut, track collection, and an
/// AliFemtoParticleCollection (which points to the output) as input. The types
/// of the tracks are determined by the template paramters, which should be
/// automatically detected by argument inspection (you don't need to specify).
///
/// The original function also specified the iterator type, but since all
/// containers are standard STL containers, it now infers the type as
/// TrackCollectionType::iterator. If the track collections change to some
/// other type, it is recommended to add TrackCollectionIterType to the
/// template list, and add the appropriate type to the function calls in
/// FillParticleCollection.
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

// This little function is used to apply ParticleCuts (TrackCuts or V0Cuts) and
// fill ParticleCollections from tacks in picoEvent. It is called from
// AliFemtoSimpleAnalysis::ProcessEvent().
//
// The actual loop implementation has been moved to the collection-generic
// DoFillParticleCollection() function
void FillHbtParticleCollection(AliFemtoParticleCut *partCut,
                               AliFemtoEvent *hbtEvent,
                               AliFemtoParticleCollection *partCollection,
                               bool performSharedDaughterCut=kFALSE)
{
  /// Fill particle collection with all particles in the event which pass
  /// the provided cut

  // determine which track collection to use based on the particle type.
  switch (partCut->Type()) {

  // cut is cutting on Tracks
  case hbtTrack:

    DoFillParticleCollection(
      (AliFemtoTrackCut*)partCut,
      hbtEvent->TrackCollection(),
      partCollection
    );

    break;

  // cut is cutting on V0s
  case hbtV0:
  {
    AliFemtoV0Cut *v0_cut = (AliFemtoV0Cut*)partCut;

    // shared daughter cut returns all passed v0s - add to particle collection
    if (performSharedDaughterCut) {
      AliFemtoV0SharedDaughterCut shared_daughter_cut;
      AliFemtoV0Collection v0_coll = shared_daughter_cut.AliFemtoV0SharedDaughterCutCollection(hbtEvent->V0Collection(), v0_cut);
      for (AliFemtoV0Iterator pIter = v0_coll.begin(); pIter != v0_coll.end(); ++pIter) {
        partCollection->push_back(new AliFemtoParticle(*pIter, v0_cut->Mass()));
      }
    } else {

      DoFillParticleCollection(
        v0_cut,
        hbtEvent->V0Collection(),
        partCollection
      );

    }

    break;
  }

  // cut is cutting on Xis
  case hbtXi:
  {
    AliFemtoXiTrackCut *xi_cut = (AliFemtoXiTrackCut*)partCut;

    // shared daughter cut returns all passed xis - add to particle collection
    if (performSharedDaughterCut)
    {
      AliFemtoXiSharedDaughterCut shared_daughter_cut;
      AliFemtoXiCollection xi_coll = shared_daughter_cut.AliFemtoXiSharedDaughterCutCollection(hbtEvent->XiCollection(), xi_cut);
      for (AliFemtoXiIterator pIter = xi_coll.begin(); pIter != xi_coll.end(); ++pIter) {
        partCollection->push_back(new AliFemtoParticle(*pIter, xi_cut->Mass()));
      }
    } 
    else
    {
      DoFillParticleCollection(
        (AliFemtoXiTrackCut*)partCut,
        hbtEvent->XiCollection(),
        partCollection
      );
    }
    break;
  }

  // cut is cutting on Kinks
  case hbtKink:

    DoFillParticleCollection(
      (AliFemtoKinkCut*)partCut,
      hbtEvent->KinkCollection(),
      partCollection
    );

    break;

  default:
    cout << "E-FillHbtParticleCollection function (in AliFemtoSimpleAnalysis.cxx): "
            "Undefined Particle Cut type!!! (" << partCut->Type() << ")\n";
  }

  partCut->FillCutMonitor(hbtEvent, partCollection);
}
//____________________________
AliFemtoSimpleAnalysis::AliFemtoSimpleAnalysis():
  fPicoEventCollectionVectorHideAway(NULL),
  fPairCut(NULL),
  fCorrFctnCollection(NULL),
  fEventCut(NULL),
  fFirstParticleCut(NULL),
  fSecondParticleCut(NULL),
  fMixingBuffer(NULL),
  fPicoEvent(NULL),
  fNumEventsToMix(0),
  fNeventsProcessed(0),
  fMinSizePartCollection(0),
  fVerbose(kTRUE),
  fPerformSharedDaughterCut(kFALSE),
  fEnablePairMonitors(kFALSE)
{
  // Default constructor
  fCorrFctnCollection = new AliFemtoCorrFctnCollection;
  fMixingBuffer = new AliFemtoPicoEventCollection;
}
//____________________________
AliFemtoSimpleAnalysis::AliFemtoSimpleAnalysis(const AliFemtoSimpleAnalysis& a):
  AliFemtoAnalysis(),
  fPicoEventCollectionVectorHideAway(NULL),
  fPairCut(NULL),
  fCorrFctnCollection(NULL),
  fEventCut(NULL),
  fFirstParticleCut(NULL),
  fSecondParticleCut(NULL),
  fMixingBuffer(NULL),
  fPicoEvent(NULL),
  fNumEventsToMix(a.fNumEventsToMix),
  fNeventsProcessed(0),
  fMinSizePartCollection(a.fMinSizePartCollection),
  fVerbose(a.fVerbose),
  fPerformSharedDaughterCut(a.fPerformSharedDaughterCut),
  fEnablePairMonitors(a.fEnablePairMonitors)
{
  /// Copy constructor

  const char msg_template[] = " AliFemtoSimpleAnalysis::AliFemtoSimpleAnalysis(const AliFemtoSimpleAnalysis& a) - %s",
            warn_template[] = " WARNING [AliFemtoSimpleAnalysis::AliFemtoSimpleAnalysis(const AliFemtoSimpleAnalysis& a)] %s";

  fCorrFctnCollection = new AliFemtoCorrFctnCollection;
  fMixingBuffer = new AliFemtoPicoEventCollection;

  // Clone the event cut
  AliFemtoEventCut *ev_cut = a.fEventCut->Clone();
  if (ev_cut) {
    SetEventCut(ev_cut);
    if (fVerbose) {
      cout << " event cut set " << endl;
    }
  } else {
    cerr << TString::Format(warn_template, "Could not clone event cut") << endl;
    // TODO: handle uncloned event cut
  }

  // Clone the pair cut
  fPairCut = a.fPairCut->Clone();
  if (fPairCut) {
    SetPairCut(fPairCut);
    if (fVerbose) {
       cout << TString::Format(msg_template, "pair cut set") << endl;
    }
  } else {
    cerr << TString::Format(warn_template, "Could not clone pair cut") << endl;
    // TODO: handle uncloned pair cut
  }


  fFirstParticleCut = a.fFirstParticleCut->Clone();
  if (fFirstParticleCut) {
    SetFirstParticleCut(fFirstParticleCut);
    if (fVerbose) {
      cout << TString::Format(msg_template, "first particle cut set") << endl;
    }
  } else {
    cerr << TString::Format(warn_template, "Could not clone first particle cut") << endl;
    // TODO: handle uncloned track cut
  }


  fSecondParticleCut = (a.fFirstParticleCut == a.fSecondParticleCut)
                     ? fFirstParticleCut
                     : a.fSecondParticleCut->Clone();
  if (fSecondParticleCut) {
    SetSecondParticleCut(fSecondParticleCut);
    if (fVerbose) {
      cout << TString::Format(msg_template, "second particle cut set") << endl;
    }
  } else {
    cerr << TString::Format(warn_template, "Could not clone second particle cut") << endl;
    // TODO: handle uncloned track cut
  }

  AliFemtoCorrFctnIterator iter;
  if (fVerbose) {
    cout << TString::Format(msg_template, "looking for correlation functions") << endl;
  }
  for (iter = a.fCorrFctnCollection->begin(); iter != a.fCorrFctnCollection->end(); ++iter) {
    AliFemtoCorrFctn* fctn = (*iter)->Clone();
    if (fctn) {
      AddCorrFctn(fctn);
    } else {
      cout << TString::Format(msg_template, "correlation function not found") << endl;
    }
  }

  if (fVerbose) {
    cout << TString::Format(msg_template, "analysis copied") << endl;
  }
}
//____________________________
AliFemtoSimpleAnalysis::~AliFemtoSimpleAnalysis()
{
  /// destructor

  if (fVerbose) {
    cout << " AliFemtoSimpleAnalysis::~AliFemtoSimpleAnalysis()" << endl;
  }

  // will not double-delete particle cut
  if (fFirstParticleCut == fSecondParticleCut) {
    fSecondParticleCut = NULL;
  }

  delete fPairCut;
  delete fEventCut;
  delete fFirstParticleCut;
  delete fSecondParticleCut;

  // delete every CorrFunction in the collection, then the collection
  if (fCorrFctnCollection) {
    for (AliFemtoCorrFctnIterator iter = fCorrFctnCollection->begin(); iter != fCorrFctnCollection->end(); iter++) {
      delete *iter;
    }
    delete fCorrFctnCollection;
  }

  // delete every PicoEvent in the EventMixingBuffer followed by the buffer
  if (fMixingBuffer) {
    for (AliFemtoPicoEventIterator piter = fMixingBuffer->begin(); piter != fMixingBuffer->end(); ++piter) {
      delete *piter;
    }
    delete fMixingBuffer;
  }
}
//______________________
AliFemtoSimpleAnalysis& AliFemtoSimpleAnalysis::operator=(const AliFemtoSimpleAnalysis& aAna)
{
  /// Assignment operator

  if (this == &aAna)
    return *this;

  // clear second particle cut to avoid double delete
  if (fFirstParticleCut == fSecondParticleCut) {
    fSecondParticleCut = NULL;
  }

  // delete current pointers
  delete fPairCut;
  delete fEventCut;
  delete fFirstParticleCut;
  delete fSecondParticleCut;

  // clear correlation functions out of fCorrFctnCollection
  if (fCorrFctnCollection) {
    for (AliFemtoCorrFctnIterator iter = fCorrFctnCollection->begin(); iter != fCorrFctnCollection->end(); ++iter) {
      delete *iter;
    }
    fCorrFctnCollection->clear();
  } else {
    cerr << " WARNING [AliFemtoSimpleAnalysis::operator=()] fCorrFctnCollection was NULL, this should not happen." << endl;
    fCorrFctnCollection = new AliFemtoCorrFctnCollection;
  }

  // clear mixing buffer
  if (fMixingBuffer) {
    for (AliFemtoPicoEventIterator piter = fMixingBuffer->begin(); piter != fMixingBuffer->end(); ++piter) {
      delete *piter;
    }
    fMixingBuffer->clear();
  } else {
    cerr << " WARNING [AliFemtoSimpleAnalysis::operator=()] fMixingBuffer was NULL, this should not happen." << endl;
    fMixingBuffer = new AliFemtoPicoEventCollection;
  }

  // clone objects
  fPairCut = aAna.fPairCut->Clone();
  fEventCut = aAna.fEventCut->Clone();
  fFirstParticleCut = aAna.fFirstParticleCut->Clone();
  fSecondParticleCut = (aAna.fFirstParticleCut == aAna.fSecondParticleCut)
                     ? fFirstParticleCut
                     : aAna.fSecondParticleCut->Clone();

  if (fPairCut) {
    SetPairCut(fPairCut);
  } else {
    cerr << " WARNING [AliFemtoSimpleAnalysis::operator=()] Could not clone pair cut." << endl;
  }

  if (fEventCut) {
    SetEventCut(fEventCut);
  } else {
    cerr << " WARNING [AliFemtoSimpleAnalysis::operator=()] Could not clone event cut." << endl;
  }

  if (fFirstParticleCut) {
    SetFirstParticleCut(fFirstParticleCut);
  } else {
    cerr << " WARNING [AliFemtoSimpleAnalysis::operator=()] Could not clone first particle cut." << endl;
  }

  if (fSecondParticleCut) {
    SetSecondParticleCut(fSecondParticleCut);
  } else {
    cerr << " WARNING [AliFemtoSimpleAnalysis::operator=()] Could not clone second particle cut." << endl;
  }

  for (AliFemtoCorrFctnIterator iter = aAna.fCorrFctnCollection->begin(); iter != aAna.fCorrFctnCollection->end(); ++iter) {
    AliFemtoCorrFctn* fctn = (*iter)->Clone();
    if (fctn) AddCorrFctn(fctn);
  }


  fNumEventsToMix = aAna.fNumEventsToMix;
  fMinSizePartCollection = aAna.fMinSizePartCollection;
  fVerbose = aAna.fVerbose;
  fPerformSharedDaughterCut = aAna.fPerformSharedDaughterCut;
  fEnablePairMonitors = aAna.fEnablePairMonitors;

  return *this;
}
//______________________
AliFemtoCorrFctn* AliFemtoSimpleAnalysis::CorrFctn(int n)
{
  /// return pointer to n-th correlation function

  if (n < 0 || n > (int)fCorrFctnCollection->size()) {
    return NULL;
  }

  AliFemtoCorrFctnIterator iter = fCorrFctnCollection->begin();
  std::advance(iter, n);
  return *iter;
}
//____________________________
AliFemtoString AliFemtoSimpleAnalysis::Report()
{
  /// Create a simple report from the analysis execution

  cout << "AliFemtoSimpleAnalysis - constructing Report..."<<endl;
  string temp = "-----------\nHbt Analysis Report:\n";
  temp += "\nEvent Cuts:\n";
  temp += fEventCut->Report();
  temp += "\nParticle Cuts - First Particle:\n";
  temp += fFirstParticleCut->Report();
  temp += "\nParticle Cuts - Second Particle:\n";
  temp += fSecondParticleCut->Report();
  temp += "\nPair Cuts:\n";
  temp += fPairCut->Report();
  temp += "\nCorrelation Functions:\n";

  if (fCorrFctnCollection->empty()) {
    cout << "AliFemtoSimpleAnalysis-Warning : no correlations functions in this analysis " << endl;
  }
  for (AliFemtoCorrFctnIterator iter = fCorrFctnCollection->begin(); iter != fCorrFctnCollection->end(); ++iter) {
    temp += (*iter)->Report();
    temp += "\n";
  }
  temp += "-------------\n";
  AliFemtoString returnThis=temp;
  return returnThis;
}
//_________________________
void AliFemtoSimpleAnalysis::ProcessEvent(const AliFemtoEvent* hbtEvent)
{
  // Add event to processed events

  // We will get a new pico event; NULL now to prevent corr fctn access to old pico event
  fPicoEvent = NULL;

  // increment number of events processed
  AddEventProcessed();

  // startup for EbyE
  EventBegin(hbtEvent);

  // event cut and event cut monitor
  bool tmpPassEvent = fEventCut->Pass(hbtEvent);

  if (!tmpPassEvent) {
    fEventCut->FillCutMonitor(hbtEvent, tmpPassEvent);
    EventEnd(hbtEvent);  // cleanup for EbyE
    return;
  }

  // Analysis likes the event -- build a pico event from it, using tracks the
  // analysis likes. This is what we will make pairs from and put in Mixing
  // Buffer.
  // No memory leak: we will delete picoevents when they come out of the
  // mixing buffer
  fPicoEvent = new AliFemtoPicoEvent;

  AliFemtoParticleCollection *collection1 = fPicoEvent->FirstParticleCollection(),
                             *collection2 = fPicoEvent->SecondParticleCollection();

  // Sanity check that both collections exist
  if (collection1 == NULL || collection2 == NULL) {
    cout << "E-AliFemtoSimpleAnalysis::ProcessEvent: new PicoEvent is missing particle collections!\n";
    EventEnd(hbtEvent);  // cleanup for EbyE
    delete fPicoEvent;
    return;
  }

  // Subroutine fills fPicoEvent'a FirstParticleCollection with tracks from
  // hbtEvent which pass fFirstParticleCut. Uses cut's "Type()" to determine
  // which track collection to pull from hbtEvent.
  FillHbtParticleCollection(fFirstParticleCut,
                            (AliFemtoEvent*)hbtEvent,
                            fPicoEvent->FirstParticleCollection(),
                            fPerformSharedDaughterCut);

  // fill second particle cut if not analyzing identical particles
  if ( !AnalyzeIdenticalParticles() ) {
      FillHbtParticleCollection(fSecondParticleCut,
                                (AliFemtoEvent*)hbtEvent,
                                fPicoEvent->SecondParticleCollection(),
                                fPerformSharedDaughterCut);
  }

  const UInt_t coll_1_size = collection1->size(),
               coll_2_size = collection2->size();

  if (fVerbose) {
    cout << "#particles in Collection 1, 2: "
         << coll_1_size << " " << coll_2_size << "\n";
  }

  // now we have created the particle collections - we fill the cut monitors
  fEventCut->FillCutMonitor(collection1, collection2); //MJ!

  const bool coll_1_size_passes = (coll_1_size >= fMinSizePartCollection),
             coll_2_size_passes = (AnalyzeIdenticalParticles() || (coll_2_size >= fMinSizePartCollection));

  // passes only if the sizes are above min
  tmpPassEvent = tmpPassEvent
              && coll_1_size_passes
              && coll_2_size_passes;

  // fill the event cut monitor
  fEventCut->FillCutMonitor(hbtEvent, tmpPassEvent);

  if (!tmpPassEvent) {
    EventEnd(hbtEvent);
    delete fPicoEvent;
    return;
  }

  //------ Make real pairs. If identical, make pairs for one collection ------//
  if (AnalyzeIdenticalParticles()) {
    collection2 = NULL;
  }

  MakePairs("real", collection1, collection2, EnablePairMonitors());

  if (fVerbose) {
    cout << "AliFemtoSimpleAnalysis::ProcessEvent() - reals done ";
  }

  //---- Make pairs for mixed events, looping over events in mixingBuffer ----//
  for (AliFemtoPicoEventIterator fPicoEventIter = MixingBuffer()->begin();
                                 fPicoEventIter != MixingBuffer()->end();
                               ++fPicoEventIter) {

    AliFemtoPicoEvent *storedEvent = *fPicoEventIter;

    // If identical - only mix the first particle collections
    if (AnalyzeIdenticalParticles()) {
      MakePairs("mixed", collection1, storedEvent->FirstParticleCollection());

    // If non-identical - mix both combinations of first and second particles
    } else {
        MakePairs("mixed", collection1,
                           storedEvent->SecondParticleCollection());

        MakePairs("mixed", storedEvent->FirstParticleCollection(),
                           collection2);
    }
  }

  if (fVerbose) {
    cout << " - mixed done   " << endl;
  }

  //--------- If mixing buffer is full, delete oldest event ---------//
  if ( MixingBufferFull() ) {
    delete MixingBuffer()->back();
    MixingBuffer()->pop_back();
  }

  //-------- Add current event (fPicoEvent) to mixing buffer --------//
  MixingBuffer()->push_front(fPicoEvent);

  EventEnd(hbtEvent);  // cleanup for EbyE
  //cout << "AliFemtoSimpleAnalysis::ProcessEvent() - return to caller ... " << endl;
}

//_________________________
void AliFemtoSimpleAnalysis::MakePairs(const char* typeIn,
                                       AliFemtoParticleCollection *partCollection1,
                                       AliFemtoParticleCollection *partCollection2,
                                       Bool_t enablePairMonitors)
{
/// Build pairs, check pair cuts, and call CFs' AddRealPair() or
/// AddMixedPair() methods. If no second particle collection is
/// specfied, make pairs within first particle collection.

  const string type = typeIn;

  //  int swpart = ((long int) partCollection1) % 2;

  // Used to swap particle 1 & 2 in identical-particle analysis
  // to avoid any implicit ordering in the event collection
  // "Seed" this here.
  bool swpart = fNeventsProcessed % 2;

  // Setup iterator ranges
  //
  // The outer loop alway starts at beginning of particle collection 1.
  // * If we are iterating over both particle collections, then the loop simply
  // runs through both from beginning to end.
  // * If we are only iterating over one particle collection, the inner loop
  // loops over all particles between the outer iterator and the end of the
  // collection. The outer loop must skip the last entry of the list.
  AliFemtoParticleConstIterator tStartOuterLoop = partCollection1->begin(),
                                tEndOuterLoop = partCollection1->end(),
                                tStartInnerLoop,
                                tEndInnerLoop;

  if (partCollection2) {                         // Two collections:
    tStartInnerLoop = partCollection2->begin();  //   Full inner & outer loops
    tEndInnerLoop   = partCollection2->end();    //
  }
  else {                                         // One collection:
    tEndOuterLoop--;                             //   Outer loop goes to next-to-last particle
    tEndInnerLoop = partCollection1->end() ;     //   Inner loop goes to last particle
  }

  // Create the pair outside the loop - only allocate once
  AliFemtoPair* tPair = new AliFemtoPair;

  // Begin the outer loop
  for (AliFemtoParticleConstIterator tPartIter1 = tStartOuterLoop;
                                     tPartIter1 != tEndOuterLoop;
                                     ++tPartIter1) {

    // If analyzing identical particles, start inner loop at the particle
    // after the current outer loop position, (loops until end)
    if (!partCollection2) {
      tStartInnerLoop = tPartIter1;
      tStartInnerLoop++;
    }

    // If we have two collections - set the first track
    if (partCollection2 != NULL) {
      tPair->SetTrack1(*tPartIter1);
    }

    // Begin the inner loop
    for (AliFemtoParticleConstIterator tPartIter2 = tStartInnerLoop;
                                       tPartIter2 != tEndInnerLoop;
                                     ++tPartIter2) {
      // If we have two collections - only set the second track
      if (partCollection2 != NULL) {
        tPair->SetTrack2(*tPartIter2);

      // Swap between first and second particles to avoid biased ordering
      } else {
        tPair->SetTrack1(swpart ? *tPartIter2 : *tPartIter1);
        tPair->SetTrack2(swpart ? *tPartIter1 : *tPartIter2);
        swpart = !swpart;
      }

      // check if the pair passes the cut
      bool tmpPassPair = fPairCut->Pass(tPair);

      // This is a condition for speed reasons
      if (enablePairMonitors) {
        fPairCut->FillCutMonitor(tPair, tmpPassPair);
      }

      // If pair passes cut, loop over CF's and add pair to real/mixed
      if (tmpPassPair) {
        for (AliFemtoCorrFctnIterator tCorrFctnIter = fCorrFctnCollection->begin();
                                      tCorrFctnIter != fCorrFctnCollection->end();
                                    ++tCorrFctnIter) {

          AliFemtoCorrFctn* tCorrFctn = *tCorrFctnIter;

          if (type == "real")
            tCorrFctn->AddRealPair(tPair);
          else if(type == "mixed")
            tCorrFctn->AddMixedPair(tPair);
          else
            cout << "Problem with pair type, type = " << type << endl;
        } // loop over corellatoin functions
      }
    }    // loop over second particle
  }      // loop over first particle

  // we are done with the pair
  delete tPair;
}
//_________________________
void AliFemtoSimpleAnalysis::EventBegin(const AliFemtoEvent* ev)
{
  /// Perform initialization operations at the beginning of the event processing
  /// cout << " AliFemtoSimpleAnalysis::EventBegin(const AliFemtoEvent* ev) " << endl;

  fFirstParticleCut->EventBegin(ev);
  fSecondParticleCut->EventBegin(ev);
  fPairCut->EventBegin(ev);
  for (AliFemtoCorrFctnIterator iter = fCorrFctnCollection->begin();
                                iter != fCorrFctnCollection->end();
                                ++iter) {
    (*iter)->EventBegin(ev);
  }
}
//_________________________
void AliFemtoSimpleAnalysis::EventEnd(const AliFemtoEvent* ev)
{
  // Finish operations at the end of event processing

  fFirstParticleCut->EventEnd(ev);
  fSecondParticleCut->EventEnd(ev);
  fPairCut->EventEnd(ev);
  for (AliFemtoCorrFctnIterator iter = fCorrFctnCollection->begin();
                                iter != fCorrFctnCollection->end();
                                ++iter) {
    (*iter)->EventEnd(ev);
  }
}
//_________________________
void AliFemtoSimpleAnalysis::Finish()
{
  // Perform finishing operations after all events are processed

  for (AliFemtoCorrFctnIterator iter = fCorrFctnCollection->begin();
                                iter != fCorrFctnCollection->end();
                                ++iter) {
    (*iter)->Finish();
  }
}
//_________________________
void AliFemtoSimpleAnalysis::AddEventProcessed()
{
  // Increase count of processed events
  fNeventsProcessed++;
}
//_________________________
TList* AliFemtoSimpleAnalysis::ListSettings()
{
  // Collect settings list

  const TString setting_prefix = "AliFemtoSimpleAnalysis.";

  TList *tListSettings = new TList();

  TList *event_cut_settings = fEventCut->ListSettings();

  if (event_cut_settings != NULL) {

    TListIter next_event_setting(event_cut_settings);
    while (TObject *obj = next_event_setting()) {
      tListSettings->Add(new TObjString(
        setting_prefix + obj->GetName()
      ));
    }

    delete event_cut_settings;
  }

  TList *p1_cut_settings = fFirstParticleCut->ListSettings();

  if (p1_cut_settings != NULL) {

    TListIter next_p1_setting(p1_cut_settings);
    while (TObject *obj = next_p1_setting()) {
      tListSettings->Add(new TObjString(
        setting_prefix + obj->GetName()
      ));
    }

    delete p1_cut_settings;
  }


  if (fSecondParticleCut != fFirstParticleCut) {
    TList *p2_cut_settings = fSecondParticleCut->ListSettings();

    if (p2_cut_settings != NULL) {

      TListIter next_p2_setting(p2_cut_settings);
      while (TObject *obj = next_p2_setting()) {
        tListSettings->Add(new TObjString(
          setting_prefix + obj->GetName()
        ));
      }

      delete p2_cut_settings;
    }
  }


  TList *pair_cut_settings = fPairCut->ListSettings();

  if (pair_cut_settings != NULL) {

    TListIter next_pair_setting(pair_cut_settings);
    while (TObject *obj = next_pair_setting()) {
      tListSettings->Add(new TObjString(
        setting_prefix + obj->GetName()
      ));
    }

    delete pair_cut_settings;
  }

  return tListSettings;
}

//_________________________
TList* AliFemtoSimpleAnalysis::GetOutputList()
{
  // Collect the list of output objects to be written

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

  TList *pairCut = fPairCut->GetOutputList();

  TIter nextpair(pairCut);
  while (TObject *obj = nextpair()) {
    tOutputList->Add(obj);
  }
  delete pairCut;

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
