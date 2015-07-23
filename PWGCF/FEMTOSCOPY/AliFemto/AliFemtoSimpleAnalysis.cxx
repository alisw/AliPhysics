///
/// \file AliFemtoSimpleAnalysis.cxx
///

#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoTrackCut.h"
#include "AliFemtoV0Cut.h"
#include "AliFemtoKinkCut.h"

#include <string>
#include <iostream>
#include <iterator>

#ifdef __ROOT__
/// \cond CLASSIMP
ClassImp(AliFemtoSimpleAnalysis)
/// \endcond
#endif

AliFemtoEventCut*    copyTheCut(AliFemtoEventCut*);
AliFemtoParticleCut* copyTheCut(AliFemtoParticleCut*);
AliFemtoPairCut*     copyTheCut(AliFemtoPairCut*);
AliFemtoCorrFctn*    copyTheCorrFctn(AliFemtoCorrFctn*);

// this little function used to apply ParticleCuts (TrackCuts or V0Cuts) and fill ParticleCollections of picoEvent
//  it is called from AliFemtoSimpleAnalysis::ProcessEvent()
void FillHbtParticleCollection(AliFemtoParticleCut *partCut,
                               AliFemtoEvent *hbtEvent,
                               AliFemtoParticleCollection *partCollection,
                               bool performSharedDaughterCut=kFALSE)
{
  /// Fill particle collections from the event
  /// by the particles that pass all the cuts

  switch (partCut->Type()) {
  case hbtTrack:       // cut is cutting on Tracks
  {
    AliFemtoTrackCut *pCut = dynamic_cast<AliFemtoTrackCut*>(partCut);
    AliFemtoTrackIterator pIter;
    AliFemtoTrackIterator startLoop = hbtEvent->TrackCollection()->begin();
    AliFemtoTrackIterator endLoop   = hbtEvent->TrackCollection()->end();
    for (pIter=startLoop; pIter!=endLoop; pIter++) {
      AliFemtoTrack *next_particle = *pIter;
      bool tmpPassParticle = pCut->Pass(next_particle);
      pCut->FillCutMonitor(next_particle, tmpPassParticle);
      if (tmpPassParticle) {
        AliFemtoParticle* particle = new AliFemtoParticle(next_particle, pCut->Mass());
        partCollection->push_back(particle);
      }
    }
    break;
  }
  case hbtV0:          // cut is cutting on V0s
    {
      AliFemtoV0Cut* pCut = (AliFemtoV0Cut*) partCut;
      AliFemtoV0* pParticle;
      AliFemtoV0Iterator pIter;

      if(performSharedDaughterCut) {
        AliFemtoV0Collection V0CorrectedCollection;
        AliFemtoV0SharedDaughterCut sharedDaughterCut;
        V0CorrectedCollection = sharedDaughterCut.AliFemtoV0SharedDaughterCutCollection(hbtEvent->V0Collection(), pCut);

        AliFemtoV0Iterator startLoop = V0CorrectedCollection.begin();
        AliFemtoV0Iterator endLoop   = V0CorrectedCollection.end();

        for (pIter=startLoop;pIter!=endLoop;pIter++) {
          pParticle = *pIter;
          AliFemtoParticle* particle = new AliFemtoParticle(pParticle,partCut->Mass());
          partCollection->push_back(particle);
        }
      }
      else { //previous, untouched loop:
        AliFemtoV0Iterator startLoop = hbtEvent->V0Collection()->begin();
        AliFemtoV0Iterator endLoop   = hbtEvent->V0Collection()->end();
        for (pIter=startLoop;pIter!=endLoop;pIter++) {
          pParticle = *pIter;
          bool tmpPassV0 = pCut->Pass(pParticle);
          pCut->FillCutMonitor(pParticle,tmpPassV0);
          if(tmpPassV0) {
            AliFemtoParticle* particle = new AliFemtoParticle(pParticle,partCut->Mass());
            partCollection->push_back(particle);
          }
        }
      }

      pCut->FillCutMonitor(hbtEvent,partCollection);// Gael 19/06/02

      break;
    }
  case hbtKink:          // cut is cutting on Kinks  -- mal 25May2001
    {
      AliFemtoKinkCut* pCut = (AliFemtoKinkCut*) partCut;
      AliFemtoKink* pParticle;
      AliFemtoKinkIterator pIter;
      AliFemtoKinkIterator startLoop = hbtEvent->KinkCollection()->begin();
      AliFemtoKinkIterator endLoop   = hbtEvent->KinkCollection()->end();
      // this following "for" loop is identical to the one above, but because of scoping, I can's see how to avoid repitition...
      for (pIter=startLoop;pIter!=endLoop;pIter++) {
	pParticle = *pIter;
	bool tmpPass = pCut->Pass(pParticle);
	pCut->FillCutMonitor(pParticle, tmpPass);
	if (tmpPass) {
	  AliFemtoParticle* particle = new AliFemtoParticle(pParticle, partCut->Mass());
	  partCollection->push_back(particle);
	}
      }
      break;
    }
  default:
    cout << "FillHbtParticleCollection function (in AliFemtoSimpleAnalysis.cxx) - undefined Particle Cut type!!! \n";
  }
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
  /// Default constructor
  ///  mControlSwitch     = 0;

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

  // now delete every CorrFunction in the Collection, and then the Collection itself
  if (fCorrFctnCollection) {
    for (AliFemtoCorrFctnIterator iter = fCorrFctnCollection->begin(); iter != fCorrFctnCollection->end(); iter++) {
      delete *iter;
    }
    delete fCorrFctnCollection;
  }

  // now delete every PicoEvent in the EventMixingBuffer and then the buffer itself
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
void AliFemtoSimpleAnalysis::ProcessEvent(const AliFemtoEvent* hbtEvent) {
  /// Add event to processed events

  fPicoEvent = NULL; // we will get a new pico event, if not prevent corr. fctn to access old pico event
  AddEventProcessed();
  // startup for EbyE
  EventBegin(hbtEvent);
  // event cut and event cut monitor
  bool tmpPassEvent = fEventCut->Pass(hbtEvent);
  if (!tmpPassEvent)
    fEventCut->FillCutMonitor(hbtEvent, tmpPassEvent);
  if (tmpPassEvent) {
    //  cout << "AliFemtoSimpleAnalysis::ProcessEvent() - Event has passed cut - build picoEvent from " <<
    //   hbtEvent->TrackCollection()->size() << " tracks in TrackCollection" << endl;
    //  cout << "Event has passed cut with " << hbtEvent->TrackCollection()->size() << " tracks" << endl;
    // OK, analysis likes the event-- build a pico event from it, using tracks the analysis likes...
    fPicoEvent = new AliFemtoPicoEvent; // this is what we will make pairs from and put in Mixing Buffer
    // no memory leak. we will delete picoevents when they come out of the mixing buffer
    FillHbtParticleCollection(fFirstParticleCut,(AliFemtoEvent*)hbtEvent,fPicoEvent->FirstParticleCollection(), fPerformSharedDaughterCut);
    if ( !(AnalyzeIdenticalParticles()) )
      FillHbtParticleCollection(fSecondParticleCut,(AliFemtoEvent*)hbtEvent,fPicoEvent->SecondParticleCollection(),
				fPerformSharedDaughterCut);
    //cout <<"AliFemtoSimpleAnalysis::ProcessEvent - #particles in First, Second Collections: " <<
//       fPicoEvent->FirstParticleCollection()->size() << " " <<
//       fPicoEvent->SecondParticleCollection()->size() << endl;


    if (fVerbose)
    cout << "#particles in Collection 1, 2: " <<
       fPicoEvent->FirstParticleCollection()->size() << " " <<
       fPicoEvent->SecondParticleCollection()->size() << endl;
    fEventCut->FillCutMonitor(fPicoEvent->FirstParticleCollection(),fPicoEvent->SecondParticleCollection()); //MJ!


    // mal - implement a switch which allows only using events with ParticleCollections containing a minimum
    // number of entries (jun2002)
    if ((fPicoEvent->FirstParticleCollection()->size() >= fMinSizePartCollection )
	&& ( AnalyzeIdenticalParticles() || (fPicoEvent->SecondParticleCollection()->size() >= fMinSizePartCollection ))) {
      fEventCut->FillCutMonitor(hbtEvent, tmpPassEvent);


//------------------------------------------------------------------------------
//   Temporary comment:
//      This whole section rewritten so that all pairs are built using the
//      same code... easier to read and manage, and MakePairs() can be called by
//      derived classes.  Also, the requirement of a full mixing buffer before
//      mixing is removed.
//                          Dan Magestro, 11/2002

      //------ Make real pairs. If identical, make pairs for one collection ------//

      if (AnalyzeIdenticalParticles()) {
        MakePairs("real", fPicoEvent->FirstParticleCollection(), 0, EnablePairMonitors());
      }
      else {
        MakePairs("real", fPicoEvent->FirstParticleCollection(),
		  fPicoEvent->SecondParticleCollection(), EnablePairMonitors() );
      }

      if (fVerbose)
       cout << "AliFemtoSimpleAnalysis::ProcessEvent() - reals done ";

      //---- Make pairs for mixed events, looping over events in mixingBuffer ----//

      AliFemtoPicoEvent* storedEvent;
      AliFemtoPicoEventIterator fPicoEventIter;
      for (fPicoEventIter=MixingBuffer()->begin();fPicoEventIter!=MixingBuffer()->end();fPicoEventIter++) {
        storedEvent = *fPicoEventIter;
        if (AnalyzeIdenticalParticles()) {
          MakePairs("mixed",fPicoEvent->FirstParticleCollection(),
                            storedEvent->FirstParticleCollection() );
        }
        else {
          MakePairs("mixed",fPicoEvent->FirstParticleCollection(),
                            storedEvent->SecondParticleCollection() );

          MakePairs("mixed",storedEvent->FirstParticleCollection(),
                            fPicoEvent->SecondParticleCollection() );
        }
      }

      if (fVerbose)
       cout << " - mixed done   " << endl;

      //--------- If mixing buffer is full, delete oldest event ---------//

      if ( MixingBufferFull() ) {
        delete MixingBuffer()->back();
        MixingBuffer()->pop_back();
      }

      //-------- Add current event (fPicoEvent) to mixing buffer --------//

      MixingBuffer()->push_front(fPicoEvent);


// Temporary comment: End of rewritten section... Dan Magestro, 11/2002
//------------------------------------------------------------------------------


    }  // if ParticleCollections are big enough (mal jun2002)
    else {
      fEventCut->FillCutMonitor(hbtEvent, !tmpPassEvent);
      delete fPicoEvent;
    }
  }   // if currentEvent is accepted by currentAnalysis
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

  string type = typeIn;

  //  int swpart = ((long int) partCollection1) % 2;
  int swpart = fNeventsProcessed % 2;

  AliFemtoPair* tPair = new AliFemtoPair;

  AliFemtoCorrFctnIterator tCorrFctnIter;

  AliFemtoParticleIterator tPartIter1, tPartIter2;

  AliFemtoParticleIterator tStartOuterLoop = partCollection1->begin();  // always
  AliFemtoParticleIterator tEndOuterLoop   = partCollection1->end();    // will be one less if identical
  AliFemtoParticleIterator tStartInnerLoop;
  AliFemtoParticleIterator tEndInnerLoop;
  if (partCollection2) {                         // Two collections:
    tStartInnerLoop = partCollection2->begin();  //   Full inner & outer loops
    tEndInnerLoop   = partCollection2->end();    //
  }
  else {                                         // One collection:
    tEndOuterLoop--;                             //   Outer loop goes to next-to-last particle
    tEndInnerLoop = partCollection1->end() ;     //   Inner loop goes to last particle
  }
  for (tPartIter1=tStartOuterLoop;tPartIter1!=tEndOuterLoop;tPartIter1++) {
    if (!partCollection2) {
      tStartInnerLoop = tPartIter1;
      tStartInnerLoop++;
    }
    tPair->SetTrack1(*tPartIter1);
    for (tPartIter2 = tStartInnerLoop; tPartIter2!=tEndInnerLoop;tPartIter2++) {
      tPair->SetTrack2(*tPartIter2);

      // The following lines have to be uncommented if you want pairCutMonitors
      // they are not in // for speed reasons
      if(enablePairMonitors) {
	bool tmpPassPair = fPairCut->Pass(tPair);
	fPairCut->FillCutMonitor(tPair, tmpPassPair);
      }
      // // if ( tmpPassPair )

      //---- If pair passes cut, loop over CF's and add pair to real/mixed ----//

      if (!partCollection2) {
	if (swpart) {
 	  tPair->SetTrack1(*tPartIter2);
 	  tPair->SetTrack2(*tPartIter1);
 	  swpart = 0;
 	}
 	else {
 	  tPair->SetTrack1(*tPartIter1);
 	  tPair->SetTrack2(*tPartIter2);
 	  swpart = 1;
	}
      }

      if (fPairCut->Pass(tPair)) {
        for (tCorrFctnIter=fCorrFctnCollection->begin(); tCorrFctnIter!=fCorrFctnCollection->end(); tCorrFctnIter++) {
          AliFemtoCorrFctn* tCorrFctn = *tCorrFctnIter;
          if (type == "real")
            tCorrFctn->AddRealPair(tPair);
	  else if(type == "mixed")
            tCorrFctn->AddMixedPair(tPair);
          else
            cout << "Problem with pair type, type = " << type.c_str() << endl;
        }
      }

    }    // loop over second particle

  }      // loop over first particle

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
  for (AliFemtoCorrFctnIterator iter=fCorrFctnCollection->begin(); iter!=fCorrFctnCollection->end();iter++) {
    (*iter)->EventBegin(ev);
  }
}
//_________________________
void AliFemtoSimpleAnalysis::EventEnd(const AliFemtoEvent* ev)
{
  /// Fiinsh operations at the end of event processing

  fFirstParticleCut->EventEnd(ev);
  fSecondParticleCut->EventEnd(ev);
  fPairCut->EventEnd(ev);
  for (AliFemtoCorrFctnIterator iter=fCorrFctnCollection->begin(); iter!=fCorrFctnCollection->end();iter++) {
    (*iter)->EventEnd(ev);
  }
}
//_________________________
void AliFemtoSimpleAnalysis::Finish()
{
  /// Perform finishing operations after all events are processed

  AliFemtoCorrFctnIterator iter;
  for (iter=fCorrFctnCollection->begin(); iter!=fCorrFctnCollection->end();iter++) {
    (*iter)->Finish();
  }
}
//_________________________
void AliFemtoSimpleAnalysis::AddEventProcessed()
{
  /// Increase count of processed events

  fNeventsProcessed++;
}
//_________________________
TList* AliFemtoSimpleAnalysis::ListSettings()
{
  /// Collect settings list

  TList *tListSettings = new TList();

  TList *p1Cut = fFirstParticleCut->ListSettings();

  TListIter nextp1(p1Cut);
  while (TObject *obj = nextp1.Next()) {
    TString cuts(obj->GetName());
    cuts.Prepend("AliFemtoSimpleAnalysis.");
    tListSettings->Add(new TObjString(cuts.Data()));
  }

  if (fSecondParticleCut != fFirstParticleCut) {
    TList *p2Cut = fSecondParticleCut->ListSettings();

    TIter nextp2(p2Cut);
    while (TObject *obj = nextp2()) {
      TString cuts(obj->GetName());
      cuts.Prepend("AliFemtoSimpleAnalysis.");
      tListSettings->Add(new TObjString(cuts.Data()));
    }
  }

  TList *pairCut = fPairCut->ListSettings();

  TIter nextpair(pairCut);
  while (TObject *obj = nextpair()) {
    TString cuts(obj->GetName());
    cuts.Prepend("AliFemtoSimpleAnalysis.");
    tListSettings->Add(new TObjString(cuts.Data()));
  }

  return tListSettings;

}

//_________________________
TList* AliFemtoSimpleAnalysis::GetOutputList()
{
  /// Collect the list of output objects
  /// to be written

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

  AliFemtoCorrFctnIterator iter;
  for (iter=fCorrFctnCollection->begin(); iter!=fCorrFctnCollection->end();iter++) {
    TList *tListCf = (*iter)->GetOutputList();

    TIter nextListCf(tListCf);
    while (TObject *obj = nextListCf()) {
      tOutputList->Add(obj);
    }
    delete tListCf;
  }

  return tOutputList;
}

