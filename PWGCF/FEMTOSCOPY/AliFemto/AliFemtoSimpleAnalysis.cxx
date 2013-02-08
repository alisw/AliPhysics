///////////////////////////////////////////////////////////////////////////
//                                                                       //
// AliFemtoSimpleAnalysis - the most basic analysis there is. All other        //
// inherit from this one. Provides basic functionality for the analysis. //
// To properly set up the analysis the following steps should be taken:  //
//                                                                       //
// - create particle cuts and add them via SetFirstParticleCut and       //
//   SetSecondParticleCut. If one analyzes identical particle            //
//   correlations, the first particle cut must be also the second        //
//   particle cut.                                                       //
//                                                                       //
// - create pair cuts and add them via SetPairCut                        //
//                                                                       //
// - create one or many correlation functions and add them via           //
//   AddCorrFctn method.                                                 //
//                                                                       //
// - specify how many events are to be strored in the mixing buffer for  //
//   background construction                                             //
//                                                                       //
// Then, when the analysis is run, for each event, the EventBegin is     //
// called before any processing is done, then the ProcessEvent is called //
// which takes care of creating real and mixed pairs and sending them    //
// to all the registered correlation functions. At the end of each event,//
// after all pairs are processed, EventEnd is called. After the whole    //
// analysis finishes (there is no more events to process) Finish() is    //
// called.                                                               //
//                                                                       //
///////////////////////////////////////////////////////////////////////////
#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoTrackCut.h"
#include "AliFemtoV0Cut.h"
#include "AliFemtoKinkCut.h"
#include <string>
#include <iostream>

// blah blah

#ifdef __ROOT__ 
ClassImp(AliFemtoSimpleAnalysis)
#endif

AliFemtoEventCut*    copyTheCut(AliFemtoEventCut*);
AliFemtoParticleCut* copyTheCut(AliFemtoParticleCut*);
AliFemtoPairCut*     copyTheCut(AliFemtoPairCut*);
AliFemtoCorrFctn*    copyTheCorrFctn(AliFemtoCorrFctn*);

// this little function used to apply ParticleCuts (TrackCuts or V0Cuts) and fill ParticleCollections of picoEvent
//  it is called from AliFemtoSimpleAnalysis::ProcessEvent()
void FillHbtParticleCollection(AliFemtoParticleCut*         partCut,
			       AliFemtoEvent*               hbtEvent,
			       AliFemtoParticleCollection*  partCollection)
{
  // Fill particle collections from the event
  // by the particles that pass all the cuts
  switch (partCut->Type()) {
  case hbtTrack:       // cut is cutting on Tracks
    {
      AliFemtoTrackCut* pCut = (AliFemtoTrackCut*) partCut;
      AliFemtoTrack* pParticle;
      AliFemtoTrackIterator pIter;
      AliFemtoTrackIterator startLoop = hbtEvent->TrackCollection()->begin();
      AliFemtoTrackIterator endLoop   = hbtEvent->TrackCollection()->end();
      for (pIter=startLoop;pIter!=endLoop;pIter++){
	pParticle = *pIter;
	bool tmpPassParticle = pCut->Pass(pParticle);
	pCut->FillCutMonitor(pParticle, tmpPassParticle);
	if (tmpPassParticle){	
	  AliFemtoParticle* particle = new AliFemtoParticle(pParticle,pCut->Mass());
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
      AliFemtoV0Iterator startLoop = hbtEvent->V0Collection()->begin();
      AliFemtoV0Iterator endLoop   = hbtEvent->V0Collection()->end();
      // this following "for" loop is identical to the one above, but because of scoping, I can's see how to avoid repitition...
      for (pIter=startLoop;pIter!=endLoop;pIter++){
	pParticle = *pIter; 
	bool tmpPassV0 = pCut->Pass(pParticle);
	pCut->FillCutMonitor(pParticle,tmpPassV0);
	if (tmpPassV0){
	  AliFemtoParticle* particle = new AliFemtoParticle(pParticle,partCut->Mass());
	  partCollection->push_back(particle);
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
      for (pIter=startLoop;pIter!=endLoop;pIter++){
	pParticle = *pIter; 
	bool tmpPass = pCut->Pass(pParticle);
	pCut->FillCutMonitor(pParticle,tmpPass);
	if (tmpPass){
	  AliFemtoParticle* particle = new AliFemtoParticle(pParticle,partCut->Mass());
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
AliFemtoSimpleAnalysis::AliFemtoSimpleAnalysis() :
  fPicoEventCollectionVectorHideAway(0), 
  fPairCut(0),            
  fCorrFctnCollection(0), 
  fEventCut(0),           
  fFirstParticleCut(0),   
  fSecondParticleCut(0),  
  fMixingBuffer(0),       
  fPicoEvent(0),          
  fNumEventsToMix(0),                     
  fNeventsProcessed(0),                   
  fMinSizePartCollection(0)
{
  // Default constructor
  //  mControlSwitch     = 0;
  fCorrFctnCollection = new AliFemtoCorrFctnCollection;
  fMixingBuffer = new AliFemtoPicoEventCollection;
}
//____________________________

AliFemtoSimpleAnalysis::AliFemtoSimpleAnalysis(const AliFemtoSimpleAnalysis& a) : 
  AliFemtoAnalysis(),
  fPicoEventCollectionVectorHideAway(0), 
  fPairCut(0),            
  fCorrFctnCollection(0), 
  fEventCut(0),           
  fFirstParticleCut(0),   
  fSecondParticleCut(0),  
  fMixingBuffer(0),       
  fPicoEvent(0),          
  fNumEventsToMix(0),                     
  fNeventsProcessed(0),                   
  fMinSizePartCollection(0)
{
  // Copy constructor
  //AliFemtoSimpleAnalysis();
  fCorrFctnCollection = new AliFemtoCorrFctnCollection;
  fMixingBuffer = new AliFemtoPicoEventCollection;

  // find the right event cut
  fEventCut = a.fEventCut->Clone();
  // find the right first particle cut
  fFirstParticleCut = a.fFirstParticleCut->Clone();
  // find the right second particle cut
  if (a.fFirstParticleCut==a.fSecondParticleCut) 
    SetSecondParticleCut(fFirstParticleCut); // identical particle hbt
  else
  fSecondParticleCut = a.fSecondParticleCut->Clone();

  fPairCut = a.fPairCut->Clone();
  
  if ( fEventCut ) {
      SetEventCut(fEventCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoSimpleAnalysis::AliFemtoSimpleAnalysis(const AliFemtoSimpleAnalysis& a) - event cut set " << endl;
  }
  if ( fFirstParticleCut ) {
      SetFirstParticleCut(fFirstParticleCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoSimpleAnalysis::AliFemtoSimpleAnalysis(const AliFemtoSimpleAnalysis& a) - first particle cut set " << endl;
  }
  if ( fSecondParticleCut ) {
      SetSecondParticleCut(fSecondParticleCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoSimpleAnalysis::AliFemtoSimpleAnalysis(const AliFemtoSimpleAnalysis& a) - second particle cut set " << endl;
  }  if ( fPairCut ) {
      SetPairCut(fPairCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoSimpleAnalysis::AliFemtoSimpleAnalysis(const AliFemtoSimpleAnalysis& a) - pair cut set " << endl;
  }

  AliFemtoCorrFctnIterator iter;
  for (iter=a.fCorrFctnCollection->begin(); iter!=a.fCorrFctnCollection->end();iter++){
    cout << " AliFemtoSimpleAnalysis::AliFemtoSimpleAnalysis(const AliFemtoSimpleAnalysis& a) - looking for correlation functions " << endl;
    AliFemtoCorrFctn* fctn = (*iter)->Clone();
    if (fctn) AddCorrFctn(fctn);
    else cout << " AliFemtoSimpleAnalysis::AliFemtoSimpleAnalysis(const AliFemtoSimpleAnalysis& a) - correlation function not found " << endl;
  }

  fNumEventsToMix = a.fNumEventsToMix;

  fMinSizePartCollection = a.fMinSizePartCollection;  // minimum # particles in ParticleCollection

  cout << " AliFemtoSimpleAnalysis::AliFemtoSimpleAnalysis(const AliFemtoSimpleAnalysis& a) - analysis copied " << endl;

}
//____________________________
AliFemtoSimpleAnalysis::~AliFemtoSimpleAnalysis(){
  // destructor
  cout << " AliFemtoSimpleAnalysis::~AliFemtoSimpleAnalysis()" << endl;
  if (fEventCut) delete fEventCut; fEventCut=0;
  if (fFirstParticleCut == fSecondParticleCut) fSecondParticleCut=0;
  if (fFirstParticleCut)  delete fFirstParticleCut; fFirstParticleCut=0;
  if (fSecondParticleCut) delete fSecondParticleCut; fSecondParticleCut=0;
  if (fPairCut) delete fPairCut; fPairCut=0;
  // now delete every CorrFunction in the Collection, and then the Collection itself
  AliFemtoCorrFctnIterator iter;
  for (iter=fCorrFctnCollection->begin(); iter!=fCorrFctnCollection->end();iter++){
    delete *iter;
  }
  delete fCorrFctnCollection;
  // now delete every PicoEvent in the EventMixingBuffer and then the Buffer itself
  if (fMixingBuffer) {
    AliFemtoPicoEventIterator piter;
    for (piter=fMixingBuffer->begin();piter!=fMixingBuffer->end();piter++){
      delete *piter;
    }
    delete fMixingBuffer;
  }
}
//______________________
AliFemtoSimpleAnalysis& AliFemtoSimpleAnalysis::operator=(const AliFemtoSimpleAnalysis& aAna) 
{
  // Assignment operator
  if (this == &aAna)
    return *this;

  if (fCorrFctnCollection) delete fCorrFctnCollection;
  fCorrFctnCollection = new AliFemtoCorrFctnCollection;
  if (fMixingBuffer) delete fMixingBuffer;
  fMixingBuffer = new AliFemtoPicoEventCollection;

  // find the right event cut
  if (fEventCut) delete fEventCut;
  fEventCut = aAna.fEventCut->Clone();
  // find the right first particle cut
  if (fFirstParticleCut) delete fFirstParticleCut;
  fFirstParticleCut = aAna.fFirstParticleCut->Clone();
  // find the right second particle cut
  if (fSecondParticleCut) delete fSecondParticleCut;
  if (aAna.fFirstParticleCut==aAna.fSecondParticleCut) 
    SetSecondParticleCut(fFirstParticleCut); // identical particle hbt
  else
    fSecondParticleCut = aAna.fSecondParticleCut->Clone();

  if (fPairCut) delete fPairCut;
  fPairCut = aAna.fPairCut->Clone();
  
  if ( fEventCut ) {
    SetEventCut(fEventCut); // this will set the myAnalysis pointer inside the cut
  }
  if ( fFirstParticleCut ) {
    SetFirstParticleCut(fFirstParticleCut); // this will set the myAnalysis pointer inside the cut
  }
  if ( fSecondParticleCut ) {
    SetSecondParticleCut(fSecondParticleCut); // this will set the myAnalysis pointer inside the cut
  }  
  if ( fPairCut ) {
    SetPairCut(fPairCut); // this will set the myAnalysis pointer inside the cut
  }

  AliFemtoCorrFctnIterator iter;
  for (iter=aAna.fCorrFctnCollection->begin(); iter!=aAna.fCorrFctnCollection->end();iter++){
    AliFemtoCorrFctn* fctn = (*iter)->Clone();
    if (fctn) AddCorrFctn(fctn);
  }

  fNumEventsToMix = aAna.fNumEventsToMix;

  fMinSizePartCollection = aAna.fMinSizePartCollection;  // minimum # particles in ParticleCollection

  return *this;
}
//______________________
AliFemtoCorrFctn* AliFemtoSimpleAnalysis::CorrFctn(int n){  
  // return pointer to n-th correlation function
  if ( n<0 || n > (int)fCorrFctnCollection->size() )
    return NULL;
  AliFemtoCorrFctnIterator iter=fCorrFctnCollection->begin();
  for (int i=0; i<n ;i++){
    iter++;
  }
  return *iter;
}
//____________________________
AliFemtoString AliFemtoSimpleAnalysis::Report()
{
  // Create a simple report from the analysis execution
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
  AliFemtoCorrFctnIterator iter;
  if ( fCorrFctnCollection->size()==0 ) {
    cout << "AliFemtoSimpleAnalysis-Warning : no correlations functions in this analysis " << endl;
  }
  for (iter=fCorrFctnCollection->begin(); iter!=fCorrFctnCollection->end();iter++){
    temp += (*iter)->Report();
    temp += "\n";
  }
  temp += "-------------\n";
  AliFemtoString returnThis=temp;
  return returnThis;
}
//_________________________
void AliFemtoSimpleAnalysis::ProcessEvent(const AliFemtoEvent* hbtEvent) {
  // Add event to processed events

  fPicoEvent=0; // we will get a new pico event, if not prevent corr. fctn to access old pico event
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
    FillHbtParticleCollection(fFirstParticleCut,(AliFemtoEvent*)hbtEvent,fPicoEvent->FirstParticleCollection());
    if ( !(AnalyzeIdenticalParticles()) )
      FillHbtParticleCollection(fSecondParticleCut,(AliFemtoEvent*)hbtEvent,fPicoEvent->SecondParticleCollection());
    //cout <<"AliFemtoSimpleAnalysis::ProcessEvent - #particles in First, Second Collections: " <<
//       fPicoEvent->FirstParticleCollection()->size() << " " <<
//       fPicoEvent->SecondParticleCollection()->size() << endl;
    
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
        MakePairs("real", fPicoEvent->FirstParticleCollection() );
      }
      else {
        MakePairs("real", fPicoEvent->FirstParticleCollection(),
                          fPicoEvent->SecondParticleCollection() );
      }
       cout << "AliFemtoSimpleAnalysis::ProcessEvent() - reals done ";

      //---- Make pairs for mixed events, looping over events in mixingBuffer ----//

      AliFemtoPicoEvent* storedEvent;
      AliFemtoPicoEventIterator fPicoEventIter;
      for (fPicoEventIter=MixingBuffer()->begin();fPicoEventIter!=MixingBuffer()->end();fPicoEventIter++){
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
    else{
      fEventCut->FillCutMonitor(hbtEvent, !tmpPassEvent);
      delete fPicoEvent;
    }
  }   // if currentEvent is accepted by currentAnalysis
  EventEnd(hbtEvent);  // cleanup for EbyE 
  //cout << "AliFemtoSimpleAnalysis::ProcessEvent() - return to caller ... " << endl;
}
//_________________________
void AliFemtoSimpleAnalysis::MakePairs(const char* typeIn, AliFemtoParticleCollection *partCollection1,
				       AliFemtoParticleCollection *partCollection2){
// Build pairs, check pair cuts, and call CFs' AddRealPair() or
// AddMixedPair() methods. If no second particle collection is
// specfied, make pairs within first particle collection.

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
  if (partCollection2) {                        // Two collections:
    tStartInnerLoop = partCollection2->begin();  //   Full inner & outer loops
    tEndInnerLoop   = partCollection2->end();    //
  }
  else {                                        // One collection:
    tEndOuterLoop--;                             //   Outer loop goes to next-to-last particle
    tEndInnerLoop = partCollection1->end() ;     //   Inner loop goes to last particle
  }
  for (tPartIter1=tStartOuterLoop;tPartIter1!=tEndOuterLoop;tPartIter1++) {
    if (!partCollection2){
      tStartInnerLoop = tPartIter1;
      tStartInnerLoop++;
    }
    tPair->SetTrack1(*tPartIter1);
    for (tPartIter2 = tStartInnerLoop; tPartIter2!=tEndInnerLoop;tPartIter2++) {
      tPair->SetTrack2(*tPartIter2);

      // The following lines have to be uncommented if you want pairCutMonitors
      // they are not in for speed reasons
      // bool tmpPassPair = fPairCut->Pass(tPair);
      // fPairCut->FillCutMonitor(tPair, tmpPassPair);
      // if ( tmpPassPair )

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

      if (fPairCut->Pass(tPair)){
        for (tCorrFctnIter=fCorrFctnCollection->begin();
             tCorrFctnIter!=fCorrFctnCollection->end();tCorrFctnIter++){
          AliFemtoCorrFctn* tCorrFctn = *tCorrFctnIter;
          if(type == "real")
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
void AliFemtoSimpleAnalysis::EventBegin(const AliFemtoEvent* ev){
  // Perform initialization operations at the beginning of the event processing
  //cout << " AliFemtoSimpleAnalysis::EventBegin(const AliFemtoEvent* ev) " << endl;
  fFirstParticleCut->EventBegin(ev);
  fSecondParticleCut->EventBegin(ev);
  fPairCut->EventBegin(ev);
  for (AliFemtoCorrFctnIterator iter=fCorrFctnCollection->begin(); iter!=fCorrFctnCollection->end();iter++){
    (*iter)->EventBegin(ev);
  }
}
//_________________________
void AliFemtoSimpleAnalysis::EventEnd(const AliFemtoEvent* ev){
  // Fiinsh operations at the end of event processing
  fFirstParticleCut->EventEnd(ev);
  fSecondParticleCut->EventEnd(ev);
  fPairCut->EventEnd(ev);
  for (AliFemtoCorrFctnIterator iter=fCorrFctnCollection->begin(); iter!=fCorrFctnCollection->end();iter++){
    (*iter)->EventEnd(ev);
  }
}
//_________________________
void AliFemtoSimpleAnalysis::Finish(){
  // Perform finishing operations after all events are processed
  AliFemtoCorrFctnIterator iter;
  for (iter=fCorrFctnCollection->begin(); iter!=fCorrFctnCollection->end();iter++){
    (*iter)->Finish();
  }
}
//_________________________
void AliFemtoSimpleAnalysis::AddEventProcessed() {
  // Increase count of processed events
  fNeventsProcessed++;
}
//_________________________
TList* AliFemtoSimpleAnalysis::ListSettings()
{
  // Collect settings list
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
  // Collect the list of output objects
  // to be written 
  TList *tOutputList = new TList();

  TList *p1Cut = fFirstParticleCut->GetOutputList();

  TListIter nextp1(p1Cut);
  while (TObject *obj = nextp1.Next()) {
    tOutputList->Add(obj);
  }

  if (fSecondParticleCut != fFirstParticleCut) {
    TList *p2Cut = fSecondParticleCut->GetOutputList();
    
    TIter nextp2(p2Cut);
    while (TObject *obj = nextp2()) {
      tOutputList->Add(obj);
    }
  }

  TList *pairCut = fPairCut->GetOutputList();

  TIter nextpair(pairCut);
  while (TObject *obj = nextpair()) {
    tOutputList->Add(obj);
  }

  TList *eventCut = fEventCut->GetOutputList();

  TIter nextevent(eventCut);
  while (TObject *obj = nextevent()) {
    tOutputList->Add(obj);
  }

  AliFemtoCorrFctnIterator iter;
  for (iter=fCorrFctnCollection->begin(); iter!=fCorrFctnCollection->end();iter++){
    TList *tListCf = (*iter)->GetOutputList();
    
    TIter nextListCf(tListCf);
    while (TObject *obj = nextListCf()) {
      tOutputList->Add(obj);
    }
  }

  return tOutputList;
  
}
