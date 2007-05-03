/***************************************************************************
 * Author: Mike Lisa, Ohio State, lisa@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *      This is the Class for Analysis objects.  Each of the simultaneous
 *      Analyses running should have one of these instantiated.  They link
 *      into the Manager in an Analysis Collection.
 *
 ***************************************************************************
 *
 *
 * Revision 1.23  2002/11/20 00:09:26  renault
 * fill a new monitor with (hbtEvent,partCollection)
 *
 * Revision 1.22  2002/11/03 16:40:31  magestro
 * Modified ProcessEvent(), added MakePairs() method, and implemented immediate event mixing
 *
 * Revision 1.21  2002/06/26 17:27:09  lisa
 * fixed small bug in AliFemtoAnalysis associated with the new feature to require ParticleCollections to have some minimum number of particles
 *
 * Revision 1.20  2002/06/22 17:53:31  lisa
 * implemented switch to allow user to require minimum number of particles in First and Second ParticleCollections - default value is zero so if user does not Set this value then behaviour is like before
 *
 * Revision 1.19  2001/11/06 20:20:53  laue
 * Order of event-mixing fixed.
 *
 * Revision 1.18  2001/05/25 23:23:59  lisa
 * Added in AliFemtoKink stuff
 *
 * Revision 1.17  2001/04/05 21:57:45  laue
 * current pico-event becomes a member of the analysis (fPicoEvent) and gets
 * an access-function (CurrentPicoEvent)
 *
 * Revision 1.15  2000/09/13 18:09:09  laue
 * Bux fix: Delete track cut only once for identical particle hbt
 *
 * Revision 1.14  2000/08/31 22:31:30  laue
 * AliFemtoAnalysis: output changed (a little bit less)
 * AliFemtoEvent: new version, members for reference mult added
 * AliFemtoIOBinary: new IO for new AliFemtoEvent version
 * AliFemtoTypes: TTree typedef to AliFemtoTTree added
 * AliFemtoVertexAnalysis: overflow and underflow added
 *
 * Revision 1.13  2000/08/11 16:35:40  rcwells
 * Added number of events processed to each HBT analysis
 *
 * Revision 1.12  2000/07/16 22:23:17  laue
 * I forgot that we moved memberfunctions out of AliFemtoBaseAnalysis.
 * So my previous check-ins didn't compile with the library.
 * Now they do.
 *
 * Revision 1.11  2000/07/16 21:38:22  laue
 * AliFemtoCoulomb.cxx AliFemtoSectoredAnalysis.cxx : updated for standalone version
 * AliFemtoV0.cc AliFemtoV0.h : some cast to prevent compiling warnings
 * AliFemtoParticle.cc AliFemtoParticle.h : pointers mTrack,mV0 initialized to 0
 * AliFemtoIOBinary.cc : some printouts in #ifdef STHBTDEBUG
 * AliFemtoEvent.cc : B-Field set to 0.25Tesla, we have to think about a better
 *                 solution
 *
 * Revision 1.10  2000/07/06 18:45:51  laue
 * Copy constructor fixed. It couldn't handle identicle particles.
 *
 * Revision 1.9  2000/04/13 20:20:22  laue
 * Event mixing corrected. Now the first collection of the current event is
 * mixed with the second collection from the mixing buffer _AND_ vice verse
 *
 * Revision 1.8  2000/03/23 23:00:01  laue
 * AliFemtoAnalysis copy constructor now uses Clone() function of cuts
 * AliFemtoTypes now has AliFemtoTF1 for fitting purposes
 *
 * Revision 1.7  2000/03/17 17:23:05  laue
 * Roberts new three particle correlations implemented.
 *
 * Revision 1.6  2000/03/16 02:07:04  laue
 * Copy constructor added to AliFemtoAnalysis (only known cuts, corrfctn).
 *
 * AliFemtoBinaryReader can now derive filename from StIOMaker and read a list
 * of files.
 *
 * AliFemtoManager now holds a collection of AliFemtoEventWriters (multiple writes
 * possible now)
 *
 * Revision 1.5  2000/02/13 17:17:12  laue
 * Calls to the EventBegin() and EventEnd() functions implemented
 * The actual analysis is moved from AliFemtoManager to AliFemtoAnalysis
 *
 * Revision 1.4  2000/01/25 17:35:16  laue
 * I. In order to run the stand alone version of the AliFemtoMaker the following
 * changes have been done:
 * a) all ClassDefs and ClassImps have been put into #ifdef __ROOT__ statements
 * b) unnecessary includes of StMaker.h have been removed
 * c) the subdirectory AliFemtoMaker/doc/Make has been created including everything
 * needed for the stand alone version
 *
 * II. To reduce the amount of compiler warning
 * a) some variables have been type casted
 * b) some destructors have been declared as virtual
 *
 * Revision 1.3  1999/10/04 15:38:53  lisa
 * include Franks new accessor methods AliFemtoAnalysis::CorrFctn and AliFemtoManager::Analysis as well as McEvent example macro
 *
 * Revision 1.2  1999/07/06 22:33:22  lisa
 * Adjusted all to work in pro and new - dev itself is broken
 *
 * Revision 1.1.1.1  1999/06/29 16:02:57  lisa
 * Installation of AliFemtoMaker
 *
 **************************************************************************/

#include "Analysis/AliFemtoAnalysis.h"
#include "Base/AliFemtoTrackCut.h"
#include "Base/AliFemtoV0Cut.h"
#include "Base/AliFemtoKinkCut.h"
#include <string>


#ifdef __ROOT__ 
ClassImp(AliFemtoAnalysis)
#endif

AliFemtoEventCut*    copyTheCut(AliFemtoEventCut*);
AliFemtoParticleCut* copyTheCut(AliFemtoParticleCut*);
AliFemtoPairCut*     copyTheCut(AliFemtoPairCut*);
AliFemtoCorrFctn*    copyTheCorrFctn(AliFemtoCorrFctn*);

// this little function used to apply ParticleCuts (TrackCuts or V0Cuts) and fill ParticleCollections of picoEvent
//  it is called from AliFemtoAnalysis::ProcessEvent()
void FillHbtParticleCollection(AliFemtoParticleCut*         partCut,
			       AliFemtoEvent*               hbtEvent,
			       AliFemtoParticleCollection*  partCollection)
{
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
    cout << "FillHbtParticleCollection function (in AliFemtoAnalysis.cxx) - undefined Particle Cut type!!! \n";
  }
}
//____________________________
AliFemtoAnalysis::AliFemtoAnalysis() :
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
  //  mControlSwitch     = 0;
  fCorrFctnCollection = new AliFemtoCorrFctnCollection;
  fMixingBuffer = new AliFemtoPicoEventCollection;
}
//____________________________

AliFemtoAnalysis::AliFemtoAnalysis(const AliFemtoAnalysis& a) : 
  AliFemtoBaseAnalysis(),
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
  //AliFemtoAnalysis();
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
      cout << " AliFemtoAnalysis::AliFemtoAnalysis(const AliFemtoAnalysis& a) - event cut set " << endl;
  }
  if ( fFirstParticleCut ) {
      SetFirstParticleCut(fFirstParticleCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoAnalysis::AliFemtoAnalysis(const AliFemtoAnalysis& a) - first particle cut set " << endl;
  }
  if ( fSecondParticleCut ) {
      SetSecondParticleCut(fSecondParticleCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoAnalysis::AliFemtoAnalysis(const AliFemtoAnalysis& a) - second particle cut set " << endl;
  }  if ( fPairCut ) {
      SetPairCut(fPairCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoAnalysis::AliFemtoAnalysis(const AliFemtoAnalysis& a) - pair cut set " << endl;
  }

  AliFemtoCorrFctnIterator iter;
  for (iter=a.fCorrFctnCollection->begin(); iter!=a.fCorrFctnCollection->end();iter++){
    cout << " AliFemtoAnalysis::AliFemtoAnalysis(const AliFemtoAnalysis& a) - looking for correlation functions " << endl;
    AliFemtoCorrFctn* fctn = (*iter)->Clone();
    if (fctn) AddCorrFctn(fctn);
    else cout << " AliFemtoAnalysis::AliFemtoAnalysis(const AliFemtoAnalysis& a) - correlation function not found " << endl;
  }

  fNumEventsToMix = a.fNumEventsToMix;

  fMinSizePartCollection = a.fMinSizePartCollection;  // minimum # particles in ParticleCollection

  cout << " AliFemtoAnalysis::AliFemtoAnalysis(const AliFemtoAnalysis& a) - analysis copied " << endl;

}
//____________________________
AliFemtoAnalysis::~AliFemtoAnalysis(){
  cout << " AliFemtoAnalysis::~AliFemtoAnalysis()" << endl;
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
AliFemtoAnalysis& AliFemtoAnalysis::operator=(const AliFemtoAnalysis& aAna) 
{
  if (this == &aAna)
    return *this;

  fCorrFctnCollection = new AliFemtoCorrFctnCollection;
  fMixingBuffer = new AliFemtoPicoEventCollection;

  // find the right event cut
  fEventCut = aAna.fEventCut->Clone();
  // find the right first particle cut
  fFirstParticleCut = aAna.fFirstParticleCut->Clone();
  // find the right second particle cut
  if (aAna.fFirstParticleCut==aAna.fSecondParticleCut) 
    SetSecondParticleCut(fFirstParticleCut); // identical particle hbt
  else
    fSecondParticleCut = aAna.fSecondParticleCut->Clone();

  fPairCut = aAna.fPairCut->Clone();
  
  if ( fEventCut ) {
      SetEventCut(fEventCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoAnalysis::AliFemtoAnalysis(const AliFemtoAnalysis& a) - event cut set " << endl;
  }
  if ( fFirstParticleCut ) {
      SetFirstParticleCut(fFirstParticleCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoAnalysis::AliFemtoAnalysis(const AliFemtoAnalysis& a) - first particle cut set " << endl;
  }
  if ( fSecondParticleCut ) {
      SetSecondParticleCut(fSecondParticleCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoAnalysis::AliFemtoAnalysis(const AliFemtoAnalysis& a) - second particle cut set " << endl;
  }  if ( fPairCut ) {
      SetPairCut(fPairCut); // this will set the myAnalysis pointer inside the cut
      cout << " AliFemtoAnalysis::AliFemtoAnalysis(const AliFemtoAnalysis& a) - pair cut set " << endl;
  }

  AliFemtoCorrFctnIterator iter;
  for (iter=aAna.fCorrFctnCollection->begin(); iter!=aAna.fCorrFctnCollection->end();iter++){
    cout << " AliFemtoAnalysis::AliFemtoAnalysis(const AliFemtoAnalysis& a) - looking for correlation functions " << endl;
    AliFemtoCorrFctn* fctn = (*iter)->Clone();
    if (fctn) AddCorrFctn(fctn);
    else cout << " AliFemtoAnalysis::AliFemtoAnalysis(const AliFemtoAnalysis& a) - correlation function not found " << endl;
  }

  fNumEventsToMix = aAna.fNumEventsToMix;

  fMinSizePartCollection = aAna.fMinSizePartCollection;  // minimum # particles in ParticleCollection

  cout << " AliFemtoAnalysis::AliFemtoAnalysis(const AliFemtoAnalysis& a) - analysis copied " << endl;

  return *this;
}
//______________________
AliFemtoCorrFctn* AliFemtoAnalysis::CorrFctn(int n){  // return pointer to n-th correlation function
  if ( n<0 || n > (int)fCorrFctnCollection->size() )
    return NULL;
  AliFemtoCorrFctnIterator iter=fCorrFctnCollection->begin();
  for (int i=0; i<n ;i++){
    iter++;
  }
  return *iter;
}
//____________________________
AliFemtoString AliFemtoAnalysis::Report()
{
  cout << "AliFemtoAnalysis - constructing Report..."<<endl;
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
    cout << "AliFemtoAnalysis-Warning : no correlations functions in this analysis " << endl;
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
void AliFemtoAnalysis::ProcessEvent(const AliFemtoEvent* hbtEvent) {
  // Add event to processed events
  fPicoEvent=0; // we will get a new pico event, if not prevent corr. fctn to access old pico event
  AddEventProcessed();
  // startup for EbyE 
  EventBegin(hbtEvent);  
  // event cut and event cut monitor
  bool tmpPassEvent = fEventCut->Pass(hbtEvent);
  fEventCut->FillCutMonitor(hbtEvent, tmpPassEvent);
  if (tmpPassEvent) {
    cout << "AliFemtoAnalysis::ProcessEvent() - Event has passed cut - build picoEvent from " <<
      hbtEvent->TrackCollection()->size() << " tracks in TrackCollection" << endl;
    // OK, analysis likes the event-- build a pico event from it, using tracks the analysis likes...
    fPicoEvent = new AliFemtoPicoEvent; // this is what we will make pairs from and put in Mixing Buffer
    // no memory leak. we will delete picoevents when they come out of the mixing buffer
    FillHbtParticleCollection(fFirstParticleCut,(AliFemtoEvent*)hbtEvent,fPicoEvent->FirstParticleCollection());
    if ( !(AnalyzeIdenticalParticles()) )
      FillHbtParticleCollection(fSecondParticleCut,(AliFemtoEvent*)hbtEvent,fPicoEvent->SecondParticleCollection());
    cout <<"AliFemtoAnalysis::ProcessEvent - #particles in First, Second Collections: " <<
      fPicoEvent->FirstParticleCollection()->size() << " " <<
      fPicoEvent->SecondParticleCollection()->size() << endl;
    
    // mal - implement a switch which allows only using events with ParticleCollections containing a minimum
    // number of entries (jun2002)
    if ((fPicoEvent->FirstParticleCollection()->size() >= fMinSizePartCollection )
	&& ( AnalyzeIdenticalParticles() || (fPicoEvent->SecondParticleCollection()->size() >= fMinSizePartCollection ))) {


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
      cout << "AliFemtoAnalysis::ProcessEvent() - reals done ";

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
      delete fPicoEvent;
    }
  }   // if currentEvent is accepted by currentAnalysis
  EventEnd(hbtEvent);  // cleanup for EbyE 
  //cout << "AliFemtoAnalysis::ProcessEvent() - return to caller ... " << endl;
}
//_________________________
void AliFemtoAnalysis::MakePairs(const char* typeIn, AliFemtoParticleCollection *partCollection1,
                                            AliFemtoParticleCollection *partCollection2){
// Build pairs, check pair cuts, and call CFs' AddRealPair() or
// AddMixedPair() methods. If no second particle collection is
// specfied, make pairs within first particle collection.

  string type = typeIn;

  AliFemtoPair* ThePair = new AliFemtoPair;

  AliFemtoCorrFctnIterator CorrFctnIter;

  AliFemtoParticleIterator PartIter1, PartIter2;

  AliFemtoParticleIterator StartOuterLoop = partCollection1->begin();  // always
  AliFemtoParticleIterator EndOuterLoop   = partCollection1->end();    // will be one less if identical
  AliFemtoParticleIterator StartInnerLoop;
  AliFemtoParticleIterator EndInnerLoop;
  if (partCollection2) {                        // Two collections:
    StartInnerLoop = partCollection2->begin();  //   Full inner & outer loops
    EndInnerLoop   = partCollection2->end();    //
  }
  else {                                        // One collection:
    EndOuterLoop--;                             //   Outer loop goes to next-to-last particle
    EndInnerLoop = partCollection1->end() ;     //   Inner loop goes to last particle
  }
  for (PartIter1=StartOuterLoop;PartIter1!=EndOuterLoop;PartIter1++) {
    if (!partCollection2){
      StartInnerLoop = PartIter1;
      StartInnerLoop++;
    }
    ThePair->SetTrack1(*PartIter1);
    for (PartIter2 = StartInnerLoop; PartIter2!=EndInnerLoop;PartIter2++) {
      ThePair->SetTrack2(*PartIter2);

      // The following lines have to be uncommented if you want pairCutMonitors
      // they are not in for speed reasons
      // bool tmpPassPair = fPairCut->Pass(ThePair);
      // fPairCut->FillCutMonitor(ThePair, tmpPassPair);
      // if ( tmpPassPair )

      //---- If pair passes cut, loop over CF's and add pair to real/mixed ----//

      if (fPairCut->Pass(ThePair)){
        for (CorrFctnIter=fCorrFctnCollection->begin();
             CorrFctnIter!=fCorrFctnCollection->end();CorrFctnIter++){
          AliFemtoCorrFctn* CorrFctn = *CorrFctnIter;
          if(type == "real")
            CorrFctn->AddRealPair(ThePair);
	  else if(type == "mixed")
            CorrFctn->AddMixedPair(ThePair);
          else
            cout << "Problem with pair type, type = " << type.c_str() << endl;
        }
      }

    }    // loop over second particle

  }      // loop over first particle

  delete ThePair;

}
//_________________________
void AliFemtoAnalysis::EventBegin(const AliFemtoEvent* ev){
  //cout << " AliFemtoAnalysis::EventBegin(const AliFemtoEvent* ev) " << endl;
  fFirstParticleCut->EventBegin(ev);
  fSecondParticleCut->EventBegin(ev);
  fPairCut->EventBegin(ev);
  for (AliFemtoCorrFctnIterator iter=fCorrFctnCollection->begin(); iter!=fCorrFctnCollection->end();iter++){
    (*iter)->EventBegin(ev);
  }
}
//_________________________
void AliFemtoAnalysis::EventEnd(const AliFemtoEvent* ev){
  fFirstParticleCut->EventEnd(ev);
  fSecondParticleCut->EventEnd(ev);
  fPairCut->EventEnd(ev);
  for (AliFemtoCorrFctnIterator iter=fCorrFctnCollection->begin(); iter!=fCorrFctnCollection->end();iter++){
    (*iter)->EventEnd(ev);
  }
}
//_________________________
void AliFemtoAnalysis::Finish(){
  AliFemtoCorrFctnIterator iter;
  for (iter=fCorrFctnCollection->begin(); iter!=fCorrFctnCollection->end();iter++){
    (*iter)->Finish();
  }
}
//_________________________
void AliFemtoAnalysis::AddEventProcessed() {
  fNeventsProcessed++;
}
