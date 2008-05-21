///////////////////////////////////////////////////////////////////////////
//                                                                       //
// This is an analysis which calculated the background from like sign    //
// pairs in the same event                                               //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "AliFemtoLikeSignAnalysis.h"
#include "AliFemtoParticleCollection.h"
#include "AliFemtoPicoEventCollectionVector.h"
#include "AliFemtoPicoEventCollectionVectorHideAway.h"

#ifdef __ROOT__ 
ClassImp(AliFemtoLikeSignAnalysis)
#endif

// this little function used to apply ParticleCuts (TrackCuts or V0Cuts) and fill ParticleCollections of picoEvent
//  it is called from AliFemtoSimpleAnalysis::ProcessEvent()


extern void FillHbtParticleCollection(AliFemtoParticleCut*         partCut,
				     AliFemtoEvent*               hbtEvent,
				     AliFemtoParticleCollection*  partCollection);

 
//____________________________
AliFemtoLikeSignAnalysis::AliFemtoLikeSignAnalysis(unsigned int bins, double min, double max) : 
  AliFemtoSimpleAnalysis(),
  fVertexBins(0),
  fOverFlow(0),  
  fUnderFlow(0)  
{
  // standard constructor
  fVertexBins = bins;
  fVertexZ[0] = min;
  fVertexZ[1] = max;
  fUnderFlow = 0; 
  fOverFlow = 0; 
  if (fMixingBuffer) delete fMixingBuffer;
  fPicoEventCollectionVectorHideAway = new AliFemtoPicoEventCollectionVectorHideAway(fVertexBins,fVertexZ[0],fVertexZ[1]);
    /* no-op */
}
//____________________________
AliFemtoLikeSignAnalysis::AliFemtoLikeSignAnalysis(const AliFemtoLikeSignAnalysis& a) : 
  AliFemtoSimpleAnalysis(a) ,
  fVertexBins(0),
  fOverFlow(0),  
  fUnderFlow(0)  
{
  // copy constructor
  fVertexBins = a.fVertexBins; 
  fVertexZ[0] = a.fVertexZ[0]; 
  fVertexZ[1] = a.fVertexZ[1];
  fUnderFlow = 0; 
  fOverFlow = 0; 
  if (fMixingBuffer) delete fMixingBuffer;
  fPicoEventCollectionVectorHideAway = new AliFemtoPicoEventCollectionVectorHideAway(fVertexBins,fVertexZ[0],fVertexZ[1]);
 }
//____________________________ 
AliFemtoLikeSignAnalysis::~AliFemtoLikeSignAnalysis(){
  // destructor
  delete fPicoEventCollectionVectorHideAway; fPicoEventCollectionVectorHideAway=0;
}
//____________________________
AliFemtoString AliFemtoLikeSignAnalysis::Report()
{  
  // prepare report
  char tCtemp[200];
  cout << "AliFemtoLikeSignAnalysis - constructing Report..."<<endl;
  AliFemtoString temp = "-----------\nHbt Analysis Report:\n";
  sprintf(tCtemp,"Events are mixed in %d bins in the range %E cm to %E cm.\n",fVertexBins,fVertexZ[0],fVertexZ[1]);
  temp += tCtemp;
  sprintf(tCtemp,"Events underflowing: %d\n",fUnderFlow);
  temp += tCtemp;
  sprintf(tCtemp,"Events overflowing: %d\n",fOverFlow);
  temp += tCtemp;
  sprintf(tCtemp,"Now adding AliFemtoSimpleAnalysis(base) Report\n");
  temp += tCtemp; 
  temp += "Adding AliFemtoSimpleAnalysis(base) Report now:\n";
  temp += AliFemtoSimpleAnalysis::Report();
  temp += "-------------\n";
  AliFemtoString returnThis=temp;
  return returnThis;
}
//_________________________
void AliFemtoLikeSignAnalysis::ProcessEvent(const AliFemtoEvent* hbtEvent) {
  // perform all the analysis tasks for a single event
  // get right mixing buffer
  double vertexZ = hbtEvent->PrimVertPos().z();
  fMixingBuffer = fPicoEventCollectionVectorHideAway->PicoEventCollection(vertexZ); 
  if (!fMixingBuffer) {
    if ( vertexZ < fVertexZ[0] ) fUnderFlow++;
    if ( vertexZ > fVertexZ[1] ) fOverFlow++;
    return;
  }

  // startup for EbyE 
  EventBegin(hbtEvent);  
  // event cut and event cut monitor
  bool tmpPassEvent = fEventCut->Pass(hbtEvent);
  fEventCut->FillCutMonitor(hbtEvent, tmpPassEvent);
  if (tmpPassEvent) {
      fNeventsProcessed++;
      cout << "AliFemtoLikeSignAnalysis::ProcessEvent() - " << hbtEvent->TrackCollection()->size();
      cout << " #track=" << hbtEvent->TrackCollection()->size();
      // OK, analysis likes the event-- build a pico event from it, using tracks the analysis likes...
      AliFemtoPicoEvent* picoEvent = new AliFemtoPicoEvent;       // this is what we will make pairs from and put in Mixing Buffer
      FillHbtParticleCollection(fFirstParticleCut,(AliFemtoEvent*)hbtEvent,picoEvent->FirstParticleCollection());
      if ( !(AnalyzeIdenticalParticles()) )
	FillHbtParticleCollection(fSecondParticleCut,(AliFemtoEvent*)hbtEvent,picoEvent->SecondParticleCollection());
      cout <<"   #particles in First, Second Collections: " <<
	picoEvent->FirstParticleCollection()->size() << " " <<
	picoEvent->SecondParticleCollection()->size() << endl;
      
      if (picoEvent->SecondParticleCollection()->size()*picoEvent->FirstParticleCollection()->size()==0) {
	delete picoEvent;
	cout << "AliFemtoLikeSignAnalysis - picoEvent deleted due to empty collection " <<endl; 
	return;
      }
      // OK, pico event is built
      // make real pairs...
      
      // Fabrice points out that we do not need to keep creating/deleting pairs all the time
      // We only ever need ONE pair, and we can just keep changing internal pointers
      // this should help speed things up
      AliFemtoPair* tThePair = new AliFemtoPair;
      
      AliFemtoParticleIterator tPartIter1;
      AliFemtoParticleIterator tPartIter2;
      AliFemtoCorrFctnIterator tCorrFctnIter;
      AliFemtoParticleIterator tStartOuterLoop = picoEvent->FirstParticleCollection()->begin();  // always
      AliFemtoParticleIterator tEndOuterLoop   = picoEvent->FirstParticleCollection()->end();    // will be one less if identical
      AliFemtoParticleIterator tStartInnerLoop;
      AliFemtoParticleIterator tEndInnerLoop;
      if (AnalyzeIdenticalParticles()) {             // only use First collection
	tEndOuterLoop--;                                               // outer loop goes to next-to-last particle in First collection
	tEndInnerLoop = picoEvent->FirstParticleCollection()->end() ;  // inner loop goes to last particle in First collection
      }
      else {                                                          // nonidentical - loop over First and Second collections
	tStartInnerLoop = picoEvent->SecondParticleCollection()->begin(); // inner loop starts at first particle in Second collection
	tEndInnerLoop   = picoEvent->SecondParticleCollection()->end() ;  // inner loop goes to last particle in Second collection
      }
      // real pairs
      for (tPartIter1=tStartOuterLoop;tPartIter1!=tEndOuterLoop;tPartIter1++){
	if (AnalyzeIdenticalParticles()){
	  tStartInnerLoop = tPartIter1;
	  tStartInnerLoop++;
	}
	tThePair->SetTrack1(*tPartIter1);
	for (tPartIter2 = tStartInnerLoop; tPartIter2!=tEndInnerLoop;tPartIter2++){
	  tThePair->SetTrack2(*tPartIter2);
	  // The following lines have to be uncommented if you want pairCutMonitors
	  // they are not in for speed reasons
	  // bool tmpPassPair = mPairCut->Pass(tThePair);
          // mPairCut->FillCutMonitor(tThePair, tmpPassPair);
	  // if ( tmpPassPair ) {
	  if (fPairCut->Pass(tThePair)){
	    for (tCorrFctnIter=fCorrFctnCollection->begin();
		 tCorrFctnIter!=fCorrFctnCollection->end();tCorrFctnIter++){
	      AliFemtoLikeSignCorrFctn* tCorrFctn = dynamic_cast<AliFemtoLikeSignCorrFctn*>(*tCorrFctnIter);
	      if (tCorrFctn) tCorrFctn->AddRealPair(tThePair);
	    }
	  }  // if passed pair cut
	}    // loop over second particle
      }      // loop over first particle
#ifdef STHBTDEBUG
      cout << "AliFemtoLikeSignAnalysis::ProcessEvent() - reals done" << endl;
#endif

      AliFemtoParticleIterator nextIter;
      AliFemtoParticleIterator prevIter;

      // like sign first partilce collection pairs
      prevIter = tEndOuterLoop;
      prevIter--;
      for (tPartIter1=tStartOuterLoop;tPartIter1!=prevIter;tPartIter1++){
	tThePair->SetTrack1(*tPartIter1);
	nextIter = tPartIter1;
	nextIter++;
	for (tPartIter2 = nextIter; tPartIter2!=tEndOuterLoop;tPartIter2++){
	  tThePair->SetTrack2(*tPartIter2);
	  // The following lines have to be uncommented if you want pairCutMonitors
	  // they are not in for speed reasons
	  // bool tmpPassPair = mPairCut->Pass(tThePair);
          // mPairCut->FillCutMonitor(tThePair, tmpPassPair);
	  // if ( tmpPassPair ) {
	  if (fPairCut->Pass(tThePair)){
	    for (tCorrFctnIter=fCorrFctnCollection->begin();
		 tCorrFctnIter!=fCorrFctnCollection->end();tCorrFctnIter++){
	      AliFemtoLikeSignCorrFctn* tCorrFctn = dynamic_cast<AliFemtoLikeSignCorrFctn*>(*tCorrFctnIter);
	      if (tCorrFctn) tCorrFctn->AddLikeSignPositivePair(tThePair);
	    }
	  }  // if passed pair cut
	}    // loop over second particle
      }      // loop over first particle
#ifdef STHBTDEBUG
      cout << "AliFemtoLikeSignAnalysis::ProcessEvent() - like sign first collection done" << endl;
#endif
      // like sign second partilce collection pairs
      prevIter = tEndInnerLoop;
      prevIter--;
      for (tPartIter1=tStartInnerLoop;tPartIter1!=prevIter;tPartIter1++){
	tThePair->SetTrack1(*tPartIter1);
	nextIter = tPartIter1;
	nextIter++;
	for (tPartIter2 = nextIter; tPartIter2!=tEndInnerLoop;tPartIter2++){
	  tThePair->SetTrack2(*tPartIter2);
	  // The following lines have to be uncommented if you want pairCutMonitors
	  // they are not in for speed reasons
	  // bool tmpPassPair = mPairCut->Pass(tThePair);
          // mPairCut->FillCutMonitor(tThePair, tmpPassPair);
	  // if ( tmpPassPair ) {
	  if (fPairCut->Pass(tThePair)){
	    for (tCorrFctnIter=fCorrFctnCollection->begin();
		 tCorrFctnIter!=fCorrFctnCollection->end();tCorrFctnIter++){
	      AliFemtoLikeSignCorrFctn* tCorrFctn = dynamic_cast<AliFemtoLikeSignCorrFctn*>(*tCorrFctnIter);
	      if (tCorrFctn) tCorrFctn->AddLikeSignNegativePair(tThePair);
	    }
	  }  // if passed pair cut
	}    // loop over second particle
      }      // loop over first particle
#ifdef STHBTDEBUG
      cout << "AliFemtoLikeSignAnalysis::ProcessEvent() - like sign second collection done" << endl;
#endif
      
      if (MixingBufferFull()){
#ifdef STHBTDEBUG
	cout << "Mixing Buffer is full - lets rock and roll" << endl;
#endif
      }
      else {
	cout << "Mixing Buffer not full -gotta wait " << MixingBuffer()->size() << endl;
      }
      if (MixingBufferFull()){
	tStartOuterLoop = picoEvent->FirstParticleCollection()->begin();
	tEndOuterLoop   = picoEvent->FirstParticleCollection()->end();
	AliFemtoPicoEvent* storedEvent;
	AliFemtoPicoEventIterator picoEventIter;
	for (picoEventIter=MixingBuffer()->begin();picoEventIter!=MixingBuffer()->end();picoEventIter++){
	  storedEvent = *picoEventIter;
	  if (AnalyzeIdenticalParticles()){
	    tStartInnerLoop = storedEvent->FirstParticleCollection()->begin();
	    tEndInnerLoop = storedEvent->FirstParticleCollection()->end();
	  }
	  else{
	    tStartInnerLoop = storedEvent->SecondParticleCollection()->begin();
	    tEndInnerLoop = storedEvent->SecondParticleCollection()->end();
	  }
	  for (tPartIter1=tStartOuterLoop;tPartIter1!=tEndOuterLoop;tPartIter1++){
	    tThePair->SetTrack1(*tPartIter1);
	    for (tPartIter2=tStartInnerLoop;tPartIter2!=tEndInnerLoop;tPartIter2++){
	      tThePair->SetTrack2(*tPartIter2);
	      // testing...	      cout << "tThePair defined... going to pair cut... ";
	      if (fPairCut->Pass(tThePair)){
		// testing...		cout << " tThePair passed PairCut... ";
		for (tCorrFctnIter=fCorrFctnCollection->begin();
		     tCorrFctnIter!=fCorrFctnCollection->end();tCorrFctnIter++){
		  AliFemtoLikeSignCorrFctn* tCorrFctn = dynamic_cast<AliFemtoLikeSignCorrFctn*>(*tCorrFctnIter);
		  if (tCorrFctn) { 
		    tCorrFctn->AddMixedPair(tThePair);
		    //cout << " tThePair has been added to MixedPair method " << endl;
		  }
		}
	      }  // if passed pair cut
	    }    // loop over second particle
	  }      // loop over first particle
	}        // loop over pico-events stored in Mixing buffer
	// Now get rid of oldest stored pico-event in buffer.
	// This means (1) delete the event from memory, (2) "pop" the pointer to it from the MixingBuffer
	delete MixingBuffer()->back();
	MixingBuffer()->pop_back();
      }  // if mixing buffer is full
      delete tThePair;
      MixingBuffer()->push_front(picoEvent);  // store the current pico-event in buffer
    }   // if currentEvent is accepted by currentAnalysis
    EventEnd(hbtEvent);  // cleanup for EbyE 
    //    cout << "AliFemtoLikeSignAnalysis::ProcessEvent() - return to caller ... " << endl;
}



