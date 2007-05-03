/***************************************************************************
 *
 * $Id$
 *
 * Author: Frank Laue, Ohio State, Laue@mps.ohio-state.edu
 ***************************************************************************
 *
 * Description: part of STAR HBT Framework: AliFemtoMaker package
 *      This is the Class for Analysis objects.  Each of the simultaneous
 *      Analyses running should have one of these instantiated.  They link
 *      into the Manager in an Analysis Collection.
 *
 ***************************************************************************/

#include "Analysis/AliFemtoLikeSignAnalysis.h"
#include "Infrastructure/AliFemtoParticleCollection.h"
#include "Base/AliFemtoTrackCut.h"
#include "Base/AliFemtoV0Cut.h"
#include "Infrastructure/AliFemtoPicoEventCollectionVector.h"
#include "Infrastructure/AliFemtoPicoEventCollectionVectorHideAway.h"

#ifdef __ROOT__ 
ClassImp(AliFemtoLikeSignAnalysis)
#endif

// this little function used to apply ParticleCuts (TrackCuts or V0Cuts) and fill ParticleCollections of picoEvent
//  it is called from AliFemtoAnalysis::ProcessEvent()


extern void FillHbtParticleCollection(AliFemtoParticleCut*         partCut,
				     AliFemtoEvent*               hbtEvent,
				     AliFemtoParticleCollection*  partCollection);

 
//____________________________
AliFemtoLikeSignAnalysis::AliFemtoLikeSignAnalysis(unsigned int bins, double min, double max) : 
  AliFemtoAnalysis(),
  fVertexBins(0),
  fOverFlow(0),  
  fUnderFlow(0)  
{
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
  AliFemtoAnalysis(a) ,
  fVertexBins(0),
  fOverFlow(0),  
  fUnderFlow(0)  
{
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
  delete fPicoEventCollectionVectorHideAway; fPicoEventCollectionVectorHideAway=0;
}
//____________________________
AliFemtoString AliFemtoLikeSignAnalysis::Report()
{  
  char Ctemp[200];
  cout << "AliFemtoLikeSignAnalysis - constructing Report..."<<endl;
  AliFemtoString temp = "-----------\nHbt Analysis Report:\n";
  sprintf(Ctemp,"Events are mixed in %d bins in the range %E cm to %E cm.\n",fVertexBins,fVertexZ[0],fVertexZ[1]);
  temp += Ctemp;
  sprintf(Ctemp,"Events underflowing: %d\n",fUnderFlow);
  temp += Ctemp;
  sprintf(Ctemp,"Events overflowing: %d\n",fOverFlow);
  temp += Ctemp;
  sprintf(Ctemp,"Now adding AliFemtoAnalysis(base) Report\n");
  temp += Ctemp; 
  temp += "Adding AliFemtoAnalysis(base) Report now:\n";
  temp += AliFemtoAnalysis::Report();
  temp += "-------------\n";
  AliFemtoString returnThis=temp;
  return returnThis;
}
//_________________________
void AliFemtoLikeSignAnalysis::ProcessEvent(const AliFemtoEvent* hbtEvent) {
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
      AliFemtoPair* ThePair = new AliFemtoPair;
      
      AliFemtoParticleIterator PartIter1;
      AliFemtoParticleIterator PartIter2;
      AliFemtoCorrFctnIterator CorrFctnIter;
      AliFemtoParticleIterator StartOuterLoop = picoEvent->FirstParticleCollection()->begin();  // always
      AliFemtoParticleIterator EndOuterLoop   = picoEvent->FirstParticleCollection()->end();    // will be one less if identical
      AliFemtoParticleIterator StartInnerLoop;
      AliFemtoParticleIterator EndInnerLoop;
      if (AnalyzeIdenticalParticles()) {             // only use First collection
	EndOuterLoop--;                                               // outer loop goes to next-to-last particle in First collection
	EndInnerLoop = picoEvent->FirstParticleCollection()->end() ;  // inner loop goes to last particle in First collection
      }
      else {                                                          // nonidentical - loop over First and Second collections
	StartInnerLoop = picoEvent->SecondParticleCollection()->begin(); // inner loop starts at first particle in Second collection
	EndInnerLoop   = picoEvent->SecondParticleCollection()->end() ;  // inner loop goes to last particle in Second collection
      }
      // real pairs
      for (PartIter1=StartOuterLoop;PartIter1!=EndOuterLoop;PartIter1++){
	if (AnalyzeIdenticalParticles()){
	  StartInnerLoop = PartIter1;
	  StartInnerLoop++;
	}
	ThePair->SetTrack1(*PartIter1);
	for (PartIter2 = StartInnerLoop; PartIter2!=EndInnerLoop;PartIter2++){
	  ThePair->SetTrack2(*PartIter2);
	  // The following lines have to be uncommented if you want pairCutMonitors
	  // they are not in for speed reasons
	  // bool tmpPassPair = mPairCut->Pass(ThePair);
          // mPairCut->FillCutMonitor(ThePair, tmpPassPair);
	  // if ( tmpPassPair ) {
	  if (fPairCut->Pass(ThePair)){
	    for (CorrFctnIter=fCorrFctnCollection->begin();
		 CorrFctnIter!=fCorrFctnCollection->end();CorrFctnIter++){
	      AliFemtoLikeSignCorrFctn* CorrFctn = dynamic_cast<AliFemtoLikeSignCorrFctn*>(*CorrFctnIter);
	      if (CorrFctn) CorrFctn->AddRealPair(ThePair);
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
      prevIter = EndOuterLoop;
      prevIter--;
      for (PartIter1=StartOuterLoop;PartIter1!=prevIter;PartIter1++){
	ThePair->SetTrack1(*PartIter1);
	nextIter = PartIter1;
	nextIter++;
	for (PartIter2 = nextIter; PartIter2!=EndOuterLoop;PartIter2++){
	  ThePair->SetTrack2(*PartIter2);
	  // The following lines have to be uncommented if you want pairCutMonitors
	  // they are not in for speed reasons
	  // bool tmpPassPair = mPairCut->Pass(ThePair);
          // mPairCut->FillCutMonitor(ThePair, tmpPassPair);
	  // if ( tmpPassPair ) {
	  if (fPairCut->Pass(ThePair)){
	    for (CorrFctnIter=fCorrFctnCollection->begin();
		 CorrFctnIter!=fCorrFctnCollection->end();CorrFctnIter++){
	      AliFemtoLikeSignCorrFctn* CorrFctn = dynamic_cast<AliFemtoLikeSignCorrFctn*>(*CorrFctnIter);
	      if (CorrFctn) CorrFctn->AddLikeSignPositivePair(ThePair);
	    }
	  }  // if passed pair cut
	}    // loop over second particle
      }      // loop over first particle
#ifdef STHBTDEBUG
      cout << "AliFemtoLikeSignAnalysis::ProcessEvent() - like sign first collection done" << endl;
#endif
      // like sign second partilce collection pairs
      prevIter = EndInnerLoop;
      prevIter--;
      for (PartIter1=StartInnerLoop;PartIter1!=prevIter;PartIter1++){
	ThePair->SetTrack1(*PartIter1);
	nextIter = PartIter1;
	nextIter++;
	for (PartIter2 = nextIter; PartIter2!=EndInnerLoop;PartIter2++){
	  ThePair->SetTrack2(*PartIter2);
	  // The following lines have to be uncommented if you want pairCutMonitors
	  // they are not in for speed reasons
	  // bool tmpPassPair = mPairCut->Pass(ThePair);
          // mPairCut->FillCutMonitor(ThePair, tmpPassPair);
	  // if ( tmpPassPair ) {
	  if (fPairCut->Pass(ThePair)){
	    for (CorrFctnIter=fCorrFctnCollection->begin();
		 CorrFctnIter!=fCorrFctnCollection->end();CorrFctnIter++){
	      AliFemtoLikeSignCorrFctn* CorrFctn = dynamic_cast<AliFemtoLikeSignCorrFctn*>(*CorrFctnIter);
	      if (CorrFctn) CorrFctn->AddLikeSignNegativePair(ThePair);
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
	StartOuterLoop = picoEvent->FirstParticleCollection()->begin();
	EndOuterLoop   = picoEvent->FirstParticleCollection()->end();
	AliFemtoPicoEvent* storedEvent;
	AliFemtoPicoEventIterator picoEventIter;
	for (picoEventIter=MixingBuffer()->begin();picoEventIter!=MixingBuffer()->end();picoEventIter++){
	  storedEvent = *picoEventIter;
	  if (AnalyzeIdenticalParticles()){
	    StartInnerLoop = storedEvent->FirstParticleCollection()->begin();
	    EndInnerLoop = storedEvent->FirstParticleCollection()->end();
	  }
	  else{
	    StartInnerLoop = storedEvent->SecondParticleCollection()->begin();
	    EndInnerLoop = storedEvent->SecondParticleCollection()->end();
	  }
	  for (PartIter1=StartOuterLoop;PartIter1!=EndOuterLoop;PartIter1++){
	    ThePair->SetTrack1(*PartIter1);
	    for (PartIter2=StartInnerLoop;PartIter2!=EndInnerLoop;PartIter2++){
	      ThePair->SetTrack2(*PartIter2);
	      // testing...	      cout << "ThePair defined... going to pair cut... ";
	      if (fPairCut->Pass(ThePair)){
		// testing...		cout << " ThePair passed PairCut... ";
		for (CorrFctnIter=fCorrFctnCollection->begin();
		     CorrFctnIter!=fCorrFctnCollection->end();CorrFctnIter++){
		  AliFemtoLikeSignCorrFctn* CorrFctn = dynamic_cast<AliFemtoLikeSignCorrFctn*>(*CorrFctnIter);
		  if (CorrFctn) { 
		    CorrFctn->AddMixedPair(ThePair);
		    //cout << " ThePair has been added to MixedPair method " << endl;
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
      delete ThePair;
      MixingBuffer()->push_front(picoEvent);  // store the current pico-event in buffer
    }   // if currentEvent is accepted by currentAnalysis
    EventEnd(hbtEvent);  // cleanup for EbyE 
    //    cout << "AliFemtoLikeSignAnalysis::ProcessEvent() - return to caller ... " << endl;
}



