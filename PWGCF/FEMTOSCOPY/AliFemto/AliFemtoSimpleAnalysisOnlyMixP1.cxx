//==============================================================================//
// dongfang.wang@cern.ch                                                        //
// AliFemtoSimpleAnalysisOnlyMixP1: Only keep update/save first particle        //
// The key is in line 200(if current event has first particle, but not have     //
// second particle, then save it)                                               //
// and line 260(Only allow first particle in mixing pool making pair            //
// with second particle in current event                                        //
//==============================================================================//

#include "AliFemtoSimpleAnalysis.h"
#include "AliFemtoSimpleAnalysisOnlyMixP1.h"
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
  ClassImp(AliFemtoSimpleAnalysisOnlyMixP1);
  /// \endcond
#endif


template <class TrackCollectionType, class TrackCutType>
void DoFillParticleCollectionOnlyMixP1(TrackCutType *cut,
                              TrackCollectionType *track_collection,
                              AliFemtoParticleCollection *output)
{
  for (const auto &track : *track_collection) {
    const Bool_t track_passes = cut->Pass(track);
    cut->FillCutMonitor(track, track_passes);
    if (track_passes) {
      output->push_back(new AliFemtoParticle(track, cut->Mass()));
    }
  }
}


// This little function is used to apply ParticleCuts (TrackCuts or V0Cuts) and
// fill ParticleCollections from tacks in picoEvent. It is called from
// AliFemtoSimpleAnalysis::ProcessEvent().
//
// The actual loop implementation has been moved to the collection-generic
// DoFillParticleCollection() function
void FillHbtParticleCollectionOnlyMixP1(AliFemtoParticleCut *partCut,
                               const AliFemtoEvent *hbtEvent,
                               AliFemtoParticleCollection *partCollection,
                               bool performSharedDaughterCut=kFALSE)
{
  /// Fill particle collection with all particles in the event which pass
  /// the provided cut

  // determine which track collection to use based on the particle type.
  switch (partCut->Type()) {

  // cut is cutting on Tracks
  case hbtTrack:
    {
      DoFillParticleCollectionOnlyMixP1(
			       (AliFemtoTrackCut*)partCut,
			       hbtEvent->TrackCollection(),
			       partCollection
			       );
    }
    break;

  // cut is cutting on V0s
  case hbtV0:
  {
    AliFemtoV0Cut *v0_cut = (AliFemtoV0Cut*)partCut;

    // shared daughter cut returns all passed v0s - add to particle collection
    if (performSharedDaughterCut) {
      AliFemtoV0SharedDaughterCut shared_daughter_cut;
      AliFemtoV0Collection v0_coll = shared_daughter_cut.AliFemtoV0SharedDaughterCutCollection(hbtEvent->V0Collection(), v0_cut);
      // for (AliFemtoV0Iterator pIter = v0_coll.begin(); pIter != v0_coll.end(); ++pIter) {
      for (auto v0 : v0_coll) {
        partCollection->push_back(new AliFemtoParticle(v0, v0_cut->Mass()));
      }
    } else {

      DoFillParticleCollectionOnlyMixP1(
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
      DoFillParticleCollectionOnlyMixP1(
        (AliFemtoXiTrackCut*)partCut,
        hbtEvent->XiCollection(),
        partCollection
      );
    }
    break;
  }

  // cut is cutting on Kinks
  case hbtKink:

    DoFillParticleCollectionOnlyMixP1(
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

// Leave this here to appease any legacy code that expected a non-const AliFemtoEvent
void FillHbtParticleCollectionOnlyMixP1(AliFemtoParticleCut *partCut,
                               AliFemtoEvent *hbtEvent,
                               AliFemtoParticleCollection *partCollection,
                               bool performSharedDaughterCut)
{
  FillHbtParticleCollectionOnlyMixP1(partCut, const_cast<const AliFemtoEvent*>(hbtEvent), partCollection, performSharedDaughterCut);
}


AliFemtoSimpleAnalysisOnlyMixP1::AliFemtoSimpleAnalysisOnlyMixP1():
    AliFemtoSimpleAnalysis()
{
	MakeCForNot = 1;
	OnlyP1Exist = 1;	
}
AliFemtoSimpleAnalysisOnlyMixP1::AliFemtoSimpleAnalysisOnlyMixP1(const AliFemtoSimpleAnalysisOnlyMixP1 &OriAnalysis):
    AliFemtoSimpleAnalysis(OriAnalysis)
{
    //copy constructor 
    MakeCForNot = OriAnalysis.MakeCForNot;
    OnlyP1Exist = OriAnalysis.OnlyP1Exist;
}
AliFemtoSimpleAnalysisOnlyMixP1::~AliFemtoSimpleAnalysisOnlyMixP1()
{
  // Destructor
}

AliFemtoSimpleAnalysisOnlyMixP1& AliFemtoSimpleAnalysisOnlyMixP1::operator=(const AliFemtoSimpleAnalysisOnlyMixP1& OriAnalysis)
{
    // assignment operator
    if (this == &OriAnalysis)
        return *this;

    AliFemtoSimpleAnalysisOnlyMixP1::operator=(OriAnalysis);
    MakeCForNot = OriAnalysis.MakeCForNot;
    OnlyP1Exist = OriAnalysis.OnlyP1Exist;
    return *this;
}

void AliFemtoSimpleAnalysisOnlyMixP1::ProcessEvent(const AliFemtoEvent* hbtEvent){
    fPicoEvent = nullptr;
    AddEventProcessed();
    EventBegin(hbtEvent);
    bool tmpPassEvent = fEventCut->Pass(hbtEvent);
    if (!tmpPassEvent) {
        //fEventCut->FillCutMonitor(hbtEvent, tmpPassEvent);
        EventEnd(hbtEvent);  // cleanup for EbyE
        return;
    }
    fPicoEvent = new AliFemtoPicoEvent;

    
     
    AliFemtoParticleCollection *collection1 = fPicoEvent->FirstParticleCollection(),
                               *collection2 = fPicoEvent->SecondParticleCollection();
    // only collection1 == nullptr, we delete this event!
    if ((collection1 == nullptr) && (OnlyP1Exist == 1)){
	fEventCut->FillCutMonitor(hbtEvent, false);
        EventEnd(hbtEvent);  // cleanup for EbyE
        delete fPicoEvent;
        return;
    }
    // dowang 2022.1.26
    // only for side band cut!
    if (OnlyP1Exist == 0){
    	if((collection1 == nullptr) && (collection2 == nullptr)){
		fEventCut->FillCutMonitor(hbtEvent, false);
		EventEnd(hbtEvent);  // cleanup for EbyE
	        delete fPicoEvent;
        	return;
	}
    }
    FillHbtParticleCollectionOnlyMixP1(fFirstParticleCut,
                            hbtEvent,
                            fPicoEvent->FirstParticleCollection(),
                            fPerformSharedDaughterCut);

    // fill second particle cut if not analyzing identical particles
    if ( !AnalyzeIdenticalParticles() ) {
        FillHbtParticleCollectionOnlyMixP1(fSecondParticleCut,
                                hbtEvent,
                                fPicoEvent->SecondParticleCollection(),
                                fPerformSharedDaughterCut);
    }

    const UInt_t coll_1_size = collection1->size(),
                 coll_2_size = collection2->size();
    fEventCut->FillCutMonitor(collection1, collection2); //MJ!
    
    const bool coll_1_size_passes = (coll_1_size >= fMinSizePartCollection),
               coll_2_size_passes = (AnalyzeIdenticalParticles() || (coll_2_size >= fMinSizePartCollection));

    tmpPassEvent = tmpPassEvent && coll_1_size_passes && coll_2_size_passes;
    // fill the event cut monitor
    fEventCut->FillCutMonitor(hbtEvent, tmpPassEvent);

    //\ dowang 1.23
    if(MakeCForNot!=1){
    	EventEnd(hbtEvent);
    	fPicoEvent = nullptr;
	return;
    }
    if (!tmpPassEvent) {
	// if collection2 == nullptr, we also update this event, but not do make pair!	
	if(coll_1_size!=0 && coll_2_size==0){
        	if ( MixingBufferFull() ) {
        		MixingBuffer()->pop_back();
		}
   	 	MixingBuffer()->push_front(fPicoEvent);
	}
	EventEnd(hbtEvent);
    	fPicoEvent = nullptr;
	return;
    }// end only

    if (AnalyzeIdenticalParticles()) {
        collection2 = nullptr;
    }

    MakePairs("real", collection1, collection2, EnablePairMonitors());
    for (auto storedEvent : *fMixingBuffer) {
	if (AnalyzeIdenticalParticles()) {
            MakePairs("mixed", collection1, storedEvent->FirstParticleCollection());
            // If non-identical - mix both combinations of first and second particles
        }else {
            // MakePairs("mixed", collection1,
            //                storedEvent->SecondParticleCollection());
            // dowang: only make first particle in pool with second particle in current event
    		MakePairs("mixed", storedEvent->FirstParticleCollection(),collection2);
	}
	
    }
    if ( MixingBufferFull() ) {
        MixingBuffer()->pop_back();
    }
    MixingBuffer()->push_front(fPicoEvent);
 EventEnd(hbtEvent);
	return;

}
void AliFemtoSimpleAnalysisOnlyMixP1::SetMakeCForNot(int fMake){
     MakeCForNot = fMake;
}
void AliFemtoSimpleAnalysisOnlyMixP1::SetOnlyP1Exist(int fOnlyP1Exist){
     OnlyP1Exist = fOnlyP1Exist;
}
