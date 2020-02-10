/*
 * AliAnalysisTaskNanoAODFemtoDreamPion.cxx
 *
 *  Created on: 16 Oct 2019
 *      Author: M. Korwieser
 */

#include "AliAnalysisTaskNanoAODFemtoDreamPion.h"
#include "AliFemtoDreamBasePart.h"
#include "AliLog.h"
#include "AliVEvent.h"
#include "AliNanoAODTrack.h"
ClassImp(AliAnalysisTaskNanoAODFemtoDreamPion)
AliAnalysisTaskNanoAODFemtoDreamPion::AliAnalysisTaskNanoAODFemtoDreamPion()
:AliAnalysisTaskSE()
,fIsMC(false)
,fOutput()
,fEvent()
,fTrack()
,fTrigger(AliVEvent::kINT7)
,fEventCuts()
,fTrackCutsPosPion()
,fTrackCutsNegPion()
,fTrackCutsProton()
,fConfig()
,fPairCleaner()
,fPartColl()
,fGTI()
,fTrackBufferSize()
{

}

AliAnalysisTaskNanoAODFemtoDreamPion::AliAnalysisTaskNanoAODFemtoDreamPion(const char *name, bool isMC)
:AliAnalysisTaskSE(name)
,fIsMC(isMC)
,fOutput()
,fEvent()
,fTrack()
,fTrigger(AliVEvent::kINT7)
,fEventCuts()
,fTrackCutsPosPion()
,fTrackCutsNegPion()
,fTrackCutsProton()
,fConfig()
,fPairCleaner()
,fPartColl()
,fGTI()
,fTrackBufferSize(2000)
{
  DefineOutput(1,TList::Class());
}

AliAnalysisTaskNanoAODFemtoDreamPion::~AliAnalysisTaskNanoAODFemtoDreamPion() {}

void AliAnalysisTaskNanoAODFemtoDreamPion::UserCreateOutputObjects() {
  fOutput = new TList();
  fOutput->SetName("Output");
  fOutput->SetOwner();

  fEvent=new AliFemtoDreamEvent(false,true,fTrigger);
  fOutput->Add(fEvent->GetEvtCutList());

  fTrack=new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(fIsMC);

  fGTI=new AliVTrack*[fTrackBufferSize];

  if (!fEventCuts) {
      AliFatal("Event Cuts not set!");
  }
    fEventCuts->InitQA();
    fOutput->Add(fEventCuts->GetHistList());

  if (!fTrackCutsPosPion) {
    AliFatal("Track Cuts for Pion+ not set!");
  }
  fTrackCutsPosPion->Init();
  fTrackCutsPosPion->SetName("Pion+");
  fOutput->Add(fTrackCutsPosPion->GetQAHists());
  if (fTrackCutsPosPion->GetIsMonteCarlo()) {
    fTrackCutsPosPion->SetMCName("MCPion+");
    fOutput->Add(fTrackCutsPosPion->GetMCQAHists());
  }

  if (!fTrackCutsNegPion) {
    AliFatal("Track Cuts for Pion- not set!");
  }
  fTrackCutsNegPion->Init();
  fTrackCutsNegPion->SetName("Pion-");
  fOutput->Add(fTrackCutsNegPion->GetQAHists());
  if (fTrackCutsNegPion->GetIsMonteCarlo()) {
    fTrackCutsNegPion->SetMCName("MCPion-");
    fOutput->Add(fTrackCutsNegPion->GetMCQAHists());
  }

  if (!fTrackCutsProton) {
    AliFatal("Track Cuts for Proton+ not set!");
  }
  fTrackCutsProton->Init();
  fTrackCutsProton->SetName("Proton+");
  fOutput->Add(fTrackCutsProton->GetQAHists());
  if (fTrackCutsProton->GetIsMonteCarlo()) {
    fTrackCutsProton->SetMCName("MCProton+");
    fOutput->Add(fTrackCutsProton->GetMCQAHists());
  }
  //Arguments for the pair cleaner as follows:
  //1. How many pairs of Tracks + Decays do you want to clean?
  //(for the purpose of this tutorial, we are going to treat the
  //second track as a decay ;)
  //2. How many decays and decays do you want to clean
  //3. Minimal booking == true means no histograms are created and filled
  //might be handy for systematic checks, in order to reduce the memory
  //usage
  fPairCleaner=new AliFemtoDreamPairCleaner(1,0,fConfig->GetMinimalBookingME());
  fOutput->Add(fPairCleaner->GetHistList());

  fPartColl=new AliFemtoDreamPartCollection(fConfig,fConfig->GetMinimalBookingME());
  fOutput->Add(fPartColl->GetHistList());
  fOutput->Add(fPartColl->GetQAList());

  PostData(1,fOutput);
}



void AliAnalysisTaskNanoAODFemtoDreamPion::UserExec(Option_t *) {
  AliVEvent *Event= fInputEvent;
  if (!Event) {
    AliWarning("No Input Event");
  } else {
    fEvent->SetEvent(Event);
    if (fEventCuts->isSelected(fEvent)) {
      ResetGlobalTrackReference();
      for(int iTrack = 0;iTrack<Event->GetNumberOfTracks();++iTrack){
        AliVTrack *track=static_cast<AliVTrack*>(Event->GetTrack(iTrack));
        if (!track) {
          AliFatal("No Standard NanoAOD");
          return;
        }
        StoreGlobalTrackReference(track);
      }
      fTrack->SetGlobalTrackInfo(fGTI,fTrackBufferSize);

      static std::vector<AliFemtoDreamBasePart> PosPions;
      PosPions.clear();
      static std::vector<AliFemtoDreamBasePart> NegPions;
      NegPions.clear();
      static std::vector<AliFemtoDreamBasePart> Protons;
      Protons.clear();

      for (int iTrack = 0;iTrack<Event->GetNumberOfTracks();++iTrack) {
        AliVTrack *track=static_cast<AliVTrack*>(Event->GetTrack(iTrack));
        if (!track) {
          AliFatal("No Standard NanoAOD");
          return;
        }
        fTrack->SetTrack(track, Event);

        if (fTrackCutsPosPion->isSelected(fTrack)) {
          PosPions.push_back(*fTrack);
        }
        if (fTrackCutsNegPion->isSelected(fTrack)) {
          NegPions.push_back(*fTrack);
        }
        if (fTrackCutsProton->isSelected(fTrack)) {
          Protons.push_back(*fTrack);
        }
      }

      //This is where the magic of selecting particles ends, and we can turn our attention to
      //calculating the results. First we need to ensure to not have any Autocorrelations by
      //selecting a track twice. Now this is hypothetical, because we are selecting opposite
      //charged particles, but imagine you want to use p+K^+ (Check for this is not yet implemented)!

      //fPairCleaner->CleanTrackAndDecay(&PosPions,&PosPions,0); // pi+ pi+
      //fPairCleaner->CleanTrackAndDecay(&PosPions,&NegPions,1); // pi+ pi-
      //fPairCleaner->CleanTrackAndDecay(&NegPions,&NegPions,2); // pi- pi-
      fPairCleaner->CleanTrackAndDecay(&PosPions, &Protons, 0);

      //The cleaner tags particles as 'bad' for use, these we don't want to give to our particle
      //pairer, that's why we call store particles, which only takes the particles marked 'good' from
      //our buffer vector.
      //First we need to reset any particles in the array!
      fPairCleaner->ResetArray();
      fPairCleaner->StoreParticle(PosPions);
      fPairCleaner->StoreParticle(NegPions);
      fPairCleaner->StoreParticle(Protons);
      //Now we can give our particlers to the particle collection where the magic happens.
      //The arguments one by one:
      //1. A vector of a vector of cleaned particles fresh from the laundromat.
      //2. The Z-Vtx of the event, this is used for event mixing, since we only want to mix events which
      //have the same acceptance in the detector to avoid acceptance effects.
      //3. Same for the multiplicity.
      //4. The centrality is really only of interest for pPb collisions (or PbPb) in order to also
      //do a binning there - kind of similar to a multiplicity binning.
      fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),fEvent->GetZVertex(),
                         fEvent->GetRefMult08(),fEvent->GetV0MCentrality());
      //For this porpouse you are done, now you only need to post the output.
      PostData(1,fOutput);
    }
  }
}

void AliAnalysisTaskNanoAODFemtoDreamPion::ResetGlobalTrackReference(){
    // for documentation see AliFemtoDreamAnalysis
    
    // Sets all the pointers to zero. To be called at
    // the beginning or end of an event
    for(UShort_t i=0;i<fTrackBufferSize;i++)
    {
        fGTI[i]=0;
    }
}
void AliAnalysisTaskNanoAODFemtoDreamPion::StoreGlobalTrackReference(AliVTrack *track){
    // for documentation see AliFemtoDreamAnalysis
    
    // Check that the id is positive
    // AliVTrack has no filtermap
    AliNanoAODTrack *nanoTrack = dynamic_cast<AliNanoAODTrack *>(track);
    const int trackID = track->GetID();
    if(trackID<0){
        return;
    }
    
    // Check id is not too big for buffer
    if(trackID>=fTrackBufferSize){
        printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n"
               ,trackID,fTrackBufferSize);
        return;
    }
    
    // Warn if we overwrite a track
    if(fGTI[trackID])
    {
        // Seems like there are FilterMap 0 tracks
        // that have zero TPCNcls, don't store these!
        if( (!nanoTrack->GetFilterMap()) && (!track->GetTPCNcls()) ){
            return;
        }
        // Imagine the other way around, the zero map zero clusters track
        // is stored and the good one wants to be added. We ommit the warning
        // and just overwrite the 'bad' track
        if( (dynamic_cast<AliNanoAODTrack *>(fGTI[trackID]))->GetFilterMap() || fGTI[trackID]->GetTPCNcls()  ){
            // If we come here, there's a problem
            printf("Warning! global track info already there!");
            printf("         TPCNcls track1 %u track2 %u",
                   ( dynamic_cast<AliNanoAODTrack *>(fGTI[trackID]))->GetTPCNcls(),track->GetTPCNcls());
            printf("         FilterMap track1 %u track2 %u\n",
                   ( dynamic_cast<AliNanoAODTrack *>(fGTI[trackID]))->GetFilterMap(),nanoTrack->GetFilterMap());
        }
    } // Two tracks same id
    
    // // There are tracks with filter bit 0,
    // // do they have TPCNcls stored?
    // if(!track->GetFilterMap()){
    //   printf("Filter map is zero, TPCNcls: %u\n"
    //     ,track->GetTPCNcls());
    // }
    
    // Assign the pointer
    (fGTI[trackID]) = track;
}
