/*
 * AliAnalysisTaskFemtoSMI.cxx
 *
 *   Created on: 2020-08-17
 *   Author: S. Schneider
 * This file is based on a copy of:
 *   AliAnalysisTaskFemtoTutorial.cxx
 * with some parts taken from:
 *   AliAnalysisTaskAODLoton.cxx
 */

#include "AliAnalysisTaskFemtoSMI.h"
#include "AliFemtoDreamBasePart.h"
#include "AliLog.h"
#include "AliVEvent.h"
ClassImp(AliAnalysisTaskFemtoSMI)
AliAnalysisTaskFemtoSMI::AliAnalysisTaskFemtoSMI()
:AliAnalysisTaskSE()
,fIsMC(false)
,fOutput()
,fEvent()
,fTrack()
,fEventCuts()
,fProton() // fProton() replaces fTrackCutsPart1()
,fLambda() // fLambda() replaces fTrackCutsPart2()

,fConfig()
,fPairCleaner()
,fPartColl()
,fGTI()

,fTrigger(AliVEvent::kINT7)

,fTrackBufferSize()
{

}

AliAnalysisTaskFemtoSMI::AliAnalysisTaskFemtoSMI(const char *name, bool isMC)
:AliAnalysisTaskSE(name)
,fIsMC(isMC)
,fOutput()
,fEvent()
,fTrack()
,fEventCuts()
,fProton()
,fLambda()

,fConfig()
,fPairCleaner()
,fPartColl()
,fGTI()

,fTrigger(AliVEvent::kINT7)

,fTrackBufferSize(2000)
{
  DefineOutput(1,TList::Class());
}

AliAnalysisTaskFemtoSMI::~AliAnalysisTaskFemtoSMI() {
  // TODO Auto-generated destructor stub
}

void AliAnalysisTaskFemtoSMI::UserCreateOutputObjects() {
  fOutput = new TList();
  fOutput->SetName("Output"); // Every output objects needs a name, be careful names can collide!
  fOutput->SetOwner();        // This tells ROOT that this list belongs to the top list / top object

  //Set up the Femto Event
  //Now not to go into the details the arguments are as follows:
  //1. for the old runs the Pile Up rejection was configured manually, this is
  //from that time, nowadays we use the AliEventCuts (offically provided by the
  // ALICE DPG), and there is no need for this
  //2. Do you want the QA from the AliEventCuts?
  //3. The trigger, if you ever switch to High Multiplicity, you need to change this
  fEvent=new AliFemtoDreamEvent(false,true,fTrigger);// used to have "AliVEvent::kINT7" as third argument
  // IC: looks different, has 4 arg's
  fOutput->Add(fEvent->GetEvtCutList());
  //Nothing special about the Femto Track, we just initialize it
  fTrack=new AliFemtoDreamTrack();
  //If this is false, the MC information is not included
  fTrack->SetUseMCInfo(fIsMC);
  //Now this is a very ugly solution, a vector would be way better. Maybe one day it will change
  //(J.C.Hrist)
  //In principal we need this to map the tracks flagged with Filterbit 128 (TPC only Tracks constrained to the
  //Primary Vertex estimated with ITS tracks) to the global tracks. This is neccessary to have the full tracking
  //information available.
  fGTI=new AliAODTrack*[fTrackBufferSize];

  if (!fProton) {
    // If the track cuts didn't arrive here, we can go home
    AliFatal("Track Cuts for Proton not set!");
  }
  //For the next step the order is really important, else you do not have histograms
   //First initialize the Histogramms in the Track Cut object. This is done to not initialize
   //the histograms already in the add task, but in the Analysis Task itself and therefore on the
   //worker node.
  fEventCuts->InitQA();
  //And add the histograms to the output list
  fOutput->Add(fEventCuts->GetHistList());
  //same as for the event cuts
  fProton->Init();
  //To avoid collision in the output list, we rename the List carrying all the histograms of this object
  fProton->SetName("Proton");
  //Now connect the output of the Track Cuts to the output
  fOutput->Add(fProton->GetQAHists());
  //If we are running over MC more histos are created and we get them seperately
  //This is done, since later you might want to use seperate output slots, so you
  //do not overload one slot.
  if (fProton->GetIsMonteCarlo()) {
    fProton->SetMCName("MCProton"); // same as above
    fOutput->Add(fProton->GetMCQAHists());
  }

  //Same game for the second track cuts object
  if (!fLambda) {
    AliFatal("Track Cuts for Particle One not set!");
  }
  fLambda->Init();
  fLambda->SetName("Lambda");
  fOutput->Add(fLambda->GetQAHists());
  if (fLambda->GetIsMonteCarlo()) {
    fLambda->SetMCName("MCLambda"); // same as above
    fOutput->Add(fLambda->GetMCQAHists());
  }

  fv0 = new AliFemtoDreamv0();
  // trying out additions for fv0, might need further additions when including MC options
  fv0->SetPDGCode(3122);
  fv0->SetPDGDaughterPos(2212);
  fv0->SetPDGDaughterNeg(211);


  //Arguments for the pair cleaner as follows:
  //1. How many pairs of Tracks + Decays do you want to clean?
  //(for the purpose of this tutorial, we are going to treat the
  //second track as a decay ;)
  //2. How many decays and decays do you want to clean
  //3. Minimal booking == true means no histograms are created and filled
  //might be handy for systematic checks, in order to reduce the memory
  //usage
  fPairCleaner=new AliFemtoDreamPairCleaner(1,1,false);
  //The output histograms have to also be added to the output
  fOutput->Add(fPairCleaner->GetHistList());
  //1. fConfig: This is the config object from your AddTask where all the things are
  //set for the calculation of the correlatin function
  //2. Minimal booking again, see above. However the Results are always filled, just
  //the QA is not existing.
  fPartColl=new AliFemtoDreamPartCollection(fConfig,false);
  //The output also has to be added to the list
  fOutput->Add(fPartColl->GetHistList());
  fOutput->Add(fPartColl->GetQAList());
  PostData(1,fOutput);
}



void AliAnalysisTaskFemtoSMI::UserExec(Option_t *) {
  AliAODEvent *Event=static_cast<AliAODEvent*>(fInputEvent);
  if (!Event) {
    AliWarning("No Input Event");
  } else {
    //Put all the event information into the Event Object
    fEvent->SetEvent(Event);
    //Use the event cut object to check if we want to use this event by
    //applying the selection criteria defined in the add task.
    if (fEventCuts->isSelected(fEvent)) {
      //if we pass the event selection criteria we can continue. First some boring
      //technical necessary stuff. We need to fill the Global Track Array (fGTI) once
      //per event with our mapping of the tracks. For beginners this is fine to ignore.
      ResetGlobalTrackReference();
      for(int iTrack = 0;iTrack<Event->GetNumberOfTracks();++iTrack){
        AliAODTrack *track=static_cast<AliAODTrack*>(Event->GetTrack(iTrack));
        if (!track) {
          AliFatal("No Standard AOD");
          return;
        }
        StoreGlobalTrackReference(track);
      }
      //Set the global track array for the track to be able to access it.
      fTrack->SetGlobalTrackInfo(fGTI,fTrackBufferSize);
      //After this is done, let the search for particles begin.
      //This is where we are going to store the particles we find for now. But hold up!
      //AliFemtoDreamBasePart?? Well, see every particle candidate type, be it track, v0 or
      //cascade inherits from this base type, which carries all the information relevant for the
      //same and mixed event distributions (e.g. momentum, angles, etc. ). E.g. the dEdx of the TPC
      //doesnt matter for this porpouse.
      static std::vector<AliFemtoDreamBasePart> Particles;
      Particles.clear();
      //Now we loop over all the tracks in the reconstructed event.
      for (int iTrack = 0;iTrack<Event->GetNumberOfTracks();++iTrack) {
        AliAODTrack *track=static_cast<AliAODTrack*>(Event->GetTrack(iTrack));
        if (!track) {
          AliFatal("No Standard AOD");
          return;
        }
        //First we fill all the track information into the track object
        fTrack->SetTrack(track);
        //then we pass it to the track cut objects, which checks if it passes our selection criteria
        //for particle 1 and if it returns true ...
        if (fProton->isSelected(fTrack)) {
          //.. we add it to our particle buffer
          Particles.push_back(*fTrack);
        }
        //same game, different name for the selection criteria 2!

		/* The following doesn't work for Lambdas anymore, using section below instead:
        if (fLambda->isSelected(fTrack)) {
          AntiParticles.push_back(*fTrack);
        }
		*/
      }
	  std::vector<AliFemtoDreamBasePart> Lambdas;
	  //std::vector<AliFemtoDreamBasePArt> AntiLambdas;
	  
  fv0->SetGlobalTrackInfo(fGTI, fTrackBufferSize);

	  for (int iv0 = 0;
      	  iv0 < static_cast<TClonesArray *>(Event->GetV0s())->GetEntriesFast();
      	  ++iv0) {
        AliAODv0* v0 = Event->GetV0(iv0);
        fv0->Setv0(Event, v0);
        if (fLambda->isSelected(fv0)) {
      	  Lambdas.push_back(*fv0);
    	}
  	  }

      //This is where the magic of selecting particles ends, and we can turn our attention to
      //calculating the results. First we need to ensure to not have any Autocorrelations by
      //selecting a track twice. Now this is hypothetical, because we are selecting opposite
      //charged particles, but imagine you want to use p+K^+ (Check for this is not yet implemented)!
      fPairCleaner->CleanTrackAndDecay(&Particles,&Lambdas,0);
      fPairCleaner->CleanDecay(&Lambdas, 0);
      //The cleaner tags particles as 'bad' for use, these we don't want to give to our particle
      //pairer, that's why we call store particles, which only takes the particles marked 'good' from
      //our buffer vector.
      //First we need to reset any particles in the array!

      fPairCleaner->ResetArray();
      fPairCleaner->StoreParticle(Particles);
      fPairCleaner->StoreParticle(Lambdas);

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

void AliAnalysisTaskFemtoSMI::ResetGlobalTrackReference(){ //IC: (almost) same
    //This method was inherited from H. Beck analysis
    
    // Sets all the pointers to zero. To be called at
    // the beginning or end of an event
    for(UShort_t i=0;i<fTrackBufferSize;i++)
    {
        fGTI[i]=0;
    }
}
void AliAnalysisTaskFemtoSMI::StoreGlobalTrackReference(AliAODTrack *track){ //IC: same
    //This method was inherited from H. Beck analysis
    
    //bhohlweg@cern.ch: We ask for the Unique Track ID that points back to the
    //ESD. Seems like global tracks have a positive ID, Tracks with Filterbit
    //128 only have negative ID, this is used to match the Tracks later to their
    //global counterparts
    
    // Stores the pointer to the global track
    
    // This was AOD073
    // // Don't use the filter bits 2 (ITS standalone) and 128 TPC only
    // // Remove this return statement and you'll see they don't have
    // // any TPC signal
    // if(track->TestFilterBit(128) || track->TestFilterBit(2))
    //   return;
    // This is AOD086
    // Another set of tracks was introduced: Global constrained.
    // We only want filter bit 1 <-- NO! we also want no
    // filter bit at all, which are the v0 tracks
    //  if(!track->TestFilterBit(1))
    //    return;
    
    // There are also tracks without any filter bit, i.e. filter map 0,
    // at the beginning of the event: they have ~id 1 to 5, 1 to 12
    // This are tracks that didn't survive the primary track filter but
    // got written cause they are V0 daughters
    
    // Check whether the track has some info
    // I don't know: there are tracks with filter bit 0
    // and no TPC signal. ITS standalone V0 daughters?
    // if(!track->GetTPCsignal()){
    //   printf("Warning: track has no TPC signal, "
    //     //    "not adding it's info! "
    //     "ID: %d FilterMap: %d\n"
    //     ,track->GetID(),track->GetFilterMap());
    //   //    return;
    // }
    
    // Check that the id is positive
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
        if( (!track->GetFilterMap()) && (!track->GetTPCNcls()) ){
            return;
        }
        // Imagine the other way around, the zero map zero clusters track
        // is stored and the good one wants to be added. We ommit the warning
        // and just overwrite the 'bad' track
        if( fGTI[trackID]->GetFilterMap() || fGTI[trackID]->GetTPCNcls()  ){
            // If we come here, there's a problem
            printf("Warning! global track info already there!");
            printf("         TPCNcls track1 %u track2 %u",
                   (fGTI[trackID])->GetTPCNcls(),track->GetTPCNcls());
            printf("         FilterMap track1 %u track2 %u\n",
                   (fGTI[trackID])->GetFilterMap(),track->GetFilterMap());
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
