/*
 * AliAnalysisTaskNanoAODFemtoDreamSigPi.cxx
 *
 *  Created on: 11 Oct 2019
 *      Author: M. Korwieser
 */

#include "AliAnalysisTaskNanoAODFemtoDreamSigPi.h"
#include "AliFemtoDreamBasePart.h"
#include "AliLog.h"
#include "AliVEvent.h"
ClassImp(AliAnalysisTaskNanoAODFemtoDreamSigPi)
AliAnalysisTaskNanoAODFemtoDreamSigPi::AliAnalysisTaskNanoAODFemtoDreamSigPi()
:AliAnalysisTaskSE()
,fIsMC(false)
,fOutput()
,fEvent()
,fTrack()
,fEventCuts()
,fTrackCutsPosPion()
,fTrackCutsNegPion()
,fTrackCutsPion0()
,fTrackCutsSigmaPlus()
,fTrackCutsSigmaMinus()
,fTrackCutsSigma0()
,fConfig()
,fPairCleaner()
,fPartColl()
,fGTI()
,fTrackBufferSize()
{

}

AliAnalysisTaskNanoAODFemtoDreamSigPi::AliAnalysisTaskNanoAODFemtoDreamSigPi(const char *name, bool isMC)
:AliAnalysisTaskSE(name)
,fIsMC(isMC)
,fOutput()
,fEvent()
,fTrack()
,fEventCuts()
,fTrackCutsPosPion()
,fTrackCutsNegPion()
,fTrackCutsPion0()
,fTrackCutsSigmaPlus()
,fTrackCutsSigmaMinus()
,fTrackCutsSigma0()
,fConfig()
,fPairCleaner()
,fPartColl()
,fGTI()
,fTrackBufferSize(2000)
{
  DefineOutput(1,TList::Class());
}

AliAnalysisTaskNanoAODFemtoDreamSigPi::~AliAnalysisTaskNanoAODFemtoDreamSigPi() {}

void AliAnalysisTaskNanoAODFemtoDreamSigPi::UserCreateOutputObjects() {
  fOutput = new TList();
  fOutput->SetName("Output");
  fOutput->SetOwner();

  fEvent=new AliFemtoDreamEvent(false,true,GetCollisionCandidates());
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

  if (!fTrackCutsPion0) {
    AliFatal("Track Cuts for Pion0 not set!");
  }
  fTrackCutsPion0->Init();
  fTrackCutsPion0->SetName("Pion0");
  fOutput->Add(fTrackCutsPion0->GetQAHists());
  if (fTrackCutsPion0->GetIsMonteCarlo()) {
    fTrackCutsPion0->SetMCName("MCPion0");
    fOutput->Add(fTrackCutsPion0->GetMCQAHists());
  }
  if (!fTrackCutsSigmaMinus) {
    AliFatal("Track Cuts for Sigma- not set!");
  }
  fTrackCutsSigmaMinus->Init();
  fTrackCutsSigmaMinus->SetName("Sigma-");
  fOutput->Add(fTrackCutsSigmaMinus->GetQAHists());
  if (fTrackCutsSigmaMinus->GetIsMonteCarlo()) {
    fTrackCutsSigmaMinus->SetMCName("MCSigma-");
    fOutput->Add(fTrackCutsSigmaMinus->GetMCQAHists());
  }

  if (!fTrackCutsSigmaPlus) {
    AliFatal("Track Cuts for Sigma+ not set!");
  }
  fTrackCutsSigmaPlus->Init();
  fTrackCutsSigmaPlus->SetName("Sigma+");
  fOutput->Add(fTrackCutsSigmaPlus->GetQAHists());
  if (fTrackCutsSigmaPlus->GetIsMonteCarlo()) {
    fTrackCutsSigmaPlus->SetMCName("MCSigma+");
    fOutput->Add(fTrackCutsSigmaPlus->GetMCQAHists());
  }

  if (!fTrackCutsSigma0) {
    AliFatal("Track Cuts for Sigm0 not set!");
  }
  fTrackCutsSigma0->Init();
  fTrackCutsSigma0->SetName("Sigm0");
  fOutput->Add(fTrackCutsSigma0->GetQAHists());
  if (fTrackCutsSigma0->GetIsMonteCarlo()) {
    fTrackCutsSigma0->SetMCName("MCSigm0");
    fOutput->Add(fTrackCutsSigma0->GetMCQAHists());
  }
  //Arguments for the pair cleaner as follows:
  //1. How many pairs of Tracks + Decays do you want to clean?
  //(for the purpose of this tutorial, we are going to treat the
  //second track as a decay ;)
  //2. How many decays and decays do you want to clean
  //3. Minimal booking == true means no histograms are created and filled
  //might be handy for systematic checks, in order to reduce the memory
  //usage
  fPairCleaner=new AliFemtoDreamPairCleaner(0,0,fConfig->GetMinimalBookingME());
  fOutput->Add(fPairCleaner->GetHistList());

  fPartColl=new AliFemtoDreamPartCollection(fConfig,fConfig->GetMinimalBookingME());
  fOutput->Add(fPartColl->GetHistList());
  fOutput->Add(fPartColl->GetQAList());

  PostData(1,fOutput);
}



void AliAnalysisTaskNanoAODFemtoDreamSigPi::UserExec(Option_t *) {
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
      static std::vector<AliFemtoDreamBasePart> Pion0s;
      Pion0s.clear();
      static std::vector<AliFemtoDreamBasePart> SigmasPlus;
      SigmasPlus.clear();
      static std::vector<AliFemtoDreamBasePart> SigmasMinus;
      SigmasMinus.clear();
      static std::vector<AliFemtoDreamBasePart> Sigmas0;
      Sigmas0.clear();


      for (int iTrack = 0;iTrack<Event->GetNumberOfTracks();++iTrack) {
        AliVTrack *track=static_cast<AliVTrack*>(Event->GetTrack(iTrack));
        if (!track) {
          AliFatal("No Standard NanoAOD");
          return;
        }
        fTrack->SetTrack(track, Event);

	int pdgcode = fTrack->GetMCPDGCode();
        //if ( pdgcode != 211 && pdgcode != -211) {
	//   printf("PDGCode value: %i \n",pdgcode);
        //}
        
        if (fTrackCutsPosPion->isSelected(fTrack) && pdgcode == 211) {
          PosPions.push_back(*fTrack);
	  //printf("PDGCode =211: %i \n",pdgcode);
        }
        if (fTrackCutsNegPion->isSelected(fTrack) && pdgcode == -211) {
          NegPions.push_back(*fTrack);
	  //printf("PDGCode =-211: %i \n",pdgcode);
        }
        if (fTrackCutsPion0->isSelected(fTrack) && pdgcode == 111) {
          Pion0s.push_back(*fTrack);
	  //printf("PDGCode =111: %i \n",pdgcode);
        }
        if (fTrackCutsSigmaPlus->isSelected(fTrack) && pdgcode == 3222) {
          SigmasPlus.push_back(*fTrack);
	  //printf("PDGCode =3222: %i \n",pdgcode);
        }
        if (fTrackCutsSigmaMinus->isSelected(fTrack) && pdgcode == 3112) {
          SigmasMinus.push_back(*fTrack);
          //printf("PDGCode =3112: %i \n",pdgcode);
        }
        if (fTrackCutsSigma0->isSelected(fTrack) && pdgcode == 3212) {
          Sigmas0.push_back(*fTrack);
	  //printf("PDGCode =3212: %i \n",pdgcode);
        }
      }
      //printf("Particle count 211: %i \n",PosPions.size());
      //printf("Particle count -211: %i \n",NegPions.size());
      //printf("Particle count 111: %i \n",Pion0s.size());
      //printf("Particle count 3222: %i \n",SigmasPlus.size());
      //printf("Particle count 3112: %i \n",SigmasMinus.size());
      //printf("Particle count 3212: %i \n",Sigmas0.size());


      //This is where the magic of selecting particles ends, and we can turn our attention to
      //calculating the results. First we need to ensure to not have any Autocorrelations by
      //selecting a track twice. Now this is hypothetical, because we are selecting opposite
      //charged particles, but imagine you want to use p+K^+ (Check for this is not yet implemented)!

      //fPairCleaner->CleanTrackAndDecay(&PosPions,&PosPions,0); // pi+ pi+
      //fPairCleaner->CleanTrackAndDecay(&PosPions,&NegPions,1); // pi+ pi-
      //fPairCleaner->CleanTrackAndDecay(&NegPions,&NegPions,2); // pi- pi-

      //The cleaner tags particles as 'bad' for use, these we don't want to give to our particle
      //pairer, that's why we call store particles, which only takes the particles marked 'good' from
      //our buffer vector.
      //First we need to reset any particles in the array!
      fPairCleaner->ResetArray();
      fPairCleaner->StoreParticle(PosPions);
      fPairCleaner->StoreParticle(NegPions);
      fPairCleaner->StoreParticle(Pion0s);
      fPairCleaner->StoreParticle(SigmasPlus);
      fPairCleaner->StoreParticle(SigmasMinus);
      fPairCleaner->StoreParticle(Sigmas0);
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

void AliAnalysisTaskNanoAODFemtoDreamSigPi::ResetGlobalTrackReference(){
    // for documentation see AliFemtoDreamAnalysis
    
    // Sets all the pointers to zero. To be called at
    // the beginning or end of an event
    for(UShort_t i=0;i<fTrackBufferSize;i++)
    {
        fGTI[i]=0;
    }
}
void AliAnalysisTaskNanoAODFemtoDreamSigPi::StoreGlobalTrackReference(AliVTrack *track){
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
