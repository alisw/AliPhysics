/*
 * AliLightNAnalysis.cxx
 *
 *  Created on: 24 Nov 2017
 *      Author: bernhardhohlweger
 */
#include "AliLightNAnalysis.h"
//#include "TClonesArray.h"
ClassImp(AliLightNAnalysis)
AliLightNAnalysis::AliLightNAnalysis()
:fMVPileUp(false)
,fEvtCutQA(false)
,fQA(0)
,fLightNTrack()
,fEvent()
,fEvtCuts()
,fTrackCutsProton()
,fTrackCutsDeuteron()
,fAntiTrackCutsProton()
,fAntiTrackCutsDeuteron()
,fTrackBufferSize(0)
,fGTI(0)
{
    
}

AliLightNAnalysis::~AliLightNAnalysis() {
    if (fLightNTrack) {
        delete fLightNTrack;
    }
}

void AliLightNAnalysis::Init() {
    fLightNTrack=new AliLightNTrack();
    fLightNTrack->SetUseMCInfo(fTrackCutsProton->GetIsMonteCarlo());
    fLightNTrack->SetUseMCInfo(fAntiTrackCutsProton->GetIsMonteCarlo());
    fLightNTrack->SetUseMCInfo(fTrackCutsDeuteron->GetIsMonteCarlo());
    fLightNTrack->SetUseMCInfo(fAntiTrackCutsDeuteron->GetIsMonteCarlo());
    fEvtCuts->InitQA();
    fTrackCutsProton->Init();
    fAntiTrackCutsProton->Init();
    fTrackCutsDeuteron->Init();
    fAntiTrackCutsDeuteron->Init();
    fGTI=new AliAODTrack*[fTrackBufferSize];
    fQA=new TList();
    fQA->SetOwner();
    fQA->SetName("QA");
    
    fEvent=new AliLightNEvent(fMVPileUp,fEvtCutQA);
    fQA->Add(fEvent->GetEvtCutList());
    
    return;
}

void AliLightNAnalysis::ResetGlobalTrackReference(){
    //This method was inherited form H. Beck analysis
    
    // Sets all the pointers to zero. To be called at
    // the beginning or end of an event
    for(UShort_t i=0;i<fTrackBufferSize;i++)
    {
        fGTI[i]=0;
    }
}
void AliLightNAnalysis::StoreGlobalTrackReference(AliAODTrack *track){
    //This method was inherited form H. Beck analysis
    
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

void AliLightNAnalysis::Make(AliAODEvent *evt) {
    if (!evt) {
        AliFatal("No Input Event");
    }
    
    Float_t lPercentile = 300;
    AliMultSelection *MultSelection = 0x0;
    MultSelection = (AliMultSelection*)evt->FindListObject("MultSelection");
    if( !MultSelection) {
        //If you get this warning (and lPercentiles 300) please check that the AliMultSelectionTask actually ran (before your task)
        AliWarning("AliMultSelection object not found!");
    }else{
        lPercentile = MultSelection->GetMultiplicityPercentile("V0M",true);
        fEvtCuts->FillV0Mlpercentile(lPercentile);
        fEvtCuts->FillV0MlpercentileHM(lPercentile);
    }
    
    fEvent->SetEvent(evt);
    if (!fEvtCuts->isSelected(fEvent)) {
        return;
    }
      //std::cout << "=============================" <<std::endl;
      //std::cout << "=============================" <<std::endl;
      //std::cout << "=========new Event===========" <<std::endl;
      //std::cout << "=============================" <<std::endl;
      //std::cout << "=============================" <<std::endl;
    
    // Loop over all MC particle that are in the MCParticles array
    if(fTrackCutsProton->GetIsMonteCarlo()){
        TClonesArray* AODMCTrackArray = dynamic_cast<TClonesArray*>(evt->FindListObject(AliAODMCParticle::StdBranchName()));
        if (AODMCTrackArray == NULL) return;
        for(Long_t i = 0; i < AODMCTrackArray->GetEntriesFast(); i++) {
            AliAODMCParticle* particle = static_cast<AliAODMCParticle*>(AODMCTrackArray->At(i));
            if (!particle) continue;
            double p = particle->P();
            //Same eta range for all particles
            if(fTrackCutsProton->GetEtaMin()<particle->Eta() && particle->Eta()<fTrackCutsProton->GetEtaMax()){
                if (particle->GetPdgCode() == 2212){
                    fTrackCutsProton->FillStackGenerated(p);
                    if(particle->IsPhysicalPrimary()){
                        fTrackCutsProton->FillStackGeneratedPrimary(p);
                    }
                }else if (particle->GetPdgCode() == -2212){
                    fAntiTrackCutsProton->FillStackGenerated(p);
                    if(particle->IsPhysicalPrimary()){
                        fAntiTrackCutsProton->FillStackGeneratedPrimary(p);
                    }
                }else if (particle->GetPdgCode() == 1000010020){
                    fTrackCutsDeuteron->FillStackGenerated(p);
                    if(particle->IsPhysicalPrimary()){
                        fTrackCutsDeuteron->FillStackGeneratedPrimary(p);
                    }
                }else if (particle->GetPdgCode() == -1000010020){
                    fAntiTrackCutsDeuteron->FillStackGenerated(p);
                    if(particle->IsPhysicalPrimary()){
                        fAntiTrackCutsDeuteron->FillStackGeneratedPrimary(p);
                    }
                }
            }
        }
    }
    ResetGlobalTrackReference();
    for(int iTrack = 0;iTrack<evt->GetNumberOfTracks();++iTrack){
        AliAODTrack *track=static_cast<AliAODTrack*>(evt->GetTrack(iTrack));
        if (!track) {
            AliFatal("No Standard AOD");
            return;
        }
        StoreGlobalTrackReference(track);
    }
    fLightNTrack->SetGlobalTrackInfo(fGTI,fTrackBufferSize);
    for (int iTrack = 0;iTrack<evt->GetNumberOfTracks();++iTrack) {
        AliAODTrack *track=static_cast<AliAODTrack*>(evt->GetTrack(iTrack));
        if (!track) {
            AliFatal("No Standard AOD");
            return;
        }
        fLightNTrack->SetTrack(track);
        fTrackCutsProton->isSelected(fLightNTrack);
        fAntiTrackCutsProton->isSelected(fLightNTrack);
        fTrackCutsDeuteron->isSelected(fLightNTrack);
        fAntiTrackCutsDeuteron->isSelected(fLightNTrack);
    }
}

