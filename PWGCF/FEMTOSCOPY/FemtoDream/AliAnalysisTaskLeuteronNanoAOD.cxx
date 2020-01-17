/*
 * AliAnalysisTaskLeuteronNanoAOD.cxx
 *
 *  Created on:	06 November 2019
 *	Author:	Michael Jung
 */

#include "AliAnalysisTaskLeuteronNanoAOD.h"
#include "AliNanoAODTrack.h"

ClassImp(AliAnalysisTaskLeuteronNanoAOD)

//  -----------------------------------------------------------------------------------------------------------------------------------------
AliAnalysisTaskLeuteronNanoAOD::AliAnalysisTaskLeuteronNanoAOD():AliAnalysisTaskSE(),
  fIsMC(false),
  fTrackBufferSize(),
  fOutput(nullptr),
  fEvent(nullptr),
  fTrack(nullptr),
  fFemtov0(nullptr),
  fEventCuts(nullptr),
  fTrackCutsPart1(nullptr),
  fTrackCutsPart2(nullptr),
  fv0CutsPart3(nullptr),
  fv0CutsPart4(nullptr),
  fConfig(nullptr),
  fPairCleaner(nullptr),
  fPartColl(nullptr),
  fGTI(nullptr)
{

}


//  -----------------------------------------------------------------------------------------------------------------------------------------
AliAnalysisTaskLeuteronNanoAOD::AliAnalysisTaskLeuteronNanoAOD(const char *name, bool isMC):AliAnalysisTaskSE(name),
  fIsMC(isMC),
  fTrackBufferSize(2000),
  fOutput(),
  fEvent(),
  fTrack(),
  fFemtov0(),
  fEventCuts(),
  fTrackCutsPart1(),
  fTrackCutsPart2(),
  fv0CutsPart3(),
  fv0CutsPart4(),
  fConfig(),
  fPairCleaner(),
  fPartColl(),
  fGTI()
{
  DefineOutput(1,TList::Class());   // define the output of the analysis
}

//  -----------------------------------------------------------------------------------------------------------------------------------------
AliAnalysisTaskLeuteronNanoAOD::~AliAnalysisTaskLeuteronNanoAOD(){	// destructor -> check if the object exists, if so delete it 

  if(fOutput){
    delete fOutput;
  }

  if(fEvent){
    delete fEvent;
  }

  if(fTrack){
    delete fTrack;
  }

  if(fFemtov0){
    delete fFemtov0;
  }

  if(fEventCuts){
    delete fEventCuts;
  }

  if(fTrackCutsPart1){
    delete fTrackCutsPart1;
  }

  if(fTrackCutsPart2){
    delete fTrackCutsPart2;
  }

  if(fv0CutsPart3){
    delete fv0CutsPart3;
  }

  if(fv0CutsPart4){
    delete fv0CutsPart4;
  }

  if(fPairCleaner){
    delete fPairCleaner;
  }

  if(fPartColl){
    delete fPartColl;
  }

}

//  -----------------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskLeuteronNanoAOD::UserCreateOutputObjects(){

  fOutput = new TList();	    // create a list were all the histograms can be added to
  fOutput->SetName("Output");	    // give the output object a name
  fOutput->SetOwner();		    // tell ROOT that this object belongs to the top list / object


  // Check if the cut objects exists, if so initialize them

  if(!fEventCuts){
    AliFatal("Event Cuts (fEventCuts) not set!\n");
  } else{
    fEventCuts->InitQA();
  }

  if(!fTrackCutsPart1){					// check if the track cuts object is set, if not, call it a day
    AliFatal("Track Cuts for Particle 1 (fTrackCutsPart1) not set!\n");
  } else{
      fTrackCutsPart1->Init();				// initialize the histograms of this object 
      fTrackCutsPart1->SetName("Particle1");		// rename the list carrying the histograms of this object (to avoid collisions in the output object)
      fOutput->Add(fTrackCutsPart1->GetQAHists());	// add the histograms of the track cuts to the output object
      if(fTrackCutsPart1->GetIsMonteCarlo()){		// check if IsMC is "true" or "false"
        fTrackCutsPart1->SetMCName("MCParticle1");	// rename the list carrying the histograms of this object (to avoid collisions in the output object)
        fOutput->Add(fTrackCutsPart1->GetMCQAHists());  // add the histograms of the Monte Carlo to the output object
      }
    }

  if(!fTrackCutsPart2){
    AliFatal("Track Cuts for Particle 2 (fTrackCutsPart2) not set!\n");
  } else{
      fTrackCutsPart2->Init();
      fTrackCutsPart2->SetName("Particle2");
      fOutput->Add(fTrackCutsPart2->GetQAHists());
      if(fTrackCutsPart2->GetIsMonteCarlo()){
        fTrackCutsPart2->SetMCName("MCParticle2");
        fOutput->Add(fTrackCutsPart2->GetMCQAHists());
      }
    }

  if(!fv0CutsPart3){
    AliFatal("V0 Cuts for Particle 3 (fv0CutsPart3) not set!\n");
  } else{
      fv0CutsPart3->Init();
      fv0CutsPart3->SetName("Particle3");
      fOutput->Add(fv0CutsPart3->GetQAHists());
      if(fv0CutsPart3->GetIsMonteCarlo()){
        fv0CutsPart3->SetMCName("MCParticle3");
        fOutput->Add(fv0CutsPart3->GetMCQAHists());
      }
    }

  if(!fv0CutsPart4){
    AliFatal("V0 Cuts for Particle 4 (fv0CutsPart4) not set!\n");
  } else{
      fv0CutsPart4->Init();
      fv0CutsPart4->SetName("Particle4");
      fOutput->Add(fv0CutsPart4->GetQAHists());
      if(fv0CutsPart4->GetIsMonteCarlo()){
        fv0CutsPart4->SetMCName("MCParticle4");
        fOutput->Add(fv0CutsPart4->GetMCQAHists());
      }
    }

  fFemtov0 = new AliFemtoDreamv0();
  fFemtov0->SetPDGCode(fv0CutsPart3->GetPDGv0());   
  fFemtov0->SetUseMCInfo(fIsMC);
  fFemtov0->SetPDGDaughterPos(fv0CutsPart3->GetPDGPosDaug());
  fFemtov0->GetPosDaughter()->SetUseMCInfo(fIsMC); 
  fFemtov0->SetPDGDaughterNeg(fv0CutsPart3->GetPDGNegDaug());
  fFemtov0->GetNegDaughter()->SetUseMCInfo(fIsMC); 

  fEvent = new AliFemtoDreamEvent(false,true,AliVEvent::kINT7);	
    // AliFemtoDreamEvent(1,2,3)
    // 1. argument (boolian) turns on the manual configuration of the pile up rejection (outdated)
    // 2. argument (boolian) provides the QA from the AliEventCuts
    // 3. argument (string) select a trigger

  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());

  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(fIsMC);

  fGTI = new AliVTrack*[fTrackBufferSize];

  fPairCleaner = new AliFemtoDreamPairCleaner(2,2,false);
    // AliFemtoDreamPairCleaner(1,2,3)
    // 1. argument (integer) number of track-decay-combinations to be cleaned (deuteron-lambda and antideuteron-antilambda)
    // 2. argument (integer) number of decay-decay-combinations to be cleaned (lambda-lambda and antilambda-antilambda)
    // 3. argument (boolian) turns on minimal booking, which means that no histograms are created and filled

  fPartColl = new AliFemtoDreamPartCollection(fConfig,false);
    // AliFemtoDreamPartCollection(1,2)
    // 1. argument (object) is the configuration object which is needed for the calculation of the correlation function
    // 2. argument (boolian) turns on minimal booking, which means the QA histograms are not created


  fOutput->Add(fEventCuts->GetHistList()); 
  fOutput->Add(fEvent->GetEvtCutList());
  fOutput->Add(fPairCleaner->GetHistList());
  fOutput->Add(fPartColl->GetHistList());
  fOutput->Add(fPartColl->GetQAList());
  PostData(1,fOutput);

}

//  -----------------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskLeuteronNanoAOD::UserExec(Option_t *){

  AliVEvent *Event = fInputEvent;

  if(!Event){
    AliFatal("No input event!\n");
  } else{
      fEvent->SetEvent(Event);						    // put all the event information into the event object
      if(fEventCuts->isSelected(fEvent)){				    // use the event cut object to check if this event fulfills the criteria
	ResetGlobalTrackReference();

  	for(int iTrack = 0;iTrack<Event->GetNumberOfTracks();++iTrack){	    // fill the fGTI once per event with the mapping of the tracks
	  AliVTrack *track=static_cast<AliVTrack*>(Event->GetTrack(iTrack));
	  if(!track){
	    AliFatal("No standard NanoAOD!\n");
	    return;
	  }
	  StoreGlobalTrackReference(track);
	}

	static std::vector<AliFemtoDreamBasePart> Particles;
	static std::vector<AliFemtoDreamBasePart> AntiParticles;
	static std::vector<AliFemtoDreamBasePart> Decays;
	static std::vector<AliFemtoDreamBasePart> AntiDecays;
		
	Particles.clear();
	AntiParticles.clear();
	Decays.clear();
	AntiDecays.clear();

	fTrack->SetGlobalTrackInfo(fGTI,fTrackBufferSize);

	for(int iTrack = 0;iTrack<Event->GetNumberOfTracks();++iTrack){	    // loop over all tracks in the event
	  AliVTrack *track=static_cast<AliVTrack*>(Event->GetTrack(iTrack));
	  if(!track){
	    AliFatal("No standard NanoAOD!\n");
	  }

	  fTrack->SetTrack(track,Event); 

	  if(fTrackCutsPart1->isSelected(fTrack)){			    // check if the track passes the selection criteria for particle 1
	    Particles.push_back(*fTrack);				    // if so, add it to the particle buffer
	  }
		
	  if(fTrackCutsPart2->isSelected(fTrack)){			    // check if the track passes the selection criteria for particle 2
	    AntiParticles.push_back(*fTrack);				    // if so, add it to the particle buffer
	  }
	}

	TClonesArray *v01 = static_cast<TClonesArray*>(dynamic_cast<AliAODEvent*>(Event)->GetV0s());
	fFemtov0->SetGlobalTrackInfo(fGTI,fTrackBufferSize);

	for(int iv0 = 0;iv0<v01->GetEntriesFast();iv0++){		    // loop over all v0 candidates
	  AliAODv0 *v0 = dynamic_cast<AliAODEvent*>(Event)->GetV0(iv0);
	  fFemtov0->Setv0(Event,v0,fEvent->GetMultiplicity()); 

	  if(fv0CutsPart3->isSelected(fFemtov0)){			    // check if the v0 candidate passes the selection criteria for particle 3
	    Decays.push_back(*fFemtov0);				    // if so, add it to the particle buffer
	  }
		    
	  if(fv0CutsPart4->isSelected(fFemtov0)){			    // check if the v0 candidate passes the selection criteria for particle 4
	    AntiDecays.push_back(*fFemtov0);				    // if so, add it to the particle buffer
	  }
	}

	fPairCleaner->CleanTrackAndDecay(&Particles,&Decays,0);		    // clean deuteron-lambda
	fPairCleaner->CleanTrackAndDecay(&AntiParticles,&AntiDecays,1);	    // clean antideuteron-antilambda

	fPairCleaner->CleanDecay(&Decays,0);				    // clean lambda-lambda
	fPairCleaner->CleanDecay(&AntiDecays,1);			    // clean antilambda-antilambda
  
	fPairCleaner->ResetArray();
	fPairCleaner->StoreParticle(Particles);
	fPairCleaner->StoreParticle(AntiParticles);
	fPairCleaner->StoreParticle(Decays);
	fPairCleaner->StoreParticle(AntiDecays);

	fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),fEvent->GetZVertex(),fEvent->GetRefMult08(),fEvent->GetV0MCentrality());
	  // SetEvent(1,2,3,4)
	  // 1. argument (vector) get vector with cleaned particles
	  // 2. argument (float) z-position of the primary vertex -> used for event mixing to avoid mixing events with different detector acceptances
	  // 3. argument (float) get mutliplicity
	  // 4. argument (float) get centrality in case of p-Pb or Pb-Pb

	PostData(1,fOutput);

	
      }
    }
}


//  -----------------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskLeuteronNanoAOD::ResetGlobalTrackReference(){

  for(int i = 0;i<fTrackBufferSize;i++){				      
    fGTI[i] = 0;							      // set all the pointers to zero
  }

}



//  -----------------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskLeuteronNanoAOD::StoreGlobalTrackReference(AliVTrack *track){

  AliNanoAODTrack *nanoTrack = dynamic_cast<AliNanoAODTrack*>(track);

  const int trackID = track->GetID();

  if(trackID<0){							      // check if the ID of the track is positive
    return;
  }

  if(trackID>=fTrackBufferSize){					      // check if ID is not too big for the buffer
    printf("Warning: track ID is too big for the buffer.\n\tID:\t %d \n\tbuffer:\t %d\n",trackID,fTrackBufferSize);
    return;
  }

  if(fGTI[trackID]){

    if( (!nanoTrack->GetFilterMap()) && (!track->GetTPCNcls()) ){
      return;
    }

    if( (dynamic_cast<AliNanoAODTrack*>(fGTI[trackID])->GetFilterMap()) || (fGTI[trackID]->GetTPCNcls()) ){
      printf("Warning: global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",(fGTI[trackID])->GetTPCNcls(),track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u",dynamic_cast<AliNanoAODTrack*>(fGTI[trackID])->GetFilterMap(),nanoTrack->GetFilterMap());
    }

  }

  (fGTI[trackID]) = track;

}

