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
  fIsHighMultV0(true),
  fBruteForceDebugging(false),
  fDeuteronSideband(false),
  fUpperSideband(false),
  fLowerSideband(false),
  fSignal(false),
  fTrackBufferSize(2000),
  fEventList(nullptr),
  fProtonList(nullptr),
  fAntiprotonList(nullptr),
  fDeuteronList(nullptr),
  fDeuteronMassSqTOF(nullptr),
  fAntideuteronList(nullptr),
  fAntideuteronMassSqTOF(nullptr),
  fLambdaList(nullptr),
  fAntilambdaList(nullptr),
  fPairCleanerList(nullptr),
  fResultsList(nullptr),
  fResultsQAList(nullptr),
  fEvent(nullptr),
  fTrack(nullptr),
  fFemtov0(nullptr),
  fEventCuts(nullptr),
  fTrackCutsPart1(nullptr),
  fTrackCutsPart2(nullptr),
  fTrackCutsPart3(nullptr),
  fTrackCutsPart4(nullptr),
  fv0CutsPart5(nullptr),
  fv0CutsPart6(nullptr),
  fConfig(nullptr),
  fPairCleaner(nullptr),
  fPartColl(nullptr),
  fGTI(nullptr)
{

}


//  -----------------------------------------------------------------------------------------------------------------------------------------
AliAnalysisTaskLeuteronNanoAOD::AliAnalysisTaskLeuteronNanoAOD(const char *name, bool isMC, bool isHighMultV0, bool BruteForceDebugging, bool DeuteronSideband, bool UpperSideband, bool LowerSideband, bool Signal):AliAnalysisTaskSE(name),
  fIsMC(isMC),
  fIsHighMultV0(isHighMultV0),
  fBruteForceDebugging(BruteForceDebugging),
  fDeuteronSideband(DeuteronSideband),
  fUpperSideband(UpperSideband),
  fLowerSideband(LowerSideband),
  fSignal(Signal),
  fTrackBufferSize(2000),
  fEventList(nullptr),
  fProtonList(nullptr),
  fAntiprotonList(nullptr),
  fDeuteronList(nullptr),
  fDeuteronMassSqTOF(nullptr),
  fAntideuteronList(nullptr),
  fAntideuteronMassSqTOF(nullptr),
  fLambdaList(nullptr),
  fAntilambdaList(nullptr),
  fPairCleanerList(nullptr),
  fResultsList(nullptr),
  fResultsQAList(nullptr),
  fEvent(nullptr),
  fTrack(nullptr),
  fFemtov0(nullptr),
  fEventCuts(nullptr),
  fTrackCutsPart1(nullptr),
  fTrackCutsPart2(nullptr),
  fTrackCutsPart3(nullptr),
  fTrackCutsPart4(nullptr),
  fv0CutsPart5(nullptr),
  fv0CutsPart6(nullptr),
  fConfig(nullptr),
  fPairCleaner(nullptr),
  fPartColl(nullptr),
  fGTI(nullptr)
{
  if(BruteForceDebugging){
    std::cout << "x-x-> AliAnalysisTaskNanoAOD: Start defining output in the named constructor" << std::endl;
  }

  DefineOutput(1,TList::Class());   // output for the event cuts
  DefineOutput(2,TList::Class());   // output for the proton cuts
  DefineOutput(3,TList::Class());   // output for the antiproton cuts
  DefineOutput(4,TList::Class());   // output for the deuteron cuts
  DefineOutput(5,TList::Class());   // output for the antideuteron cuts
  DefineOutput(6,TList::Class());   // output for the lambda cuts
  DefineOutput(7,TList::Class());   // output for the antilambda cuts
  DefineOutput(8,TList::Class());   // output for the pair cleaner
  DefineOutput(9,TList::Class());   // output for the results
  DefineOutput(10,TList::Class());  // output for the results QA
  
  if(BruteForceDebugging){
    std::cout << "x-x-> AliAnalysisTaskNanoAOD: Output in the named constructor defined" << std::endl;
  }
}

//  -----------------------------------------------------------------------------------------------------------------------------------------
AliAnalysisTaskLeuteronNanoAOD::~AliAnalysisTaskLeuteronNanoAOD(){	// destructor -> check if the object exists, if so delete it 

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

  if(fTrackCutsPart3){
    delete fTrackCutsPart3;
  }

  if(fTrackCutsPart4){
    delete fTrackCutsPart4;
  }

  if(fv0CutsPart5){
    delete fv0CutsPart5;
  }

  if(fv0CutsPart6){
    delete fv0CutsPart6;
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


  if(fBruteForceDebugging){
    std::cout << "x-x-> AliAnalysisTaskNanoAOD: Begin of UserCreateOutputObjects" << std::endl;
  }

  fResultsQAList = new TList();
  fResultsQAList->SetName("ResultsQA");
  fResultsQAList->SetOwner();

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
      fProtonList = fTrackCutsPart1->GetQAHists();	// add the histograms of the track cuts to the output object

      if(fBruteForceDebugging){
	std::cout << "x-x-> AliAnalysisTaskNanoAOD: fTrackCutsPart1->GetQAHists() done" << std::endl;
      }

      if(fTrackCutsPart1->GetIsMonteCarlo()){		// check if IsMC is "true" or "false"
        fTrackCutsPart1->SetMCName("MCParticle1");	// rename the list carrying the histograms of this object (to avoid collisions in the output object)
        fProtonList->Add(fTrackCutsPart1->GetMCQAHists());  // add the histograms of the Monte Carlo to the output object
      }
    }
  
  if(fBruteForceDebugging){
    std::cout << "x-x-> AliAnalysisTaskNanoAOD: fTrackCutsPart1 (Proton) initialized" << std::endl;
  }

  if(!fTrackCutsPart2){
    AliFatal("Track Cuts for Particle 2 (fTrackCutsPart2) not set!\n");
  } else{
      fTrackCutsPart2->Init();
      fTrackCutsPart2->SetName("Particle2");
      fAntiprotonList = fTrackCutsPart2->GetQAHists();

      if(fBruteForceDebugging){
	std::cout << "x-x-> AliAnalysisTaskNanoAOD: fTrackCutsPart2->GetQAHists() done" << std::endl;
      }

      if(fTrackCutsPart2->GetIsMonteCarlo()){
        fTrackCutsPart2->SetMCName("MCParticle2");
        fAntiprotonList->Add(fTrackCutsPart2->GetMCQAHists());
      }
    }

  if(fBruteForceDebugging){
    std::cout << "x-x-> AliAnalysisTaskNanoAOD: fTrackCutsPart2 (Antiproton) initialized" << std::endl;
  }

  if(!fTrackCutsPart3){
    AliFatal("Track Cuts for Particle 3 (fTrackCutsPart3) not set!\n");
  } else{
      fTrackCutsPart3->Init();
      fTrackCutsPart3->SetName("Particle3");
      fDeuteronList = fTrackCutsPart3->GetQAHists();
      
      if(fBruteForceDebugging){
	std::cout << "x-x-> AliAnalysisTaskNanoAOD: fTrackCutsPart3->GetQAHists() done" << std::endl;
      }

      if(fTrackCutsPart3->GetIsMonteCarlo()){
        fTrackCutsPart3->SetMCName("MCParticle3");
        fDeuteronList->Add(fTrackCutsPart3->GetMCQAHists());
      }
    }

  if(fBruteForceDebugging){
    std::cout << "x-x-> AliAnalysisTaskNanoAOD: fTrackCutsPart3 (Deuteron) initialized" << std::endl;
  }

  if(!fTrackCutsPart4){
    AliFatal("Track Cuts for Particle 4 (fTrackCutsPart4) not set!\n");
  } else{
      fTrackCutsPart4->Init();
      fTrackCutsPart4->SetName("Particle4");
      fAntideuteronList = fTrackCutsPart4->GetQAHists();

      if(fBruteForceDebugging){
	std::cout << "x-x-> AliAnalysisTaskNanoAOD: fTrackCutsPart4->GetQAHists() done" << std::endl;
      }

      if(fTrackCutsPart4->GetIsMonteCarlo()){
        fTrackCutsPart4->SetMCName("MCParticle4");
        fAntideuteronList->Add(fTrackCutsPart4->GetMCQAHists());
      }
    }

  if(fBruteForceDebugging){
    std::cout << "x-x-> AliAnalysisTaskNanoAOD: fTrackCutsPart4 (Antideuteron) initialized" << std::endl;
  }

  if(!fv0CutsPart5){
    AliFatal("V0 Cuts for Particle 5 (fv0CutsPart5) not set!\n");
  } else{
      fv0CutsPart5->Init();
      fv0CutsPart5->SetName("Particle5");
      fLambdaList = fv0CutsPart5->GetQAHists();

      if(fBruteForceDebugging){
	std::cout << "x-x-> AliAnalysisTaskNanoAOD: fv0CutsPart5->GetQAHists() done" << std::endl;
      }

      if(fv0CutsPart5->GetIsMonteCarlo()){
        fv0CutsPart5->SetMCName("MCParticle5");
        fLambdaList->Add(fv0CutsPart5->GetMCQAHists());
      }
    }

  if(fBruteForceDebugging){
    std::cout << "x-x-> AliAnalysisTaskNanoAOD: fv0CutsPart5 (Lambda) initialized" << std::endl;
  }

  if(!fv0CutsPart6){
    AliFatal("V0 Cuts for Particle 6 (fv0CutsPart6) not set!\n");
  } else{
      fv0CutsPart6->Init();
      fv0CutsPart6->SetName("Particle6");
      fAntilambdaList = fv0CutsPart6->GetQAHists();

      if(fBruteForceDebugging){
	std::cout << "x-x-> AliAnalysisTaskNanoAOD: fv0CutsPart6->GetQAHists() done" << std::endl;
      }

      if(fv0CutsPart6->GetIsMonteCarlo()){
        fv0CutsPart6->SetMCName("MCParticle6");
        fAntilambdaList->Add(fv0CutsPart6->GetMCQAHists());
      }
    }

  if(fBruteForceDebugging){
    std::cout << "x-x-> AliAnalysisTaskNanoAOD: fv0CutsPart6 (Antilambda) initialized" << std::endl;
  }

  fFemtov0 = new AliFemtoDreamv0();
  fFemtov0->SetPDGCode(fv0CutsPart5->GetPDGv0());   
  fFemtov0->SetUseMCInfo(fIsMC);
  fFemtov0->SetPDGDaughterPos(fv0CutsPart5->GetPDGPosDaug());
  fFemtov0->GetPosDaughter()->SetUseMCInfo(fIsMC); 
  fFemtov0->SetPDGDaughterNeg(fv0CutsPart5->GetPDGNegDaug());
  fFemtov0->GetNegDaughter()->SetUseMCInfo(fIsMC); 

  if(fIsHighMultV0){

    fEvent = new AliFemtoDreamEvent(false,true,AliVEvent::kHighMultV0);	
    // AliFemtoDreamEvent(1,2,3)
    // 1. argument (boolean) turns on the manual configuration of the pile up rejection (outdated)
    // 2. argument (boolean) provides the QA from the AliEventCuts
    // 3. argument (string) select a trigger

  } else{

    fEvent = new AliFemtoDreamEvent(false,true,AliVEvent::kINT7);	

  }

  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());

  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(fIsMC);

  fGTI = new AliVTrack*[fTrackBufferSize];

  fPairCleaner = new AliFemtoDreamPairCleaner(2,2,false);
    // AliFemtoDreamPairCleaner(1,2,3)
    // 1. argument (integer) number of track-decay-combinations to be cleaned (deuteron-lambda and antideuteron-antilambda)
    // 2. argument (integer) number of decay-decay-combinations to be cleaned (lambda-lambda and antilambda-antilambda)
    // 3. argument (boolean) turns on minimal booking, which means that no histograms are created and filled

  fPartColl = new AliFemtoDreamPartCollection(fConfig,false);
    // AliFemtoDreamPartCollection(1,2)
    // 1. argument (object) is the configuration object which is needed for the calculation of the correlation function
    // 2. argument (boolean) turns on minimal booking, which means the QA histograms are not created


  // Create and fill the deuteron and antideuteron mass2 histograms
  fDeuteronMassSqTOF = new TH2F("fDeuteronMassSqTOF","Deuterons",50,0.0,5.0,400,0.0,8.0);
  fDeuteronMassSqTOF->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fDeuteronMassSqTOF->GetYaxis()->SetTitle("m^{2} (GeV^{2}/c^{4})");
  fDeuteronList->Add(fDeuteronMassSqTOF);

  fAntideuteronMassSqTOF = new TH2F("fAntideuteronMassSqTOF","Antideuterons",50,0.0,5.0,400,0.0,8.0);
  fAntideuteronMassSqTOF->GetXaxis()->SetTitle("p_{T} (GeV/c)");
  fAntideuteronMassSqTOF->GetYaxis()->SetTitle("m^{2} (GeV^{2}/c^{4})");
  fAntideuteronList->Add(fAntideuteronMassSqTOF);

  if(!fEventCuts->GetMinimalBooking()){
    fEventList = fEventCuts->GetHistList();
  } else{
      fEventList = new TList();		    // create a list were all the histograms can be added to
      fEventList->SetName("EventCuts");	    // give the output object name
      fEventList->SetOwner();		    // tell ROOT that this object belongs to the top list / object
    }

  fResultsList = fPartColl->GetHistList();
  fResultsQAList->Add(fPartColl->GetQAList());
  fEventList->Add(fEvent->GetEvtCutList());
  fPairCleanerList = fPairCleaner->GetHistList();

  PostData(1,fEventList);
  PostData(2,fProtonList);
  PostData(3,fAntiprotonList);
  PostData(4,fDeuteronList);
  PostData(5,fAntideuteronList);
  PostData(6,fLambdaList);
  PostData(7,fAntilambdaList);
  PostData(8,fPairCleanerList);
  PostData(9,fResultsList);
  PostData(10,fResultsQAList);

  if(fBruteForceDebugging){
    std::cout << "x-x-> AliAnalysisTaskNanoAOD: End of UserCreateOutputObjects" << std::endl;
  }
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

	double mass2 = 0.0;
	double pT = 0.0;
	double mean = 0.0;
	double limit1 = 0.0;	
	double limit2 = 0.0;

	static std::vector<AliFemtoDreamBasePart> ProtonParticles;
	static std::vector<AliFemtoDreamBasePart> AntiprotonParticles;
	static std::vector<AliFemtoDreamBasePart> DeuteronParticles;
	static std::vector<AliFemtoDreamBasePart> AntideuteronParticles;
	static std::vector<AliFemtoDreamBasePart> Decays;
	static std::vector<AliFemtoDreamBasePart> AntiDecays;
		
	ProtonParticles.clear();
	AntiprotonParticles.clear();
	DeuteronParticles.clear();
	AntideuteronParticles.clear();
	Decays.clear();
	AntiDecays.clear();

	fTrack->SetGlobalTrackInfo(fGTI,fTrackBufferSize);

	for(int iTrack = 0;iTrack<Event->GetNumberOfTracks();++iTrack){	    // loop over all tracks in the event
	  AliVTrack *track=static_cast<AliVTrack*>(Event->GetTrack(iTrack));
	  if(!track){
	    AliFatal("No standard NanoAOD!\n");
	  }

	  fTrack->SetTrack(track,Event); 
	  
	  // protons
	  if(fTrackCutsPart1->isSelected(fTrack)){			    // check if the track passes the selection criteria for particle 1
	    ProtonParticles.push_back(*fTrack);				    // if so, add it to the particle buffer
	  }
	  
	  // antiprotons
	  if(fTrackCutsPart2->isSelected(fTrack)){
	    AntiprotonParticles.push_back(*fTrack);
	  }
		
	  // deuterons
	  if(fTrackCutsPart3->isSelected(fTrack)){

	    mass2 = CalculateMassSqTOF(fTrack);
	    pT = fTrack->GetPt();
	    mean = GetDeuteronMass2Mean_pp(pT);

	    if(fDeuteronSideband){

              // upper sideband
              if(fUpperSideband){

                limit1 = GetLimit(pT,mean,+1,0.24);
                limit2 = GetLimit(pT,limit1,+1,0.24);

              }   

              // lower sideband
              if(fLowerSideband){

                limit2 = GetLimit(pT,mean,-1,0.30);
                limit1 = GetLimit(pT,limit2,-1,0.30);

              }   

              // signal
              if(fSignal){

                limit1 = GetLimit(pT,mean,-1,0.30);
                limit2 = GetLimit(pT,mean,+1,0.24);

              }   

              if(mass2 >= limit1 && mass2 <= limit2){

                DeuteronParticles.push_back(*fTrack);
                fDeuteronMassSqTOF->Fill(fTrack->GetPt(),mass2); 

              }


	    }else{
	      DeuteronParticles.push_back(*fTrack);
	      fDeuteronMassSqTOF->Fill(fTrack->GetPt(),mass2);
	    }
	  }

	  // antideuterons
	  if(fTrackCutsPart4->isSelected(fTrack)){

	    mass2 = CalculateMassSqTOF(fTrack);
	    pT = fTrack->GetPt();
	    mean = GetAntideuteronMass2Mean_pp(pT);

	    if(fDeuteronSideband){

              // upper sideband
              if(fUpperSideband){

                limit1 = GetLimit(pT,mean,+1,0.24);
                limit2 = GetLimit(pT,limit1,+1,0.24);

              }

              // lower sideband
              if(fLowerSideband){

                limit2 = GetLimit(pT,mean,-1,0.30);
                limit1 = GetLimit(pT,limit2,-1,0.30);

              }

              // signal
              if(fSignal){

                limit1 = GetLimit(pT,mean,-1,0.30);
                limit2 = GetLimit(pT,mean,+1,0.24);

              }

              if(mass2 >= limit1 && mass2 <= limit2){

                AntideuteronParticles.push_back(*fTrack);
                fAntideuteronMassSqTOF->Fill(fTrack->GetPt(),mass2); 

              }


 
	    }else{
	      AntideuteronParticles.push_back(*fTrack);
	      fAntideuteronMassSqTOF->Fill(fTrack->GetPt(),CalculateMassSqTOF(fTrack));
	    }
	  }
	}

	TClonesArray *v01 = static_cast<TClonesArray*>(dynamic_cast<AliAODEvent*>(Event)->GetV0s());
	fFemtov0->SetGlobalTrackInfo(fGTI,fTrackBufferSize);

	for(int iv0 = 0;iv0<v01->GetEntriesFast();iv0++){		    // loop over all v0 candidates
	  AliAODv0 *v0 = dynamic_cast<AliAODEvent*>(Event)->GetV0(iv0);
	  fFemtov0->Setv0(Event,v0); 

	  if(fv0CutsPart5->isSelected(fFemtov0)){			    // check if the v0 candidate passes the selection criteria for particle 3
	    Decays.push_back(*fFemtov0);				    // if so, add it to the particle buffer
	  }
		    
	  if(fv0CutsPart6->isSelected(fFemtov0)){			    // check if the v0 candidate passes the selection criteria for particle 4
	    AntiDecays.push_back(*fFemtov0);				    // if so, add it to the particle buffer
	  }
	}
	
	fPairCleaner->CleanTrackAndDecay(&DeuteronParticles,&Decays,0);		    // clean deuteron-lambda
	fPairCleaner->CleanTrackAndDecay(&AntideuteronParticles,&AntiDecays,1);	    // clean antideuteron-antilambda

	fPairCleaner->CleanDecay(&Decays,0);					    // clean lambda-lambda
	fPairCleaner->CleanDecay(&AntiDecays,1);				    // clean antilambda-antilambda
  
	fPairCleaner->ResetArray();
	fPairCleaner->StoreParticle(ProtonParticles);
	fPairCleaner->StoreParticle(AntiprotonParticles);
	fPairCleaner->StoreParticle(DeuteronParticles);
	fPairCleaner->StoreParticle(AntideuteronParticles);
	fPairCleaner->StoreParticle(Decays);
	fPairCleaner->StoreParticle(AntiDecays);

	fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),fEvent->GetZVertex(),fEvent->GetRefMult08(),fEvent->GetV0MCentrality());
	  // SetEvent(1,2,3,4)
	  // 1. argument (vector) get vector with cleaned particles
	  // 2. argument (float) z-position of the primary vertex -> used for event mixing to avoid mixing events with different detector acceptances
	  // 3. argument (float) get mutliplicity
	  // 4. argument (float) get centrality in case of p-Pb or Pb-Pb

	PostData(1,fEventList);
        PostData(2,fProtonList);
        PostData(3,fAntiprotonList);
        PostData(4,fDeuteronList);
        PostData(5,fAntideuteronList);
        PostData(6,fLambdaList);
        PostData(7,fAntilambdaList);
        PostData(8,fPairCleanerList);
        PostData(9,fResultsList);
        PostData(10,fResultsQAList);

	
      }
    }
}


//  -----------------------------------------------------------------------------------------------------------------------------------------
Float_t AliAnalysisTaskLeuteronNanoAOD::CalculateMassSqTOF(AliFemtoDreamTrack *track){

  Float_t p = track->GetP();
  Float_t beta = track->GetbetaTOF();
  Float_t massSq = -999;

  if(beta > 0.0){
    massSq = ((1/(beta*beta))-1) * (p*p);
  }

  return massSq;
}

//  -----------------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskLeuteronNanoAOD::ResetGlobalTrackReference(){

  for(int i = 0;i<fTrackBufferSize;i++){				      
    fGTI[i] = 0;							      // set all the pointers to zero
  }

}

//  -----------------------------------------------------------------------------------------------------------------------------------------
Double_t AliAnalysisTaskLeuteronNanoAOD::GetDeuteronMass2Mean_pp(float pT){

// These values were obtained by fitting the mean values of the deuteron mass2 projections calculated with the NanoAOD dataset (FemtoTUM_LHC16_17_18_HM)
  TF1 *fit = new TF1("fit","[0]+[1]*pow((1-([2]/(x))),[3])",0.6,4.0);
  fit->FixParameter(0,3.53009e+00);
  fit->FixParameter(1,7.09944e-14);
  fit->FixParameter(2,-4.51303e+04);
  fit->FixParameter(3,2.68906e+00);

  Double_t value = fit->Eval(pT);
  fit->Delete();
  return value;

}

//  -----------------------------------------------------------------------------------------------------------------------------------------

Double_t AliAnalysisTaskLeuteronNanoAOD::GetLimit(float pT, double mean, double sign, double offset){

  TF1 *fit = new TF1("fit","[0] + ([1] *([2] + ([3]*(x)) + ([4]*(x)*(x)) + ([5]*(x)*(x)*(x))))",1.0,4.0);
  fit->FixParameter(0,mean);
  fit->FixParameter(1,sign);
  fit->FixParameter(2,offset);
  fit->FixParameter(3,0.003);
  fit->FixParameter(4,0.002);
  fit->FixParameter(5,0.009);

  Double_t value = fit->Eval(pT);
  fit->Delete();
  return value;

}


//  -----------------------------------------------------------------------------------------------------------------------------------------
Double_t AliAnalysisTaskLeuteronNanoAOD::GetAntideuteronMass2Mean_pp(float pT){

// These values were obtained by fitting the mean values of the antideuteron mass2 projections calculated with the NanoAOD dataset (FemtoTUM_LHC16_17_18_HM)
  TF1 *fit = new TF1("fit","[0]+[1]*pow((1-([2]/(x))),[3])",0.6,4.0);
  fit->FixParameter(0,3.57317e+00);
  fit->FixParameter(1,6.63106e-15);
  fit->FixParameter(2,-6.62276e+03);
  fit->FixParameter(3,3.51473e+00);

  Double_t value = fit->Eval(pT);
  fit->Delete();
  return value;

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

