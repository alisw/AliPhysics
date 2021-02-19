/*
 * AliAnalysisTaskLeuteronAOD.cxx
 *
 *  Created on:	19 December 2019
 *	Author:	Michael Jung
 */

#include "AliAnalysisTaskLeuteronAOD.h"
#include "AliAODTrack.h"

ClassImp(AliAnalysisTaskLeuteronAOD)

//  -----------------------------------------------------------------------------------------------------------------------------------------
AliAnalysisTaskLeuteronAOD::AliAnalysisTaskLeuteronAOD():AliAnalysisTaskSE(),
  fIsMC(false),
  fIsHighMultV0(true),
  fBruteForceDebugging(false),
  fDeuteronSideband(false),
  fDeuteronSigmaLeft(0.0),
  fDeuteronSigmaRight(0.0),
  fAntideuteronSigmaLeft(0.0),
  fAntideuteronSigmaRight(0.0),
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
AliAnalysisTaskLeuteronAOD::AliAnalysisTaskLeuteronAOD(const char *name, bool isMC, bool isHighMultV0, bool BruteForceDebugging,bool DeuteronSideband, double DeuteronSigmaLeft, double DeuteronSigmaRight, double AntideuteronSigmaLeft, double AntideuteronSigmaRight):AliAnalysisTaskSE(name),
  fIsMC(isMC),
  fIsHighMultV0(isHighMultV0),
  fBruteForceDebugging(BruteForceDebugging),
  fDeuteronSideband(DeuteronSideband),
  fDeuteronSigmaLeft(DeuteronSigmaLeft),
  fDeuteronSigmaRight(DeuteronSigmaRight),
  fAntideuteronSigmaLeft(DeuteronSigmaLeft),
  fAntideuteronSigmaRight(DeuteronSigmaRight),
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
    std::cout << "x-x-> AliAnalysisTaskAOD: Start defining output in the named constructor" << std::endl;
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
    std::cout << "x-x-> AliAnalysisTaskAOD: Output in the named constructor defined" << std::endl;
  }

}

//  -----------------------------------------------------------------------------------------------------------------------------------------
AliAnalysisTaskLeuteronAOD::~AliAnalysisTaskLeuteronAOD(){	// destructor -> check if the object exists, if so delete it 

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
void AliAnalysisTaskLeuteronAOD::UserCreateOutputObjects(){

  if(fBruteForceDebugging){
    std::cout << "x-x-> AliAnalysisTaskAOD: Begin of UserCreateOutputObjects" << std::endl;
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

  if(!fTrackCutsPart1){						      // check if the track cuts object is set, if not, call it a day
    AliFatal("Track Cuts for Protons (fTrackCutsPart1) not set!\n");
  } else{
      fTrackCutsPart1->Init();					      // initialize the histograms of this object 
      fTrackCutsPart1->SetName("Proton");			      // rename the list carrying the histograms of this object (avoid collisions in the output object)
      fProtonList = fTrackCutsPart1->GetQAHists();		      // add the histograms of the track cuts to the output object

      if(fBruteForceDebugging){
	std::cout << "x-x-> AliAnalysisTaskAOD: fTrackCutsPart1->GetQAHists() done" << std::endl;
      }

      if(fTrackCutsPart1->GetIsMonteCarlo()){			      // check if IsMC is "true" or "false"
        fTrackCutsPart1->SetMCName("MCProton");			      // rename the list carrying the histograms of this object (avoid collisions in the output object)
        fProtonList->Add(fTrackCutsPart1->GetMCQAHists());	      // add the histograms of the Monte Carlo to the output object
      }
  }

  if(fBruteForceDebugging){
    std::cout << "x-x-> AliAnalysisTaskAOD: fTrackCutsPart1 (Proton) initialized" << std::endl;
  }

  if(!fTrackCutsPart2){
    AliFatal("Track Cuts for Antiprotons (fTrackCutsPart2) not set!\n");
  } else{
      fTrackCutsPart2->Init();
      fTrackCutsPart2->SetName("Antiprotons");
      fAntiprotonList = fTrackCutsPart2->GetQAHists();
      
      if(fBruteForceDebugging){
	std::cout << "x-x-> AliAnalysisTaskAOD: fTrackCutsPart2->GetQAHists() done" << std::endl;
      }

      if(fTrackCutsPart2->GetIsMonteCarlo()){
        fTrackCutsPart2->SetMCName("MCAntiprotons");
        fAntiprotonList->Add(fTrackCutsPart2->GetMCQAHists());
      }
  }
  
  if(fBruteForceDebugging){
    std::cout << "x-x-> AliAnalysisTaskAOD: fTrackCutsPart2 (Antiproton) initialized" << std::endl;
  }

  if(!fTrackCutsPart3){
    AliFatal("Track Cuts for Deuterons (fTrackCutsPart3) not set!\n");
  } else{
      fTrackCutsPart3->Init();
      fTrackCutsPart3->SetName("Deuteron");
      fDeuteronList = fTrackCutsPart3->GetQAHists();
      
      if(fBruteForceDebugging){
	std::cout << "x-x-> AliAnalysisTaskAOD: fTrackCutsPart3->GetQAHists() done" << std::endl;
      }

      if(fTrackCutsPart3->GetIsMonteCarlo()){
        fTrackCutsPart3->SetMCName("MCDeuteron");
        fDeuteronList->Add(fTrackCutsPart3->GetMCQAHists());
      }
    }
  
  if(fBruteForceDebugging){
    std::cout << "x-x-> AliAnalysisTaskAOD: fTrackCutsPart3 (Deuteron) initialized" << std::endl;
  }

  if(!fTrackCutsPart4){
    AliFatal("Track Cuts for Antideuterons (fTrackCutsPart4) not set!\n");
  } else{
      fTrackCutsPart4->Init();
      fTrackCutsPart4->SetName("Antideuteron");
      fAntideuteronList = fTrackCutsPart4->GetQAHists();
      
      if(fBruteForceDebugging){
	std::cout << "x-x-> AliAnalysisTaskAOD: fTrackCutsPart4->GetQAHists() done" << std::endl;
      }

      if(fTrackCutsPart4->GetIsMonteCarlo()){
        fTrackCutsPart4->SetMCName("MCAntideuteron");
        fAntideuteronList->Add(fTrackCutsPart4->GetMCQAHists());
      }
  }

  if(fBruteForceDebugging){
    std::cout << "x-x-> AliAnalysisTaskAOD: fTrackCutsPart4 (Antideuteron) initialized" << std::endl;
  }

  if(!fv0CutsPart5){
    AliFatal("V0 Cuts for Lambdas (fv0CutsPart5) not set!\n");
  } else{
      fv0CutsPart5->Init();
      fv0CutsPart5->SetName("Lambda");
      fLambdaList = fv0CutsPart5->GetQAHists();
      if(fv0CutsPart5->GetIsMonteCarlo()){
        fv0CutsPart5->SetMCName("MCLambda");
        fLambdaList = fv0CutsPart5->GetMCQAHists();
      }
    }
  
  if(fBruteForceDebugging){
    std::cout << "x-x-> AliAnalysisTaskAOD: fv0CutsPart5 (Lambda) initialized" << std::endl;
  }

  if(!fv0CutsPart6){
    AliFatal("V0 Cuts for Antilambdas (fv0CutsPart6) not set!\n");
  } else{
      fv0CutsPart6->Init();
      fv0CutsPart6->SetName("Antilambda");
      fAntilambdaList = fv0CutsPart6->GetQAHists();
      if(fv0CutsPart6->GetIsMonteCarlo()){
        fv0CutsPart6->SetMCName("MCAntilambda");
        fAntilambdaList = fv0CutsPart6->GetMCQAHists();
      }
    }

  if(fBruteForceDebugging){
    std::cout << "x-x-> AliAnalysisTaskAOD: fv0CutsPart6 (Antilambda) initialized" << std::endl;
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

  fGTI = new AliAODTrack*[fTrackBufferSize];

  fPairCleaner = new AliFemtoDreamPairCleaner(2,2,false);
    // AliFemtoDreamPairCleaner(1,2,3)
    // 1. argument (integer) number of track-decay-combinations to be cleaned (proton-lambda, antiproton-antilambda, deuteron-lambda and antideuteron-antilambda)
    // 2. argument (integer) number of decay-decay-combinations to be cleaned (lambda-lambda and antilambda-antilambda)
    // 3. argument (boolean) turns on minimal booking, which means that no histograms are created and filled

  fPartColl = new AliFemtoDreamPartCollection(fConfig,false);
    // AliFemtoDreamPartCollection(1,2)
    // 1. argument (object) is the configuration object which is needed for the calculation of the correlation function
    // 2. argument (boolean) turns on minimal booking, which means the QA histograms are not created

  if(!fEventCuts->GetMinimalBooking()){
    fEventList = fEventCuts->GetHistList();
  } else{
      fEventList = new TList();		    // create a list were all the histograms can be added to
      fEventList->SetName("EventCuts");	    // give the output object name
      fEventList->SetOwner();		    // tell ROOT that this object belongs to the top list / object
    }

  // Create and fill the deuteron and antideuteron mass2 histograms
  fDeuteronMassSqTOF = new TH2F("fDeuteronMassSqTOF","Deuterons",50,0.0,5.0,400,0.0,8.0);
  fDeuteronMassSqTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fDeuteronMassSqTOF->GetYaxis()->SetTitle("#it{m}^{2} (GeV^{2}/#it{c}^{4})");
  fDeuteronList->Add(fDeuteronMassSqTOF);

  fAntideuteronMassSqTOF = new TH2F("fAntideuteronMassSqTOF","Antideuterons",50,0.0,5.0,400,0.0,8.0);
  fAntideuteronMassSqTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fAntideuteronMassSqTOF->GetYaxis()->SetTitle("#it{m}^{2} (GeV^{2}/#it{c}^{4})");
  fAntideuteronList->Add(fAntideuteronMassSqTOF);

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
    std::cout << "x-x-> AliAnalysisTaskAOD: End of UserCreateOutputObjects" << std::endl;
  }


}

//  -----------------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskLeuteronAOD::UserExec(Option_t *){

  AliAODEvent *Event=static_cast<AliAODEvent*>(fInputEvent);

  if(!Event){
    AliFatal("No input event!\n");
  } else{
      fEvent->SetEvent(Event);						    // put all the event information into the event object

      if(fEventCuts->isSelected(fEvent)){				    // use the event cut object to check if this event fulfills the criteria
	ResetGlobalTrackReference();

  	for(int iTrack = 0;iTrack<Event->GetNumberOfTracks();++iTrack){	    // fill the fGTI once per event with the mapping of the tracks
	  AliAODTrack *track=static_cast<AliAODTrack*>(Event->GetTrack(iTrack));
	  if(!track){
	    AliFatal("No standard AOD!\n");
	    return;
	  }
	  StoreGlobalTrackReference(track);
	}

	double mass2;
	double mean;
	double sigma;

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
	  AliAODTrack *track=static_cast<AliAODTrack*>(Event->GetTrack(iTrack));
	  if(!track){
	    AliFatal("No standard AOD!\n");
	  }

	  fTrack->SetTrack(track); 

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

	    // deuterons (sideband only)
	    if(fDeuteronSideband){

	      mass2 = CalculateMassSqTOF(fTrack); 
	      mean = GetDeuteronMass2Mean_pp(fTrack->GetPt());
	      sigma = GetDeuteronMass2Sigma_pp(fTrack->GetPt());

	      if((mass2 >= mean + (fDeuteronSigmaLeft*sigma)) && (mass2 <= (mean + fDeuteronSigmaRight*sigma))){
                DeuteronParticles.push_back(*fTrack);
		fDeuteronMassSqTOF->Fill(fTrack->GetPt(),mass2); 
	      }
	    // deuterons (all particles)
	    }else{
	      mass2 = CalculateMassSqTOF(fTrack); 
              DeuteronParticles.push_back(*fTrack);
              fDeuteronMassSqTOF->Fill(fTrack->GetPt(),mass2);
            }   
          }   
	
	  // antideuterons
	  if(fTrackCutsPart4->isSelected(fTrack)){

	    // antideuterons (sideband only)
	    if(fDeuteronSideband){
            
	      mass2 = CalculateMassSqTOF(fTrack);
	      mean = GetAntideuteronMass2Mean_pp(fTrack->GetPt());
	      sigma = GetAntideuteronMass2Sigma_pp(fTrack->GetPt());
              
		if((mass2 >= mean + (fAntideuteronSigmaLeft*sigma)) && (mass2 <= mean + (fAntideuteronSigmaRight*sigma))){
		  AntideuteronParticles.push_back(*fTrack);
		  fAntideuteronMassSqTOF->Fill(fTrack->GetPt(),mass2);
		}
	    // antideuterons (all particles)
	    }else{
	      mass2 = CalculateMassSqTOF(fTrack);
              AntideuteronParticles.push_back(*fTrack);
              fAntideuteronMassSqTOF->Fill(fTrack->GetPt(),CalculateMassSqTOF(fTrack));
            }
          }

	}

	TClonesArray *v01 = static_cast<TClonesArray*>(Event->GetV0s());
	fFemtov0->SetGlobalTrackInfo(fGTI,fTrackBufferSize);

	for(int iv0 = 0;iv0<v01->GetEntriesFast();iv0++){		    // loop over all v0 candidates
	  AliAODv0 *v0 = Event->GetV0(iv0);
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
Float_t AliAnalysisTaskLeuteronAOD::CalculateMassSqTOF(AliFemtoDreamTrack *track){

  Float_t p = track->GetP();
  Float_t beta = track->GetbetaTOF();
  Float_t massSq = -999;

  if(beta > 0.0){
    massSq = ((1/(beta*beta))-1) * (p*p);
  }
                                                                                                                                                                      
  return massSq;
}




//  -----------------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskLeuteronAOD::ResetGlobalTrackReference(){

  for(int i = 0;i<fTrackBufferSize;i++){				      
    fGTI[i] = 0;							      // set all the pointers to zero
  }

}

//  -----------------------------------------------------------------------------------------------------------------------------------------
Double_t AliAnalysisTaskLeuteronAOD::GetDeuteronMass2Mean_pp(float pT){

// These values were obtained by fitting the mean values of the deuteron mass2 projections calculated with the AOD dataset (2016,2017,2018)
  TF1 *fit = new TF1("fit","[0]+[1]*pow((1-([2]/(x))),[3])",0.6,4.0);
  fit->FixParameter(0,3.52478e+00);
  fit->FixParameter(1,1.35169e-13);
  fit->FixParameter(2,-7.78509e+04);
  fit->FixParameter(3,2.50751e+00);

  Double_t value = fit->Eval(pT);
  fit->Delete();
  return value;

}

//  -----------------------------------------------------------------------------------------------------------------------------------------
Double_t AliAnalysisTaskLeuteronAOD::GetDeuteronMass2Sigma_pp(float pT){

// These values were obtained by fitting the sigma values of the deuteron mass2 projections calculated with the AOD dataset (2016,2017,2018)
  TF1 *fit = new TF1("fit","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x",0.6,4.0);
  fit->FixParameter(0,1.03688e+00);
  fit->FixParameter(1,-2.14030e+00);
  fit->FixParameter(2,1.86493e+00);
  fit->FixParameter(3,-7.79532e-01);
  fit->FixParameter(4,1.58122e-01);
  fit->FixParameter(5,-1.23721e-02);

  Double_t value = fit->Eval(pT);
  fit->Delete();
  return value;

}

//  -----------------------------------------------------------------------------------------------------------------------------------------
Double_t AliAnalysisTaskLeuteronAOD::GetAntideuteronMass2Mean_pp(float pT){

// These values were obtained by fitting the mean values of the antideuteron mass2 projections calculated with the AOD dataset (2016,2017,2018)
  TF1 *fit = new TF1("fit","[0]+[1]*pow((1-([2]/(x))),[3])",0.6,4.0);
  fit->FixParameter(0,3.54817e+00);
  fit->FixParameter(1,1.97729e-14);
  fit->FixParameter(2,-1.35810e+04);
  fit->FixParameter(3,3.15295e+00);

  Double_t value = fit->Eval(pT);
  fit->Delete();
  return value;

}

//  -----------------------------------------------------------------------------------------------------------------------------------------
Double_t AliAnalysisTaskLeuteronAOD::GetAntideuteronMass2Sigma_pp(float pT){

// These values were obtained by fitting the sigma values of the antideuteron mass2 projections calculated with the AOD dataset (2016,2017,2018)
  TF1 *fit = new TF1("fit","[0]+[1]*x+[2]*x*x+[3]*x*x*x+[4]*x*x*x*x+[5]*x*x*x*x*x",0.6,4.0);
  fit->FixParameter(0,1.29360e+00);
  fit->FixParameter(1,-2.78238e+00);
  fit->FixParameter(2,2.45586e+00);
  fit->FixParameter(3,-1.02786e+00);
  fit->FixParameter(4,2.06204e-01);
  fit->FixParameter(5,-1.57982e-02);

  Double_t value = fit->Eval(pT);
  fit->Delete();
  return value;

}





//  -----------------------------------------------------------------------------------------------------------------------------------------
void AliAnalysisTaskLeuteronAOD::StoreGlobalTrackReference(AliAODTrack *track){

  const int trackID = track->GetID();

  if(trackID<0){							      // check if the ID of the track is positive
    return;
  }

  if(trackID>=fTrackBufferSize){					      // check if ID is not too big for the buffer
    printf("Warning: track ID is too big for the buffer.\n\tID:\t %d \n\tbuffer:\t %d\n",trackID,fTrackBufferSize);
    return;
  }

  if(fGTI[trackID]){

    if( (!track->GetFilterMap()) && (!track->GetTPCNcls()) ){
      return;
    }

    if( (fGTI[trackID])->GetFilterMap() || (fGTI[trackID])->GetTPCNcls() ){
      printf("Warning: global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",(fGTI[trackID])->GetTPCNcls(),track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u",fGTI[trackID]->GetFilterMap(),track->GetFilterMap());
    }

  }

  (fGTI[trackID]) = track;	  // assign the pointer

}
