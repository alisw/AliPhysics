/*
 * AliAnalysisTaskLeuteronAOD.cxx
 *
 *  Created on:	19 December 2019
 *	Author:	Michael Jung
 */

#include "AliAnalysisTaskLeuteronAOD.h"
#include "AliAODTrack.h"
#include "AliMultSelection.h"

ClassImp(AliAnalysisTaskLeuteronAOD)

//  -----------------------------------------------------------------------------------------------------------------------------------------
AliAnalysisTaskLeuteronAOD::AliAnalysisTaskLeuteronAOD():AliAnalysisTaskSE(),
  fIsMC(false),
  fIsHighMultV0(true),
  fBruteForceDebugging(false),
  fisSidebandSignal(false),
  fisUpperSideband(false),
  fisLowerSideband(false),
  fisPbPb(false),
  fisCentral(false),
  fTrackBufferSize(50000),
  fEventList(nullptr),
  fProtonList(nullptr),
  fAntiprotonList(nullptr),
  fDeuteronList(nullptr),
  fDeuteronMassSqTOF(nullptr),
  fDeuteronMassSqTOFFullPt(nullptr),
  fDeuteronTPCnSigma(nullptr),
  fAntideuteronList(nullptr),
  fAntideuteronMassSqTOF(nullptr),
  fAntideuteronMassSqTOFFullPt(nullptr),
  fAntideuteronTPCnSigma(nullptr),
  fLambdaList(nullptr),
  fAntilambdaList(nullptr),
  fPairCleanerList(nullptr),
  fResultsList(nullptr),
  fResultsQAList(nullptr),
  fSimpleEventCounter(nullptr),
  fEventCentrality(nullptr),
  fEvent(nullptr),
  fTrack(nullptr),
  fFemtov0(nullptr),
  fEventCuts(nullptr),
  fTrackCutsPart1(nullptr),
  fTrackCutsPart2(nullptr),
  fTrackCutsPart3(nullptr),
  fTrackCutsPart3Mass(nullptr),
  fTrackCutsPart3Sigma(nullptr),
  fTrackCutsPart4(nullptr),
  fTrackCutsPart4Mass(nullptr),
  fTrackCutsPart4Sigma(nullptr),
  fv0CutsPart5(nullptr),
  fv0CutsPart6(nullptr),
  fConfig(nullptr),
  fEnableEventQAPlots(false),
  fEnableResultsQAPlots(false),
  fPairCleaner(nullptr),
  fPartColl(nullptr),
  fGTI(nullptr)
{
 
}


//  -----------------------------------------------------------------------------------------------------------------------------------------
AliAnalysisTaskLeuteronAOD::AliAnalysisTaskLeuteronAOD(const char *name, bool isMC, bool isHighMultV0, bool BruteForceDebugging,bool isSidebandSignal, bool isUpperSideband, bool isLowerSideband,bool isPbPb,bool doEventQAPlots,bool doResultsQAPlots,bool isCentral):AliAnalysisTaskSE(name),
  fIsMC(isMC),
  fIsHighMultV0(isHighMultV0),
  fBruteForceDebugging(BruteForceDebugging),
  fisSidebandSignal(isSidebandSignal),
  fisUpperSideband(isUpperSideband),
  fisLowerSideband(isLowerSideband),
  fisPbPb(isPbPb),
  fisCentral(isCentral),
  fTrackBufferSize(50000),
  fEventList(nullptr),
  fProtonList(nullptr),
  fAntiprotonList(nullptr),
  fDeuteronList(nullptr),
  fDeuteronMassSqTOF(nullptr),
  fDeuteronMassSqTOFFullPt(nullptr),
  fDeuteronTPCnSigma(nullptr),
  fAntideuteronList(nullptr),
  fAntideuteronMassSqTOF(nullptr),
  fAntideuteronMassSqTOFFullPt(nullptr),
  fAntideuteronTPCnSigma(nullptr),
  fLambdaList(nullptr),
  fAntilambdaList(nullptr),
  fPairCleanerList(nullptr),
  fResultsList(nullptr),
  fResultsQAList(nullptr),
  fSimpleEventCounter(nullptr),
  fEventCentrality(nullptr),
  fEvent(nullptr),
  fTrack(nullptr),
  fFemtov0(nullptr),
  fEventCuts(nullptr),
  fTrackCutsPart1(nullptr),
  fTrackCutsPart2(nullptr),
  fTrackCutsPart3(nullptr),
  fTrackCutsPart3Mass(nullptr),
  fTrackCutsPart3Sigma(nullptr),
  fTrackCutsPart4(nullptr),
  fTrackCutsPart4Mass(nullptr),
  fTrackCutsPart4Sigma(nullptr),
  fv0CutsPart5(nullptr),
  fv0CutsPart6(nullptr),
  fConfig(nullptr),
  fEnableEventQAPlots(doEventQAPlots),
  fEnableResultsQAPlots(doResultsQAPlots),
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

  if(fTrackCutsPart3Mass){
    delete fTrackCutsPart3Mass;
  }

  if(fTrackCutsPart3Sigma){
    delete fTrackCutsPart3Sigma;
  }

  if(fTrackCutsPart4){
    delete fTrackCutsPart4;
  }

  if(fTrackCutsPart4Mass){
    delete fTrackCutsPart4Mass;
  }

  if(fTrackCutsPart4Sigma){
    delete fTrackCutsPart4Sigma;
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


  if(!fTrackCutsPart3Mass){
    AliFatal("Track Cuts for DeuteronsMass (fTrackCutsPart3Mass) not set!\n");
  } else{
      fTrackCutsPart3Mass->Init();
      fTrackCutsPart3Mass->SetName("DeuteronMass");

      if(fTrackCutsPart3Mass->GetIsMonteCarlo()){
        fTrackCutsPart3Mass->SetMCName("MCDeuteronMass");
      }
   }


  if(!fTrackCutsPart3Sigma){
    AliFatal("Track Cuts for DeuteronsSigma (fTrackCutsPart3Sigma) not set!\n");
  } else{
      fTrackCutsPart3Sigma->Init();
      fTrackCutsPart3Sigma->SetName("DeuteronSigma");

      if(fTrackCutsPart3Sigma->GetIsMonteCarlo()){
        fTrackCutsPart3Sigma->SetMCName("MCDeuteronSigma");
      }
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

  if(!fTrackCutsPart4Mass){
    AliFatal("Track Cuts for DeuteronsMass (fTrackCutsPart4Mass) not set!\n");
  } else{
      fTrackCutsPart4Mass->Init();
      fTrackCutsPart4Mass->SetName("AntiDeuteronMass");

      if(fTrackCutsPart4Mass->GetIsMonteCarlo()){
        fTrackCutsPart4Mass->SetMCName("MCAntiDeuteronMass");
      }
   }


  if(!fTrackCutsPart4Sigma){
    AliFatal("Track Cuts for AntiDeuteronsSigma (fTrackCutsPart4Sigma) not set!\n");
  } else{
      fTrackCutsPart4Sigma->Init();
      fTrackCutsPart4Sigma->SetName("DeuteronSigma");

      if(fTrackCutsPart4Sigma->GetIsMonteCarlo()){
        fTrackCutsPart4Sigma->SetMCName("MCAntiDeuteronSigma");
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
    fEvent = new AliFemtoDreamEvent(false,fEnableEventQAPlots,AliVEvent::kHighMultV0);	
    // AliFemtoDreamEvent(1,2,3)
    // 1. argument (boolean) turns on the manual configuration of the pile up rejection (outdated)
    // 2. argument (boolean) provides the QA from the AliEventCuts
    // 3. argument (string) select a trigger
  } else{
    fEvent = new AliFemtoDreamEvent(false,fEnableEventQAPlots,AliVEvent::kINT7);	
  }

  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());

  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(fIsMC);

  fGTI = new AliAODTrack*[fTrackBufferSize];

  fPairCleaner = new AliFemtoDreamPairCleaner(12,2,false);
    // AliFemtoDreamPairCleaner(1,2,3)
    // 1. argument (integer) number of track-decay-combinations to be cleaned (proton-lambda, antiproton-antilambda, deuteron-lambda and antideuteron-antilambda)
    // 2. argument (integer) number of decay-decay-combinations to be cleaned (lambda-lambda and antilambda-antilambda)
    // 3. argument (boolean) turns on minimal booking, which means that no histograms are created and filled

  fPartColl = new AliFemtoDreamPartCollection(fConfig,!fEnableResultsQAPlots);
    // AliFemtoDreamPartCollection(1,2)
    // 1. argument (object) is the configuration object which is needed for the calculation of the correlation function
    // 2. argument (boolean) turns on minimal booking, which means the ResultsQA histograms are not created

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

  fDeuteronMassSqTOFFullPt = new TH2F("fDeuteronMassSqTOFFullPt","Deuterons",50,0.0,5.0,400,0.0,8.0);
  fDeuteronMassSqTOFFullPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fDeuteronMassSqTOFFullPt->GetYaxis()->SetTitle("#it{m}^{2} (GeV^{2}/#it{c}^{4})");
  fDeuteronList->Add(fDeuteronMassSqTOFFullPt);

  fDeuteronTPCnSigma = new TH2F("fDeuteronTPCnSigma","Deuterons",120,0.0,6.0,1200,-60.0,60.0);
  fDeuteronTPCnSigma->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fDeuteronTPCnSigma->GetYaxis()->SetTitle("#it{n#sigma}_{TPC}");
  fDeuteronList->Add(fDeuteronTPCnSigma);

  fAntideuteronMassSqTOF = new TH2F("fAntideuteronMassSqTOF","Antideuterons",50,0.0,5.0,400,0.0,8.0);
  fAntideuteronMassSqTOF->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fAntideuteronMassSqTOF->GetYaxis()->SetTitle("#it{m}^{2} (GeV^{2}/#it{c}^{4})");
  fAntideuteronList->Add(fAntideuteronMassSqTOF);

  fAntideuteronMassSqTOFFullPt = new TH2F("fAntideuteronMassSqTOFFullPt","Antideuterons",50,0.0,5.0,400,0.0,8.0);
  fAntideuteronMassSqTOFFullPt->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
  fAntideuteronMassSqTOFFullPt->GetYaxis()->SetTitle("#it{m}^{2} (GeV^{2}/#it{c}^{4})");
  fAntideuteronList->Add(fAntideuteronMassSqTOFFullPt);

  fAntideuteronTPCnSigma = new TH2F("fDeuteronTPCnSigma","Deuterons",120,0.0,6.0,1200,-60.0,60.0);
  fAntideuteronTPCnSigma->GetXaxis()->SetTitle("#it{p}_{TPC} (GeV/#it{c})");
  fAntideuteronTPCnSigma->GetYaxis()->SetTitle("#it{n#sigma}_{TPC}");
  fAntideuteronList->Add(fAntideuteronTPCnSigma);

  fSimpleEventCounter = new TH1F("fSimpleEventCounter","Number of events",1,0.0,1.0);
  fSimpleEventCounter->GetYaxis()->SetTitle("number of events");

  fEventCentrality = new TH1F("fEventCentrality","Centrality",100,0.0,100.0);
  fEventCentrality->GetXaxis()->SetTitle("centrality (%)");
  fEventCentrality->GetYaxis()->SetTitle("counts");

  fEventList->Add(fSimpleEventCounter);
  fEventList->Add(fEventCentrality);
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

    double centrality_min = 0.0;
    double centrality_max = 100.0;
    double centrality = -999.0;

    if(fisCentral == true){

      centrality_min = 0.0;
      centrality_max = 10.0;

    }else{

      centrality_min = 30.0;
      centrality_max = 50.0;

    }


    AliMultSelection *MultSelection = (AliMultSelection*) Event->FindListObject("MultSelection");
    centrality = MultSelection->GetMultiplicityPercentile("V0M");

    if(fisPbPb){

      if((centrality < centrality_min) || (centrality > centrality_max)) return;

    }

    // Fill centrality histo
    fEventCentrality->Fill(centrality);

    // fill the simple event counter
    fSimpleEventCounter->Fill(0.5);

	double mass2 = 0.0;
	double pT = 0.0;
	//double mean = 0.0;
	//double limit1 = 0.0;
	//double limit2 = 0.0;
	//double buffer = 0.0;
	//double SidebandWidth = 0.2;

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
	  pT = fTrack->GetPt();

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

	  if(pT > 1.0){
	    bool isDeuteron = CheckDeuteronMassSquarePID(fTrack,3);
	    if(isDeuteron == false) continue;
	  }else{
  
	    bool isDeuteron = CheckTPCDeuteronPID(fTrack,3);
	    if(isDeuteron == false) continue;
	  }

	  mass2 = CalculateMassSqTOF(fTrack); 
          DeuteronParticles.push_back(*fTrack);
          fDeuteronMassSqTOF->Fill(fTrack->GetPt(),mass2);

	  }

	 /* 
	    // deuterons (sideband analysis)
	    if((fisSidebandSignal == true) || (fisLowerSideband == true) || (fisUpperSideband == true)){

	      mass2 = CalculateMassSqTOF(fTrack); 
	      pT = fTrack->GetPt();
	      mean = GetDeuteronMass2Mean_pp(pT);

	      // upper sideband
	      if(fisUpperSideband){

		limit1 = GetLimit(pT,mean,+1,0.30,0.009);
		buffer = GetLimit(pT,mean,+1,0.24,0.0065);
		limit2 = GetLimit(pT,buffer,+1,SidebandWidth,0.0065);

	      }

	      // lower sideband
	      if(fisLowerSideband){

		limit2 = GetLimit(pT,mean,-1,0.36,0.009);
		buffer = GetLimit(pT,mean,-1,0.30,0.0065);
		limit1 = GetLimit(pT,buffer,-1,SidebandWidth,0.0065);

	      }

	      // signal
	      if(fisSidebandSignal){

		limit1 = GetLimit(pT,mean,-1,0.30,0.009);
		limit2 = GetLimit(pT,mean,+1,0.24,0.009);

	      }

	      if(mass2 >= limit1 && mass2 <= limit2){

                DeuteronParticles.push_back(*fTrack);
		fDeuteronMassSqTOF->Fill(fTrack->GetPt(),mass2); 

	      }

	    // deuterons (all particles)
	    }else{

	      mass2 = CalculateMassSqTOF(fTrack); 
              DeuteronParticles.push_back(*fTrack);
              fDeuteronMassSqTOF->Fill(fTrack->GetPt(),mass2);

            }   
	  

	  }	*/

	  
	  // deuteron mass
	  if(fTrackCutsPart3Mass->isSelected(fTrack)){

	      mass2 = CalculateMassSqTOF(fTrack); 
              fDeuteronMassSqTOFFullPt->Fill(fTrack->GetPt(),mass2);

	  }

	  // deuteron sigma
	  if(fTrackCutsPart3Sigma->isSelected(fTrack)){

              fDeuteronTPCnSigma->Fill(fTrack->GetMomTPC(),fTrack->GetnSigmaTPC((int) (AliPID::kDeuteron)));

	  }
  
	
	  // antideuterons
	  if(fTrackCutsPart4->isSelected(fTrack)){

	  if(pT > 1.0){
	    bool isAntiDeuteron = CheckAntiDeuteronMassSquarePID(fTrack,3);
	    if(isAntiDeuteron == false) continue;
	  }else{
  	    bool isAntiDeuteron = CheckTPCDeuteronPID(fTrack,3);
	    if(isAntiDeuteron == false) continue;
	  }

	    mass2 = CalculateMassSqTOF(fTrack);
	    AntideuteronParticles.push_back(*fTrack);
	    fAntideuteronMassSqTOF->Fill(fTrack->GetPt(),CalculateMassSqTOF(fTrack));


	  }


	  /*
	    // antideuterons (sideband only)
	    if((fisSidebandSignal == true) || (fisLowerSideband == true) || (fisUpperSideband == true)){

	      mass2 = CalculateMassSqTOF(fTrack); 
	      pT = fTrack->GetPt();
	      mean = GetAntideuteronMass2Mean_pp(pT);

	      // upper sideband
	      if(fisUpperSideband){

		limit1 = GetLimit(pT,mean,+1,0.30,0.009);
		buffer = GetLimit(pT,mean,+1,0.24,0.0065);
		limit2 = GetLimit(pT,buffer,+1,SidebandWidth,0.0065);

	      }

	      // lower sideband
	      if(fisLowerSideband){

		limit2 = GetLimit(pT,mean,-1,0.36,0.009);
		buffer = GetLimit(pT,mean,-1,0.30,0.0065);
		limit1 = GetLimit(pT,buffer,-1,SidebandWidth,0.0065);

	      }

	      // signal
	      if(fisSidebandSignal){

		limit1 = GetLimit(pT,mean,-1,0.30,0.009);
		limit2 = GetLimit(pT,mean,+1,0.24,0.009);

	      }

	      if(mass2 >= limit1 && mass2 <= limit2){

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
	      */

	  // antideuteron mass
	  if(fTrackCutsPart4Mass->isSelected(fTrack)){

	      mass2 = CalculateMassSqTOF(fTrack); 
              fAntideuteronMassSqTOFFullPt->Fill(fTrack->GetPt(),mass2);

	  }

	  // antideuteron sigma
	  if(fTrackCutsPart4Sigma->isSelected(fTrack)){

              fAntideuteronTPCnSigma->Fill(fTrack->GetMomTPC(),fTrack->GetnSigmaTPC((int) (AliPID::kDeuteron)));

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

	fPairCleaner->CleanTrackAndDecay(&DeuteronParticles,&Decays,0);			    // clean deuteron-lambda
	fPairCleaner->CleanTrackAndDecay(&AntideuteronParticles,&AntiDecays,1);		    // clean antideuteron-antilambda

	fPairCleaner->CleanTrackAndDecay(&DeuteronParticles,&AntiDecays,2);		    // clean deuteron-antilambda
	fPairCleaner->CleanTrackAndDecay(&AntideuteronParticles,&Decays,3);		    // clean antideuteron-lambda

	fPairCleaner->CleanTrackAndDecay(&ProtonParticles,&Decays,4);			    // clean proton-lambda
	fPairCleaner->CleanTrackAndDecay(&AntiprotonParticles,&AntiDecays,5);		    // clean antiproton-antilambda

	fPairCleaner->CleanTrackAndDecay(&ProtonParticles,&AntiDecays,6);		    // clean proton-antilambda
	fPairCleaner->CleanTrackAndDecay(&AntiprotonParticles,&Decays,7);		    // clean antiproton-lambda

	fPairCleaner->CleanTrackAndDecay(&DeuteronParticles,&ProtonParticles,8);	    // clean deuteron-proton
	fPairCleaner->CleanTrackAndDecay(&AntideuteronParticles,&AntiprotonParticles,9);    // clean antideuteron-antiproton
      
	fPairCleaner->CleanTrackAndDecay(&DeuteronParticles,&AntideuteronParticles,10);	    // clean deuteron-antideuteron
	fPairCleaner->CleanTrackAndDecay(&ProtonParticles,&AntiprotonParticles,11);	    // clean proton-antiproton

	fPairCleaner->CleanDecay(&Decays,0);						    // clean lambda-lambda
	fPairCleaner->CleanDecay(&AntiDecays,1);					    // clean antilambda-antilambda
 
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



bool AliAnalysisTaskLeuteronAOD::CheckTPCDeuteronPID(AliFemtoDreamTrack *track, double nSigma){

  bool isDeuteron = false;
  double sigma = track->GetnSigmaTPC(AliPID::kDeuteron);

  if(TMath::Abs(sigma) > nSigma) isDeuteron = false;

  return isDeuteron;

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
  TF1 *fit = new TF1("fit","[0]+[1]*pow((1-([2]/(x))),[3])",1.0,4.0);
  fit->FixParameter(0,3.532e+00);
  fit->FixParameter(1,3.899e-14);
  fit->FixParameter(2,-2.217e+04);
  fit->FixParameter(3,2.957e+00);

  Double_t value = fit->Eval(pT);
  fit->Delete();
  return value;

}

//  -----------------------------------------------------------------------------------------------------------------------------------------
Double_t AliAnalysisTaskLeuteronAOD::GetLimit(float pT, double mean, double sign, double offset, double lastpar){

  TF1 *fit = new TF1("fit","[0] + ([1] *([2] + ([3]*(x)) + ([4]*(x)*(x)) + ([5]*(x)*(x)*(x))))",1.0,4.0);
  fit->FixParameter(0,mean);
  fit->FixParameter(1,sign);
  fit->FixParameter(2,offset);
  fit->FixParameter(3,0.003);
  fit->FixParameter(4,0.002);
  fit->FixParameter(5,lastpar);

  Double_t value = fit->Eval(pT);
  fit->Delete();
  return value;

}


//  -----------------------------------------------------------------------------------------------------------------------------------------
bool AliAnalysisTaskLeuteronAOD::CheckDeuteronMassSquarePID(AliFemtoDreamTrack *track,double nSigma){

  bool isDeuteron = false;
  double pT = track->GetPt();
  double massSq = CalculateMassSqTOF(track);
  double parameter[35][2] = {
{4.17615,0.298251},
{4.08097,0.297396},
{3.95197,0.240426},
{3.82068,0.149232},
{3.74164,0.118375},
{3.69264,0.11247},
{3.65695,0.112983},
{3.63148,0.111041},
{3.61447,0.111148},
{3.60399,0.112195},
{3.59677,0.113264},
{3.59217,0.114995},
{3.59126,0.118023},
{3.59045,0.120715},
{3.58861,0.122859},
{3.58585,0.126467},
{3.58287,0.12996},
{3.5818,0.134162},
{3.58132,0.138716},
{3.58191,0.142358},
{3.58124,0.147029},
{3.5818,0.152781},
{3.58151,0.156665},
{3.58148,0.160992},
{3.58288,0.16659},
{3.58375,0.171199},
{3.58526,0.1769},
{3.58636,0.182282},
{3.58728,0.188404},
{3.58966,0.194788},
{3.59245,0.199839},
{3.59187,0.206191},
{3.59449,0.213257},
{3.59563,0.216757}
}; // end of array definition

  int row = 0; // pt < 0.5

  if(pT > 0.6) row = 1;
  if(pT > 0.7) row = 2;
  if(pT > 0.8) row = 3;
  if(pT > 0.9) row = 4;
  if(pT > 1.0) row = 5;
  if(pT > 1.1) row = 6;
  if(pT > 1.2) row = 7;
  if(pT > 1.3) row = 8;
  if(pT > 1.4) row = 9;
  if(pT > 1.5) row = 10;
  if(pT > 1.6) row = 11;
  if(pT > 1.7) row = 12;
  if(pT > 1.8) row = 13;
  if(pT > 1.9) row = 14;
  if(pT > 2.0) row = 15;
  if(pT > 2.1) row = 16;
  if(pT > 2.2) row = 17;
  if(pT > 2.3) row = 18;
  if(pT > 2.4) row = 19;
  if(pT > 2.5) row = 20;
  if(pT > 2.6) row = 21;
  if(pT > 2.7) row = 22;
  if(pT > 2.8) row = 23;
  if(pT > 2.9) row = 24;
  if(pT > 3.0) row = 25;
  if(pT > 3.1) row = 26;
  if(pT > 3.2) row = 27;
  if(pT > 3.3) row = 28;
  if(pT > 3.4) row = 29;
  if(pT > 3.5) row = 30;
  if(pT > 3.6) row = 31;
  if(pT > 3.7) row = 32;
  if(pT > 3.8) row = 33;
  if(pT > 3.9) row = 34;

  double mean = parameter[row][0];
  double sigma = parameter[row][1];

  double left = mean - nSigma*sigma;
  double right = mean + nSigma*sigma;

  if(massSq > left && massSq < right) isDeuteron = true;

//  std::cout << "left: " << left << "\t m²: " << massSq << "\t right: " << right << "\t pT: " << pT << "\t isDeuteron: " << isDeuteron << std::endl;

  return isDeuteron;

}



//  -----------------------------------------------------------------------------------------------------------------------------------------
bool AliAnalysisTaskLeuteronAOD::CheckAntiDeuteronMassSquarePID(AliFemtoDreamTrack *track,double nSigma){

  bool isAntiDeuteron = false;
  double pT = track->GetPt();
  double massSq = CalculateMassSqTOF(track);
  double parameter[35][2] = {
{4.34701,0.293789},
{4.1721,0.270784},
{3.97304,0.210971},
{3.83358,0.16},
{3.74727,0.130793},
{3.69333,0.119528},
{3.65658,0.115693},
{3.63029,0.113118},
{3.6127,0.112189},
{3.60345,0.113028},
{3.59819,0.115984},
{3.59424,0.116708},
{3.59396,0.120272},
{3.59399,0.123617},
{3.59312,0.12665},
{3.59111,0.129108},
{3.58911,0.13153},
{3.58753,0.135724},
{3.58664,0.140496},
{3.58885,0.146367},
{3.58733,0.150429},
{3.58957,0.152747},
{3.58825,0.156932},
{3.59141,0.162311},
{3.59057,0.166413},
{3.59221,0.172248},
{3.59319,0.178037},
{3.59445,0.184472},
{3.5954,0.190762},
{3.59904,0.198087},
{3.6033,0.20339},
{3.60191,0.207736},
{3.60263,0.213331},
{3.60588,0.214263},
}; // end of array definition

  int row = 0; // pt < 0.5

  if(pT > 0.6) row = 1;
  if(pT > 0.7) row = 2;
  if(pT > 0.8) row = 3;
  if(pT > 0.9) row = 4;
  if(pT > 1.0) row = 5;
  if(pT > 1.1) row = 6;
  if(pT > 1.2) row = 7;
  if(pT > 1.3) row = 8;
  if(pT > 1.4) row = 9;
  if(pT > 1.5) row = 10;
  if(pT > 1.6) row = 11;
  if(pT > 1.7) row = 12;
  if(pT > 1.8) row = 13;
  if(pT > 1.9) row = 14;
  if(pT > 2.0) row = 15;
  if(pT > 2.1) row = 16;
  if(pT > 2.2) row = 17;
  if(pT > 2.3) row = 18;
  if(pT > 2.4) row = 19;
  if(pT > 2.5) row = 20;
  if(pT > 2.6) row = 21;
  if(pT > 2.7) row = 22;
  if(pT > 2.8) row = 23;
  if(pT > 2.9) row = 24;
  if(pT > 3.0) row = 25;
  if(pT > 3.1) row = 26;
  if(pT > 3.2) row = 27;
  if(pT > 3.3) row = 28;
  if(pT > 3.4) row = 29;
  if(pT > 3.5) row = 30;
  if(pT > 3.6) row = 31;
  if(pT > 3.7) row = 32;
  if(pT > 3.8) row = 33;
  if(pT > 3.9) row = 34;

  double mean = parameter[row][0];
  double sigma = parameter[row][1];

  double left = mean - nSigma*sigma;
  double right = mean + nSigma*sigma;

  if(massSq > left && massSq < right) isAntiDeuteron = true;

  std::cout << "left: " << left << "\t m²: " << massSq << "\t right: " << right << "\t pT: " << pT << "\t isAntiDeuteron: " << isAntiDeuteron << std::endl;

  return isAntiDeuteron;

}




//  -----------------------------------------------------------------------------------------------------------------------------------------
Double_t AliAnalysisTaskLeuteronAOD::GetAntideuteronMass2Mean_pp(float pT){

// These values were obtained by fitting the mean values of the antideuteron mass2 projections calculated with the AOD dataset (2016,2017,2018)
  TF1 *fit = new TF1("fit","[0]+[1]*pow((1-([2]/(x))),[3])",1.0,4.0);
  fit->FixParameter(0,3.549e+00);
  fit->FixParameter(1,3.146e-14);
  fit->FixParameter(2,-2.134e+04);
  fit->FixParameter(3,2.948e+00);

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
