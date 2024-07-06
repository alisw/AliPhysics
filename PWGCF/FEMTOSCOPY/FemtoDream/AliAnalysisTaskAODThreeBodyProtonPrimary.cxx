/*
 * AliAnalysisTaskAODThreeBodyProtonPrimary.cxx
 *
 *  Created on: May 13, 2019
 *      Authors: Raffaele Del Grande, Marcel Lesch
 *      Based on AliAnalysisTaskThreeBodyFemtoAOD.cxx from Laura Serksnyte
 */

#include "AliAnalysisTaskAODThreeBodyProtonPrimary.h"
#include "AliFemtoDreamHigherPairMath.h"
#include "AliNanoAODTrack.h"
#include "Riostream.h"
#include "AliAODInputHandler.h"
#include "AliAnalysisManager.h"


ClassImp(AliAnalysisTaskAODThreeBodyProtonPrimary)
AliAnalysisTaskAODThreeBodyProtonPrimary::AliAnalysisTaskAODThreeBodyProtonPrimary()
    : AliAnalysisTaskSE(),
      fisLightWeight(false),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fTrack(nullptr),
      fProton(nullptr),
      fProtonList(nullptr),
      fProtonMCList(nullptr),
      fAntiProton(nullptr),
      fAntiProtonList(nullptr),
      fAntiProtonMCList(nullptr),
      fPrimary(nullptr),
      fPrimaryList(nullptr),
      fPrimaryMCList(nullptr),
      fAntiPrimary(nullptr),
      fAntiPrimaryList(nullptr),
      fAntiPrimaryMCList(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsThreeBody(nullptr),
      fSameEvent(nullptr),
      fMixedEvent(nullptr),
      fSameEventMult(nullptr),
      fMixedEventMult(nullptr),
      fSameEventmT(nullptr),
      fMixedEventmT(nullptr),
      fSameEventPhiTheta_SamePair(nullptr),
      fMixedEventPhiTheta_SamePair(nullptr),
      fSameEventPhiTheta_DifferentPair(nullptr),
      fMixedEventPhiTheta_DifferentPair(nullptr),
      fInvMassList(nullptr),
      fP1Histos(nullptr),
      fKHistos(nullptr),
      fInvMassVSmTList(nullptr),

      fRunThreeBody(true),
      fRunPlotInvMass(false),
      fRunPlotPhiTheta(true),
      fRunPlotMult(true),
      fRunPairMultThreeBody(false),
      fRunOfficialTwoBody(false),
      fRunPlotInvMassVSmT(false),

      fClosePairRejectionForAll(false),
      fturnoffClosePairRejectionCompletely(false),
      fQ3LimitForDeltaPhiDeltaEta(1.),
      fDeltaPhiMaxPP(0.017),
      fDeltaEtaMaxPP(0.017),
      fDeltaPhiMaxPPrim(0.04),
      fDeltaEtaMaxPPrim(0.012),
      fDeltaPhiMaxPAPrim(0.04),
      fDeltaEtaMaxPAPrim(0.012),
      fQ3cutValue(1.),
      fQ3MinValue(0.),
      fDoOnlyThreeBody(true),
      fStandardMixing(true),
      fRunmTPlots(false), 
      fIsMC(false),
      fRemoveMCResonances(false),
      fRunProjector(false), 
      fGetMomentumResolution(false),

      fSameEventTripletArray(nullptr),
      fSameEventTripletMultArray(nullptr),
      fSameEventTripletMultArray12(nullptr),
      fSameEventTripletMultArray23(nullptr),
      fSameEventTripletMultArray31(nullptr),
      fSameEventTripletmTArray12(nullptr),
      fSameEventTripletmTArray23(nullptr),
      fSameEventTripletmTArray31(nullptr),
      fSameEventTripletPhiThetaArray_SamePair(nullptr), 
      fSameEventTripletPhiThetaArray_DifferentPair(nullptr),

      fSameEventPairArray_TwoBody(nullptr),
      fSameEventMultArray_TwoBody(nullptr),
      fSameEventDphiArray_TwoBody(nullptr),
      fSameEventPairMultArray_TwoBody(nullptr),
      fSameEventPairPhiThetaArray_TwoBody(nullptr),

      fPairTranverseMass_TwoBody(nullptr), 
      fPairTranverseMassVSkstar_TwoBody(nullptr),
      fMixingConfig(nullptr),

      fPartContainer(0),
      fPartContainerPPP(0),
      fPartContainerPPPrim(0),
      fPartContainerPPAPrim(0),
      fPartContainerPP(0),
      fPartContainerPPrim(0),
      fPartContainerPAPrim(0),
      fVectPartContainers(0),

      fMixedEventTripletArray(nullptr),
      fMixedEventTripletMultArray(nullptr),
      fMixedEventTripletMultArray12(nullptr),
      fMixedEventTripletMultArray23(nullptr),
      fMixedEventTripletMultArray31(nullptr),
      fMixedEventTripletmTArray12(nullptr),
      fMixedEventTripletmTArray23(nullptr),
      fMixedEventTripletmTArray31(nullptr),
      fMixedEventTripletPhiThetaArray_SamePair(nullptr), 
      fMixedEventTripletPhiThetaArray_DifferentPair(nullptr),
      
      fMixedEventPairArray_TwoBody(nullptr),
      fMixedEventPairMultArray_TwoBody(nullptr),
      fMixedEventMultArray_TwoBody(nullptr),
      fMixedEventDphiArray_TwoBody(nullptr),
      fMixedEventPairPhiThetaArray_TwoBody(nullptr),

      fInvMass12(nullptr),
      fInvMass23(nullptr),
      fInvMass31(nullptr),
      fInvMassVSmT12(nullptr),
      fInvMassVSmT23(nullptr),
      fInvMassVSmT31(nullptr),
      fInvMassVSmT12_MixedEvent(nullptr),
      fInvMassVSmT23_MixedEvent(nullptr),
      fInvMassVSmT31_MixedEvent(nullptr),
      fProjectorData(nullptr),
      fSameEventTripletResolutionList(nullptr), 
      fSameEventTripletResolution(nullptr), 
      fSameEventTripletResolutionAll(nullptr), 
      fMixedEventTripletResolutionList(nullptr),  
      fMixedEventTripletResolution(nullptr), 
      fMixedEventTripletResolutionAll(nullptr), 

      fResultsQA(nullptr),
      fSample(nullptr),
      fResultsSample(nullptr),
      fResultsSampleQA(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
}

AliAnalysisTaskAODThreeBodyProtonPrimary::AliAnalysisTaskAODThreeBodyProtonPrimary(const char* name, bool isMC)
    : AliAnalysisTaskSE(name),
      fisLightWeight(false),
      fEvent(nullptr),
      fEventCuts(nullptr),
      fEvtList(nullptr),
      fTrack(nullptr),
      fProton(nullptr),
      fProtonList(nullptr),
      fProtonMCList(nullptr),
      fAntiProton(nullptr),
      fAntiProtonList(nullptr),
      fAntiProtonMCList(nullptr),
      fPrimary(nullptr),
      fPrimaryList(nullptr),
      fPrimaryMCList(nullptr),
      fAntiPrimary(nullptr),
      fAntiPrimaryList(nullptr),
      fAntiPrimaryMCList(nullptr),
      fConfig(nullptr),
      fPairCleaner(nullptr),
      fPartColl(nullptr),
      fResults(nullptr),
      fResultsThreeBody(nullptr),
      fSameEvent(nullptr),
      fMixedEvent(nullptr),
      fSameEventMult(nullptr),
      fMixedEventMult(nullptr),
      fSameEventmT(nullptr),
      fMixedEventmT(nullptr),
      fSameEventPhiTheta_SamePair(nullptr),
      fMixedEventPhiTheta_SamePair(nullptr),
      fSameEventPhiTheta_DifferentPair(nullptr),
      fMixedEventPhiTheta_DifferentPair(nullptr),
      fInvMassList(nullptr),
      fP1Histos(nullptr),
      fKHistos(nullptr),
      fInvMassVSmTList(nullptr),

      fRunThreeBody(true),
      fRunPlotInvMass(false),
      fRunPlotPhiTheta(true),
      fRunPlotMult(true),
      fRunPairMultThreeBody(false),
      fRunOfficialTwoBody(false),
      fRunPlotInvMassVSmT(false),

      fClosePairRejectionForAll(false),
      fturnoffClosePairRejectionCompletely(false),
      fQ3LimitForDeltaPhiDeltaEta(1.),
      fDeltaPhiMaxPP(0.017),
      fDeltaEtaMaxPP(0.017),
      fDeltaPhiMaxPPrim(0.04),
      fDeltaEtaMaxPPrim(0.012),
      fDeltaPhiMaxPAPrim(0.04),
      fDeltaEtaMaxPAPrim(0.012),
      fQ3cutValue(1.),
      fQ3MinValue(0.),
      fDoOnlyThreeBody(true),
      fStandardMixing(true),
      fRunmTPlots(false), 
      fIsMC(false),
      fRemoveMCResonances(false),
      fRunProjector(false), 
      fGetMomentumResolution(false),

      fSameEventTripletArray(nullptr),
      fSameEventTripletMultArray(nullptr),
      fSameEventTripletMultArray12(nullptr),
      fSameEventTripletMultArray23(nullptr),
      fSameEventTripletMultArray31(nullptr),
      fSameEventTripletmTArray12(nullptr),
      fSameEventTripletmTArray23(nullptr),
      fSameEventTripletmTArray31(nullptr),
      fSameEventTripletPhiThetaArray_SamePair(nullptr), 
      fSameEventTripletPhiThetaArray_DifferentPair(nullptr),

      fSameEventPairArray_TwoBody(nullptr),
      fSameEventMultArray_TwoBody(nullptr),
      fSameEventDphiArray_TwoBody(nullptr),
      fSameEventPairMultArray_TwoBody(nullptr),
      fSameEventPairPhiThetaArray_TwoBody(nullptr),

      fPairTranverseMass_TwoBody(nullptr),
      fPairTranverseMassVSkstar_TwoBody(nullptr),
      fMixingConfig(nullptr),

      fPartContainer(0),
      fPartContainerPPP(0),
      fPartContainerPPPrim(0),
      fPartContainerPPAPrim(0),
      fPartContainerPP(0),
      fPartContainerPPrim(0),
      fPartContainerPAPrim(0),
      fVectPartContainers(0),

      fMixedEventTripletArray(nullptr),
      fMixedEventTripletMultArray(nullptr),
      fMixedEventTripletMultArray12(nullptr),
      fMixedEventTripletMultArray23(nullptr),
      fMixedEventTripletMultArray31(nullptr),
      fMixedEventTripletmTArray12(nullptr),
      fMixedEventTripletmTArray23(nullptr),
      fMixedEventTripletmTArray31(nullptr),
      fMixedEventTripletPhiThetaArray_SamePair(nullptr), 
      fMixedEventTripletPhiThetaArray_DifferentPair(nullptr),
      
      fMixedEventPairArray_TwoBody(nullptr),
      fMixedEventPairMultArray_TwoBody(nullptr),
      fMixedEventMultArray_TwoBody(nullptr),
      fMixedEventDphiArray_TwoBody(nullptr),
      fMixedEventPairPhiThetaArray_TwoBody(nullptr),

      fInvMass12(nullptr),
      fInvMass23(nullptr),
      fInvMass31(nullptr),
      fInvMassVSmT12(nullptr),
      fInvMassVSmT23(nullptr),
      fInvMassVSmT31(nullptr),
      fInvMassVSmT12_MixedEvent(nullptr),
      fInvMassVSmT23_MixedEvent(nullptr),
      fInvMassVSmT31_MixedEvent(nullptr),
      fProjectorData(nullptr),
      fSameEventTripletResolutionList(nullptr), 
      fSameEventTripletResolution(nullptr), 
      fSameEventTripletResolutionAll(nullptr), 
      fMixedEventTripletResolutionList(nullptr),  
      fMixedEventTripletResolution(nullptr), 
      fMixedEventTripletResolutionAll(nullptr), 

      fResultsQA(nullptr),
      fSample(nullptr),
      fResultsSample(nullptr),
      fResultsSampleQA(nullptr),
      fTrackBufferSize(2000),
      fGTI(nullptr) {
        DefineOutput(1, TList::Class());  //Output for the Event Cuts
        DefineOutput(2, TList::Class());  //Output for the Proton Cuts
        DefineOutput(3, TList::Class());  //Output for the AntiProton Cuts
        DefineOutput(4, TList::Class());  //Output for the Primary Cuts
        DefineOutput(5, TList::Class());  //Output for the AntiPrimary Cuts
        DefineOutput(6, TList::Class());  //Output for the Results
        DefineOutput(7, TList::Class());  //Output for the Results QA
        DefineOutput(8, TList::Class());  //Output for the Results
        DefineOutput(9, TList::Class());  //Output for the Results QA
        DefineOutput(10, TList::Class());  //Output for the Results Three body
        if (isMC) {
          DefineOutput(11, TList::Class());  //Output for the Track MC
          DefineOutput(12, TList::Class());  //Output for the Anti Track MC
          DefineOutput(13, TList::Class());  //Output for the Primary MC
          DefineOutput(14, TList::Class());  //Output for the Anti Primary MC
        }
      }


//==================================================================================================================================================

AliAnalysisTaskAODThreeBodyProtonPrimary::~AliAnalysisTaskAODThreeBodyProtonPrimary() {
  if (fEvent) {
    delete fEvent;
  }
  if (fEventCuts) {
    delete fEventCuts;
  }
  if (fTrack) {
    delete fTrack;
  }
  if (fProton) {
    delete fProton;
  }
  if (fAntiProton) {
    delete fAntiProton;
  }
  if (fPrimary) {
    delete fPrimary;
  }
  if (fAntiPrimary) {
    delete fAntiPrimary;
  }
  if (fPairCleaner) {
    delete fPairCleaner;
  }
  if (fPartColl) {
    delete fPartColl;
  }
  if (fSample) {
    delete fSample;
  }
} //AliAnalysisTaskAODThreeBodyProtonPrimary::~AliAnalysisTaskAODThreeBodyProtonPrimary()

//==================================================================================================================================================

void AliAnalysisTaskAODThreeBodyProtonPrimary::UserCreateOutputObjects() {
  fGTI = new AliAODTrack*[fTrackBufferSize];

  if (!fEventCuts) {
    AliError("No Event cuts \n");
  } else {
    fEventCuts->InitQA();
  }
  if (!fProton) {
    AliError("No Proton cuts \n");
  } else {
    fProton->Init();
  }
  if (!fAntiProton) {
    AliError("No AntiProton cuts \n");
  } else {
    fAntiProton->Init();
  }
  if (!fPrimary) {
    AliError("No Primary cuts \n");
  } else {
    fPrimary->Init();
  }
  if (!fAntiPrimary) {
    AliError("No AntiPrimary cuts \n");
  } else {
    fAntiPrimary->Init();
  }
  if (!fConfig) {
    AliError("No Correlation Config \n");
  } else {
    fPartColl = new AliFemtoDreamPartCollection(fConfig,
                                                fConfig->GetMinimalBookingME());
    fPairCleaner = new AliFemtoDreamPairCleaner(2, 2,
                                                fConfig->GetMinimalBookingME());
    if (fConfig->GetUsePhiSpinning()) {
      fSample = new AliFemtoDreamControlSample(fConfig);
    }
  }

  //fEvent = new AliFemtoDreamEvent(true, !fisLightWeight,
                                //  GetCollisionCandidates(), false);
  //fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());
  //fEvent->SetCalcSpherocity(fEventCuts->GetDoSpherocityCuts());

  fEvent = new AliFemtoDreamEvent(true, !fisLightWeight,GetCollisionCandidates(), true);
  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());

  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(
      fProton->GetIsMonteCarlo() || fAntiProton->GetIsMonteCarlo()||fPrimary->GetIsMonteCarlo() || fAntiPrimary->GetIsMonteCarlo());


  if (!fEventCuts->GetMinimalBooking()) {
    fEvtList = fEventCuts->GetHistList();
  } else {
    fEvtList = new TList();
    fEvtList->SetName("EventCuts");
    fEvtList->SetOwner();
  }

  fProtonList = fProton->GetQAHists();
  fAntiProtonList = fAntiProton->GetQAHists();
  fPrimaryList = fPrimary->GetQAHists();
  fAntiPrimaryList = fAntiPrimary->GetQAHists();

  fResultsQA = new TList();
  fResultsQA->SetOwner();
  fResultsQA->SetName("ResultsQA");

  if (fConfig->GetUseEventMixing()) {
    fResults = fPartColl->GetHistList();
    if (!fConfig->GetMinimalBookingME()) {
      fResultsQA->Add(fPartColl->GetQAList());
      fResultsQA->Add(fPairCleaner->GetHistList());
    }
  } else {
    fResults = new TList();
    fResults->SetOwner();
    fResults->SetName("Results");
  }


  //=================================================================================================================================

  if (fRunThreeBody){

    fResultsThreeBody = new TList();
    fResultsThreeBody->SetOwner();
    fResultsThreeBody->SetName("ResultsThreeBody");

    fMixingConfig = new TH1F("fMixingConfig", "fMixingConfig", 3, 0.,3.); 
    fMixingConfig->GetXaxis()->SetBinLabel(1,"Standard Mixing"); 
    fMixingConfig->GetXaxis()->SetBinLabel(2,"Triplet Mixing");
    fMixingConfig->GetXaxis()->SetBinLabel(3,"Mixing Depth");

    if(fStandardMixing){
      fMixingConfig->SetBinContent(1,1);
    } else {
      fMixingConfig->SetBinContent(2,1);
    }
    fMixingConfig->SetBinContent(3,fConfig->GetMixingDepth());
    fResultsThreeBody->Add(fMixingConfig); 
    //...............................................................................................
    // Same event distributions 1D
    fSameEvent = new TList();
    fSameEvent->SetOwner();
    fSameEvent->SetName("SameEvent");

    //Q3 Distribution for (xxy) and (xx)(y) events
    fSameEventTripletArray = new TH1F*[16];
    
    if(fDoOnlyThreeBody){
	    for (int i = 0; i < 16; ++i) {
       TString histTitle; 
       if(i<6){
         histTitle = "sameEventDistribution"+fParticleNames[fTripletCombinations[i][0]]+fParticleNames[fTripletCombinations[i][1]]+fParticleNames[fTripletCombinations[i][2]];
       } else {
         histTitle = "sameEventDistribution"+fParticleNames[fSameMixedCominations[i-6][0]]+fParticleNames[fSameMixedCominations[i-6][1]]+"Same"+fParticleNames[fSameMixedCominations[i-6][2]]+"Mixed";
       }
	     fSameEventTripletArray[i] =  new TH1F(histTitle,histTitle, 8000, 0, 8);
       fSameEvent->Add(fSameEventTripletArray[i]);
	   }
    }

    //k* Distribution for (xy)
    fSameEventPairArray_TwoBody = new TH1F*[6];

    //mT Distribution for (xy)
    fPairTranverseMass_TwoBody = new TH1F*[6]; 
    
    if(!fDoOnlyThreeBody){
	   for (int i = 0; i < 6; ++i) {
       TString histTitle = "sameEventDistribution"+fParticleNames[fPairCombinations[i][0]]+fParticleNames[fPairCombinations[i][1]];
	     fSameEventPairArray_TwoBody[i] =  new TH1F(histTitle,histTitle, 8000, 0, 8);
	     fSameEvent->Add(fSameEventPairArray_TwoBody[i]);

      if(fRunOfficialTwoBody){
        TString histTitle2 = "pairTranverseMassDistribution"+fParticleNames[fPairCombinations[i][0]]+fParticleNames[fPairCombinations[i][1]];
        fPairTranverseMass_TwoBody[i] =  new TH1F(histTitle2,histTitle2, 8000, 0, 8);  
        fSameEvent->Add(fPairTranverseMass_TwoBody[i]);  
      }
	   }
    }

    //...............................................................................................
    // Same event multiplicity dist
    fSameEventMult = new TList();
    fSameEventMult->SetOwner();
    fSameEventMult->SetName("SameEventMult");

    fSameEventmT = new TList();
    fSameEventmT->SetOwner();
    fSameEventmT->SetName("fSameEventmT");

    //Q3 Distribution vs mult bin for (xxy) and (xx)(y) events
    fSameEventTripletMultArray = new TH2F*[16];

    if(fDoOnlyThreeBody){
	    for (int i = 0; i < 16; ++i) {
        TString histTitle; 
        if(i<6){
          histTitle = "sameEventDistributionMult"+fParticleNames[fTripletCombinations[i][0]]+fParticleNames[fTripletCombinations[i][1]]+fParticleNames[fTripletCombinations[i][2]];
        } else {
          histTitle = "sameEventDistributionMult"+fParticleNames[fSameMixedCominations[i-6][0]]+fParticleNames[fSameMixedCominations[i-6][1]]+"Same"+fParticleNames[fSameMixedCominations[i-6][2]]+"Mixed";
        }
	      fSameEventTripletMultArray[i] =  new TH2F(histTitle,histTitle, 8000, 0, 8.,26,1,27);
        if(fRunPlotMult){fSameEventMult->Add(fSameEventTripletMultArray[i]);}
	     }
    }

    //k* of pair (ij) vs mult bin 
    fSameEventTripletMultArray12 = new TH2F*[16]; //here we also have the case of 2same1mixed: pair (ij), i,j = 1,2 same event, 3 mixed event
    fSameEventTripletMultArray23 = new TH2F*[6];
    fSameEventTripletMultArray31 = new TH2F*[6];
    
    //Q3 vs mT of a pair. In case of 2same1mixed: pair (ij), i,j = 1,2 same event, 3 mixed event
    fSameEventTripletmTArray12 = new TH2F*[16];
    fSameEventTripletmTArray23 = new TH2F*[16];
    fSameEventTripletmTArray31 = new TH2F*[16];
    
    if(fDoOnlyThreeBody){
      for (int i = 0; i < 16; ++i) { 
       
        TString histTitle; 
        if(i<6){
          histTitle = "sameEventDistributionMult"+fParticleNames[fTripletCombinations[i][0]]+fParticleNames[fTripletCombinations[i][1]]+fParticleNames[fTripletCombinations[i][2]];

           fSameEventTripletMultArray23[i] = new TH2F(histTitle+"23",histTitle+"23", 8000, 0, 8,26,1,27);
          if(fRunPairMultThreeBody){fSameEventMult->Add(fSameEventTripletMultArray23[i]);}
           fSameEventTripletMultArray31[i] = new TH2F(histTitle+"31",histTitle+"31", 8000, 0, 8,26,1,27);
          if(fRunPairMultThreeBody){fSameEventMult->Add(fSameEventTripletMultArray31[i]);}

        } else {
          histTitle = "sameEventDistributionMult"+fParticleNames[fSameMixedCominations[i-6][0]]+fParticleNames[fSameMixedCominations[i-6][1]]+"Same"+fParticleNames[fSameMixedCominations[i-6][2]]+"Mixed";
        }
        
        fSameEventTripletMultArray12[i] = new TH2F(histTitle+"12",histTitle+"12", 8000, 0, 8,26,1,27);
          if(fRunPairMultThreeBody){fSameEventMult->Add(fSameEventTripletMultArray12[i]);}  
      }

      for (int i = 0; i < 16; ++i) {
        TString histTitle; 
        if(i<6){
          histTitle  = "sameEventDistributionTransMass"+fParticleNames[fTripletCombinations[i][0]]+fParticleNames[fTripletCombinations[i][1]]+fParticleNames[fTripletCombinations[i][2]];
        } else {
          histTitle = "sameEventDistributionTransMass"+fParticleNames[fSameMixedCominations[i-6][0]]+fParticleNames[fSameMixedCominations[i-6][1]]+"Same"+fParticleNames[fSameMixedCominations[i-6][2]]+"Mixed";
        }

        fSameEventTripletmTArray12[i] = new TH2F(histTitle+"12",histTitle+"12", 4000, 0, 8, 100, 0., 5.); 
          if(fRunmTPlots){fSameEventmT->Add(fSameEventTripletmTArray12[i]);}
        fSameEventTripletmTArray23[i] = new TH2F(histTitle+"23",histTitle+"23", 4000, 0, 8, 100, 0., 5.);
          if(fRunmTPlots){fSameEventmT->Add(fSameEventTripletmTArray23[i]);}
        fSameEventTripletmTArray31[i] = new TH2F(histTitle+"31",histTitle+"31", 4000, 0, 8, 100, 0., 5.);
          if(fRunmTPlots){fSameEventmT->Add(fSameEventTripletmTArray31[i]);}
       }
    }


    //Momentum Resolution Triplet SE
    fSameEventTripletResolutionList = new TList();
    fSameEventTripletResolutionList->SetOwner();
    fSameEventTripletResolutionList->SetName("SameEventTripletResolution");

    fSameEventTripletResolution = new TH2F*[6];
    fSameEventTripletResolutionAll = new TH2F*[6];
    TString histTitleSameResolution; 

    if(fDoOnlyThreeBody && fIsMC && fGetMomentumResolution){
      for (int i = 0; i < 6; ++i) {
        histTitleSameResolution =  "sameEventResolution"+fParticleNames[fTripletCombinations[i][0]]+fParticleNames[fTripletCombinations[i][1]]+fParticleNames[fTripletCombinations[i][2]];
        fSameEventTripletResolution[i] = new TH2F(histTitleSameResolution,histTitleSameResolution, 1500, 0, 1.5, 1500, 0, 1.5);

        histTitleSameResolution =  "sameEventResolutionAll"+fParticleNames[fTripletCombinations[i][0]]+fParticleNames[fTripletCombinations[i][1]]+fParticleNames[fTripletCombinations[i][2]];
        fSameEventTripletResolutionAll[i] = new TH2F(histTitleSameResolution,histTitleSameResolution, 1500, 0, 1.5, 1500, 0, 1.5);
      
        fSameEventTripletResolutionList->Add(fSameEventTripletResolution[i]);
        fSameEventTripletResolutionList->Add(fSameEventTripletResolutionAll[i]);
       }
       fResultsThreeBody->Add(fSameEventTripletResolutionList);
    }

    //Two-Body .................................................
    fSameEventPairMultArray_TwoBody = new TH2F*[6];
    fSameEventMultArray_TwoBody = new TH2F*[6];
    fSameEventDphiArray_TwoBody = new TH2F*[6];
    fPairTranverseMassVSkstar_TwoBody = new TH2F*[6]; 
    
    if(!fDoOnlyThreeBody){
	   for (int i = 0; i < 6; ++i) {
        TString histTitle = "sameEventDistributionMult"+fParticleNames[fPairCombinations[i][0]]+fParticleNames[fPairCombinations[i][1]];
	      fSameEventPairMultArray_TwoBody[i] =  new TH2F(histTitle,histTitle,8000, 0, 8,26,1,27);

        TString histTitle2 = "sameEventDistributionMult_TwoBody"+fParticleNames[fPairCombinations[i][0]]+fParticleNames[fPairCombinations[i][1]];
        fSameEventMultArray_TwoBody[i] =  new TH2F(histTitle2,histTitle2,8000, 0, 8,26,1,27);

        TString histTitle3 = "sameEventDistributionDphi"+fParticleNames[fPairCombinations[i][0]]+fParticleNames[fPairCombinations[i][1]];
        fSameEventDphiArray_TwoBody[i] =  new TH2F(histTitle3,histTitle3,8000, 0, 8,200,-10,10);

        if(fRunPlotMult) fSameEventMult->Add(fSameEventPairMultArray_TwoBody[i]);

	      if(fRunOfficialTwoBody){
             TString histTitle4 = "PairTranverseMassVSkstarDistribution"+fParticleNames[fPairCombinations[i][0]]+fParticleNames[fPairCombinations[i][1]];
             fPairTranverseMassVSkstar_TwoBody[i] =  new TH2F(histTitle4,histTitle4, 800, 0, 8, 800, 0, 8); 
             fSameEvent->Add(fPairTranverseMassVSkstar_TwoBody[i]); 
        }
      }
    }

    //...............................................................................................
    // Mixed event 1D
    fMixedEvent = new TList();
    fMixedEvent->SetOwner();
    fMixedEvent->SetName("MixedEvent");

    //Q3 Distribution for (x)(x)(y)
    fMixedEventTripletArray = new TH1F*[6];
    
    if(fDoOnlyThreeBody){
	    for (int i = 0; i < 6; ++i) {
        TString histTitle = "mixedEventDistribution"+fParticleNames[fTripletCombinations[i][0]]+fParticleNames[fTripletCombinations[i][1]]+fParticleNames[fTripletCombinations[i][2]];
       
	      fMixedEventTripletArray[i] = new TH1F(histTitle,histTitle, 8000, 0, 8);
	      fMixedEvent->Add(fMixedEventTripletArray[i]);
	     }
    }

   
    fMixedEventPairArray_TwoBody = new TH1F*[6];
    
    if(!fDoOnlyThreeBody){
	     for (int i = 0; i < 6; ++i) {
        TString histTitle = "mixedEventDistribution"+fParticleNames[fPairCombinations[i][0]]+fParticleNames[fPairCombinations[i][1]];
	      fMixedEventPairArray_TwoBody[i] =  new TH1F(histTitle,histTitle,  8000, 0, 8);
	      fMixedEvent->Add(fMixedEventPairArray_TwoBody[i]);
	     }
    }

    //...............................................................................................
    // Mixed event multiplicity dist 
    fMixedEventMult = new TList();
    fMixedEventMult->SetOwner();
    fMixedEventMult->SetName("MixedEventMult");

    fMixedEventmT = new TList();
    fMixedEventmT->SetOwner();
    fMixedEventmT->SetName("fMixedEventmT");

    //Q3 Distribution vs mult bin for (x)(x)(y)
    fMixedEventTripletMultArray = new TH2F*[6];

    //k* of pair (ij) vs mult bin (mixed event)
    fMixedEventTripletMultArray12 = new TH2F*[6];
    fMixedEventTripletMultArray23 = new TH2F*[6];
    fMixedEventTripletMultArray31 = new TH2F*[6];
   
    //Mixed Event Q3 vs mT Pair
    fMixedEventTripletmTArray12 = new TH2F*[6];
    fMixedEventTripletmTArray23 = new TH2F*[6];
    fMixedEventTripletmTArray31 = new TH2F*[6];
    
    if(fDoOnlyThreeBody){
	    for (int i = 0; i < 6; ++i) {
        TString histTitle = "mixedEventDistributionMult"+fParticleNames[fTripletCombinations[i][0]]+fParticleNames[fTripletCombinations[i][1]]+fParticleNames[fTripletCombinations[i][2]];

	      fMixedEventTripletMultArray[i] = new TH2F(histTitle,histTitle, 8000, 0, 8,26,1,27);
          if(fRunPlotMult){fMixedEventMult->Add(fMixedEventTripletMultArray[i]);}
        fMixedEventTripletMultArray12[i] = new TH2F(histTitle+"12",histTitle+"12", 8000, 0, 8,26,1,27);
          if(fRunPairMultThreeBody){fMixedEventMult->Add(fMixedEventTripletMultArray12[i]);}
        fMixedEventTripletMultArray23[i] = new TH2F(histTitle+"23",histTitle+"23", 8000, 0, 8,26,1,27);
          if(fRunPairMultThreeBody){fMixedEventMult->Add(fMixedEventTripletMultArray23[i]);}
        fMixedEventTripletMultArray31[i] = new TH2F(histTitle+"31",histTitle+"31", 8000, 0, 8,26,1,27);
          if(fRunPairMultThreeBody){fMixedEventMult->Add(fMixedEventTripletMultArray31[i]);}

        TString histTitle2 = "mixedEventDistribution"+fParticleNames[fTripletCombinations[i][0]]+fParticleNames[fTripletCombinations[i][1]]+fParticleNames[fTripletCombinations[i][2]];
        fMixedEventTripletmTArray12[i] = new TH2F(histTitle2+"12",histTitle2+"12", 4000, 0, 8, 100, 0., 5.);
          if(fRunmTPlots && i<6){fMixedEventmT->Add(fMixedEventTripletmTArray12[i]);} 
        fMixedEventTripletmTArray23[i] = new TH2F(histTitle2+"23",histTitle2+"23", 4000, 0, 8, 100, 0., 5.);
          if(fRunmTPlots && i<6){fMixedEventmT->Add(fMixedEventTripletmTArray23[i]);}
        fMixedEventTripletmTArray31[i] = new TH2F(histTitle2+"31",histTitle2+"31", 4000, 0, 8, 100, 0., 5.);
          if(fRunmTPlots && i<6){fMixedEventmT->Add(fMixedEventTripletmTArray31[i]);}
	     }

       fProjectorData = new TH2F*[6];
       if(fRunProjector && fDoOnlyThreeBody){
       for (int i = 0; i < 6; ++i) {
        if(i == 2){
          TString histTitle = "Projector"+fParticleNames[fTripletCombinations[i][0]]+fParticleNames[fTripletCombinations[i][1]]+fParticleNames[fTripletCombinations[i][2]];
	        fProjectorData[i] = new TH2F(histTitle,histTitle, 1000,0.,1., 1000, 0, 1.);
	        fMixedEventMult->Add(fProjectorData[i]);
        } else {
          fProjectorData[i] = nullptr; 
        }

	     }
    }

    }

    //Momentum Resolution Triplet ME
    fMixedEventTripletResolutionList = new TList();
    fMixedEventTripletResolutionList->SetOwner();
    fMixedEventTripletResolutionList->SetName("MixedEventTripletResolution");

    fMixedEventTripletResolution = new TH2F*[6];
    fMixedEventTripletResolutionAll = new TH2F*[6];
    TString histTitleMixedResolution; 

    if(fDoOnlyThreeBody && fIsMC && fGetMomentumResolution){
      for (int i = 0; i < 6; ++i) {
        histTitleMixedResolution =  "mixedEventResolution"+fParticleNames[fTripletCombinations[i][0]]+fParticleNames[fTripletCombinations[i][1]]+fParticleNames[fTripletCombinations[i][2]];
        fMixedEventTripletResolution[i] = new TH2F(histTitleMixedResolution,histTitleMixedResolution, 1500, 0, 1.5, 1500, 0, 1.5);

        histTitleMixedResolution =  "mixedEventResolutionAll"+fParticleNames[fTripletCombinations[i][0]]+fParticleNames[fTripletCombinations[i][1]]+fParticleNames[fTripletCombinations[i][2]];
        fMixedEventTripletResolutionAll[i] = new TH2F(histTitleMixedResolution,histTitleMixedResolution, 1500, 0, 1.5, 1500, 0, 1.5);
      
        fMixedEventTripletResolutionList->Add(fMixedEventTripletResolution[i]);
        fMixedEventTripletResolutionList->Add(fMixedEventTripletResolutionAll[i]);
       }
       fResultsThreeBody->Add(fMixedEventTripletResolutionList);
    }

    //Two-Body ......................................
    fMixedEventPairMultArray_TwoBody = new TH2F*[6];
    fMixedEventMultArray_TwoBody = new TH2F*[6];
    fMixedEventDphiArray_TwoBody = new TH2F*[6];
  
    if(!fDoOnlyThreeBody){
	     for (int i = 0; i < 6; ++i) {
        TString histTitle = "mixedEventDistributionMult"+fParticleNames[fPairCombinations[i][0]]+fParticleNames[fPairCombinations[i][1]];
	      fMixedEventPairMultArray_TwoBody[i] =  new TH2F(histTitle,histTitle, 8000, 0, 8,26,1,27);
        fMixedEventMultArray_TwoBody[i] =  new TH2F(histTitle,histTitle, 8000, 0, 8,26,1,27);

        TString histTitle2 = "mixedEventDistributionDphi"+fParticleNames[fPairCombinations[i][0]]+fParticleNames[fPairCombinations[i][1]];
        fMixedEventDphiArray_TwoBody[i] =  new TH2F(histTitle2,histTitle2, 8000, 0, 8,200,-10,10);
             if(fRunPlotMult){ fMixedEventMult->Add(fMixedEventPairMultArray_TwoBody[i]); }
	     }
    }

    if(fRunPlotMult){
      fResultsThreeBody->Add(fSameEventMult);
      fResultsThreeBody->Add(fMixedEventMult);
    }else{
      fResultsThreeBody->Add(fSameEvent);
      fResultsThreeBody->Add(fMixedEvent);
    }
    if(fRunmTPlots){
      fResultsThreeBody->Add(fSameEventmT);
      fResultsThreeBody->Add(fMixedEventmT);
    }

    //...............................................................................................
    // Same event phi theta distribution Same Pair
    fSameEventPhiTheta_SamePair = new TList();
    fSameEventPhiTheta_SamePair->SetOwner();
    fSameEventPhiTheta_SamePair->SetName("SameEventPhiTheta");

    fSameEventTripletPhiThetaArray_SamePair = new TH2F*[32];

    if(fDoOnlyThreeBody){
	    for(int i=0;i<16;i++){
        TString histTitle; 
        if(i<6){
          histTitle = "sameEventPhiEtaMult"+fParticleNames[fTripletCombinations[i][0]]+fParticleNames[fTripletCombinations[i][1]]+fParticleNames[fTripletCombinations[i][2]];
        } else {
          histTitle = "sameEventPhiEtaMult"+fParticleNames[fSameMixedCominations[i-6][0]]+fParticleNames[fSameMixedCominations[i-6][1]]+"Same"+fParticleNames[fSameMixedCominations[i-6][2]]+"Mixed";
        }

	      fSameEventTripletPhiThetaArray_SamePair[i] = new TH2F(histTitle+"Before",histTitle+"Before", 500, -0.15,0.15,500,-0.15,0.15);
	      fSameEventTripletPhiThetaArray_SamePair[16+i] = new TH2F(histTitle+"After",histTitle+"After", 500, -0.15,0.15,500,-0.15,0.15);

        if(fRunPlotPhiTheta){
	        fSameEventPhiTheta_SamePair->Add(fSameEventTripletPhiThetaArray_SamePair[i]);
	        fSameEventPhiTheta_SamePair->Add(fSameEventTripletPhiThetaArray_SamePair[16+i]);
        }
	    }
    }

    fSameEventPairPhiThetaArray_TwoBody = new TH2F*[12];
    
    //Two-Body ........................
    if(!fDoOnlyThreeBody){
	    for(int i=0;i<6;i++){
        TString histTitle = "sameEventPhiEtaMult"+fParticleNames[fPairCombinations[i][0]]+fParticleNames[fPairCombinations[i][1]];
	      fSameEventPairPhiThetaArray_TwoBody[i] = new TH2F(histTitle+"Before",histTitle+"Before", 500, -0.15,0.15,500,-0.15,0.15);
	      fSameEventPairPhiThetaArray_TwoBody[6+i] = new TH2F(histTitle+"After",histTitle+"After", 500, -0.15,0.15,500,-0.15,0.15);

         if(fRunPlotPhiTheta){
      	    fSameEventPhiTheta_SamePair->Add(fSameEventPairPhiThetaArray_TwoBody[i]);
      	    fSameEventPhiTheta_SamePair->Add(fSameEventPairPhiThetaArray_TwoBody[6+i]);
        }
	    }
    }

    //...............................................................................................
    // Same event phi theta distribution Different Pair
    fSameEventPhiTheta_DifferentPair = new TList();
    fSameEventPhiTheta_DifferentPair->SetOwner();
    fSameEventPhiTheta_DifferentPair->SetName("SameEventPhiTheta_DifferentPair");

    fSameEventTripletPhiThetaArray_DifferentPair = new TH2F*[32];

    if(fDoOnlyThreeBody){
      for(int i=0;i<16;i++){
         TString histTitle; 
        if(i<6){
          histTitle = "sameEventPhiEtaMultDifferentPair"+fParticleNames[fTripletCombinations[i][0]]+fParticleNames[fTripletCombinations[i][1]]+fParticleNames[fTripletCombinations[i][2]];
        } else {
          histTitle = "sameEventPhiEtaMultDifferentPair"+fParticleNames[fSameMixedCominations[i-6][0]]+fParticleNames[fSameMixedCominations[i-6][1]]+"Same"+fParticleNames[fSameMixedCominations[i-6][2]]+"Mixed";
        }

        fSameEventTripletPhiThetaArray_DifferentPair[i] = new TH2F(histTitle+"Before",histTitle+"Before", 500, -0.15,0.15,500,-0.15,0.15);
        fSameEventTripletPhiThetaArray_DifferentPair[16+i] = new TH2F(histTitle+"After",histTitle+"After", 500, -0.15,0.15,500,-0.15,0.15);

        if(fRunPlotPhiTheta){
          fSameEventPhiTheta_DifferentPair->Add(fSameEventTripletPhiThetaArray_DifferentPair[i]);
          fSameEventPhiTheta_DifferentPair->Add(fSameEventTripletPhiThetaArray_DifferentPair[16+i]);
        }
      }
    }

    //...............................................................................................
    // Mixed event phi theta distribution
    fMixedEventPhiTheta_SamePair = new TList();
    fMixedEventPhiTheta_SamePair->SetOwner();
    fMixedEventPhiTheta_SamePair->SetName("MixedEventPhiTheta");

    fMixedEventTripletPhiThetaArray_SamePair = new TH2F*[12];

    if(fDoOnlyThreeBody){
	    for(int i=0;i<6;i++){
        TString histTitle = "mixedEventPhiEta"+fParticleNames[fTripletCombinations[i][0]]+fParticleNames[fTripletCombinations[i][1]]+fParticleNames[fTripletCombinations[i][2]];

	      fMixedEventTripletPhiThetaArray_SamePair[i] = new TH2F(histTitle+"Before",histTitle+"Before", 500, -0.15,0.15,500,-0.15,0.15);
	      fMixedEventTripletPhiThetaArray_SamePair[6+i] = new TH2F(histTitle+"After",histTitle+"After", 500, -0.15,0.15,500,-0.15,0.15);

             if(fRunPlotPhiTheta){
      	        fMixedEventPhiTheta_SamePair->Add(fMixedEventTripletPhiThetaArray_SamePair[i]);
      	        fMixedEventPhiTheta_SamePair->Add(fMixedEventTripletPhiThetaArray_SamePair[6+i]);
             }
	    }
    }

    fMixedEventPairPhiThetaArray_TwoBody = new TH2F*[12];
    
    if(!fDoOnlyThreeBody){
	    for(int i=0;i<6;i++){

        TString histTitle = "mixeEventPhiEtaMult"+fParticleNames[fPairCombinations[i][0]]+fParticleNames[fPairCombinations[i][1]];

        
	      fMixedEventPairPhiThetaArray_TwoBody[i] = new TH2F(histTitle+"Before",histTitle+"Before", 500, -0.15,0.15,500,-0.15,0.15);
	      fMixedEventPairPhiThetaArray_TwoBody[6+i] = new TH2F(histTitle+"After",histTitle+"After", 500, -0.15,0.15,500,-0.15,0.15);

              if(fRunPlotPhiTheta){
                fMixedEventPhiTheta_SamePair->Add(fMixedEventPairPhiThetaArray_TwoBody[i]);
	              fMixedEventPhiTheta_SamePair->Add(fMixedEventPairPhiThetaArray_TwoBody[6+i]);
              }
	    }
    }

    //...............................................................................................
    // Mixed event phi theta distribution Different Pair
    fMixedEventPhiTheta_DifferentPair = new TList();
    fMixedEventPhiTheta_DifferentPair->SetOwner();
    fMixedEventPhiTheta_DifferentPair->SetName("MixedEventPhiTheta_DifferentPair");

    fMixedEventTripletPhiThetaArray_DifferentPair = new TH2F*[12];

    if(fDoOnlyThreeBody){
      for(int i=0;i<6;i++){

        TString histTitle = "mixedEventPhiEtaDifferentPair"+fParticleNames[fTripletCombinations[i][0]]+fParticleNames[fTripletCombinations[i][1]]+fParticleNames[fTripletCombinations[i][2]];

        fMixedEventTripletPhiThetaArray_DifferentPair[i] = new TH2F(histTitle+"Before",histTitle+"Before", 500, -0.15,0.15,500,-0.15,0.15);
        fMixedEventTripletPhiThetaArray_DifferentPair[6+i] = new TH2F(histTitle+"After",histTitle+"After", 500, -0.15,0.15,500,-0.15,0.15);

        if(fRunPlotPhiTheta){
          fMixedEventPhiTheta_DifferentPair->Add(fMixedEventTripletPhiThetaArray_DifferentPair[i]);
          fMixedEventPhiTheta_DifferentPair->Add(fMixedEventTripletPhiThetaArray_DifferentPair[6+i]);
        }
      }
    }

    if(fRunPlotPhiTheta){
      fResultsThreeBody->Add(fSameEventPhiTheta_SamePair);
      fResultsThreeBody->Add(fMixedEventPhiTheta_SamePair);
      fResultsThreeBody->Add(fSameEventPhiTheta_DifferentPair);
      fResultsThreeBody->Add(fMixedEventPhiTheta_DifferentPair);
     }

    //...............................................................................................
    //invariant Mass of pair (ij) vs Q3 of the triplet for (xxy) and (xx)(y) events. In case of  (xx)(y): Array(ij), i,j = 1,2 same event, 3 mixed event
    fInvMass12 = new TH2F*[16];
    fInvMass23 = new TH2F*[16];
    fInvMass31 = new TH2F*[16];
    
    if (fRunPlotInvMass){
      fInvMassList = new TList();
      fInvMassList->SetOwner();
      fInvMassList->SetName("InvMass");
      for(int i=0;i<16;i++){

        TString histTitle; 
        if(i<6){
          histTitle = "InvMass"+fParticleNames[fTripletCombinations[i][0]]+fParticleNames[fTripletCombinations[i][1]]+fParticleNames[fTripletCombinations[i][2]];
        } else {
          histTitle = "InvMass"+fParticleNames[fSameMixedCominations[i-6][0]]+fParticleNames[fSameMixedCominations[i-6][1]]+"Same"+fParticleNames[fSameMixedCominations[i-6][2]]+"Mixed";
        }

        fInvMass12[i] = new TH2F(histTitle+"_12",histTitle+"_12", 500, 0., 5., 12, 0., 1.2);
        fInvMassList->Add(fInvMass12[i]);
        fInvMass23[i] = new TH2F(histTitle+"_23",histTitle+"_23", 500, 0., 5., 12, 0., 1.2);
        fInvMassList->Add(fInvMass23[i]);
        fInvMass31[i] = new TH2F(histTitle+"_31",histTitle+"_31", 500, 0., 5., 12, 0., 1.2);
        fInvMassList->Add(fInvMass31[i]);
      }

      fResultsThreeBody->Add(fInvMassList);  
    }

    //invariant Mass of pair (ij) vs mT of pair (ij) for (xxy) and (xx)(y) events. In case of  (xx)(y): Array(ij), i,j = 1,2 same event, 3 mixed event
    fInvMassVSmT12 = new TH2F*[16];
    fInvMassVSmT23 = new TH2F*[16];
    fInvMassVSmT31 = new TH2F*[16];
    
    //invariant Mass of pair (ij) vs mT of pair (ij) for (x)(x)(y) events
    fInvMassVSmT12_MixedEvent = new TH2F*[6];
    fInvMassVSmT23_MixedEvent = new TH2F*[6];
    fInvMassVSmT31_MixedEvent = new TH2F*[6];
    
    if (fRunPlotInvMassVSmT){
      fInvMassVSmTList = new TList();
      fInvMassVSmTList->SetOwner();
      fInvMassVSmTList->SetName("InvMassVSmT");

      for(int i=0;i<16;i++){
        TString histTitle; 
        if(i<6){
          histTitle = "InvMassVSmT"+fParticleNames[fTripletCombinations[i][0]]+fParticleNames[fTripletCombinations[i][1]]+fParticleNames[fTripletCombinations[i][2]];
        } else {
          histTitle = "InvMassVSmT"+fParticleNames[fSameMixedCominations[i-6][0]]+fParticleNames[fSameMixedCominations[i-6][1]]+"Same"+fParticleNames[fSameMixedCominations[i-6][2]]+"Mixed";
        }

        fInvMassVSmT12[i] = new TH2F(histTitle+"_12",histTitle+"_12", 500, 0., 5., 100, 0., 5.);
        fInvMassVSmTList->Add(fInvMassVSmT12[i]);
        fInvMassVSmT23[i] = new TH2F(histTitle+"_23",histTitle+"_23", 500, 0., 5., 100, 0., 5.);
        fInvMassVSmTList->Add(fInvMassVSmT23[i]);
        fInvMassVSmT31[i] = new TH2F(histTitle+"_31",histTitle+"_31", 500, 0., 5., 100, 0., 5.);
        fInvMassVSmTList->Add(fInvMassVSmT31[i]);
        if(i<6){

          TString histTitle2 = "InvMassVSmTMixedEvent"+fParticleNames[fTripletCombinations[i][0]]+fParticleNames[fTripletCombinations[i][1]]+fParticleNames[fTripletCombinations[i][2]];
          fInvMassVSmT12_MixedEvent[i] = new TH2F(histTitle2+"_12",histTitle2+"_12", 500, 0., 5., 100, 0., 5.);
          fInvMassVSmTList->Add(fInvMassVSmT12_MixedEvent[i]);
          fInvMassVSmT23_MixedEvent[i] = new TH2F(histTitle2+"_23",histTitle2+"_23", 500, 0., 5., 100, 0., 5.);
          fInvMassVSmTList->Add(fInvMassVSmT23_MixedEvent[i]);
          fInvMassVSmT31_MixedEvent[i] = new TH2F(histTitle2+"_31",histTitle2+"_31", 500, 0., 5., 100, 0., 5.);
          fInvMassVSmTList->Add(fInvMassVSmT31_MixedEvent[i]);
        }
      }

      fResultsThreeBody->Add(fInvMassVSmTList);
    }


   }//if (fRunThreeBody)


  //...............................................................................................
  fResultsSampleQA = new TList();
  fResultsSampleQA->SetOwner();
  fResultsSampleQA->SetName("ResultsSampleQA");

  if (fConfig->GetUsePhiSpinning()) {
    fResultsSample = fSample->GetHistList();

    if (!fConfig->GetMinimalBookingSample()) {
      fResultsSampleQA->Add(fSample->GetQAList());
      fResultsQA->Add(fPairCleaner->GetHistList());
    }
  } else {
    fResultsSample = new TList();
    fResultsSample->SetOwner();
    fResultsSample->SetName("ResultsSample");
  }

  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fPrimaryList);
  PostData(5, fAntiPrimaryList);

  PostData(6, fResults);
  PostData(7, fResultsQA);
  PostData(8, fResultsSample);
  PostData(9, fResultsSampleQA);
  PostData(10, fResultsThreeBody);
  if (fProton->GetIsMonteCarlo()) {
    if (!fProton->GetMinimalBooking()) {
      fProtonMCList = fProton->GetMCQAHists();
    } else {
      fProtonMCList = new TList();
      fProtonMCList->SetName("MCTrkCuts");
      fProtonMCList->SetOwner();
    }
    PostData(11, fProtonMCList);
  }
  if (fAntiProton->GetIsMonteCarlo()) {
    if (!fAntiProton->GetMinimalBooking()) {
      fAntiProtonMCList = fAntiProton->GetMCQAHists();
    } else {
      fAntiProtonMCList = new TList();
      fAntiProtonMCList->SetName("MCAntiTrkCuts");
      fAntiProtonMCList->SetOwner();
    }
    PostData(12, fAntiProtonMCList);
  }

  if (fPrimary->GetIsMonteCarlo()) {
    if (!fPrimary->GetMinimalBooking()) {
      fPrimaryMCList = fPrimary->GetMCQAHists();
    } else {
      fPrimaryMCList = new TList();
      fPrimaryMCList->SetName("MCPrimaryTrkCuts");
      fPrimaryMCList->SetOwner();
    }
    PostData(13, fPrimaryMCList);
  }
  if (fAntiPrimary->GetIsMonteCarlo()) {
    if (!fAntiPrimary->GetMinimalBooking()) {
      fAntiPrimaryMCList = fAntiPrimary->GetMCQAHists();
    } else {
      fAntiPrimaryMCList = new TList();
      fAntiPrimaryMCList->SetName("MCAntiPrimaryTrkCuts");
      fAntiPrimaryMCList->SetOwner();
    }
    PostData(14, fAntiPrimaryMCList);
  }

 // Mixed event distribution ------------------------------------------------------------------------------
      // Take care of the mixing PartContainer for three particles
    auto ZVtxBinsSize = fConfig->GetNZVtxBins();
    auto MultBinsSize = fConfig->GetNMultBins();

    static std::vector<int> PDGCodes = fConfig->GetPDGCodes();

    for(int iZVtx = 0; iZVtx<ZVtxBinsSize; iZVtx++){
      std::vector<std::vector<AliFemtoDreamPartContainer>> MultContainer;

      std::vector<std::vector<AliFemtoDreamPartContainer>> MultContainerPPP;
      std::vector<std::vector<AliFemtoDreamPartContainer>> MultContainerPPPrim;
      std::vector<std::vector<AliFemtoDreamPartContainer>> MultContainerPPAPrim;

      for(int iMult = 0; iMult<MultBinsSize; iMult++){
        std::vector<AliFemtoDreamPartContainer> AllUsedParticles;

        std::vector<AliFemtoDreamPartContainer> AllUsedParticlesPPP;
        std::vector<AliFemtoDreamPartContainer> AllUsedParticlesPPPrim;
        std::vector<AliFemtoDreamPartContainer> AllUsedParticlesPPAPrim;

        for(unsigned int iSpecies = 0; iSpecies<PDGCodes.size(); iSpecies++){
          auto tempPartContainer = new AliFemtoDreamPartContainer(fConfig->GetMixingDepth());
          AllUsedParticles.push_back(*tempPartContainer);

          auto tempPartContainerPPP = new AliFemtoDreamPartContainer(fConfig->GetMixingDepth());
          auto tempPartContainerPPPrim = new AliFemtoDreamPartContainer(fConfig->GetMixingDepth());
          auto tempPartContainerPPAPrim = new AliFemtoDreamPartContainer(fConfig->GetMixingDepth());

          AllUsedParticlesPPP.push_back(*tempPartContainerPPP);
          AllUsedParticlesPPPrim.push_back(*tempPartContainerPPPrim);
          AllUsedParticlesPPAPrim.push_back(*tempPartContainerPPAPrim);
        }

        MultContainer.push_back(AllUsedParticles);

        MultContainerPPP.push_back(AllUsedParticlesPPP);
        MultContainerPPPrim.push_back(AllUsedParticlesPPPrim);
        MultContainerPPAPrim.push_back(AllUsedParticlesPPAPrim);
      }

      if(fStandardMixing){
         fPartContainer.push_back(MultContainer);
      } else {
        if(fDoOnlyThreeBody){
          fPartContainerPPP.push_back(MultContainerPPP);
          fPartContainerPPPrim.push_back(MultContainerPPPrim);
          fPartContainerPPAPrim.push_back(MultContainerPPAPrim);
        } else {
          fPartContainerPP.push_back(MultContainer);
          fPartContainerPPrim.push_back(MultContainer);
          fPartContainerPAPrim.push_back(MultContainer);
        }
      }
    }

    if(fStandardMixing){
        fVectPartContainers.push_back(fPartContainer);
    }else{
      if(fDoOnlyThreeBody){
        fVectPartContainers.push_back(fPartContainerPPP);
        fVectPartContainers.push_back(fPartContainerPPPrim);
        fVectPartContainers.push_back(fPartContainerPPAPrim);
      } else {
        fVectPartContainers.push_back(fPartContainerPP);
        fVectPartContainers.push_back(fPartContainerPPrim);
        fVectPartContainers.push_back(fPartContainerPAPrim);

      }
    }

} //void AliAnalysisTaskAODThreeBodyProtonPrimary::UserCreateOutputObjects()

//==================================================================================================================================================

void AliAnalysisTaskAODThreeBodyProtonPrimary::UserExec(Option_t *option) {
//  AliVEvent *fInputEvent = InputEvent();

//List of actions:
//	a) Particle Selections
//	   a.0) Proton
//	   a.1) General Primary (e.g. K+, K-, Pi+, Pi-)
//      b) Optional pair cleaning
//      c) Start Three Body Calculus
//	   c.0) Same event distribution
//	   c.1) Mixed event distribution
//		c.1.0) Same 2, Mixed 1
//		c.1.1) Normal mixed
//	   c.2)

  AliAODEvent *Event = static_cast<AliAODEvent*>(InputEvent());

  if (!Event) {
    AliWarning("No Input Event");
    return;
  } 

  fEvent->SetEvent(Event);
  if (!fEventCuts->isSelected(fEvent)) {
       return;
  }

  //a) Particle Selection +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
    AliAODTrack *track = static_cast<AliAODTrack*>(Event->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return;
    }
    StoreGlobalTrackReference(track);
  }//for (int iTrack = 0; iTrack < fInputEvent->GetNumberOfTracks(); ++iTrack)

  //a.0) Proton Selection ------------------------------
  std::vector<AliFemtoDreamBasePart> Protons;
  std::vector<AliFemtoDreamBasePart> AntiProtons;
  std::vector<AliFemtoDreamBasePart> Primaries;
  std::vector<AliFemtoDreamBasePart> AntiPrimaries;

  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);


  for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) { // TO DO: think about double track selection
    AliAODTrack *track = static_cast<AliAODTrack*>(Event->GetTrack(iTrack));
    fTrack->SetTrack(track);

    if (fIsMC && fRemoveMCResonances) { 
        TClonesArray *mcarray = dynamic_cast<TClonesArray *>(Event->FindListObject(AliAODMCParticle::StdBranchName()));
        if (!mcarray) {
          AliError("SPTrack: MC Array not found");
        }
        if (fTrack->GetID() >= 0) {
          AliAODMCParticle *mcPart = (AliAODMCParticle *)mcarray->At(fTrack->GetID());
          if (!(mcPart)) {
            continue;
          }
          if(IsResonance(mcPart->GetPdgCode())){
             continue; 
          }
          int motherID = mcPart->GetMother();
          int lastMother = motherID;
          AliAODMCParticle *mcMother = nullptr;
          /*bool RemoveTrack = false;
          while (motherID != -1) {
            lastMother = motherID;
            mcMother = (AliAODMCParticle *)mcarray->At(motherID);
            motherID = mcMother->GetMother();
            if(IsResonance(mcMother->GetPdgCode())){
               fTrack->SetMotherPDG(mcMother->GetPdgCode()); //Change the PDG of the mother so it is set to the resonance. The Mother ID keeps set to the original parton
               RemoveTrack = true;
            }
          }*/ 
          if ((lastMother != -1)) {
            mcMother = (AliAODMCParticle *)mcarray->At(lastMother);
          }
          if (mcMother) {
            int motherPDG = mcMother->GetPdgCode(); 
            if(IsResonance(motherPDG)){
              fTrack->SetMotherPDG(motherPDG); //Change the PDG of the mother so it is set to the resonance. The Mother ID keeps set to the original parton
              //RemoveTrack = true;
            }
          }
          //if (RemoveTrack && fRemoveMCResonanceDaughters){
          //  continue; 
          //}
        } else {
          continue;  // if we don't have MC Information, don't use that track
        }
      } //if (fIsMC && fRemoveMCResonances)

      if (fProton->isSelected(fTrack)) {
        Protons.push_back(*fTrack);
      }
      if (fAntiProton->isSelected(fTrack)) {
        AntiProtons.push_back(*fTrack);
      }
      if (fPrimary->isSelected(fTrack)) {
        Primaries.push_back(*fTrack);
      }
      if (fAntiPrimary->isSelected(fTrack)) {
        AntiPrimaries.push_back(*fTrack);
      }
  }

  //b) Optional pair cleaning +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  fPairCleaner->ResetArray();
  fPairCleaner->StoreParticle(Protons);
  fPairCleaner->StoreParticle(AntiProtons);
  fPairCleaner->StoreParticle(Primaries);
  fPairCleaner->StoreParticle(AntiPrimaries);

  //c) Start Three Body Calculus +++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  if(fRunThreeBody){
    static std::vector<int> PDGCodes = fConfig->GetPDGCodes();
    int bins[2] = { 0, 0 };
    float ZVtx = fEvent->GetZVertex();
    float Mult = fEvent->GetMultiplicity();
    fPartColl->FindBin(ZVtx, Mult, bins);

    // c.0) Same event distribution --------------------------------
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector = fPairCleaner->GetCleanParticles();

    if(fDoOnlyThreeBody){
       for(int iComb=0; iComb<6; iComb++){
         FillTripletDistribution( ParticleVector, fTripletCombinations[iComb][0], fTripletCombinations[iComb][1], fTripletCombinations[iComb][2], 
                                 fSameEventTripletArray[iComb],PDGCodes, 
                                 bins[1], fSameEventTripletMultArray[iComb], fSameEventTripletMultArray12[iComb] ,fSameEventTripletMultArray23[iComb], fSameEventTripletMultArray31[iComb], 
                                 fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair,iComb, 
                                 *fConfig,
                                 fInvMass12[iComb], fInvMass23[iComb], fInvMass31[iComb], 
                                 fSameEventTripletmTArray12[iComb], fSameEventTripletmTArray23[iComb], fSameEventTripletmTArray31[iComb], 
                                 fInvMassVSmT12[iComb], fInvMassVSmT23[iComb], fInvMassVSmT31[iComb], 
                                 fSameEventTripletResolution[iComb], fSameEventTripletResolutionAll[iComb]);
       }      
    }//if(fDoOnlyThreeBody)
    
    else {
      //Two Body Analyses...........
      for(int iComb=0; iComb<6; iComb++){
         FillPairDistribution( ParticleVector, fPairCombinations[iComb][0], fPairCombinations[iComb][1], 
                              fSameEventPairArray_TwoBody[iComb],PDGCodes, 
                              bins[1],fSameEventPairMultArray_TwoBody[iComb], fSameEventDphiArray_TwoBody[iComb], fSameEventMultArray_TwoBody[iComb], 
                              fSameEventPairPhiThetaArray_TwoBody,iComb, *fConfig);

       if(fRunOfficialTwoBody){
         FillPairTransverseMass(ParticleVector, fPairCombinations[iComb][0], fPairCombinations[iComb][1], fPairTranverseMass_TwoBody[iComb], PDGCodes, fPairTranverseMassVSkstar_TwoBody[iComb]);
       }
      }
    }//else

    // c.1) Mixed event distribution --------------------------------

    // TAKE CARE OF MULT AND ZVtx!
    if (!(bins[0] == -99 || bins[1] == -99)) {

      std::vector<std::vector<AliFemtoDreamPartContainer>*> VectItMult;

      if(fStandardMixing){
        auto itZVtx = fVectPartContainers[0].begin()+ bins[0];
        auto itMult = itZVtx->begin() + bins[1];
        VectItMult.push_back(&(*itMult));
      } else {
        for(int i=0; i<3; i++){
          auto itZVtx = fVectPartContainers[i].begin()+ bins[0];
          auto itMult = itZVtx->begin() + bins[1];
          VectItMult.push_back(&(*itMult));
        }

      }

      //c.1.0) Same 2, Mixed 1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      int ContainerId = 0; 

      if(fDoOnlyThreeBody){

        for(int iComb=0; iComb<10; iComb++){
           if(fStandardMixing){
             ContainerId = 0;
           }else{
             ContainerId = fSameMixedContainerID[iComb];
           } 

           FillTripletDistributionSE2ME1(ParticleVector, *VectItMult[ContainerId], fSameMixedCominations[iComb][0], fSameMixedCominations[iComb][1], fSameMixedCominations[iComb][2], 
                                        fSameEventTripletArray[iComb+6], PDGCodes,
                                        bins[1], fSameEventTripletMultArray[iComb+6], fSameEventTripletMultArray12[iComb+6],
                                        fSameEventTripletPhiThetaArray_SamePair, fSameEventTripletPhiThetaArray_DifferentPair, iComb+6, 
                                        *fConfig, 
                                        fInvMass12[iComb+6], fInvMass23[iComb+6], fInvMass31[iComb+6], 
                                        fSameEventTripletmTArray12[iComb+6], fSameEventTripletmTArray23[iComb+6], fSameEventTripletmTArray31[iComb+6], 
                                        fInvMassVSmT12[iComb+6], fInvMassVSmT23[iComb+6], fInvMassVSmT31[iComb+6]);

        }
      }//if(fDoOnlyThreeBody)

      //c.1.1) Normal mixing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      if(fDoOnlyThreeBody){
        for(int iComb=0; iComb<6; iComb++){
           if(fStandardMixing){
             ContainerId = 0;
           }else{
             ContainerId = fTripletContainerID[iComb];
           } 

          FillTripletDistributionME(ParticleVector, *VectItMult[ContainerId], fTripletCombinations[iComb][0], fTripletCombinations[iComb][1], fTripletCombinations[iComb][2],
                                   fMixedEventTripletArray[iComb], PDGCodes, 
                                   bins[1], fMixedEventTripletMultArray[iComb], fMixedEventTripletMultArray12[iComb], fMixedEventTripletMultArray23[iComb], fMixedEventTripletMultArray31[iComb], 
                                   fMixedEventTripletPhiThetaArray_SamePair, fMixedEventTripletPhiThetaArray_DifferentPair, iComb, *fConfig, 
                                   fMixedEventTripletmTArray12[iComb], fMixedEventTripletmTArray23[iComb], fMixedEventTripletmTArray31[iComb], 
                                   fInvMassVSmT12_MixedEvent[iComb], fInvMassVSmT23_MixedEvent[iComb], fInvMassVSmT31_MixedEvent[iComb], fProjectorData[iComb], 
                                   fMixedEventTripletResolution[iComb], fMixedEventTripletResolutionAll[iComb]);

        }  
      } else {
        //Two Body Analyses...........
        for(int iComb=0; iComb<6; iComb++){
           if(fStandardMixing){
             ContainerId = 0;
           }else{
             ContainerId = fPairContainerID[iComb];
           } 
          
         FillPairDistributionME( ParticleVector, *VectItMult[ContainerId], fPairCombinations[iComb][0], fPairCombinations[iComb][1],
                                fMixedEventPairArray_TwoBody[iComb],PDGCodes, 
                                bins[1],fMixedEventPairMultArray_TwoBody[iComb], fMixedEventDphiArray_TwoBody[iComb], 
                                fMixedEventMultArray_TwoBody[iComb], fMixedEventPairPhiThetaArray_TwoBody,iComb, *fConfig);
        
        }
      } //else 

      // Update the particle container with current event
      if(fStandardMixing) {
        SetMixedEvent(ParticleVector, VectItMult[0]);
      }else{
        if(fDoOnlyThreeBody){
          SetMixedEventPPP(ParticleVector, VectItMult[0]);
          SetMixedEventPPPrim(ParticleVector, VectItMult[1]);
          SetMixedEventPPAPrim(ParticleVector, VectItMult[2]);
        } else {
          SetMixedEventPP(ParticleVector, VectItMult[0]);
          SetMixedEventPPrim(ParticleVector, VectItMult[1]);
          SetMixedEventPAPrim(ParticleVector, VectItMult[2]);
        }
      }
    }//if (!(bins[0] == -99 || bins[1] == -99))

  }//if(fRunThreeBody)


  //----- OFFICIAL CODE for the 2-body -----------------
  if(fRunOfficialTwoBody){  
    if (fPairCleaner->GetCounter() > 0) {
      if (fConfig->GetUseEventMixing()) {
        fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),
                            fEvent);
      }
      if (fConfig->GetUsePhiSpinning()) {
        fSample->SetEvent(fPairCleaner->GetCleanParticles(), fEvent);
      }
    }
  }
  //-----------------------------------------------------

  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fPrimaryList);
  PostData(5, fAntiPrimaryList);

  if(fRunOfficialTwoBody){
    PostData(6, fResults);
    PostData(7, fResultsQA);
    PostData(8, fResultsSample);
    PostData(9, fResultsSampleQA);
  }
  PostData(10, fResultsThreeBody);
  if (fProton->GetIsMonteCarlo()) {
    PostData(11, fProtonMCList);
  }
  if (fAntiProton->GetIsMonteCarlo()) {
    PostData(12, fAntiProtonMCList);
  }
  if (fPrimary->GetIsMonteCarlo()) {
    PostData(13, fPrimaryMCList);
  }
  if (fAntiPrimary->GetIsMonteCarlo()) {
    PostData(14, fAntiPrimaryMCList);
  }

}

//==================================================================================================================================================

void AliAnalysisTaskAODThreeBodyProtonPrimary::ResetGlobalTrackReference() {

  for (UShort_t i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//==================================================================================================================================================

void AliAnalysisTaskAODThreeBodyProtonPrimary::StoreGlobalTrackReference(AliAODTrack *track){
  // see AliFemtoDreamAnalysis for details
  const int trackID = track->GetID();
  if (trackID < 0) {
    return;
  }
  if (trackID >= fTrackBufferSize) {
    printf("Warning: track ID too big for buffer: ID: %d, buffer %d\n", trackID,
           fTrackBufferSize);
    return;
  }

  if (fGTI[trackID]) {
    if ((!track->GetFilterMap()) && (!track->GetTPCNcls())) {
      return;
    }
    if ((fGTI[trackID])->GetFilterMap() || fGTI[trackID]->GetTPCNcls()) {
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             (fGTI[trackID])->GetFilterMap(), track->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}

//==================================================================================================================================================

TLorentzVector AliAnalysisTaskAODThreeBodyProtonPrimary::RelativePairMomentum(
    TLorentzVector &PartOne, TLorentzVector &PartTwo) {
  // The q12 components can be calculated as:
  //             q^mu = (p1-p2)^mu /2 - ((p1-p2)*P/(2P^2))*P^mu
  //             where P = p1+p2
  // Reference: https://www.annualreviews.org/doi/pdf/10.1146/annurev.nucl.55.090704.151533
  // In the following code the above written equation will be expressed as:
  //             q = trackDifference/2 -  scaling * trackSum
  // where scaling is a float number:
  //             scaling = trackDifference*trackSum/(2*trackSum^2) = ((p1-p2)*P/(2P^2))
  // PR on 15-06-2020 -> don't use the reduced vector - no division by 2
  TLorentzVector trackSum = PartOne + PartTwo;
  TLorentzVector trackDifference = PartOne - PartTwo;
  float scaling = trackDifference*trackSum/(trackSum*trackSum);
  TLorentzVector qPartOnePartTwo = trackDifference -  scaling * trackSum;

  return qPartOnePartTwo;
}

//==================================================================================================================================================

float AliAnalysisTaskAODThreeBodyProtonPrimary::BoostOneParticle(
    TLorentzVector &PartOne, TLorentzVector &PartTwo, TLorentzVector &PartThree) {
    // The code boost PartThree in the rest frame of PartOne and PartTwo
  TLorentzVector trackSum = PartOne + PartTwo;

  float beta = trackSum.Beta();
  float betax = beta * cos(trackSum.Phi()) * sin(trackSum.Theta());
  float betay = beta * sin(trackSum.Phi()) * sin(trackSum.Theta());
  float betaz = beta * cos(trackSum.Theta());

  TLorentzVector PartThreeCMS = PartThree;

  PartThreeCMS.Boost(-betax, -betay, -betaz);

  return PartThreeCMS.P();
}

//==================================================================================================================================================

void AliAnalysisTaskAODThreeBodyProtonPrimary::FillTripletDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F* hist2d12, TH2F* hist2d23, TH2F* hist2d31, TH2F **fEventTripletPhiThetaArray_SamePair, TH2F **fEventTripletPhiThetaArray_DifferentPair, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* InvMass12, TH2F* InvMass23, TH2F* InvMass31, TH2F* histmTQ312, TH2F* histmTQ323, TH2F* histmTQ331, TH2F* InvMassVsmT12, TH2F* InvMassVsmT23, TH2F* InvMassVsmT31, TH2F* Res, TH2F* ResAll){
  // This function creates a triplet distribution in Q3 bins (defined lower).
  // It requires the particle vector from PairCleaner() and the three indices of particles of interest. So
  // if you want to get distribution for particles that are saved in particle vector as 1 2 3 element, just
  // call the function with firstSpecies=1,secondSpecies=2,thirdSpecies=3

  auto Particle1Vector = ParticleVector.begin()+firstSpecies;
  auto Particle2Vector = ParticleVector.begin()+secondSpecies;
  auto Particle3Vector = ParticleVector.begin()+thirdSpecies;

  // Get the PID codes std::vector<int>
  auto itPDGPar1 = PDGCodes.begin()+firstSpecies;
  auto itPDGPar2 = PDGCodes.begin()+secondSpecies;
  auto itPDGPar3 = PDGCodes.begin()+thirdSpecies;
  // Get particle masses
  auto massparticle1 = TDatabasePDG::Instance()->GetParticle(*itPDGPar1)->Mass();
  auto massparticle2 = TDatabasePDG::Instance()->GetParticle(*itPDGPar2)->Mass();
  auto massparticle3 = TDatabasePDG::Instance()->GetParticle(*itPDGPar3)->Mass();
  unsigned int DaughterPart1 = 0;
  unsigned int DaughterPart2 = 0;
  unsigned int DaughterPart3 = 0;
  if(abs(*itPDGPar1)==2212) DaughterPart1 = 1;
  if(abs(*itPDGPar1)==211) DaughterPart1 = 1;
  if(abs(*itPDGPar2)==2212) DaughterPart2 = 1;
  if(abs(*itPDGPar2)==211) DaughterPart2 = 1;
  if(abs(*itPDGPar3)==2212) DaughterPart3 = 1;
  if(abs(*itPDGPar3)==211) DaughterPart3 = 1;
  unsigned int DoThisPair12 = DaughterPart1*10+DaughterPart2;
  unsigned int DoThisPair23 = DaughterPart2*10+DaughterPart3;
  unsigned int DoThisPair31 = DaughterPart3*10+DaughterPart1;

  // loop over first particle
  for (auto iPart1 = Particle1Vector->begin(); iPart1 != Particle1Vector->end(); ++iPart1) {
    // if second particle species is different than first - start with the first particle in the vector
    auto iPart2 = Particle2Vector->begin();
    // if second particle  and first are the species, start second loop from the next particle (to not double count)
    if (firstSpecies==secondSpecies) iPart2 = iPart1+1;
    // loop over second particle ...
    for (; iPart2 != Particle2Vector->end(); ++iPart2) {
      auto iPart3 = Particle3Vector->begin();
      if (firstSpecies==thirdSpecies) iPart3 = iPart1+1;
      if (secondSpecies==thirdSpecies) iPart3 = iPart2+1;
      for ( ; iPart3 != Particle3Vector->end(); ++iPart3) {


        // Now we have the three particles, lets create their Lorentz vectors
        TLorentzVector part1_LorVec, part2_LorVec, part3_LorVec;
        part1_LorVec.SetPxPyPzE(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(),
        iPart1->GetMomentum().Z(), sqrt(pow(iPart1->GetP(),2)+pow(massparticle1,2)));
        part2_LorVec.SetPxPyPzE(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(),
        iPart2->GetMomentum().Z(), sqrt(pow(iPart2->GetP(),2)+pow(massparticle2,2)));
        part3_LorVec.SetPxPyPzE(iPart3->GetMomentum().X(), iPart3->GetMomentum().Y(),
        iPart3->GetMomentum().Z(), sqrt(pow(iPart3->GetP(),2)+pow(massparticle3,2)));

        // Now when we have the lorentz vectors, we can calculate the Lorentz invariant relative momenta q12, q23, q31
        TLorentzVector q12 = AliAnalysisTaskAODThreeBodyProtonPrimary::RelativePairMomentum(part1_LorVec,part2_LorVec);
        TLorentzVector q23 = AliAnalysisTaskAODThreeBodyProtonPrimary::RelativePairMomentum(part2_LorVec,part3_LorVec);
        TLorentzVector q31 = AliAnalysisTaskAODThreeBodyProtonPrimary::RelativePairMomentum(part3_LorVec,part1_LorVec);

        float RelativeMomentum12 = AliFemtoDreamHigherPairMath::RelativePairMomentum(part1_LorVec,part2_LorVec);
        float RelativeMomentum23 = AliFemtoDreamHigherPairMath::RelativePairMomentum(part2_LorVec,part3_LorVec);
        float RelativeMomentum31 = AliFemtoDreamHigherPairMath::RelativePairMomentum(part3_LorVec,part1_LorVec);

        // The particles in current methodology are put in bins of:
        //                 Q3=sqrt(q12^2+q23^2+q31^2)
        float Q32 = q12*q12+q23*q23+q31*q31;
        // From 3 pion paper, the q must be multiplied by -1 before taking quare root
        float Q3 = sqrt(-Q32); // the minus from pion paper



        bool Pair12 = true;
        bool Pair23 = true;
        bool Pair31 = true;

        if(!fturnoffClosePairRejectionCompletely){

          if(fClosePairRejectionForAll){
            if(abs(*itPDGPar1)==abs(*itPDGPar2)){
              Pair12 = DeltaEtaDeltaPhi(firstSpecies, secondSpecies, *iPart1,*iPart2, *itPDGPar1, *itPDGPar2, true,  DoThisPair12, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[16+phiEtaHistNo],Config, Q3);
            }else{
              Pair12 = DeltaEtaDeltaPhi(firstSpecies, secondSpecies,*iPart1,*iPart2, *itPDGPar1, *itPDGPar2, true,  DoThisPair12, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[16+phiEtaHistNo],Config, Q3);
            }
            if(abs(*itPDGPar2)==abs(*itPDGPar3)){
              Pair23 = DeltaEtaDeltaPhi(secondSpecies, thirdSpecies, *iPart2,*iPart3, *itPDGPar2, *itPDGPar3, true,  DoThisPair23, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[16+phiEtaHistNo],Config, Q3);
            }else{
              Pair23 = DeltaEtaDeltaPhi(secondSpecies, thirdSpecies, *iPart2,*iPart3, *itPDGPar2, *itPDGPar3, true,  DoThisPair23, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[16+phiEtaHistNo],Config, Q3);
            }
            if(abs(*itPDGPar3)==abs(*itPDGPar1)){
              Pair31 = DeltaEtaDeltaPhi(thirdSpecies, firstSpecies, *iPart3,*iPart1, *itPDGPar3, *itPDGPar1, true,  DoThisPair31, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[16+phiEtaHistNo],Config, Q3);
            }else{
              Pair31 = DeltaEtaDeltaPhi(thirdSpecies, firstSpecies, *iPart3,*iPart1, *itPDGPar3, *itPDGPar1, true,  DoThisPair31, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[16+phiEtaHistNo],Config, Q3);
            }
          }
          else {
           
              if(DoThisPair12==11){
                if(abs(*itPDGPar1)==abs(*itPDGPar2)){
                  Pair12 = DeltaEtaDeltaPhi(firstSpecies, secondSpecies, *iPart1,*iPart2, *itPDGPar1, *itPDGPar2, true,  DoThisPair12, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[16+phiEtaHistNo],Config, Q3);
                }else{
                  Pair12 = DeltaEtaDeltaPhi(firstSpecies, secondSpecies, *iPart1,*iPart2, *itPDGPar1, *itPDGPar2, true,  DoThisPair12, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[16+phiEtaHistNo],Config, Q3);
                }
              }
              if(DoThisPair23==11){
                if(abs(*itPDGPar2)==abs(*itPDGPar3)){
                  Pair23 = DeltaEtaDeltaPhi(secondSpecies, thirdSpecies, *iPart2,*iPart3, *itPDGPar2, *itPDGPar3, true,  DoThisPair23, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[16+phiEtaHistNo],Config, Q3);
                }else{
                  Pair23 = DeltaEtaDeltaPhi(secondSpecies, thirdSpecies, *iPart2,*iPart3, *itPDGPar2, *itPDGPar3, true,  DoThisPair23, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[16+phiEtaHistNo],Config, Q3);
                }
              }
              if(DoThisPair31==11){
                if(abs(*itPDGPar3)==abs(*itPDGPar1)){
                  Pair31 = DeltaEtaDeltaPhi(thirdSpecies, firstSpecies, *iPart3,*iPart1, *itPDGPar3, *itPDGPar1, true,  DoThisPair31, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[16+phiEtaHistNo],Config, Q3);
                }else{
                  Pair31 = DeltaEtaDeltaPhi(thirdSpecies, firstSpecies, *iPart3,*iPart1, *itPDGPar3, *itPDGPar1, true,  DoThisPair31, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[16+phiEtaHistNo],Config, Q3);
                }
              }
          }

        }


        if(!Pair12||!Pair23||!Pair31) {continue;}

        if(fIsMC && fRemoveMCResonances){
          bool CommonMother12 = CommonMotherResonance(&(*iPart1), &(*iPart2));
          bool CommonMother13 = CommonMotherResonance(&(*iPart1), &(*iPart3)); 
          bool CommonMother23 = CommonMotherResonance(&(*iPart2), &(*iPart3)); 

          if(CommonMother12 || CommonMother13 || CommonMother23){
            continue; 
          }
        }

        hist->Fill(Q3);
        hist2d->Fill(Q3,mult+1);


        if(fRunPairMultThreeBody && fQ3MinValue <= Q3 && Q3<fQ3cutValue){
          hist2d12->Fill(RelativeMomentum12,mult+1);
          hist2d23->Fill(RelativeMomentum23,mult+1);
          hist2d31->Fill(RelativeMomentum31,mult+1);
        }
   
        if(fRunPlotInvMass)
        {
           TLorentzVector Sum12 = part1_LorVec + part2_LorVec;
           TLorentzVector Sum23 = part2_LorVec + part3_LorVec;
           TLorentzVector Sum31 = part3_LorVec + part1_LorVec;
           InvMass12->Fill(Sum12.Mag(),Q3);
           InvMass23->Fill(Sum23.Mag(),Q3);
           InvMass31->Fill(Sum31.Mag(),Q3);

        }//if(fRunPlotInvMass)

        if(fRunmTPlots){
          float mT12 = GetmT(part1_LorVec, massparticle1, part2_LorVec, massparticle2);
          float mT23 = GetmT(part2_LorVec, massparticle2, part3_LorVec, massparticle3);
          float mT31 = GetmT(part3_LorVec, massparticle3, part1_LorVec, massparticle1);

          histmTQ312->Fill(Q3,mT12);
          histmTQ323->Fill(Q3,mT23);
          histmTQ331->Fill(Q3,mT31);
        }

        if(fRunPlotInvMassVSmT){
          TLorentzVector Sum12 = part1_LorVec + part2_LorVec;
          TLorentzVector Sum23 = part2_LorVec + part3_LorVec;
          TLorentzVector Sum31 = part3_LorVec + part1_LorVec;
          float mT12 = GetmT(part1_LorVec, massparticle1, part2_LorVec, massparticle2);
          float mT23 = GetmT(part2_LorVec, massparticle2, part3_LorVec, massparticle3);
          float mT31 = GetmT(part3_LorVec, massparticle3, part1_LorVec, massparticle1);
          InvMassVsmT12->Fill(Sum12.Mag(),mT12);
          InvMassVsmT23->Fill(Sum23.Mag(),mT23);
          InvMassVsmT31->Fill(Sum31.Mag(),mT31);

        }

         if(fIsMC && fGetMomentumResolution){
           MomentumResolution(ResAll, Res, *iPart1, *itPDGPar1, massparticle1,*iPart2, *itPDGPar2, massparticle2, *iPart3,*itPDGPar3, massparticle3, Q3) ;
         }
      }
    }
  }
}

//==================================================================================================================================================
void AliAnalysisTaskAODThreeBodyProtonPrimary::FillPairDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F* histDeltaPhi, TH2F* histLowDeltaPhi, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){//, TH2F* InvMassSame){
  // Two Body Same Event

  auto Particle1Vector = ParticleVector.begin()+firstSpecies;
  auto Particle2Vector = ParticleVector.begin()+secondSpecies;

  // Get the PID codes std::vector<int>
  auto itPDGPar1 = PDGCodes.begin()+firstSpecies;
  auto itPDGPar2 = PDGCodes.begin()+secondSpecies;

  // Get particle masses
  auto massparticle1 = TDatabasePDG::Instance()->GetParticle(*itPDGPar1)->Mass();
  auto massparticle2 = TDatabasePDG::Instance()->GetParticle(*itPDGPar2)->Mass();

  unsigned int DaughterPart1 = 0;
  unsigned int DaughterPart2 = 0;

  if(abs(*itPDGPar1)==2212) DaughterPart1 = 1;
  if(abs(*itPDGPar1)==211) DaughterPart1 = 1;
  if(abs(*itPDGPar2)==2212) DaughterPart2 = 1;
  if(abs(*itPDGPar2)==211) DaughterPart2 = 1;


  unsigned int DoThisPair12 = DaughterPart1*10+DaughterPart2;

  // loop over first particle
  for (auto iPart1 = Particle1Vector->begin(); iPart1 != Particle1Vector->end(); ++iPart1) {
    // if second particle species is different than first - start with the first particle in the vector
    auto iPart2 = Particle2Vector->begin();
    // if second particle  and first are the species, start second loop from the next particle (to not double count)
    if (firstSpecies==secondSpecies) iPart2 = iPart1+1;
    // loop over second particle ...
    for (; iPart2 != Particle2Vector->end(); ++iPart2) {

        // Now we have the three particles, lets create their Lorentz vectors
        TLorentzVector Particle1_LV, Particle2_LV;
        Particle1_LV.SetXYZM(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(),iPart1->GetMomentum().Z(), massparticle1);
        Particle2_LV.SetXYZM(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(),iPart2->GetMomentum().Z(), massparticle2);
        // Get momentum
        float RelativeMomentum = AliFemtoDreamHigherPairMath::RelativePairMomentum(Particle1_LV, Particle2_LV);

        bool Pair12 = true;

        if(!fturnoffClosePairRejectionCompletely){

          if(fClosePairRejectionForAll){
             Pair12 = DeltaEtaDeltaPhi(firstSpecies, secondSpecies, *iPart1,*iPart2, *itPDGPar1, *itPDGPar2, true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[6+phiEtaHistNo],Config); 
          }
          else{

             if(DoThisPair12==11){
               Pair12 = DeltaEtaDeltaPhi(firstSpecies, secondSpecies, *iPart1,*iPart2, *itPDGPar1, *itPDGPar2, true,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[6+phiEtaHistNo],Config); 
             }
          } 
        }

        if(!Pair12) {continue;}

        hist->Fill(RelativeMomentum);
        hist2d->Fill(RelativeMomentum,mult+1);

        double deltaPhi = (iPart2->GetPhi().at(0) - iPart1->GetPhi().at(0));

        histDeltaPhi->Fill(RelativeMomentum,deltaPhi);
        if(abs(deltaPhi)>0.5&&abs(deltaPhi)<2.5){
          histLowDeltaPhi->Fill(RelativeMomentum,mult+1);
        }

    }
  }
} //void AliAnalysisTaskAODThreeBodyProtonPrimary::FillPairDistribution(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config)

//==================================================================================================================================================

//void AliAnalysisTaskAODThreeBodyProtonPrimary::FillTripletDistributionPPL(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies,int thirdSpecies, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* InvMassSame, TH2F* InvMassDET,TH2F* InvMassPDG)

//==================================================================================================================================================

void AliAnalysisTaskAODThreeBodyProtonPrimary::SetMixedEvent(
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  for(unsigned int iSpecies = 0; iSpecies<ParticleVector.size(); iSpecies++){
    if ((ParticleVector.begin()+iSpecies)->size() > 0) {
      (fPartContainer->begin()+iSpecies)->SetEvent(*(ParticleVector.begin()+iSpecies));
    }
  }
}

//==================================================================================================================================================

void AliAnalysisTaskAODThreeBodyProtonPrimary::SetMixedEventPP(
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  // THIS WORKS ONLY IF 0 and 2 is proton and lambda, 1 and 3 is antiproton antilambda.
  // Fill the particles only if both lambda and proton are present in the event, so later on for mixing
  // one would be able to know what the lambda and proton are not from the same event
  if ((ParticleVector.begin())->size() > 1) {
    (fPartContainer->begin())->SetEvent(*(ParticleVector.begin()));
  }
  if ((ParticleVector.begin()+1)->size() > 1) {
    (fPartContainer->begin()+1)->SetEvent(*(ParticleVector.begin()+1));
  }
}

//==================================================================================================================================================

void AliAnalysisTaskAODThreeBodyProtonPrimary::SetMixedEventPPrim(
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  // THIS WORKS ONLY IF 0 and 2 is proton and lambda, 1 and 3 is antiproton antilambda.
  // Fill the particles only if both lambda and proton are present in the event, so later on for mixing
  // one would be able to know what the lambda and proton are not from the same event
  if ((ParticleVector.begin())->size() > 0 && (ParticleVector.begin()+2)->size() > 0) {
    (fPartContainer->begin())->SetEvent(*(ParticleVector.begin()));
    (fPartContainer->begin()+2)->SetEvent(*(ParticleVector.begin()+2));
  }
  if ((ParticleVector.begin()+1)->size() > 0 && (ParticleVector.begin()+3)->size() > 0) {
    (fPartContainer->begin()+1)->SetEvent(*(ParticleVector.begin()+1));
    (fPartContainer->begin()+3)->SetEvent(*(ParticleVector.begin()+3));
  }
}

//==================================================================================================================================================

void AliAnalysisTaskAODThreeBodyProtonPrimary::SetMixedEventPAPrim(  
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  // THIS WORKS ONLY IF 0 and 2 is proton and lambda, 1 and 3 is antiproton antilambda.
  // Fill the particles only if both lambda and proton are present in the event, so later on for mixing
  // one would be able to know what the lambda and proton are not from the same event
  if ((ParticleVector.begin())->size() > 0 && (ParticleVector.begin()+3)->size() > 0) {
    (fPartContainer->begin())->SetEvent(*(ParticleVector.begin()));
    (fPartContainer->begin()+3)->SetEvent(*(ParticleVector.begin()+3));
  }
  if ((ParticleVector.begin()+1)->size() > 0 && (ParticleVector.begin()+2)->size() > 0) {
    (fPartContainer->begin()+1)->SetEvent(*(ParticleVector.begin()+1));
    (fPartContainer->begin()+2)->SetEvent(*(ParticleVector.begin()+2));
  }
}

//==================================================================================================================================================

void AliAnalysisTaskAODThreeBodyProtonPrimary::SetMixedEventPPP(
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  // THIS WORKS ONLY IF 0 and 2 is proton and primary, 1 and 3 is antiproton antiprimary.
  // Fill the particles only if both lambda and proton are present in the event, so later on for mixing
  // one would be able to know what the lambda and proton are not from the same event
  if ((ParticleVector.begin())->size() > 2) {
    (fPartContainer->begin())->SetEvent(*(ParticleVector.begin()));
  }
  if ((ParticleVector.begin()+1)->size() > 2) {
    (fPartContainer->begin()+1)->SetEvent(*(ParticleVector.begin()+1));
  }
}
//==================================================================================================================================================

void AliAnalysisTaskAODThreeBodyProtonPrimary::SetMixedEventPPPrim(
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  // THIS WORKS ONLY IF 0 and 2 is proton and primary, 1 and 3 is antiproton antiprimary.
  // Fill the particles only if at least 2 protons and a primary particle are present in the event
  if ((ParticleVector.begin())->size() > 1 && (ParticleVector.begin()+2)->size() > 0) {
    (fPartContainer->begin())->SetEvent(*(ParticleVector.begin()));
    (fPartContainer->begin()+2)->SetEvent(*(ParticleVector.begin()+2));
  }
  if ((ParticleVector.begin()+1)->size() > 1 && (ParticleVector.begin()+3)->size() > 0) {
    (fPartContainer->begin()+1)->SetEvent(*(ParticleVector.begin()+1));
    (fPartContainer->begin()+3)->SetEvent(*(ParticleVector.begin()+3));
  }
}

//==================================================================================================================================================

void AliAnalysisTaskAODThreeBodyProtonPrimary::SetMixedEventPPAPrim(
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *fPartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  // THIS WORKS ONLY IF 0 and 2 is proton and primary, 1 and 3 is antiproton antiprimary.
  // Fill the particles only if at least 2 prontons and a primary particle are present in the event
  if ((ParticleVector.begin())->size() > 1 && (ParticleVector.begin()+3)->size() > 0) {
    (fPartContainer->begin())->SetEvent(*(ParticleVector.begin()));
    (fPartContainer->begin()+3)->SetEvent(*(ParticleVector.begin()+3));
  }
  if ((ParticleVector.begin()+1)->size() > 1 && (ParticleVector.begin()+2)->size() > 0) {
    (fPartContainer->begin()+1)->SetEvent(*(ParticleVector.begin()+1));
    (fPartContainer->begin()+2)->SetEvent(*(ParticleVector.begin()+2));
  }
}

//==================================================================================================================================================

void AliAnalysisTaskAODThreeBodyProtonPrimary::FillTripletDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, int speciesME2, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F* hist2d12, TH2F* hist2d23, TH2F* hist2d31, TH2F **fEventTripletPhiThetaArray_SamePair, TH2F **fEventTripletPhiThetaArray_DifferentPair, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* histmTQ312, TH2F* histmTQ323, TH2F* histmTQ331, TH2F* InvMassVsmT12, TH2F* InvMassVsmT23, TH2F* InvMassVsmT31, TH2F* Projector, TH2F* Res, TH2F* ResAll){//, TH2F* InvMassMixed, TH2F* Q3VskDistribution12Mixed, TH2F*  Q3VskDistribution23Mixed){
  // Description of function given in AliAnalysisTaskAODThreeBodyProtonPrimary::FillTripletDistribution
  // In this function, only one particle is used from current event, and the other two - from other two events

  // Current behavior with the mixed events:
  //              1) implemented ONLY IF ME1 and ME2 are the same species!!!!!!!!!!!
  //              2) the first of the two ME particles: takes Nth event from ME, takes every particle in this event
  //              2) the second of the two ME particles: takes (N+1)th event from ME, takes every particle in this event


  if(speciesME1!=speciesME2) {
    AliError("You chose different species ME1 and ME2 for mixing. This is not yet implemented! \n");
  }
  auto ParticleSE = ParticleVector.begin()+speciesSE;
  auto MixedEvent1Container = fPartContainer.begin()+speciesME1;
  auto MixedEvent2Container = fPartContainer.begin()+speciesME2;

  // Get the PID codes std::vector<int>
  auto itPDGParSE = PDGCodes.begin()+speciesSE;
  auto itPDGParME1 = PDGCodes.begin()+speciesME1;
  auto itPDGParME2 = PDGCodes.begin()+speciesME2;

  // Get particle masses
  auto massParticleSE = TDatabasePDG::Instance()->GetParticle(*itPDGParSE)->Mass();
  auto massParticleME1 = TDatabasePDG::Instance()->GetParticle(*itPDGParME1)->Mass();
  auto massParticleME2 = TDatabasePDG::Instance()->GetParticle(*itPDGParME2)->Mass();

  unsigned int DaughterPart1 = 0;
  unsigned int DaughterPart2 = 0;
  unsigned int DaughterPart3 = 0;

  if(abs(*itPDGParSE)==2212) DaughterPart1 = 1;
  if(abs(*itPDGParSE)==211) DaughterPart1 = 1;
  if(abs(*itPDGParME1)==2212) DaughterPart2 = 1;
  if(abs(*itPDGParME1)==211) DaughterPart2 = 1;
  if(abs(*itPDGParME2)==2212) DaughterPart3 = 1;
  if(abs(*itPDGParME2)==211) DaughterPart3 = 1;
  unsigned int DoThisPair12 = DaughterPart1*10+DaughterPart2;
  unsigned int DoThisPair23 = DaughterPart2*10+DaughterPart3;
  unsigned int DoThisPair31 = DaughterPart3*10+DaughterPart1;

  // loop over first particle
  for (auto iPart1 = ParticleSE->begin(); iPart1 != ParticleSE->end(); ++iPart1) {
    // loop over second particle ...
    for (int iDepth1 = 0; iDepth1 < (int) MixedEvent1Container->GetMixingDepth(); ++iDepth1) {
      std::vector<AliFemtoDreamBasePart> iEvent2 = MixedEvent1Container->GetEvent(iDepth1);
      for ( auto iPart2 = iEvent2.begin(); iPart2 != iEvent2.end(); ++iPart2) {
        int iDepth2 = 0;
        if(speciesME1==speciesME2) iDepth2 = iDepth1+1;
        for ( ; iDepth2 < (int) MixedEvent2Container->GetMixingDepth(); ++iDepth2) {
          std::vector<AliFemtoDreamBasePart> iEvent3 = MixedEvent2Container->GetEvent(iDepth2);
          for ( auto iPart3 = iEvent3.begin(); iPart3 != iEvent3.end(); ++iPart3) {


            TLorentzVector part1_LorVec, part2_LorVec, part3_LorVec;
            part1_LorVec.SetPxPyPzE(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(),
            iPart1->GetMomentum().Z(), sqrt(pow(iPart1->GetP(),2)+pow(massParticleSE,2)));
            part2_LorVec.SetPxPyPzE(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(),
            iPart2->GetMomentum().Z(), sqrt(pow(iPart2->GetP(),2)+pow(massParticleME1,2)));
            part3_LorVec.SetPxPyPzE(iPart3->GetMomentum().X(), iPart3->GetMomentum().Y(),
            iPart3->GetMomentum().Z(), sqrt(pow(iPart3->GetP(),2)+pow(massParticleME2,2)));
            // Now when we have the lorentz vectors, we can calculate the Lorentz invariant relative momenta q12, q23, q31
            TLorentzVector q12 = AliAnalysisTaskAODThreeBodyProtonPrimary::RelativePairMomentum(part1_LorVec,part2_LorVec);
            TLorentzVector q23 = AliAnalysisTaskAODThreeBodyProtonPrimary::RelativePairMomentum(part2_LorVec,part3_LorVec);
            TLorentzVector q31 = AliAnalysisTaskAODThreeBodyProtonPrimary::RelativePairMomentum(part3_LorVec,part1_LorVec);

            float RelativeMomentum12 = AliFemtoDreamHigherPairMath::RelativePairMomentum(part1_LorVec,part2_LorVec);
            float RelativeMomentum23 = AliFemtoDreamHigherPairMath::RelativePairMomentum(part2_LorVec,part3_LorVec);
            float RelativeMomentum31 = AliFemtoDreamHigherPairMath::RelativePairMomentum(part3_LorVec,part1_LorVec);

            // The particles in current methodology are put in bins of:
            //                 Q3=sqrt(q12^2+q23^2+q31^2)
            float Q32 = q12*q12+q23*q23+q31*q31;
            // From 3 pion paper, the q must be multiplied by -1 before taking quare root
            float Q3 = sqrt(-Q32); // the minus from pion paper


            bool Pair12 = true;
            bool Pair23 = true;
            bool Pair31 = true;


            if(!fturnoffClosePairRejectionCompletely){

              if(fClosePairRejectionForAll){
                if(abs(*itPDGParSE)==abs(*itPDGParME1)){
                  Pair12 = DeltaEtaDeltaPhi(speciesSE, speciesME1, *iPart1,*iPart2, *itPDGParSE, *itPDGParME1, true,  DoThisPair12, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[6+phiEtaHistNo],Config, Q3);
                }else{
                  Pair12 = DeltaEtaDeltaPhi(speciesSE, speciesME1, *iPart1,*iPart2, *itPDGParSE, *itPDGParME1, true,  DoThisPair12, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[6+phiEtaHistNo],Config, Q3);
                }
                if(abs(*itPDGParME1)==abs(*itPDGParME2)){
                  Pair23 = DeltaEtaDeltaPhi(speciesME1, speciesME2, *iPart2,*iPart3, *itPDGParME1, *itPDGParME2, true,  DoThisPair23, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[6+phiEtaHistNo],Config, Q3);
                }else{
                  Pair23 = DeltaEtaDeltaPhi(speciesME1, speciesME2, *iPart2,*iPart3, *itPDGParME1, *itPDGParME2, true,  DoThisPair23, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[6+phiEtaHistNo],Config, Q3);
                }
                if(abs(*itPDGParME2)==abs(*itPDGParSE)){
                  Pair31 = DeltaEtaDeltaPhi(speciesME2, speciesSE, *iPart3,*iPart1, *itPDGParME2, *itPDGParSE, true,  DoThisPair31, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[6+phiEtaHistNo],Config, Q3);
                }else{
                  Pair31 = DeltaEtaDeltaPhi(speciesME2, speciesSE, *iPart3,*iPart1, *itPDGParME2, *itPDGParSE, true,  DoThisPair31, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[6+phiEtaHistNo],Config, Q3);
                }
              }
              if(!fClosePairRejectionForAll){
               
                  if(DoThisPair12==11){
                    if(abs(*itPDGParSE)==abs(*itPDGParME1)){
                      Pair12 = DeltaEtaDeltaPhi(speciesSE, speciesME1, *iPart1,*iPart2, *itPDGParSE, *itPDGParME1, true,  DoThisPair12, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[6+phiEtaHistNo],Config, Q3);
                    }else{
                      Pair12 = DeltaEtaDeltaPhi(speciesSE, speciesME1, *iPart1,*iPart2, *itPDGParSE, *itPDGParME1, true,  DoThisPair12, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[6+phiEtaHistNo],Config, Q3);
                    }
                  }
                  if(DoThisPair23==11){
                    if(abs(*itPDGParME1)==abs(*itPDGParME2)){
                      Pair23 = DeltaEtaDeltaPhi(speciesME1, speciesME2, *iPart2,*iPart3, *itPDGParME1, *itPDGParME2, true,  DoThisPair23, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[6+phiEtaHistNo],Config, Q3);
                    }else{
                      Pair23 = DeltaEtaDeltaPhi(speciesME1, speciesME2, *iPart2,*iPart3, *itPDGParME1, *itPDGParME2, true,  DoThisPair23, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[6+phiEtaHistNo],Config, Q3);
                    }
                  }
                  if(DoThisPair31==11){
                    if(abs(*itPDGParME2)==abs(*itPDGParSE)){
                      Pair31 = DeltaEtaDeltaPhi(speciesME2, speciesSE, *iPart3,*iPart1, *itPDGParME2, *itPDGParSE, true,  DoThisPair31, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[6+phiEtaHistNo],Config, Q3);
                    }else{
                      Pair31 = DeltaEtaDeltaPhi(speciesME2, speciesSE, *iPart3,*iPart1, *itPDGParME2, *itPDGParSE, true,  DoThisPair31, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[6+phiEtaHistNo],Config, Q3);
                    }
                  }  

              }

            }

            if(!Pair12||!Pair23||!Pair31) {continue;}


            // Now we have the three particles, lets create their Lorentz vectors

            hist->Fill(Q3);
            hist2d->Fill(Q3,mult+1);
    
            if(fRunPairMultThreeBody && fQ3MinValue <= Q3 && Q3<fQ3cutValue){
              hist2d12->Fill(RelativeMomentum12,mult+1);
              hist2d23->Fill(RelativeMomentum23,mult+1);
              hist2d31->Fill(RelativeMomentum31,mult+1);
            }

           if(fRunmTPlots){
             float mT12 = GetmT(part1_LorVec, massParticleSE, part2_LorVec, massParticleME1);
             float mT23 = GetmT(part2_LorVec, massParticleME1, part3_LorVec, massParticleME2);
             float mT31 = GetmT(part3_LorVec, massParticleME2, part1_LorVec, massParticleSE);

             histmTQ312->Fill(Q3,mT12);
             histmTQ323->Fill(Q3,mT23);
             histmTQ331->Fill(Q3,mT31);
           }


        if(fRunPlotInvMassVSmT){
          TLorentzVector Sum12 = part1_LorVec + part2_LorVec;
          TLorentzVector Sum23 = part2_LorVec + part3_LorVec;
          TLorentzVector Sum31 = part3_LorVec + part1_LorVec;
          float mT12 = GetmT(part1_LorVec, massParticleSE, part2_LorVec, massParticleME1);
          float mT23 = GetmT(part2_LorVec, massParticleME1, part3_LorVec, massParticleME2);
          float mT31 = GetmT(part3_LorVec, massParticleME2, part1_LorVec, massParticleSE);
          InvMassVsmT12->Fill(Sum12.Mag(),mT12);
          InvMassVsmT23->Fill(Sum23.Mag(),mT23);
          InvMassVsmT31->Fill(Sum31.Mag(),mT31);

        }

        if(fRunProjector && phiEtaHistNo == 2){ //only run it for  (p)-(p)-Prim
          Projector->Fill(Q3, RelativeMomentum12); 
        }

        if(fIsMC && fGetMomentumResolution){
          MomentumResolution(ResAll, Res, *iPart1, *itPDGParSE, massParticleSE,*iPart2, *itPDGParME1, massParticleME1, *iPart3,*itPDGParME2, massParticleME2, Q3) ;
        }

          }
        }
      }
    }
  }
}


//==================================================================================================================================================

void AliAnalysisTaskAODThreeBodyProtonPrimary::FillPairDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME1, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F* histDeltaPhi, TH2F* histLowDeltaPhi, TH2F **fEventTripletPhiThetaArray, int phiEtaHistNo, AliFemtoDreamCollConfig Config){ 

  //  Check if reproduces the FemtoDream framework result
  auto ParticleSE = ParticleVector.begin()+speciesSE;
  auto MixedEvent1Container = fPartContainer.begin()+speciesME1;

  // Get the PID codes std::vector<int>
  auto itPDGParSE = PDGCodes.begin()+speciesSE;
  auto itPDGParME1 = PDGCodes.begin()+speciesME1;

  unsigned int DaughterPart1 = 0;
  unsigned int DaughterPart2 = 0;

  if(abs(*itPDGParSE)==2212) DaughterPart1 = 1;
  if(abs(*itPDGParSE)==211) DaughterPart1 = 1;
  if(abs(*itPDGParME1)==2212) DaughterPart2 = 1;
  if(abs(*itPDGParME1)==211) DaughterPart2 = 1;

  unsigned int DoThisPair12 = DaughterPart1*10+DaughterPart2;

  // Get particle masses
  auto massParticleSE = TDatabasePDG::Instance()->GetParticle(*itPDGParSE)->Mass();
  auto massParticleME1 = TDatabasePDG::Instance()->GetParticle(*itPDGParME1)->Mass();

  // loop over first particle
  for (auto iPart1 = ParticleSE->begin(); iPart1 != ParticleSE->end(); ++iPart1) {
    // loop over second particle ...
    for (int iDepth1 = 0; iDepth1 < (int) MixedEvent1Container->GetMixingDepth(); ++iDepth1) {
      std::vector<AliFemtoDreamBasePart> iEvent2 = MixedEvent1Container->GetEvent(iDepth1); //GANESHA Here check if two species are the same
      for ( auto iPart2 = iEvent2.begin(); iPart2 != iEvent2.end(); ++iPart2) {

        
        bool Pair12 = true;

         if(!fturnoffClosePairRejectionCompletely){

          if(fClosePairRejectionForAll){
             Pair12 = DeltaEtaDeltaPhi(speciesSE, speciesME1, *iPart1,*iPart2, *itPDGParSE, *itPDGParME1, false,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[6+phiEtaHistNo],Config); 
          }
          else{

             if(DoThisPair12==11){
               Pair12 = DeltaEtaDeltaPhi(speciesSE, speciesME1,*iPart1,*iPart2, *itPDGParSE, *itPDGParME1, false,  DoThisPair12, fEventTripletPhiThetaArray[phiEtaHistNo],fEventTripletPhiThetaArray[6+phiEtaHistNo],Config); 
             }
          } 
        }

        if(!Pair12) {continue;}

        // Now we have the three particles, lets create their Lorentz vectors
        TLorentzVector part1_LorVec, part2_LorVec;
        part1_LorVec.SetXYZM(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(),iPart1->GetMomentum().Z(),massParticleSE);
        part2_LorVec.SetXYZM(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(),iPart2->GetMomentum().Z(),massParticleME1);
        // Get momentum
        float RelativeK = AliFemtoDreamHigherPairMath::RelativePairMomentum(part1_LorVec, part2_LorVec);
        // No need to check pair selection because p lambda
        hist->Fill(RelativeK);
        hist2d->Fill(RelativeK,mult+1);


        double deltaPhi = (iPart2->GetPhi().at(0) - iPart1->GetPhi().at(0));

        histDeltaPhi->Fill(RelativeK,deltaPhi);
        if(abs(deltaPhi)>0.5&&abs(deltaPhi)<2.5){
          histLowDeltaPhi->Fill(RelativeK,mult+1);
        }



      }
    }
  }
}

//==================================================================================================================================================

void AliAnalysisTaskAODThreeBodyProtonPrimary::FillTripletDistributionSE2ME1(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> &fPartContainer, int speciesSE1, int speciesSE2, int speciesME, TH1F* hist, std::vector<int> PDGCodes, int mult, TH2F* hist2d, TH2F* hist2d12, TH2F **fEventTripletPhiThetaArray_SamePair, TH2F **fEventTripletPhiThetaArray_DifferentPair, int phiEtaHistNo, AliFemtoDreamCollConfig Config, TH2F* InvMass12, TH2F* InvMass23, TH2F* InvMass31, TH2F* histmTQ312, TH2F* histmTQ323, TH2F* histmTQ331, TH2F* InvMassVsmT12, TH2F* InvMassVsmT23, TH2F* InvMassVsmT31){
  // Description of function given in AliAnalysisTaskAODThreeBodyProtonPrimary::FillTripletDistribution
  // In this function, two particles are used from current event, and one - from other event
  auto ParticleSE1 = ParticleVector.begin()+speciesSE1;
  auto ParticleSE2 = ParticleVector.begin()+speciesSE2;
  auto MixedEventContainer = fPartContainer.begin()+speciesME;


  if(!fStandardMixing){
    if((speciesSE1==0&&speciesSE2==0&&speciesME==0)||(speciesSE1==1&&speciesSE2==1&&speciesME==1)){ // (p-p)-p
      if(ParticleVector[speciesSE1].size()<3) return;
    }
    if((speciesSE1==0&&speciesSE2==0&&speciesME==2)||(speciesSE1==1&&speciesSE2==1&&speciesME==3)){ // (p-p)-Prim
      if(ParticleVector[speciesSE1].size()<2||ParticleVector[speciesME].size()<1) return;
    }
    if((speciesSE1==0&&speciesSE2==2&&speciesME==0)||(speciesSE1==1&&speciesSE2==3&&speciesME==1)){ // (p-Prim)-p
      if(ParticleVector[speciesSE1].size()<2||ParticleVector[speciesSE2].size()<1) return;
    }
    if((speciesSE1==0&&speciesSE2==0&&speciesME==3)||(speciesSE1==1&&speciesSE2==1&&speciesME==2)){ // (p-p)-APrim
      if(ParticleVector[speciesSE1].size()<2||ParticleVector[speciesME].size()<1) return;
    }
    if((speciesSE1==0&&speciesSE2==3&&speciesME==0)||(speciesSE1==1&&speciesSE2==2&&speciesME==1)){ // (p-APrim)-p
      if(ParticleVector[speciesSE1].size()<2||ParticleVector[speciesSE2].size()<1) return;
    }
  }

  // Get the PID codes std::vector<int>
  auto itPDGParSE1 = PDGCodes.begin()+speciesSE1;
  auto itPDGParSE2 = PDGCodes.begin()+speciesSE2;
  auto itPDGParME = PDGCodes.begin()+speciesME;

  // Get particle masses
  auto massParticleSE1 = TDatabasePDG::Instance()->GetParticle(*itPDGParSE1)->Mass();
  auto massParticleSE2 = TDatabasePDG::Instance()->GetParticle(*itPDGParSE2)->Mass();
  auto massParticleME = TDatabasePDG::Instance()->GetParticle(*itPDGParME)->Mass();
  // get info for close pair rejection
  unsigned int DaughterPart1 = 0;
  unsigned int DaughterPart2 = 0;
  unsigned int DaughterPart3 = 0;

  if(abs(*itPDGParSE1)==2212) DaughterPart1 = 1;
  if(abs(*itPDGParSE1)==211) DaughterPart1 = 1;
  if(abs(*itPDGParSE2)==2212) DaughterPart2 = 1;
  if(abs(*itPDGParSE2)==211) DaughterPart2 = 1;
  if(abs(*itPDGParME)==2212) DaughterPart3 = 1;
  if(abs(*itPDGParME)==211) DaughterPart3 = 1;
  unsigned int DoThisPair12 = DaughterPart1*10+DaughterPart2;
  unsigned int DoThisPair23 = DaughterPart2*10+DaughterPart3;
  unsigned int DoThisPair31 = DaughterPart3*10+DaughterPart1;
  // loop over first particle
  for (auto iPart1 = ParticleSE1->begin(); iPart1 != ParticleSE1->end(); ++iPart1) {
    auto iPart2 = ParticleSE2->begin();
    if(speciesSE1==speciesSE2) iPart2 = iPart1+1;
    for(;iPart2 != ParticleSE2->end(); ++iPart2)
    {
      for ( int iDepth1 = 0; iDepth1 < (int) MixedEventContainer->GetMixingDepth(); ++iDepth1) {
        std::vector<AliFemtoDreamBasePart> iEvent3 = MixedEventContainer->GetEvent(iDepth1);
        for ( auto iPart3 = iEvent3.begin(); iPart3 != iEvent3.end(); ++iPart3) {
          // Now we have the three particles, lets create their Lorentz vectors


          TLorentzVector part1_LorVec, part2_LorVec, part3_LorVec;
          part1_LorVec.SetPxPyPzE(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(),
          iPart1->GetMomentum().Z(), sqrt(pow(iPart1->GetP(),2)+pow(massParticleSE1,2)));
          part2_LorVec.SetPxPyPzE(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(),
          iPart2->GetMomentum().Z(), sqrt(pow(iPart2->GetP(),2)+pow(massParticleSE2,2)));
          part3_LorVec.SetPxPyPzE(iPart3->GetMomentum().X(), iPart3->GetMomentum().Y(),
          iPart3->GetMomentum().Z(), sqrt(pow(iPart3->GetP(),2)+pow(massParticleME,2)));
          // Now when we have the lorentz vectors, we can calculate the Lorentz invariant relative momenta q12, q23, q31
          TLorentzVector q12 = AliAnalysisTaskAODThreeBodyProtonPrimary::RelativePairMomentum(part1_LorVec,part2_LorVec);
          TLorentzVector q23 = AliAnalysisTaskAODThreeBodyProtonPrimary::RelativePairMomentum(part2_LorVec,part3_LorVec);
          TLorentzVector q31 = AliAnalysisTaskAODThreeBodyProtonPrimary::RelativePairMomentum(part3_LorVec,part1_LorVec);
          // The particles in current methodology are put in bins of:
          //                 Q3=sqrt(q12^2+q23^2+q31^2)
          float Q32 = q12*q12+q23*q23+q31*q31;
          // From 3 pion paper, the q must be multiplied by -1 before taking quare root
          float Q3 = sqrt(-Q32); // the minus from pion paper


          bool Pair12 = true;
          bool Pair23 = true;
          bool Pair31 = true;
          if(!fturnoffClosePairRejectionCompletely){

            if(fClosePairRejectionForAll){
              if(abs(*itPDGParSE1)==abs(*itPDGParSE2)){
                Pair12 = DeltaEtaDeltaPhi(speciesSE1, speciesSE2, *iPart1,*iPart2, *itPDGParSE1, *itPDGParSE2, true,  DoThisPair12, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[16+phiEtaHistNo],Config, Q3);
              }else{
                Pair12 = DeltaEtaDeltaPhi(speciesSE1, speciesSE2, *iPart1,*iPart2, *itPDGParSE1, *itPDGParSE2, true,  DoThisPair12, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[16+phiEtaHistNo],Config, Q3);
              }
              if(abs(*itPDGParSE2)==abs(*itPDGParME)){
                Pair23 = DeltaEtaDeltaPhi(speciesSE2, speciesME, *iPart2,*iPart3, *itPDGParSE2, *itPDGParME, true,  DoThisPair23, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[16+phiEtaHistNo],Config, Q3);
              }else{
                Pair23 = DeltaEtaDeltaPhi(speciesSE2, speciesME, *iPart2,*iPart3, *itPDGParSE2, *itPDGParME, true,  DoThisPair23, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[16+phiEtaHistNo],Config, Q3);
              }
              if(abs(*itPDGParSE1)==abs(*itPDGParME)){
                Pair31 = DeltaEtaDeltaPhi(speciesME, speciesSE1, *iPart3,*iPart1, *itPDGParME, *itPDGParSE1, true,  DoThisPair31, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[16+phiEtaHistNo],Config, Q3);
              }else{
                Pair31 = DeltaEtaDeltaPhi(speciesME, speciesSE1, *iPart3,*iPart1, *itPDGParME, *itPDGParSE1, true,  DoThisPair31, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[16+phiEtaHistNo],Config, Q3);
              }
            }
            if(!fClosePairRejectionForAll){
          
                if(DoThisPair12==11){
                  if(abs(*itPDGParSE1)==abs(*itPDGParSE2)){
                    Pair12 = DeltaEtaDeltaPhi(speciesSE1, speciesSE2, *iPart1,*iPart2, *itPDGParSE1, *itPDGParSE2, true,  DoThisPair12, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[16+phiEtaHistNo],Config, Q3);
                  }else{
                    Pair12 = DeltaEtaDeltaPhi(speciesSE1, speciesSE2, *iPart1,*iPart2, *itPDGParSE1, *itPDGParSE2, true,  DoThisPair12, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[16+phiEtaHistNo],Config, Q3);
                  }
                }
                if(DoThisPair23==11){
                  if(abs(*itPDGParSE2)==abs(*itPDGParME)){
                    Pair23 = DeltaEtaDeltaPhi(speciesSE2, speciesME, *iPart2,*iPart3, *itPDGParSE2, *itPDGParME, true,  DoThisPair23, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[16+phiEtaHistNo],Config, Q3);
                  }else{
                    Pair23 = DeltaEtaDeltaPhi(speciesSE2, speciesME, *iPart2,*iPart3, *itPDGParSE2, *itPDGParME, true,  DoThisPair23, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[16+phiEtaHistNo],Config, Q3);
                  }
                }
                if(DoThisPair31==11){
                  if(abs(*itPDGParSE1)==abs(*itPDGParME)){
                    Pair31 = DeltaEtaDeltaPhi(speciesME, speciesSE1, *iPart3,*iPart1, *itPDGParME, *itPDGParSE1, true,  DoThisPair31, fEventTripletPhiThetaArray_SamePair[phiEtaHistNo],fEventTripletPhiThetaArray_SamePair[16+phiEtaHistNo],Config, Q3);
                  }else{
                    Pair31 = DeltaEtaDeltaPhi(speciesME, speciesSE1, *iPart3,*iPart1, *itPDGParME, *itPDGParSE1, true,  DoThisPair31, fEventTripletPhiThetaArray_DifferentPair[phiEtaHistNo],fEventTripletPhiThetaArray_DifferentPair[16+phiEtaHistNo],Config, Q3);
                  }
                }
              
            }

          }

          if(!Pair12||!Pair23||!Pair31) {continue;}

          if(fIsMC && fRemoveMCResonances){
            bool CommonMother12 = CommonMotherResonance(&(*iPart1), &(*iPart2)); 
            bool CommonMother13 = CommonMotherResonance(&(*iPart1), &(*iPart3)); //this should always be false by definition
            bool CommonMother23 = CommonMotherResonance(&(*iPart2), &(*iPart3)); //this should always be false by definition

            if(CommonMother12 || CommonMother13 || CommonMother23){
              continue; 
            }
          }


          hist->Fill(Q3);
          hist2d->Fill(Q3,mult+1);

        if(fRunPairMultThreeBody && fQ3MinValue <= Q3 && Q3<fQ3cutValue){
          float RelativeMomentum12 = AliFemtoDreamHigherPairMath::RelativePairMomentum(part1_LorVec,part2_LorVec);
          hist2d12->Fill(RelativeMomentum12,mult+1);
        }

        if(fRunmTPlots){
            float mT12 = GetmT(part1_LorVec, massParticleSE1, part2_LorVec, massParticleSE2);
            float mT23 = GetmT(part2_LorVec, massParticleSE2, part3_LorVec, massParticleME);
            float mT31 = GetmT(part3_LorVec, massParticleME, part1_LorVec, massParticleSE1);

            histmTQ312->Fill(Q3,mT12);
            histmTQ323->Fill(Q3,mT23);
            histmTQ331->Fill(Q3,mT31);
        }

        if(fRunPlotInvMass){
           TLorentzVector Sum12 = part1_LorVec + part2_LorVec;
           TLorentzVector Sum23 = part2_LorVec + part3_LorVec;
           TLorentzVector Sum31 = part3_LorVec + part1_LorVec;
           InvMass12->Fill(Sum12.Mag(),Q3);
           InvMass23->Fill(Sum23.Mag(),Q3);
           InvMass31->Fill(Sum31.Mag(),Q3);
        }//if(fRunPlotInvMass)

        if(fRunPlotInvMassVSmT){
          TLorentzVector Sum12 = part1_LorVec + part2_LorVec;
          TLorentzVector Sum23 = part2_LorVec + part3_LorVec;
          TLorentzVector Sum31 = part3_LorVec + part1_LorVec;
          float mT12 = GetmT(part1_LorVec, massParticleSE1, part2_LorVec, massParticleSE2);
          float mT23 = GetmT(part2_LorVec, massParticleSE2, part3_LorVec, massParticleME);
          float mT31 = GetmT(part3_LorVec, massParticleME, part1_LorVec, massParticleSE1);
          InvMassVsmT12->Fill(Sum12.Mag(),mT12);
          InvMassVsmT23->Fill(Sum23.Mag(),mT23);
          InvMassVsmT31->Fill(Sum31.Mag(),mT31);
        }

        }
      }
    }
  }
}//AliAnalysisTaskAODThreeBodyProtonPrimary::FillTripletDistributionSE2ME1

//==================================================================================================================================================

bool AliAnalysisTaskAODThreeBodyProtonPrimary::DeltaEtaDeltaPhi(int species1, int species2,
                                                   AliFemtoDreamBasePart &part1,
                                                   AliFemtoDreamBasePart &part2,
                                                   int part1PDGcode,
                                                   int part2PDGcode,
                                                   bool SEorME,  unsigned int DoThisPair, TH2F* beforeHist,TH2F* afterHist,
                                                   AliFemtoDreamCollConfig Config) {
  // DoThisPair = ij where i is the number of daughters for first particle, j for the second

  static const float piHi = TMath::Pi();
  //auto fDeltaPhiSqMax = Config.GetDeltaPhiMax() * Config.GetDeltaPhiMax();
  //auto fDeltaEtaSqMax = Config.GetDeltaEtaMax() * Config.GetDeltaEtaMax() ;

  double DeltaPhiSqMaxValue = 0;
  double DeltaEtaSqMaxValue = 0;


  //cout<<part1.GetPDGCode()<<endl;

  
  if ((species1==2&&species2==0)||(species1==0&&species2==2)||
      (species1==3&&species2==1)||(species1==1&&species2==3))
  {
    //fDeltaPhiSqMax = 0.04*0.04;
    //fDeltaEtaSqMax = 0.012*0.012;
    DeltaPhiSqMaxValue = fDeltaPhiMaxPPrim*fDeltaPhiMaxPPrim;
    DeltaEtaSqMaxValue = fDeltaEtaMaxPPrim*fDeltaEtaMaxPPrim;

  }
  else if ((species1==0&&species2==0)||(species1==1&&species2==1))
  {
    //fDeltaPhiSqMax = 0.017*0.017;
    //fDeltaEtaSqMax = 0.017*0.017;
    DeltaPhiSqMaxValue = fDeltaPhiMaxPP*fDeltaPhiMaxPP;
    DeltaEtaSqMaxValue = fDeltaEtaMaxPP*fDeltaEtaMaxPP;

  }
  else if ((species1==3&&species2==0)||(species1==0&&species2==3)||
      (species1==2&&species2==1)||(species1==1&&species2==2))
  {
    //fDeltaPhiSqMax = 0.04*0.04;
    //fDeltaEtaSqMax = 0.012*0.012;
    DeltaPhiSqMaxValue = fDeltaPhiMaxPAPrim*fDeltaPhiMaxPAPrim;
    DeltaEtaSqMaxValue = fDeltaEtaMaxPAPrim*fDeltaEtaMaxPAPrim;

  }

  bool pass = true;
  // if nDaug == 1 => Single Track, else decay
  unsigned int nDaug1 = (unsigned int) DoThisPair / 10;
  if (nDaug1 > 9) {
    AliWarning("you are doing something wrong \n");
  }
  if (nDaug1 > part1.GetPhiAtRaidius().size()) {
    TString outMessage =
        TString::Format(
            "For pair number %u your number of Daughters 1 (%u) and Radii 1 (%u) do not correspond \n",
            DoThisPair, nDaug1, (unsigned int)part1.GetPhiAtRaidius().size());
    AliWarning(outMessage.Data());
  }
  unsigned int nDaug2 = (unsigned int) DoThisPair % 10;

  if (nDaug2 > part2.GetPhiAtRaidius().size()) {
    TString outMessage =
        TString::Format(
            "For pair number %u your number of Daughters 2 (%u) and Radii 2 (%u) do not correspond \n",
            DoThisPair, nDaug2, (unsigned int)part2.GetPhiAtRaidius().size());
    AliWarning(outMessage.Data());
  }
  std::vector<float> eta1 = part1.GetEta();
  std::vector<float> eta2 = part2.GetEta();

  for (unsigned int iDaug1 = 0; iDaug1 < nDaug1; ++iDaug1) {
    std::vector<float> PhiAtRad1 = part1.GetPhiAtRaidius().at(iDaug1);
    float etaPar1;
    if (nDaug1 == 1) {
      etaPar1 = eta1.at(0);
    } else {
      etaPar1 = eta1.at(iDaug1 + 1);
    }
    for (unsigned int iDaug2 = 0; iDaug2 < nDaug2; ++iDaug2) {
      std::vector<float> phiAtRad2 = part2.GetPhiAtRaidius().at(iDaug2);
      float etaPar2;
      if (nDaug2 == 1) {
        etaPar2 = eta2.at(0);
      } else {
        etaPar2 = eta2.at(iDaug2 + 1);
      }
      float deta = etaPar1 - etaPar2;
      const int size =
          (PhiAtRad1.size() > phiAtRad2.size()) ?
              phiAtRad2.size() : PhiAtRad1.size();
      float dphiAvg = 0;
      for (int iRad = 0; iRad < size; ++iRad) {
        float dphi = PhiAtRad1.at(iRad) - phiAtRad2.at(iRad);
        if (dphi > piHi) {
          dphi += -piHi * 2;
        } else if (dphi < -piHi) {
          dphi += piHi * 2;
        }
        dphi = TVector2::Phi_mpi_pi(dphi);

        dphiAvg += dphi;
      }
      if(fRunPlotPhiTheta){
        beforeHist->Fill(dphiAvg/ (float) size, deta);
      }
      if (pass) {
        if ((dphiAvg / (float) size) * (dphiAvg / (float) size) / DeltaPhiSqMaxValue
            + deta * deta / DeltaEtaSqMaxValue < 1.) {
          pass = false;
        }
        else{
          if(fRunPlotPhiTheta){
            afterHist->Fill(dphiAvg/ (float) size, deta);
          }
        }
      }
    }
  }
  return pass;
}

//==================================================================================================================================================

bool AliAnalysisTaskAODThreeBodyProtonPrimary::DeltaEtaDeltaPhi(int species1, int species2,
                                                   AliFemtoDreamBasePart &part1,
                                                   AliFemtoDreamBasePart &part2,
                                                   int part1PDGcode,
                                                   int part2PDGcode,
                                                   bool SEorME,  unsigned int DoThisPair, TH2F* beforeHist,TH2F* afterHist,
                                                   AliFemtoDreamCollConfig Config, double Q3) {
  // DoThisPair = ij where i is the number of daughters for first particle, j for the second

  static const float piHi = TMath::Pi();
  //auto fDeltaPhiSqMax = Config.GetDeltaPhiMax() * Config.GetDeltaPhiMax();
  //auto fDeltaEtaSqMax = Config.GetDeltaEtaMax() * Config.GetDeltaEtaMax() ;

  double DeltaPhiSqMaxValue = 0;
  double DeltaEtaSqMaxValue = 0;

  //cout<<part1PDGcode<<endl;

 //bool test = false;

 if ((species1==2&&species2==0)||(species1==0&&species2==2)||
      (species1==3&&species2==1)||(species1==1&&species2==3))
  {
    //fDeltaPhiSqMax = 0.04*0.04;
    //fDeltaEtaSqMax = 0.012*0.012;

    DeltaPhiSqMaxValue = fDeltaPhiMaxPPrim*fDeltaPhiMaxPPrim;
    DeltaEtaSqMaxValue = fDeltaEtaMaxPPrim*fDeltaEtaMaxPPrim;
  //  test = true;
  }
  else if ((species1==0&&species2==0)||(species1==1&&species2==1))
  {
    //fDeltaPhiSqMax = 0.017*0.017;
    //fDeltaEtaSqMax = 0.017*0.017;

    DeltaPhiSqMaxValue = fDeltaPhiMaxPP*fDeltaPhiMaxPP;
    DeltaEtaSqMaxValue = fDeltaEtaMaxPP*fDeltaEtaMaxPP;

  //  test = true;
  }
  else if ((species1==3&&species2==0)||(species1==0&&species2==3)||
      (species1==2&&species2==1)||(species1==1&&species2==2))
  {
    //fDeltaPhiSqMax = 0.04*0.04;
    //fDeltaEtaSqMax = 0.012*0.012;
    DeltaPhiSqMaxValue = fDeltaPhiMaxPAPrim*fDeltaPhiMaxPAPrim;
    DeltaEtaSqMaxValue = fDeltaEtaMaxPAPrim*fDeltaEtaMaxPAPrim;

  }

  //cout<<species1<<"   "<<species2<<"  "<<test<<endl;

  bool pass = true;
  // if nDaug == 1 => Single Track, else decay
  unsigned int nDaug1 = (unsigned int) DoThisPair / 10;
  if (nDaug1 > 9) {
    AliWarning("you are doing something wrong \n");
  }
  if (nDaug1 > part1.GetPhiAtRaidius().size()) {
    TString outMessage =
        TString::Format(
            "For pair number %u your number of Daughters 1 (%u) and Radii 1 (%u) do not correspond \n",
            DoThisPair, nDaug1, (unsigned int)part1.GetPhiAtRaidius().size());
    AliWarning(outMessage.Data());
  }
  unsigned int nDaug2 = (unsigned int) DoThisPair % 10;

  if (nDaug2 > part2.GetPhiAtRaidius().size()) {
    TString outMessage =
        TString::Format(
            "For pair number %u your number of Daughters 2 (%u) and Radii 2 (%u) do not correspond \n",
            DoThisPair, nDaug2, (unsigned int)part2.GetPhiAtRaidius().size());
    AliWarning(outMessage.Data());
  }
  std::vector<float> eta1 = part1.GetEta();
  std::vector<float> eta2 = part2.GetEta();

  for (unsigned int iDaug1 = 0; iDaug1 < nDaug1; ++iDaug1) {
    std::vector<float> PhiAtRad1 = part1.GetPhiAtRaidius().at(iDaug1);
    float etaPar1;
    if (nDaug1 == 1) {
      etaPar1 = eta1.at(0);
    } else {
      etaPar1 = eta1.at(iDaug1 + 1);
    }
    for (unsigned int iDaug2 = 0; iDaug2 < nDaug2; ++iDaug2) {
      std::vector<float> phiAtRad2 = part2.GetPhiAtRaidius().at(iDaug2);
      float etaPar2;
      if (nDaug2 == 1) {
        etaPar2 = eta2.at(0);
      } else {
        etaPar2 = eta2.at(iDaug2 + 1);
      }
      float deta = etaPar1 - etaPar2;
      const int size =
          (PhiAtRad1.size() > phiAtRad2.size()) ?
              phiAtRad2.size() : PhiAtRad1.size();
      float dphiAvg = 0;
      for (int iRad = 0; iRad < size; ++iRad) {
        float dphi = PhiAtRad1.at(iRad) - phiAtRad2.at(iRad);
        if (dphi > piHi) {
          dphi += -piHi * 2;
        } else if (dphi < -piHi) {
          dphi += piHi * 2;
        }
        dphi = TVector2::Phi_mpi_pi(dphi);

        dphiAvg += dphi;
      }
      if(fRunPlotPhiTheta){
        if(Q3<fQ3LimitForDeltaPhiDeltaEta){
          beforeHist->Fill(dphiAvg/ (float) size, deta);
        }
      }
      if (pass) {
        if ((dphiAvg / (float) size) * (dphiAvg / (float) size) / DeltaPhiSqMaxValue
            + deta * deta / DeltaEtaSqMaxValue < 1.) {
          pass = false;
        }
        else{
          if(fRunPlotPhiTheta){
            if(Q3<fQ3LimitForDeltaPhiDeltaEta){
              afterHist->Fill(dphiAvg/ (float) size, deta);
            }
          }
        }
      }
    }
  }
  return pass;
}

//==================================================================================================================================================

void AliAnalysisTaskAODThreeBodyProtonPrimary::FillPairInvMass( AliFemtoDreamBasePart &part1,
                        AliFemtoDreamBasePart &part2,
                        AliFemtoDreamBasePart &part3,
                        TH2F* hist, float Q3) {
  TVector3 momPart1 = part1.GetMomentum();
  TVector3 momPart2 = part2.GetMomentum();
  TVector3 momPart3 = part3.GetMomentum();
  TLorentzVector track1, track2, track3;
  track1.SetXYZM(momPart1.Px(), momPart1.Py(), momPart1.Pz(),
                   part1.GetInvMass());
  track2.SetXYZM(momPart2.Px(), momPart2.Py(), momPart2.Pz(),
                   part2.GetInvMass());
  track3.SetXYZM(momPart3.Px(), momPart3.Py(), momPart3.Pz(),
                   part3.GetInvMass());
  TLorentzVector trackSum = track1 + track2 + track3;

  hist->Fill(Q3, trackSum.M());

}

//==================================================================================================================================================

void AliAnalysisTaskAODThreeBodyProtonPrimary::FillPDGPairInvMass( AliFemtoDreamBasePart &part1, float massPart1,
                        AliFemtoDreamBasePart &part2, float massPart2,
                        AliFemtoDreamBasePart &part3, float massPart3,
                        TH2F* hist, float Q3) {

  TVector3 momPart1 = part1.GetMomentum();
  TVector3 momPart2 = part2.GetMomentum();
  TVector3 momPart3 = part3.GetMomentum();
  TLorentzVector track1, track2, track3;
  track1.SetXYZM(momPart1.Px(), momPart1.Py(), momPart1.Pz(),
                   massPart1);
  track2.SetXYZM(momPart2.Px(), momPart2.Py(), momPart2.Pz(),
                   massPart2);
  track3.SetXYZM(momPart3.Px(), momPart3.Py(), momPart3.Pz(),
                   massPart3);
  TLorentzVector trackSum = track1 + track2 + track3;

  hist->Fill(Q3, trackSum.M());

}

void AliAnalysisTaskAODThreeBodyProtonPrimary::FillPairTransverseMass(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector,
                        int firstSpecies, int secondSpecies, TH1F* hist1, std::vector<int> PDGCodes, TH2F* hist2) {  

  // Two Body Same Event

  auto Particle1Vector = ParticleVector.begin()+firstSpecies;
  auto Particle2Vector = ParticleVector.begin()+secondSpecies;

  // Get the PID codes std::vector<int>
  auto itPDGPar1 = PDGCodes.begin()+firstSpecies;
  auto itPDGPar2 = PDGCodes.begin()+secondSpecies;

  // Get particle masses
  auto massPart1 = TDatabasePDG::Instance()->GetParticle(*itPDGPar1)->Mass();
  auto massPart2 = TDatabasePDG::Instance()->GetParticle(*itPDGPar2)->Mass();


//  TVector3 momPart1 = part1.GetMomentum();
//  TVector3 momPart2 = part2.GetMomentum();
//  TLorentzVector track1, track2;


  unsigned int DaughterPart1 = 0;
  unsigned int DaughterPart2 = 0;

  if(abs(*itPDGPar1)==2212) DaughterPart1 = 1;
  if(abs(*itPDGPar1)==211) DaughterPart1 = 1;
  if(abs(*itPDGPar2)==2212) DaughterPart2 = 1;
  if(abs(*itPDGPar2)==211) DaughterPart2 = 1;

  // loop over first particle
  for (auto iPart1 = Particle1Vector->begin(); iPart1 != Particle1Vector->end(); ++iPart1) {
    // if second particle species is different than first - start with the first particle in the vector
    auto iPart2 = Particle2Vector->begin();
    // if second particle  and first are the species, start second loop from the next particle (to not double count)
    if (firstSpecies==secondSpecies) iPart2 = iPart1+1;
    // loop over second particle ...
    for (; iPart2 != Particle2Vector->end(); ++iPart2) {

        // Now we have the three particles, lets create their Lorentz vectors
        TLorentzVector Particle1_LV, Particle2_LV;
        Particle1_LV.SetXYZM(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(),iPart1->GetMomentum().Z(), massPart1);
        Particle2_LV.SetXYZM(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(),iPart2->GetMomentum().Z(), massPart2);
        // Get momentum
        float RelativeMomentum = AliFemtoDreamHigherPairMath::RelativePairMomentum(Particle1_LV, Particle2_LV);

        TLorentzVector trackSum_LV = Particle1_LV + Particle2_LV;

        Double_t pairMT = TMath::Sqrt(pow(0.5*trackSum_LV.Pt(),2.) + pow(0.5*(massPart1+massPart2),2.));

        hist2->Fill(RelativeMomentum, pairMT);
        hist1->Fill(pairMT);

    }
  }

}

float AliAnalysisTaskAODThreeBodyProtonPrimary::GetmT(TLorentzVector &Part1, float mass1, TLorentzVector &Part2, float mass2){
  float mT = 0.; 

  TLorentzVector Sum;
  Sum = Part1 + Part2;
  float kT = 0.5 * Sum.Pt();

  float averageMass = 0.5 * (mass1 + mass2);

  mT = TMath::Sqrt(pow(kT, 2.) + pow(averageMass, 2.));
  return mT;
}

bool AliAnalysisTaskAODThreeBodyProtonPrimary::CommonAncestors(AliFemtoDreamBasePart& part1, AliFemtoDreamBasePart& part2) {
    bool IsCommon = false;
    if(part1.GetMotherID() == part2.GetMotherID()){
      IsCommon = true;
    }else if(part1.GetMotherID() != part2.GetMotherID()){
      IsCommon = false;
    }
    return IsCommon;
}

bool AliAnalysisTaskAODThreeBodyProtonPrimary::CommonMotherResonance(AliFemtoDreamBasePart* part1, AliFemtoDreamBasePart* part2) {

  if(part1->GetMotherID() != part2->GetMotherID()) {
    return false; //MotherID is different -> no common resonance
  }

  if(part1->GetMotherPDG() != part2->GetMotherPDG()) { //the ID is the same, but the PDG different -> Two tracks from same hard scattering but different resonances.
    return false; 
  }

  bool HasCommonMotherResonance = true;
  HasCommonMotherResonance = IsResonance(part1->GetMotherPDG()); //the resonance is of the type that should be removed 
  return HasCommonMotherResonance; 
}

bool AliAnalysisTaskAODThreeBodyProtonPrimary::IsResonance(int PDG) {

  int ProtonAntiPion[33] = {2114, 12112, 1214, 22112, 32114, 1212, 32112, 2116, 12116, 12114, 42112, 21214, 31214, 11212, 9902114, 1216, 9902112, 9912112, 21212, 22114, 9912114, 2118, 11216, 9902116, 9922112, 9922114, 1218, 9901218, 99021110, 99121110, 99012112, 99021112, 3122};
  int ProtonPion[12] = {2224, 32224, 2222, 12224, 12222, 2226, 22222, 22224, 2228, 12226, 9902228, 99022212};

  // When the element is not found, std::find returns the end of the range
  if ( std::find(std::begin(ProtonAntiPion), std::end(ProtonAntiPion), abs(PDG)) != std::end(ProtonAntiPion) ) {
    return true;
  } else if ( std::find(std::begin(ProtonPion), std::end(ProtonPion), abs(PDG)) != std::end(ProtonPion) ) {
    return true;
  } else {
    return false;
  }
}

//Taken from https://github.com/alisw/AliPhysics/blob/master/PWGCF/FEMTOSCOPY/FemtoDream/AliAnalysisTaskThreeBodyFemto.cxx
void AliAnalysisTaskAODThreeBodyProtonPrimary::MomentumResolution( TH2F* histAll, TH2F* hist,
    AliFemtoDreamBasePart& part1, int PDGPart1, float mass1,
    AliFemtoDreamBasePart& part2, int PDGPart2, float mass2,
    AliFemtoDreamBasePart& part3, int PDGPart3, float mass3,
    float Q3Reconstructed) {

  TLorentzVector part1_LorVec, part2_LorVec, part3_LorVec;

  float mom1 = sqrt(part1.GetMCMomentum().X()*part1.GetMCMomentum().X()+part1.GetMCMomentum().Y()*part1.GetMCMomentum().Y()+part1.GetMCMomentum().Z()*part1.GetMCMomentum().Z());
  float mom2 = sqrt(part2.GetMCMomentum().X()*part2.GetMCMomentum().X()+part2.GetMCMomentum().Y()*part2.GetMCMomentum().Y()+part2.GetMCMomentum().Z()*part2.GetMCMomentum().Z());
  float mom3 = sqrt(part3.GetMCMomentum().X()*part3.GetMCMomentum().X()+part3.GetMCMomentum().Y()*part3.GetMCMomentum().Y()+part3.GetMCMomentum().Z()*part3.GetMCMomentum().Z());
  
  part1_LorVec.SetPxPyPzE(part1.GetMCMomentum().X(), part1.GetMCMomentum().Y(), 
  part1.GetMCMomentum().Z(), sqrt(pow(mom1,2)+pow(mass1,2)));
  part2_LorVec.SetPxPyPzE(part2.GetMCMomentum().X(), part2.GetMCMomentum().Y(), 
  part2.GetMCMomentum().Z(), sqrt(pow(mom2,2)+pow(mass2,2)));
  part3_LorVec.SetPxPyPzE(part3.GetMCMomentum().X(), part3.GetMCMomentum().Y(), 
  part3.GetMCMomentum().Z(), sqrt(pow(mom3,2)+pow(mass3,2)));

  TLorentzVector q12 = AliAnalysisTaskAODThreeBodyProtonPrimary::RelativePairMomentum(part1_LorVec,part2_LorVec);
  TLorentzVector q23 = AliAnalysisTaskAODThreeBodyProtonPrimary::RelativePairMomentum(part2_LorVec,part3_LorVec);
  TLorentzVector q31 = AliAnalysisTaskAODThreeBodyProtonPrimary::RelativePairMomentum(part3_LorVec,part1_LorVec);
  float Q32 = q12*q12+q23*q23+q31*q31;
  float Q3Real = sqrt(-Q32); 
          

  // check if particles are identified correctly

  histAll->Fill(Q3Real, Q3Reconstructed);
  if ((TMath::Abs(PDGPart1) == TMath::Abs(part1.GetMCPDGCode()))
      && ((TMath::Abs(PDGPart2) == TMath::Abs(part2.GetMCPDGCode())))
      && ((TMath::Abs(PDGPart3) == TMath::Abs(part3.GetMCPDGCode())))) {
      hist->Fill(Q3Real, Q3Reconstructed);
  }
}

