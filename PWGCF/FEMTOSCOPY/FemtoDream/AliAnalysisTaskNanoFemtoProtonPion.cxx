/*
 * AliAnalysisTaskNanoFemtoProtonPion.cxx
 *
 *  Created on: 11 Mar 2022
 *  Author: Lesch Marcel
 */

#include "AliAnalysisTaskNanoFemtoProtonPion.h"
#include "AliAnalysisTaskNanoPt.h"
#include "AliAnalysisManager.h"
#include "AliInputEventHandler.h"
#include "AliAODInputHandler.h"
#include "AliMCEvent.h"
#include "AliMultSelection.h"
#include "AliNanoAODTrack.h"
#include <sstream>

ClassImp(AliAnalysisTaskNanoFemtoProtonPion)
AliAnalysisTaskNanoFemtoProtonPion::AliAnalysisTaskNanoFemtoProtonPion()
  : AliAnalysisTaskSE(),
    fisLightWeight(false),
    fTrackBufferSize(),
    fIsMC(false),
    fDoPairCleaning(false),
    fCombinationInput(""),
    fNameTagInput(""),
    fDoOfficialFemto(false),
    fDoThreeDFemto(false),
    fRunPlotMult(false),
    fRunPlotPhiTheta(false),
    fDoClosePairRejection(false),
    fEvent(nullptr),
    fTrack(nullptr),
    fEventCuts(nullptr),
    fTrackCutsPion(nullptr),
    fTrackCutsAntiPion(nullptr),
    fTrackCutsProton(nullptr),
    fTrackCutsAntiProton(nullptr),
    fConfig(nullptr),
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fGTI(nullptr),
    fEvtList(nullptr),
    fProtonList(nullptr),
    fProtonMCList(nullptr),
    fAntiProtonList(nullptr),
    fAntiProtonMCList(nullptr),
    fPionList(nullptr),
    fPionMCList(nullptr),
    fAntiPionList(nullptr),
    fAntiPionMCList(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr),
    fResultsThreeDFemto(nullptr),
    fPartContainer(0),
    fSameEvent_List_OneDimensional(nullptr),
    fSameEventDeltaEtaDeltaPhi_List_OneDimensional(nullptr),
    fSameEvent_OneDimensional(nullptr),
    fSameEventMult_OneDimensional(nullptr),
    fSameEventPhiTheta_OneDimensional(nullptr),
    fMixedEvent_List_OneDimensional(nullptr),
    fMixedEventDeltaEtaDeltaPhi_List_OneDimensional(nullptr),
    fMixedEvent_OneDimensional(nullptr),
    fMixedEventMult_OneDimensional(nullptr),
    fMixedEventTripletPhiTheta_OneDimensional(nullptr){
}

AliAnalysisTaskNanoFemtoProtonPion::AliAnalysisTaskNanoFemtoProtonPion(
  const char *name, bool isMC)
  : AliAnalysisTaskSE(name),
    fisLightWeight(false),
    fTrackBufferSize(2000),
    fIsMC(isMC),
    fDoPairCleaning(false),
    fCombinationInput(""),
    fNameTagInput(""),
    fDoOfficialFemto(false),
    fDoThreeDFemto(false),
    fRunPlotMult(false),
    fRunPlotPhiTheta(false),
    fDoClosePairRejection(true),
    fEvent(nullptr),
    fTrack(nullptr),
    fEventCuts(nullptr),
    fTrackCutsPion(nullptr),
    fTrackCutsAntiPion(nullptr),
    fTrackCutsProton(nullptr),
    fTrackCutsAntiProton(nullptr),
    fConfig(nullptr),
    fPairCleaner(nullptr),
    fPartColl(nullptr),
    fGTI(nullptr),
    fEvtList(nullptr),
    fProtonList(nullptr),
    fProtonMCList(nullptr),
    fAntiProtonList(nullptr),
    fAntiProtonMCList(nullptr),
    fPionList(nullptr),
    fPionMCList(nullptr),
    fAntiPionList(nullptr),
    fAntiPionMCList(nullptr),
    fResults(nullptr),
    fResultsQA(nullptr),
    fResultsThreeDFemto(nullptr),
    fPartContainer(0),
    fSameEvent_List_OneDimensional(nullptr),
    fSameEventDeltaEtaDeltaPhi_List_OneDimensional(nullptr),
    fSameEvent_OneDimensional(nullptr),
    fSameEventMult_OneDimensional(nullptr),
    fSameEventPhiTheta_OneDimensional(nullptr),
    fMixedEvent_List_OneDimensional(nullptr),
    fMixedEventDeltaEtaDeltaPhi_List_OneDimensional(nullptr),
    fMixedEvent_OneDimensional(nullptr),
    fMixedEventMult_OneDimensional(nullptr),
    fMixedEventTripletPhiTheta_OneDimensional(nullptr){
  DefineOutput(1, TList::Class());  //Output for the Event Cuts
  DefineOutput(2, TList::Class());  //Output for the Proton Cuts
  DefineOutput(3, TList::Class());  //Output for the AntiProton Cuts
  DefineOutput(4, TList::Class());  //Output for the Pion Cuts
  DefineOutput(5, TList::Class());  //Output for the AntiPion Cuts
  DefineOutput(6, TList::Class());  //Output for the Results
  DefineOutput(7, TList::Class());  //Output for the Results QA 
  DefineOutput(8, TList::Class());  //Output for the Results 3D Femto
  if (fIsMC) {
    DefineOutput(9, TList::Class());  //Output for the Proton MC
    DefineOutput(10, TList::Class());  //Output for the AntiProton MC
    DefineOutput(11, TList::Class());  //Output for the Pion MC
    DefineOutput(12, TList::Class());  //Output for the AntiPion MC
  }
}

AliAnalysisTaskNanoFemtoProtonPion::~AliAnalysisTaskNanoFemtoProtonPion() {
  delete fEvent;
  delete fTrack;
  delete fTrackCutsPion; 
  delete fTrackCutsAntiPion;
  delete fTrackCutsProton;
  delete fTrackCutsAntiProton;
  delete fPairCleaner;
  delete fPartColl;
}

//==================================================================================================================================================

void AliAnalysisTaskNanoFemtoProtonPion::UserCreateOutputObjects() {

  fGTI = new AliVTrack*[fTrackBufferSize];

  if (!fEventCuts) {
    AliError("No Event cuts \n");
  } else {
    fEventCuts->InitQA();
  }

  if (!fTrackCutsProton) {
    AliError("No Proton cuts \n");
  } else {
    fTrackCutsProton->Init();
    fProtonList = fTrackCutsProton->GetQAHists();
    if (fIsMC) {
    fProtonMCList = fTrackCutsProton->GetMCQAHists();
    }
  }

  if (!fTrackCutsAntiProton) {
    AliError("No AntiProton cuts \n");
  } else {
    fTrackCutsAntiProton->Init();
    fAntiProtonList = fTrackCutsAntiProton->GetQAHists();
    if (fIsMC) {
      fAntiProtonMCList = fTrackCutsAntiProton->GetMCQAHists();
    }
  }

  if (!fTrackCutsPion) {
    AliError("No Proton cuts \n");
  } else {
    fTrackCutsPion->Init();
    fPionList = fTrackCutsPion->GetQAHists();
    if (fIsMC) {
      fPionMCList = fTrackCutsPion->GetMCQAHists();
    }
  }

  if (!fTrackCutsAntiPion) {
    AliError("No Proton cuts \n");
  } else {
    fTrackCutsAntiPion->Init();
    fAntiPionList = fTrackCutsAntiPion->GetQAHists();
    if (fIsMC) {
      fAntiPionMCList = fTrackCutsAntiPion->GetMCQAHists();
    }
  }
  //////////////////////////////////////////////////////////////////////////
  if (!fConfig) {
    AliError("No Correlation Config \n");
  } else {
    fPartColl = new AliFemtoDreamPartCollection(fConfig,
        fConfig->GetMinimalBookingME());
    fPairCleaner = new AliFemtoDreamPairCleaner(2, 0,
        fConfig->GetMinimalBookingME()); 
  }

  fEvent = new AliFemtoDreamEvent(true, !fisLightWeight,GetCollisionCandidates(), true);
  fEvent->SetMultiplicityEstimator(fConfig->GetMultiplicityEstimator());
  fTrack = new AliFemtoDreamTrack();
  fTrack->SetUseMCInfo(fIsMC);

  if (!fEventCuts->GetMinimalBooking()) {
    fEvtList = fEventCuts->GetHistList();
  } else {
    fEvtList = new TList();
    fEvtList->SetName("EventCuts");
    fEvtList->SetOwner();
  }

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


  //////////////////////////////////////////////////////////////////////
  //Specific histograms from 3D femtoscopy 

  this->InitializeArrays();

  int CounterPassedCombinations = 0; 
  for(int i=0; i<10; i++)
  {
     if(fCombinations[i][0]<0)
	break;
     CounterPassedCombinations++; 
  }

  const int PassedCombinations = CounterPassedCombinations;
  if(PassedCombinations>10)
	AliFatal("More than 10 combinations passed");

   if(fDoThreeDFemto)
   {
    fResultsThreeDFemto = new TList();
    fResultsThreeDFemto->SetOwner();
    fResultsThreeDFemto->SetName("ResultsThreeDFemto");

    //-------------------------------------------------------------------
    //1D Same Event
    fSameEvent_List_OneDimensional = new TList();
    fSameEvent_List_OneDimensional->SetOwner();
    fSameEvent_List_OneDimensional->SetName("SameEventOneDimensional");

    fSameEvent_OneDimensional = new TH1F*[10];
    for (int i = 0; i < PassedCombinations; i++) {
        std::string title = "SameEvent_OneDimensional_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
	fSameEvent_OneDimensional[i] =  new TH1F(title.data(),title.data(), 8000, 0, 8);
	fSameEvent_List_OneDimensional->Add(fSameEvent_OneDimensional[i]);
    }

    fSameEventMult_OneDimensional = new TH2F*[10];
    for (int i = 0; i < PassedCombinations; ++i) {
        std::string title = "SameEventMult_OneDimensional_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
	fSameEventMult_OneDimensional[i] =  new TH2F(title.data(),title.data(), 8000, 0, 8,26,1,27);
        if(fRunPlotMult){fSameEvent_List_OneDimensional->Add(fSameEventMult_OneDimensional[i]);}
    }

    fSameEventDeltaEtaDeltaPhi_List_OneDimensional = new TList();
    fSameEventDeltaEtaDeltaPhi_List_OneDimensional->SetOwner();
    fSameEventDeltaEtaDeltaPhi_List_OneDimensional->SetName("SameEventDeltaEtaDeltaPhiOneDimensional");
  
    fSameEventPhiTheta_OneDimensional = new TH2F*[20]; 
    for (int i = 0; i < PassedCombinations; ++i) {
    	std::string titlebefore = "SameEventDeltaEtaDeltaPhi_OneDimensional_Before_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
    	std::string titleafter = "SameEventDeltaEtaDeltaPhi_OneDimensional_After_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
        fSameEventPhiTheta_OneDimensional[i] = new TH2F(titlebefore.data(),titlebefore.data(), 500, -0.15,0.15,500,-0.15,0.15);
        fSameEventPhiTheta_OneDimensional[10+i] = new TH2F(titleafter.data(),titleafter.data(), 500, -0.15,0.15,500,-0.15,0.15);
	if(fRunPlotPhiTheta){
	  fSameEventDeltaEtaDeltaPhi_List_OneDimensional->Add(fSameEventPhiTheta_OneDimensional[i]);
	  fSameEventDeltaEtaDeltaPhi_List_OneDimensional->Add(fSameEventPhiTheta_OneDimensional[10+i]);
	}
    }

    //-------------------------------------------------------------------
    //1D Mixed Event
    fMixedEvent_List_OneDimensional = new TList();
    fMixedEvent_List_OneDimensional->SetOwner();
    fMixedEvent_List_OneDimensional->SetName("MixedEventOneDimensional");

    fMixedEvent_OneDimensional = new TH1F*[10];
    for (int i = 0; i < PassedCombinations; i++) {
	std::string title = "MixedEvent_OneDimensional_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
	fMixedEvent_OneDimensional[i] =  new TH1F(title.data(),title.data(), 8000, 0, 8);
	fMixedEvent_List_OneDimensional->Add(fMixedEvent_OneDimensional[i]);
    }

    fMixedEventMult_OneDimensional = new TH2F*[10];
    for (int i = 0; i < PassedCombinations; ++i) {
	std::string title = "MixedEventMult_OneDimensional_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
	fMixedEventMult_OneDimensional[i] =  new TH2F(title.data(),title.data(), 8000, 0, 8,26,1,27);
        if(fRunPlotMult){fMixedEvent_List_OneDimensional->Add(fMixedEventMult_OneDimensional[i]);}
    }

    fMixedEventDeltaEtaDeltaPhi_List_OneDimensional = new TList();
    fMixedEventDeltaEtaDeltaPhi_List_OneDimensional->SetOwner();
    fMixedEventDeltaEtaDeltaPhi_List_OneDimensional->SetName("MixedEventDeltaEtaDeltaPhiOneDimensional"); 

    fMixedEventTripletPhiTheta_OneDimensional = new TH2F*[20]; 
    for (int i = 0; i < PassedCombinations; ++i) {
    	std::string titlebefore = "MixedEventDeltaEtaDeltaPhi_OneDimensional_Before_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
    	std::string titleafter = "MixedEventDeltaEtaDeltaPhi_OneDimensional_After_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
        fMixedEventTripletPhiTheta_OneDimensional[i] = new TH2F(titlebefore.data(),titlebefore.data(), 500, -0.15,0.15,500,-0.15,0.15);
        fMixedEventTripletPhiTheta_OneDimensional[10+i] = new TH2F(titleafter.data(),titleafter.data(), 500, -0.15,0.15,500,-0.15,0.15);
	if(fRunPlotPhiTheta){
	  fMixedEventDeltaEtaDeltaPhi_List_OneDimensional->Add(fMixedEventTripletPhiTheta_OneDimensional[i]);
	  fMixedEventDeltaEtaDeltaPhi_List_OneDimensional->Add(fMixedEventTripletPhiTheta_OneDimensional[10+i]);
        }
    }

    //-------------------------------------------------------------------
    fResultsThreeDFemto->Add(fSameEvent_List_OneDimensional); 
    if(fRunPlotPhiTheta){ fResultsThreeDFemto->Add(fSameEventDeltaEtaDeltaPhi_List_OneDimensional); }
    fResultsThreeDFemto->Add(fMixedEvent_List_OneDimensional); 
    if(fRunPlotPhiTheta){ fResultsThreeDFemto->Add(fMixedEventDeltaEtaDeltaPhi_List_OneDimensional); }
 

  } //if(fDoThreeDFemto)

  //////////////////////////////////////////////////////////////////////

  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fPionList);
  PostData(5, fAntiPionList);
  PostData(6, fResults);
  PostData(7, fResultsQA); 
  PostData(8, fResultsThreeDFemto);

  if (fTrackCutsProton->GetIsMonteCarlo()) {
    PostData(9, fProtonMCList);
  }
  if (fTrackCutsAntiProton->GetIsMonteCarlo()) {
    PostData(10, fAntiProtonMCList);
  }

  if (fTrackCutsPion->GetIsMonteCarlo()) {
    PostData(11, fPionMCList);
  }
  if (fTrackCutsAntiPion->GetIsMonteCarlo()) {
    PostData(12, fAntiPionMCList);
  }

 
  // Mixed event distribution ------------------------------------------------------------------------------
  // Take care of the mixing PartContainer
  auto ZVtxBinsSize = fConfig->GetNZVtxBins();
  auto MultBinsSize = fConfig->GetNMultBins();

  static std::vector<int> PDGCodes = fConfig->GetPDGCodes();

  for(int iZVtx = 0; iZVtx<ZVtxBinsSize; iZVtx++){
    std::vector<std::vector<AliFemtoDreamPartContainer>> MultContainer;
    for(int iMult = 0; iMult<MultBinsSize; iMult++){
      std::vector<AliFemtoDreamPartContainer> AllUsedParticles;
      for(unsigned int iSpecies = 0; iSpecies<PDGCodes.size(); iSpecies++){
        auto tempPartContainer = new AliFemtoDreamPartContainer(fConfig->GetMixingDepth());
        AllUsedParticles.push_back(*tempPartContainer);
      }
      MultContainer.push_back(AllUsedParticles);
    }
    fPartContainer.push_back(MultContainer);
  } 
 
} //void AliAnalysisTaskNanoFemtoProtonPion::UserCreateOutputObjects()

//==================================================================================================================================================

void AliAnalysisTaskNanoFemtoProtonPion::UserExec(Option_t*) {
  AliVEvent *Event= fInputEvent;

  if (!Event) {
    AliWarning("No Input Event");
  } else {

  fEvent->SetEvent(Event);
  if (fEventCuts->isSelected(fEvent)) {
       return;
  }

  ResetGlobalTrackReference();
  for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack*>(Event->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard NanoAOD");
      return;
    }
    StoreGlobalTrackReference(track);
  }

  fTrack->SetGlobalTrackInfo(fGTI, fTrackBufferSize);
  static std::vector<AliFemtoDreamBasePart> SelectedPions; 
  SelectedPions.clear();
  static std::vector<AliFemtoDreamBasePart> SelectedAntiPions; 
  SelectedAntiPions.clear();
  static std::vector<AliFemtoDreamBasePart> SelectedProtons;
  SelectedProtons.clear();
  static std::vector<AliFemtoDreamBasePart> SelectedAntiProtons;
  SelectedAntiProtons.clear();

  //Now we loop over all the tracks in the reconstructed event.
  for (int iTrack = 0; iTrack < Event->GetNumberOfTracks(); ++iTrack) {
    AliVTrack *track = static_cast<AliVTrack*>(Event->GetTrack(iTrack));
    if (!track) {
      AliFatal("No Standard AOD");
      return;
    }

    fTrack->SetTrack(track, fInputEvent);
    if (fTrackCutsProton->isSelected(fTrack)) {
      SelectedProtons.push_back(*fTrack);
    }
    if (fTrackCutsAntiProton->isSelected(fTrack)) {
      SelectedAntiProtons.push_back(*fTrack);
    }
    if (fTrackCutsPion->isSelected(fTrack)){ 
      SelectedPions.push_back(*fTrack);
    }
    if (fTrackCutsAntiPion->isSelected(fTrack)){
      SelectedAntiPions.push_back(*fTrack);
    }
  }

  //loop once over the MC stack to calculate Efficiency/Purity
  if (fIsMC) {
  AliAODInputHandler *eventHandler = dynamic_cast<AliAODInputHandler*>(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
  AliMCEvent* fMC = eventHandler->MCEvent();

  for (int iPart = 0; iPart < (fMC->GetNumberOfTracks()); iPart++) {
      AliAODMCParticle *mcPart = (AliAODMCParticle*) fMC->GetTrack(iPart);
      if (mcPart->IsPhysicalPrimary()) {
        if (mcPart->GetPdgCode() == fTrackCutsProton->GetPDGCode()) {
          fTrackCutsProton->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fTrackCutsAntiProton->GetPDGCode()) {
          fTrackCutsAntiProton->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fTrackCutsPion->GetPDGCode()) {
          fTrackCutsPion->FillGenerated(mcPart->Pt());
        } else if (mcPart->GetPdgCode() == fTrackCutsAntiPion->GetPDGCode()) {
          fTrackCutsAntiPion->FillGenerated(mcPart->Pt());
        }
      }
    }
  }

  if(fDoPairCleaning){
    fPairCleaner->CleanTrackAndDecay(&SelectedProtons, &SelectedPions, 0); 
    fPairCleaner->CleanTrackAndDecay(&SelectedAntiProtons, &SelectedAntiPions, 0); 
  }

  fPairCleaner->ResetArray();
  fPairCleaner->StoreParticle(SelectedProtons);
  fPairCleaner->StoreParticle(SelectedAntiProtons);
  fPairCleaner->StoreParticle(SelectedPions);
  fPairCleaner->StoreParticle(SelectedAntiPions);

  //Official FemtoDream Two-Body Calculations
  if(fDoOfficialFemto){ fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),fEvent->GetZVertex(), fEvent->GetRefMult08(), fEvent->GetV0MCentrality()); }

  //Three Dimensional Two-Body Calculations
 if(fDoThreeDFemto)
  {
    int CounterPassedCombinations = 0; 
    for(int i=0; i<10; i++)
    {
      if(fCombinations[i][0]<0)
	 break;
      CounterPassedCombinations++; 
    }
    const int PassedCombinations = CounterPassedCombinations;

    static std::vector<int> PDGCodes = fConfig->GetPDGCodes();
    int bins[2] = { 0, 0 };
    float ZVtx = fEvent->GetZVertex();
    float Mult = fEvent->GetMultiplicity();
    fPartColl->FindBin(ZVtx, Mult, bins); //is this needed? 

    // Same event distribution ------------------------------------------
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector = fPairCleaner->GetCleanParticles();

    for (int i = 0; i < PassedCombinations; i++) {
    	FillPairDistributionSE(ParticleVector,fCombinations[i][0],fCombinations[i][1],PDGCodes,bins[1],fSameEvent_OneDimensional[i],fSameEventMult_OneDimensional[i], fSameEventPhiTheta_OneDimensional,i,*fConfig);
    }

    // Mixed event distribution -----------------------------------------
    if (!(bins[0] == -99 || bins[1] == -99)) {
      auto itZVtx = fPartContainer.begin()+ bins[0];
      auto itMult = itZVtx->begin() + bins[1];

      for (int i = 0; i < PassedCombinations; i++) {
    	 FillPairDistributionME(ParticleVector,*itMult,fCombinations[i][0],fCombinations[i][1],PDGCodes,bins[1],fMixedEvent_OneDimensional[i],fMixedEventMult_OneDimensional[i], fMixedEventTripletPhiTheta_OneDimensional,i,*fConfig);

	 if(fCombinations[i][0] != fCombinations[i][1]) //if the two particles are not the same species, we can mix a second time
	 {
	     FillPairDistributionME(ParticleVector,*itMult,fCombinations[i][1],fCombinations[i][0],PDGCodes,bins[1],fMixedEvent_OneDimensional[i],fMixedEventMult_OneDimensional[i], fMixedEventTripletPhiTheta_OneDimensional,i,*fConfig);
	 }	 
      } //for (int i = 0; i < PassedCombinations; i++)

     SetMixedEvent(ParticleVector, &(*itMult));

    }//if (!(bins[0] == -99 || bins[1] == -99))
   
  }//if(fDoThreeDFemto) 
 
  } //else 

  PostData(1, fEvtList);
  PostData(2, fProtonList);
  PostData(3, fAntiProtonList);
  PostData(4, fPionList);
  PostData(5, fAntiPionList);
  PostData(6, fResults);
  PostData(7, fResultsQA);
  PostData(8, fResultsThreeDFemto);
//-----------------------------------------MCTracksStorage------------------------------------------------------------------------------
  if (fTrackCutsProton->GetIsMonteCarlo()) {
    PostData(9, fProtonMCList);
  }
  if (fTrackCutsAntiProton->GetIsMonteCarlo()) {
    PostData(10, fAntiProtonMCList);
  }
  if (fTrackCutsPion->GetIsMonteCarlo()) {
    PostData(11, fPionMCList);
  }
  if (fTrackCutsAntiPion->GetIsMonteCarlo()) {
    PostData(12, fAntiPionMCList);
  }
}//void AliAnalysisTaskNanoFemtoProtonPion::UserExec(Option_t*)

//==================================================================================================================================================

void AliAnalysisTaskNanoFemtoProtonPion::ResetGlobalTrackReference() {

  for (UShort_t i = 0; i < fTrackBufferSize; i++) {
    fGTI[i] = 0;
  }
}

//==================================================================================================================================================

void AliAnalysisTaskNanoFemtoProtonPion::StoreGlobalTrackReference(AliVTrack *track) {
  AliNanoAODTrack *nanoTrack = dynamic_cast<AliNanoAODTrack*>(track);
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
    if ((!nanoTrack->GetFilterMap()) && (!track->GetTPCNcls())) {
      return;
    }
    if (dynamic_cast<AliNanoAODTrack*>(fGTI[trackID])->GetFilterMap()
        || fGTI[trackID]->GetTPCNcls()) {
      printf("Warning! global track info already there!");
      printf("         TPCNcls track1 %u track2 %u",
             (fGTI[trackID])->GetTPCNcls(), track->GetTPCNcls());
      printf("         FilterMap track1 %u track2 %u\n",
             dynamic_cast<AliNanoAODTrack*>(fGTI[trackID])->GetFilterMap(),
             nanoTrack->GetFilterMap());
    }
  }
  (fGTI[trackID]) = track;
}//void AliAnalysisTaskNanoFemtoProtonPion::StoreGlobalTrackReference(AliVTrack *track)

//==================================================================================================================================================

void AliAnalysisTaskNanoFemtoProtonPion::FillPairDistributionSE(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies, std::vector<int> PDGCodes, int mult, TH1F* hist, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int CombinationNumber, AliFemtoDreamCollConfig Config){

  auto Particle1Vector = ParticleVector.begin()+firstSpecies;
  auto Particle2Vector = ParticleVector.begin()+secondSpecies;

  // Get the PID codes std::vector<int>
  auto itPDGPar1 = PDGCodes.begin()+firstSpecies;
  auto itPDGPar2 = PDGCodes.begin()+secondSpecies;

  // Get particle masses
  auto massParticle1 = TDatabasePDG::Instance()->GetParticle(*itPDGPar1)->Mass();
  auto massParticle2 = TDatabasePDG::Instance()->GetParticle(*itPDGPar2)->Mass();

  unsigned int NumberDaughtersParticle1 = 0;
  unsigned int NumberDaughtersParticle2 = 0;

  //Number of daughter particles 
  if(abs(*itPDGPar1)==211) NumberDaughtersParticle1 = 1;  
  if(abs(*itPDGPar1)==2212) NumberDaughtersParticle1 = 1;
  if(abs(*itPDGPar2)==211) NumberDaughtersParticle2 = 1;
  if(abs(*itPDGPar2)==2212) NumberDaughtersParticle2 = 1;

  unsigned int PairDaughterIdentifier = NumberDaughtersParticle1*10+NumberDaughtersParticle2;

  // Loop over first particle
  for (auto iPart1 = Particle1Vector->begin(); iPart1 != Particle1Vector->end(); ++iPart1) {
    // If second particle species is different than first - start with the first particle in the vector
    auto iPart2 = Particle2Vector->begin();
    // If second particle  and first are the species, start second loop from the next particle (to not double count)
    if (firstSpecies==secondSpecies) iPart2 = iPart1+1;
    // Loop over second particle ...
    for (; iPart2 != Particle2Vector->end(); ++iPart2) {

        // Now we have the three particles, lets create their Lorentz vectors
        TLorentzVector Particle1_LV, Particle2_LV;
        Particle1_LV.SetXYZM(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(),iPart1->GetMomentum().Z(), massParticle1);
        Particle2_LV.SetXYZM(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(),iPart2->GetMomentum().Z(), massParticle2);
        // Get momentum
        float RelativeMomentum = AliFemtoDreamHigherPairMath::RelativePairMomentum(Particle1_LV, Particle2_LV);

        bool PassedClosePairRejection = true;
	if(fDoClosePairRejection)
	{
	  PassedClosePairRejection =  DeltaEtaDeltaPhi(firstSpecies, secondSpecies, *iPart1,*iPart2, *itPDGPar1, *itPDGPar2, PairDaughterIdentifier, fSameEventPhiTheta_OneDimensional[CombinationNumber],fSameEventPhiTheta_OneDimensional[10+CombinationNumber],Config,RelativeMomentum); 
	}
        if(!PassedClosePairRejection) {continue;}

        hist->Fill(RelativeMomentum);
        hist2d->Fill(RelativeMomentum,mult+1);
    }
  }
} //AliAnalysisTaskNanoFemtoProtonPion::FillPairDistributionSE

//==================================================================================================================================================


void AliAnalysisTaskNanoFemtoProtonPion::FillPairDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME, std::vector<int> PDGCodes, int mult, TH1F* hist, TH2F* hist2d, TH2F **fEventTripletPhiThetaArray, int CombinationNumber, AliFemtoDreamCollConfig Config){ 

  auto ParticleSE = ParticleVector.begin()+speciesSE;
  auto MixedEventContainer = fPartContainer.begin()+speciesME;

  // Get the PID codes std::vector<int>
  auto itPDGParSE = PDGCodes.begin()+speciesSE;
  auto itPDGParME = PDGCodes.begin()+speciesME;

  unsigned int NumberDaughtersParticle1 = 0;
  unsigned int NumberDaughtersParticle2 = 0;

  //Number of daughter particles 
  if(abs(*itPDGParSE)==211) NumberDaughtersParticle1 = 1;  
  if(abs(*itPDGParSE)==2212) NumberDaughtersParticle1 = 1;
  if(abs(*itPDGParME)==211) NumberDaughtersParticle2 = 1;
  if(abs(*itPDGParME)==2212) NumberDaughtersParticle2 = 1;

  unsigned int PairDaughterIdentifier = NumberDaughtersParticle1*10+NumberDaughtersParticle2;

  // Get particle masses
  auto massParticleSE = TDatabasePDG::Instance()->GetParticle(*itPDGParSE)->Mass();
  auto massParticleME = TDatabasePDG::Instance()->GetParticle(*itPDGParME)->Mass();

  // Loop over first particle
  for (auto iPart1 = ParticleSE->begin(); iPart1 != ParticleSE->end(); ++iPart1) {
    // Loop over second particle ...
    for (int iDepth1 = 0; iDepth1 < (int) MixedEventContainer->GetMixingDepth(); ++iDepth1) {
      std::vector<AliFemtoDreamBasePart> iEvent2 = MixedEventContainer->GetEvent(iDepth1); 
      for ( auto iPart2 = iEvent2.begin(); iPart2 != iEvent2.end(); ++iPart2) {

        // Now we have the three particles, lets create their Lorentz vectors
        TLorentzVector part1_LorVec, part2_LorVec;
        part1_LorVec.SetXYZM(iPart1->GetMomentum().X(), iPart1->GetMomentum().Y(),iPart1->GetMomentum().Z(),massParticleSE);
        part2_LorVec.SetXYZM(iPart2->GetMomentum().X(), iPart2->GetMomentum().Y(),iPart2->GetMomentum().Z(),massParticleME);
        // Get momentum
        float RelativeMomentum = AliFemtoDreamHigherPairMath::RelativePairMomentum(part1_LorVec, part2_LorVec);

        bool PassedClosePairRejection = true;
	if(fDoClosePairRejection)
	{
	  PassedClosePairRejection =  DeltaEtaDeltaPhi(speciesSE, speciesME, *iPart1,*iPart2, *itPDGParSE, *itPDGParME, PairDaughterIdentifier, fMixedEventTripletPhiTheta_OneDimensional[CombinationNumber],fMixedEventTripletPhiTheta_OneDimensional[10+CombinationNumber],Config, RelativeMomentum); 
	}
        if(!PassedClosePairRejection) {continue;}
   
        hist->Fill(RelativeMomentum);
        hist2d->Fill(RelativeMomentum,mult+1);
      }
    }
  }
} //void AliAnalysisTaskNanoFemtoProtonPion::FillPairDistributionME
 
//==================================================================================================================================================

bool AliAnalysisTaskNanoFemtoProtonPion::DeltaEtaDeltaPhi(int species1, int species2,
                                                   AliFemtoDreamBasePart &part1,
                                                   AliFemtoDreamBasePart &part2,
                                                   int part1PDGcode,
                                                   int part2PDGcode, unsigned int PairDaughterIdentifier, TH2F* beforeHist,TH2F* afterHist,
                                                   AliFemtoDreamCollConfig Config, double RelativeMomentum) {

  // PairDaughterIdentifier = ij where i is the number of daughters for first particle, j for the second

  static const float piHi = TMath::Pi();
  //auto fDeltaPhiSqMax = Config.GetDeltaPhiMax() * Config.GetDeltaPhiMax();
  //auto fDeltaEtaSqMax = Config.GetDeltaEtaMax() * Config.GetDeltaEtaMax() ;

  double DeltaPhiSqMaxValue = 0;
  double DeltaEtaSqMaxValue = 0;

  DeltaPhiSqMaxValue = Config.GetDeltaPhiMax() * Config.GetDeltaPhiMax();
  DeltaEtaSqMaxValue = Config.GetDeltaEtaMax() * Config.GetDeltaEtaMax() ;

  bool pass = true;
  // if nDaug == 1 => Single Track, else decay
  unsigned int nDaug1 = (unsigned int) PairDaughterIdentifier / 10;
  if (nDaug1 > 9) {
    AliWarning("you are doing something wrong \n");
  }
  if (nDaug1 > part1.GetPhiAtRaidius().size()) {
    TString outMessage =
        TString::Format(
            "For pair number %u your number of Daughters 1 (%u) and Radii 1 (%u) do not correspond \n",
            PairDaughterIdentifier, nDaug1, (unsigned int)part1.GetPhiAtRaidius().size());
    AliWarning(outMessage.Data());
  }
  unsigned int nDaug2 = (unsigned int) PairDaughterIdentifier % 10;

  if (nDaug2 > part2.GetPhiAtRaidius().size()) {
    TString outMessage =
        TString::Format(
            "For pair number %u your number of Daughters 2 (%u) and Radii 2 (%u) do not correspond \n",
            PairDaughterIdentifier, nDaug2, (unsigned int)part2.GetPhiAtRaidius().size());
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
        if(RelativeMomentum<2.){
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
            if(RelativeMomentum<2.){
              afterHist->Fill(dphiAvg/ (float) size, deta);
            }
          }
        }
      }
    }
  }
  return pass;
}//bool AliAnalysisTaskNanoFemtoProtonPion::DeltaEtaDeltaPhi

//==================================================================================================================================================

void AliAnalysisTaskNanoFemtoProtonPion::SetMixedEvent(
    std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer> *PartContainer) {
  // Feed this function with GetCleanParticles output and fill the mixed events for different particles
  for(unsigned int iSpecies = 0; iSpecies<ParticleVector.size(); iSpecies++){
    if ((ParticleVector.begin()+iSpecies)->size() > 0) {
      (PartContainer->begin()+iSpecies)->SetEvent(*(ParticleVector.begin()+iSpecies));
    }
  }
}

//==================================================================================================================================================

void AliAnalysisTaskNanoFemtoProtonPion::InitializeArrays()
{
  std::istringstream issCombination(fCombinationInput);
  int Counter = 0;

  for(int i=0; i<10; i++)
  {
    fCombinations[i][0] = -1; 
    fCombinations[i][1] = -1;   
  }

  do {
    if(Counter>9){
	AliFatal("Max. 10 entries of combinations allowed"); 
	break; 
    }

    std::string subs;
    issCombination >> subs;

    if(subs.size() == 0)
	continue;

    int FirstEntry = ((int)subs[0]) - ((int)'0');
    int SecondEntry = ((int)subs[1]) - ((int)'0');

    if(FirstEntry < 0 || SecondEntry < 0) 
	AliFatal("No number below 0 allowed in passed combinations");
    if(FirstEntry > 3 || SecondEntry > 3) 
	AliFatal("No number above 3 allowed in passed combinations");

    fCombinations[Counter][0] = FirstEntry; 
    fCombinations[Counter][1] = SecondEntry;   

    Counter++; 
   } while (issCombination);

  std::istringstream issNameTag(fNameTagInput);
  Counter = 0;

  while (issNameTag) {
    if(Counter>3)
	break;

    std::string subs;
    issNameTag >> subs;

    if(subs.size() == 0)
    	continue; 

    fNameTags[Counter] = subs;
    Counter++; 
   } 

   if(Counter!=4){
      AliFatal("Passed NameTag has to have exaclty 4 entries!");
   }

} // void AliAnalysisTaskNanoFemtoProtonPion::InitializeArrays()

//==================================================================================================================================================
