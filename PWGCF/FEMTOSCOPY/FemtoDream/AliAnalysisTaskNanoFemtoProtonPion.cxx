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
    fClosePairRejectionInput(""),
    fDoOfficialFemto(false),
    fDoOwnFemto(false),
    fDoThreeDFemto(false),
    fRunPlotMult(false),
    fRunPlotPhiTheta(false),
    fDoAncestors(false),
    fRemoveMCResonances(true),
    fRemoveMCResonanceDaughters(true),
    fDoInvMassPlot(false), 
    fDoResonanceLorentzFactor(true),
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
    fSameEvent_OneDimensional(nullptr),
    fSameEventMult_OneDimensional(nullptr),
    fSameEvent_OneDimensional_Ancestors(nullptr),
    fSameEventMult_OneDimensional_Ancestors(nullptr),
    fSameEvent_InvMass(nullptr),
    fSameEvent_InvMass_MCResonance(nullptr),
    fMixedEvent_List_OneDimensional(nullptr),
    fMixedEvent_OneDimensional(nullptr),
    fMixedEventMult_OneDimensional(nullptr),
    fMixedEvent_InvMass(nullptr),
    fSameEvent_List_ThreeDimensional(nullptr),
    fSameEvent_ThreeDimensional(nullptr),
    fMixedEvent_List_ThreeDimensional(nullptr),
    fMixedEvent_ThreeDimensional(nullptr),
    fSameEventDeltaEtaDeltaPhi_List(nullptr),
    fSameEventPhiTheta(nullptr),
    fSameEventPhiTheta_Ancestors(nullptr),
    fMixedEventDeltaEtaDeltaPhi_List(nullptr),
    fMixedEventPhiTheta(nullptr),
    fResonanceLorentzFactor(nullptr),
    fInvMassResonancesMCTruth(nullptr){
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
    fClosePairRejectionInput(""),
    fDoOfficialFemto(false),
    fDoOwnFemto(false),
    fDoThreeDFemto(false),
    fRunPlotMult(false),
    fRunPlotPhiTheta(false),
    fDoAncestors(false),
    fRemoveMCResonances(true),
    fRemoveMCResonanceDaughters(true),
    fDoInvMassPlot(false), 
    fDoResonanceLorentzFactor(true),
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
    fSameEvent_OneDimensional(nullptr),
    fSameEventMult_OneDimensional(nullptr),
    fSameEvent_OneDimensional_Ancestors(nullptr),
    fSameEventMult_OneDimensional_Ancestors(nullptr),
    fSameEvent_InvMass(nullptr),
    fSameEvent_InvMass_MCResonance(nullptr),
    fMixedEvent_List_OneDimensional(nullptr),
    fMixedEvent_OneDimensional(nullptr),
    fMixedEventMult_OneDimensional(nullptr),
    fMixedEvent_InvMass(nullptr),
    fSameEvent_List_ThreeDimensional(nullptr),
    fSameEvent_ThreeDimensional(nullptr),
    fMixedEvent_List_ThreeDimensional(nullptr),
    fMixedEvent_ThreeDimensional(nullptr),
    fSameEventDeltaEtaDeltaPhi_List(nullptr),
    fSameEventPhiTheta(nullptr),
    fSameEventPhiTheta_Ancestors(nullptr),
    fMixedEventDeltaEtaDeltaPhi_List(nullptr),
    fMixedEventPhiTheta(nullptr),
    fResonanceLorentzFactor(nullptr),
    fInvMassResonancesMCTruth(nullptr){
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
  //fEvent->SetCalcSpherocity(fEventCuts->GetDoSpherocityCuts()); 
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

  if(fDoAncestors && !fIsMC)
  { AliFatal("Cannot do Ancestor study for non MC study");  }
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

   if(fDoOwnFemto)
   {
    fResultsThreeDFemto = new TList();
    fResultsThreeDFemto->SetOwner();
    fResultsThreeDFemto->SetName("ResultsThreeDFemto");

    //-------------------------------------------------------------------
    //1D Same Event
    fSameEvent_List_OneDimensional = new TList();
    fSameEvent_List_OneDimensional->SetOwner();
    fSameEvent_List_OneDimensional->SetName("SameEventOneDimensional");

    //Close Pair Rejection Plots Same Event 
    fSameEventDeltaEtaDeltaPhi_List = new TList();
    fSameEventDeltaEtaDeltaPhi_List->SetOwner();
    fSameEventDeltaEtaDeltaPhi_List->SetName("SameEventDeltaEtaDeltaPhiOneDimensional");

    //1D SE Objects for Data and total MC ~~~~~~~~~~~~~~~~~~~
    fSameEvent_OneDimensional = new TH1F*[10];
    fSameEventMult_OneDimensional = new TH2F*[10];
    fSameEvent_InvMass = new TH1F*[10];
    fSameEvent_InvMass_MCResonance = new TH1F*[10];

    if(!fDoThreeDFemto){
       for (int i = 0; i < PassedCombinations; i++) {
          std::string title = "SameEvent_OneDimensional_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
          fSameEvent_OneDimensional[i] =  new TH1F(title.data(),title.data(), 3000, 0, 3);
          fSameEvent_List_OneDimensional->Add(fSameEvent_OneDimensional[i]);
       }

       for (int i = 0; i < PassedCombinations; ++i) {
          std::string title = "SameEventMult_OneDimensional_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
          fSameEventMult_OneDimensional[i] =  new TH2F(title.data(),title.data(), 3000, 0, 3,26,1,27);
          if(fRunPlotMult){fSameEvent_List_OneDimensional->Add(fSameEventMult_OneDimensional[i]);}
       }

       fSameEventPhiTheta = new TH2F*[20]; 
        for (int i = 0; i < PassedCombinations; ++i) {
          std::string titlebefore = "SameEventDeltaEtaDeltaPhi_OneDimensional_Before_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
          std::string titleafter = "SameEventDeltaEtaDeltaPhi_OneDimensional_After_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
          fSameEventPhiTheta[i] = new TH2F(titlebefore.data(),titlebefore.data(), 500, -0.15,0.15,500,-0.15,0.15);
          fSameEventPhiTheta[10+i] = new TH2F(titleafter.data(),titleafter.data(), 500, -0.15,0.15,500,-0.15,0.15);
          if(fRunPlotPhiTheta){
            fSameEventDeltaEtaDeltaPhi_List->Add(fSameEventPhiTheta[i]);
            fSameEventDeltaEtaDeltaPhi_List->Add(fSameEventPhiTheta[10+i]);
          }
        }
       if(fDoInvMassPlot){
         for (int i = 0; i < PassedCombinations; i++) {
           std::string title = "SameEvent_InvMass_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
           fSameEvent_InvMass[i] =  new TH1F(title.data(),title.data(), 3000, 0, 3.);
           fSameEvent_List_OneDimensional->Add(fSameEvent_InvMass[i]);
         }

         if(fIsMC){ 
          for (int i = 0; i < PassedCombinations; i++) {
            std::string title = "SameEvent_InvMass_MCTruth_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
            fSameEvent_InvMass_MCResonance[i] =  new TH1F(title.data(),title.data(), 3000, 0, 3.);
            fSameEvent_List_OneDimensional->Add(fSameEvent_InvMass_MCResonance[i]);
          }
         }
       }
    }//if(!fDoThreeDFemto)


    //SE Objects for ancestor studies in MC ~~~~~~~~~~~~~~~~~~~
    fSameEvent_OneDimensional_Ancestors = new TH1F*[20]; //0-9 common ancestors, 10-19 non common
    fSameEventMult_OneDimensional_Ancestors = new TH2F*[20]; //0-9 common ancestors, 10-19 non common
    fSameEventPhiTheta_Ancestors = new TH2F*[40];  //0-9 common ancestors before, 10-19 common after, 20-29 non common before, 30-39 non common after

    if(fDoAncestors && !fDoThreeDFemto)
    {
       for (int i = 0; i < PassedCombinations; i++) {
          std::string title_common = "SameEvent_OneDimensional_Common_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
          std::string title_noncommon = "SameEvent_OneDimensional_NonCommon_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
          fSameEvent_OneDimensional_Ancestors[i] =  new TH1F(title_common.data(),title_common.data(), 3000, 0, 5);
          fSameEvent_OneDimensional_Ancestors[10+i] =  new TH1F(title_noncommon.data(),title_noncommon.data(), 3000, 0, 5);
          fSameEvent_List_OneDimensional->Add(fSameEvent_OneDimensional_Ancestors[i]);
          fSameEvent_List_OneDimensional->Add(fSameEvent_OneDimensional_Ancestors[10+i]);
       }
   
       for (int i = 0; i < PassedCombinations; ++i) {
          std::string title_common = "SameEventMult_OneDimensional_Common_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
          std::string title_noncommon = "SameEventMult_OneDimensional_NonCommon_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
          fSameEventMult_OneDimensional_Ancestors[i] =  new TH2F(title_common.data(),title_common.data(), 3000, 0, 3,26,1,27);
          fSameEventMult_OneDimensional_Ancestors[10+i] =  new TH2F(title_noncommon.data(),title_noncommon.data(), 3000, 0, 3,26,1,27);
          if(fRunPlotMult){
            fSameEvent_List_OneDimensional->Add(fSameEventMult_OneDimensional_Ancestors[i]);
            fSameEvent_List_OneDimensional->Add(fSameEventMult_OneDimensional_Ancestors[10+i]);
          }
       }
 
       for (int i = 0; i < PassedCombinations; ++i) {
        std::string titlebefore_common = "SameEventDeltaEtaDeltaPhi_OneDimensional_Common_Before_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
        std::string titleafter_common = "SameEventDeltaEtaDeltaPhi_OneDimensional_Common_After_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
        std::string titlebefore_noncommon = "SameEventDeltaEtaDeltaPhi_OneDimensional_NonCommon_Before_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
        std::string titleafter_noncommon = "SameEventDeltaEtaDeltaPhi_OneDimensional_NonCommon_After_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
        
        fSameEventPhiTheta_Ancestors[i] = new TH2F(titlebefore_common.data(),titlebefore_common.data(), 500, -0.15,0.15,500,-0.15,0.15);
        fSameEventPhiTheta_Ancestors[10+i] = new TH2F(titleafter_common.data(),titleafter_common.data(), 500, -0.15,0.15,500,-0.15,0.15);
        fSameEventPhiTheta_Ancestors[20+i] = new TH2F(titlebefore_noncommon.data(),titlebefore_noncommon.data(), 500, -0.15,0.15,500,-0.15,0.15);
        fSameEventPhiTheta_Ancestors[30+i] = new TH2F(titleafter_noncommon.data(),titleafter_noncommon.data(), 500, -0.15,0.15,500,-0.15,0.15);
        if(fRunPlotPhiTheta){
           fSameEventDeltaEtaDeltaPhi_List->Add(fSameEventPhiTheta_Ancestors[i]);
           fSameEventDeltaEtaDeltaPhi_List->Add(fSameEventPhiTheta_Ancestors[10+i]);
           fSameEventDeltaEtaDeltaPhi_List->Add(fSameEventPhiTheta_Ancestors[20+i]);
           fSameEventDeltaEtaDeltaPhi_List->Add(fSameEventPhiTheta_Ancestors[30+i]);
         }
       }
       
    }//if(fDoAncestors)


    //3D Same Event ~~~~~~~~~~~~~~~~~~~
    fSameEvent_List_ThreeDimensional = new TList();
    fSameEvent_List_ThreeDimensional->SetOwner();
    fSameEvent_List_ThreeDimensional->SetName("SameEventThreeDimensional");

    fSameEvent_ThreeDimensional = new TH2F*[30]; //0-9: out, 10-19: side, 20-29: long
    if(fDoThreeDFemto){
       for (int i = 0; i < PassedCombinations; ++i) {
         std::string title = "SameEvent_ThreeDimensional_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
         fSameEvent_ThreeDimensional[i] =  new TH2F(title.data(),title.data(), 3000, 0, 3,26,1,27);
         fSameEvent_List_ThreeDimensional->Add(fSameEvent_ThreeDimensional[i]);
       }
    }//if(fDoThreeDFemto)

    //-------------------------------------------------------------------
    //1D Mixed Event ~~~~~~~~~~~~~~~~~~~
    fMixedEvent_List_OneDimensional = new TList();
    fMixedEvent_List_OneDimensional->SetOwner();
    fMixedEvent_List_OneDimensional->SetName("MixedEventOneDimensional");

    fMixedEvent_OneDimensional = new TH1F*[10];
    fMixedEventMult_OneDimensional = new TH2F*[10];
    fMixedEvent_InvMass = new TH1F*[10];

    if(!fDoThreeDFemto){
      for (int i = 0; i < PassedCombinations; i++) {
        std::string title = "MixedEvent_OneDimensional_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
        fMixedEvent_OneDimensional[i] =  new TH1F(title.data(),title.data(), 3000, 0, 3);
        fMixedEvent_List_OneDimensional->Add(fMixedEvent_OneDimensional[i]);
      }

      for (int i = 0; i < PassedCombinations; ++i) {
        std::string title = "MixedEventMult_OneDimensional_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
        fMixedEventMult_OneDimensional[i] =  new TH2F(title.data(),title.data(), 3000, 0, 3,26,1,27);
        if(fRunPlotMult){fMixedEvent_List_OneDimensional->Add(fMixedEventMult_OneDimensional[i]);}
      }
  
      if(fDoInvMassPlot){
        for (int i = 0; i < PassedCombinations; i++) {
          std::string title = "MixedEvent_InvMass_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
          fMixedEvent_InvMass[i] =  new TH1F(title.data(),title.data(), 3000, 0, 3.);
          fMixedEvent_List_OneDimensional->Add(fMixedEvent_InvMass[i]);
        }
      }
    }//if(!fDoThreeDFemto)

    //3D Mixed Event ~~~~~~~~~~~~~~~~~~~
    fMixedEvent_List_ThreeDimensional = new TList();
    fMixedEvent_List_ThreeDimensional->SetOwner();
    fMixedEvent_List_ThreeDimensional->SetName("MixedEventThreeDimensional");

    fMixedEvent_ThreeDimensional = new TH2F*[30]; //0-9: out, 10-19: side, 20-29: long
    if(fDoThreeDFemto){
      for (int i = 0; i < PassedCombinations; ++i) {
        std::string title = "MixedEvent_ThreeDimensional_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
        fMixedEvent_ThreeDimensional[i] =  new TH2F(title.data(),title.data(), 3000, 0, 3,26,1,27);
        fMixedEvent_List_ThreeDimensional->Add(fMixedEvent_ThreeDimensional[i]);
      }
    }//if(fDoThreeDFemto)

    //Close Pair Rejecton Plots Mixed Event ~~~~~~~~~~~~~~~~~~~
    fMixedEventDeltaEtaDeltaPhi_List = new TList();
    fMixedEventDeltaEtaDeltaPhi_List->SetOwner();
    fMixedEventDeltaEtaDeltaPhi_List->SetName("MixedEventDeltaEtaDeltaPhiOneDimensional"); 

    fMixedEventPhiTheta = new TH2F*[20]; 
    for (int i = 0; i < PassedCombinations; ++i) {
      std::string titlebefore = "MixedEventDeltaEtaDeltaPhi_OneDimensional_Before_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
      std::string titleafter = "MixedEventDeltaEtaDeltaPhi_OneDimensional_After_"+fNameTags[fCombinations[i][0]]+fNameTags[fCombinations[i][1]];
      fMixedEventPhiTheta[i] = new TH2F(titlebefore.data(),titlebefore.data(), 500, -0.15,0.15,500,-0.15,0.15);
      fMixedEventPhiTheta[10+i] = new TH2F(titleafter.data(),titleafter.data(), 500, -0.15,0.15,500,-0.15,0.15);
      if(fRunPlotPhiTheta){
        fMixedEventDeltaEtaDeltaPhi_List->Add(fMixedEventPhiTheta[i]);
        fMixedEventDeltaEtaDeltaPhi_List->Add(fMixedEventPhiTheta[10+i]);
      }
    }

    //-------------------------------------------------------------------
   
    if(!fDoThreeDFemto){  
      fResultsThreeDFemto->Add(fSameEvent_List_OneDimensional); 
      fResultsThreeDFemto->Add(fMixedEvent_List_OneDimensional); 
    }
    if(fRunPlotPhiTheta){  
     fResultsThreeDFemto->Add(fSameEventDeltaEtaDeltaPhi_List);
     fResultsThreeDFemto->Add(fMixedEventDeltaEtaDeltaPhi_List); 
    }
 
  } //if(fDoOwnFemto)

  if(fIsMC){
   
    int ProtonAntiPion[33] = {2114, 12112, 1214, 22112, 32114, 1212, 32112, 2116, 12116, 12114, 42112, 21214, 31214, 11212, 9902114, 1216, 9902112, 9912112, 21212, 22114, 9912114, 2118, 11216, 9902116, 9922112, 9922114, 1218, 9901218, 99021110, 99121110, 99012112, 99021112, 3122};
    int ProtonPion[12] = {2224, 32224, 2222, 12224, 12222, 2226, 22222, 22224, 2228, 12226, 9902228, 99022212};
    
    if(fDoResonanceLorentzFactor){
      fResonanceLorentzFactor = new TH2F("fResonanceLorentzFactor","fResonanceLorentzFactor", 1990, 1.,200., 45,0.,45.);
      for(int i=0; i<33; i++){
        fResonanceLorentzFactor->GetYaxis()->SetBinLabel(1+i,Form("%d",ProtonAntiPion[i]));
      }
      for(int i=0; i<12; i++){
        fResonanceLorentzFactor->GetYaxis()->SetBinLabel(34+i,Form("%d",ProtonPion[i]));
      }
      fResults->Add(fResonanceLorentzFactor); 
    }
    if(fDoInvMassPlot){
      fInvMassResonancesMCTruth = new TH2F("fInvMassResonancesMCTruth","fInvMassResonancesMCTruth", 3000, 0.,3., 45,0.,45.);
      for(int i=0; i<33; i++){
        fInvMassResonancesMCTruth->GetYaxis()->SetBinLabel(1+i,Form("%d",ProtonAntiPion[i]));
      }
      for(int i=0; i<12; i++){
        fInvMassResonancesMCTruth->GetYaxis()->SetBinLabel(34+i,Form("%d",ProtonPion[i]));
      }
      fResults->Add(fInvMassResonancesMCTruth); 
    }
  } 

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
    return;
  } 

  fEvent->SetEvent(Event);
  if (!fEventCuts->isSelected(fEvent)) {
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
      continue;
    }

    fTrack->SetTrack(track, fInputEvent);

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
        bool RemoveTrack = false;
        while (motherID != -1) {
          lastMother = motherID;
          mcMother = (AliAODMCParticle *)mcarray->At(motherID);
          motherID = mcMother->GetMother();
          if(IsResonance(mcMother->GetPdgCode())){
             fTrack->SetMotherPDG(mcMother->GetPdgCode()); //Change the PDG of the mother so it is set to the resonance. The Mother ID keeps set to the original parton
             RemoveTrack = true;
          }
        }
        if ((lastMother != -1)) {
          mcMother = (AliAODMCParticle *)mcarray->At(lastMother);
        }
        if (mcMother) {
          int motherPDG = mcMother->GetPdgCode(); 
          if(IsResonance(motherPDG)){
             fTrack->SetMotherPDG(motherPDG); //Change the PDG of the mother so it is set to the resonance. The Mother ID keeps set to the original parton
             RemoveTrack = true;
          }
        }
        if (RemoveTrack && fRemoveMCResonanceDaughters){
           continue; 
        }
      } else {
        continue;  // if we don't have MC Information, don't use that track
      }
    } 

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
      if(fDoResonanceLorentzFactor){
        if(IsResonance(mcPart->GetPdgCode())){
          fResonanceLorentzFactor->Fill(mcPart->E()/mcPart->M(), Form("%d",abs(mcPart->GetPdgCode())), 1.);
        }
      }
      if(fDoInvMassPlot){ 
        if(IsResonance(mcPart->GetPdgCode())){
          fInvMassResonancesMCTruth->Fill(mcPart->M(), Form("%d",abs(mcPart->GetPdgCode())), 1.);
        }
      }
    }
  }

  if(fDoPairCleaning){
    fPairCleaner->CleanTrackAndDecay(&SelectedProtons, &SelectedPions, 0); 
    fPairCleaner->CleanTrackAndDecay(&SelectedAntiProtons, &SelectedAntiPions, 1); 
  }

  fPairCleaner->ResetArray();
  fPairCleaner->StoreParticle(SelectedProtons);
  fPairCleaner->StoreParticle(SelectedAntiProtons);
  fPairCleaner->StoreParticle(SelectedPions);
  fPairCleaner->StoreParticle(SelectedAntiPions);

  //Official FemtoDream Two-Body Calculations
  if(fDoOfficialFemto){ fPartColl->SetEvent(fPairCleaner->GetCleanParticles(),fEvent->GetZVertex(), fEvent->GetRefMult08(), fEvent->GetV0MCentrality()); }

  //Three Dimensional Two-Body Calculations
  if(fDoOwnFemto)
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
        if(!fDoAncestors){
         FillPairDistributionSE(ParticleVector,fCombinations[i][0],fCombinations[i][1],PDGCodes,bins[1],fClosePairRejection[i],fSameEvent_OneDimensional[i],fSameEventMult_OneDimensional[i],fSameEvent_InvMass[i],fSameEvent_InvMass_MCResonance[i],fSameEventPhiTheta,i,*fConfig);
      } else {
         FillPairDistributionSEAncestors(ParticleVector,fCombinations[i][0],fCombinations[i][1],PDGCodes,bins[1],fClosePairRejection[i],fSameEvent_OneDimensional[i],fSameEventMult_OneDimensional[i],fSameEvent_InvMass[i],fSameEventPhiTheta,fSameEvent_OneDimensional_Ancestors,fSameEventMult_OneDimensional_Ancestors,fSameEventPhiTheta_Ancestors,i,*fConfig);
      }
    }

    // Mixed event distribution -----------------------------------------
    if (!(bins[0] == -99 || bins[1] == -99)) {
      auto itZVtx = fPartContainer.begin()+ bins[0];
      auto itMult = itZVtx->begin() + bins[1];

      for (int i = 0; i < PassedCombinations; i++) {
       FillPairDistributionME(ParticleVector,*itMult,fCombinations[i][0],fCombinations[i][1],PDGCodes,bins[1],fClosePairRejection[i],fMixedEvent_OneDimensional[i],fMixedEventMult_OneDimensional[i],fMixedEvent_InvMass[i],fMixedEventPhiTheta,i,*fConfig);

        if(fCombinations[i][0] != fCombinations[i][1]) //if the two particles are not the same species, we can mix a second time
        {
            FillPairDistributionME(ParticleVector,*itMult,fCombinations[i][1],fCombinations[i][0],PDGCodes,bins[1],fClosePairRejection[i],fMixedEvent_OneDimensional[i],fMixedEventMult_OneDimensional[i],fMixedEvent_InvMass[i],fMixedEventPhiTheta,i,*fConfig);
        }   
      } //for (int i = 0; i < PassedCombinations; i++)

     SetMixedEvent(ParticleVector, &(*itMult));

    }//if (!(bins[0] == -99 || bins[1] == -99))
   
  }//if(fDoOwnFemto) 
 

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

void AliAnalysisTaskNanoFemtoProtonPion::FillPairDistributionSE(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies, std::vector<int> PDGCodes, int mult, bool DoClosePairRejection, TH1F* hist, TH2F* hist2d, TH1F* HistInvMass, TH1F* HistInvMassMCResonance, TH2F **SameEventPhiTheta_OneDimensional, int CombinationNumber, AliFemtoDreamCollConfig Config){

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
        if(DoClosePairRejection)
        {
          PassedClosePairRejection =  DeltaEtaDeltaPhi(firstSpecies, secondSpecies, *iPart1,*iPart2, *itPDGPar1, *itPDGPar2, PairDaughterIdentifier, SameEventPhiTheta_OneDimensional[CombinationNumber],SameEventPhiTheta_OneDimensional[10+CombinationNumber],Config,RelativeMomentum); 
        }
        if(!PassedClosePairRejection) {continue;}

        hist->Fill(RelativeMomentum);
        hist2d->Fill(RelativeMomentum,mult+1);

        if(fDoInvMassPlot){
          TLorentzVector Sum = Particle1_LV + Particle2_LV; 
          if(Sum.M() >= 0.){
              HistInvMass->Fill(Sum.M()); 

              if(fIsMC){
                bool HasCommonAncestor = CommonAncestors(*iPart1, *iPart2);
                if(HasCommonAncestor){
                  bool HasCommonMotherResonance = CommonMotherResonance(*iPart1, *iPart2);
                  if(HasCommonMotherResonance){
                    HistInvMassMCResonance->Fill(Sum.M());
                  }
                } 
              }//if(fIsMC)
          }
        }//if(fDoInvMassPlot)
    }
  }
} //AliAnalysisTaskNanoFemtoProtonPion::FillPairDistributionSE

//==================================================================================================================================================

void AliAnalysisTaskNanoFemtoProtonPion::FillPairDistributionSEAncestors(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, int firstSpecies,int secondSpecies, std::vector<int> PDGCodes, int mult, bool DoClosePairRejection, TH1F* hist, TH2F* hist2d, TH1F* HistInvMass, TH2F **SameEventPhiTheta_OneDimensional, TH1F **histAncestor, TH2F **hist2dAncestor, TH2F **SameEventPhiTheta_OneDimensionalAncestor, int CombinationNumber, AliFemtoDreamCollConfig Config){ 

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

        bool HasCommonAncestor = CommonAncestors(*iPart1, *iPart2);

        bool HasCommonMotherResonance = true; 
        if(HasCommonAncestor && fRemoveMCResonances){
          HasCommonMotherResonance = CommonMotherResonance(*iPart1, *iPart2);
          if(HasCommonMotherResonance){
            continue;
          }
        }

        bool PassedClosePairRejection = true; //Close pair rejection for any MC type particle
        bool PassedClosePairRejection_Ancestor = true; //Close pair rejection for any common or non common ancestor MC particles
        
        if(DoClosePairRejection)
        {
          PassedClosePairRejection =  DeltaEtaDeltaPhi(firstSpecies, secondSpecies, *iPart1,*iPart2, *itPDGPar1, *itPDGPar2, PairDaughterIdentifier, SameEventPhiTheta_OneDimensional[CombinationNumber],SameEventPhiTheta_OneDimensional[10+CombinationNumber],Config,RelativeMomentum); 
          
          if(HasCommonAncestor){
            PassedClosePairRejection_Ancestor =  DeltaEtaDeltaPhi(firstSpecies, secondSpecies, *iPart1,*iPart2, *itPDGPar1, *itPDGPar2, PairDaughterIdentifier, SameEventPhiTheta_OneDimensionalAncestor[CombinationNumber],SameEventPhiTheta_OneDimensionalAncestor[10+CombinationNumber],Config,RelativeMomentum);
          } else {
            PassedClosePairRejection_Ancestor =  DeltaEtaDeltaPhi(firstSpecies, secondSpecies, *iPart1,*iPart2, *itPDGPar1, *itPDGPar2, PairDaughterIdentifier, SameEventPhiTheta_OneDimensionalAncestor[20+CombinationNumber],SameEventPhiTheta_OneDimensionalAncestor[30+CombinationNumber],Config,RelativeMomentum);
          }
        }
  
        if(PassedClosePairRejection){
          hist->Fill(RelativeMomentum); 
          hist2d->Fill(RelativeMomentum,mult+1); 
          if(fDoInvMassPlot){
            TLorentzVector Sum = Particle1_LV + Particle2_LV; 
            if(Sum.M() >= 0.){
                HistInvMass->Fill(Sum.M()); 
            }
          }
        }
        if(PassedClosePairRejection_Ancestor){
          if(HasCommonAncestor){
            histAncestor[CombinationNumber]->Fill(RelativeMomentum); 
            hist2dAncestor[CombinationNumber]->Fill(RelativeMomentum,mult+1); 
          } else {
            histAncestor[10+CombinationNumber]->Fill(RelativeMomentum); 
            hist2dAncestor[10+CombinationNumber]->Fill(RelativeMomentum,mult+1); 
          }
       }
        
    }
  }
} //AliAnalysisTaskNanoFemtoProtonPion::FillPairDistributionSEAncestors

//==================================================================================================================================================


void AliAnalysisTaskNanoFemtoProtonPion::FillPairDistributionME(std::vector<std::vector<AliFemtoDreamBasePart>> &ParticleVector, std::vector<AliFemtoDreamPartContainer>  &fPartContainer, int speciesSE, int speciesME, std::vector<int> PDGCodes, int mult, bool DoClosePairRejection, TH1F* hist, TH2F* hist2d, TH1F* HistInvMass, TH2F **EventPhiThetaArray, int CombinationNumber, AliFemtoDreamCollConfig Config){ 

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
        if(DoClosePairRejection)
        {
          PassedClosePairRejection =  DeltaEtaDeltaPhi(speciesSE, speciesME, *iPart1,*iPart2, *itPDGParSE, *itPDGParME, PairDaughterIdentifier, EventPhiThetaArray[CombinationNumber],EventPhiThetaArray[10+CombinationNumber],Config, RelativeMomentum); 
        }
        if(!PassedClosePairRejection) {continue;}
   
        hist->Fill(RelativeMomentum);
        hist2d->Fill(RelativeMomentum,mult+1);

        if(fDoInvMassPlot){
          TLorentzVector Sum = part1_LorVec + part2_LorVec; 
          if(Sum.M() >= 0.){
              HistInvMass->Fill(Sum.M()); 
          }
        }

      }
    }
  }
} //void AliAnalysisTaskNanoFemtoProtonPion::FillPairDistributionME
 
//==================================================================================================================================================

double AliAnalysisTaskNanoFemtoProtonPion::GetQOutLCMS(const TLorentzVector Particle1, const TLorentzVector Particle2)
{  
  const double Px = Particle1.Px() + Particle2.Px();
  const double Py = Particle1.Py() + Particle2.Py();
  const double Pt = sqrt(Px*Px + Py*Py);

  const double dPx = Particle1.Px() - Particle2.Px();
  const double dPy = Particle1.Py() - Particle2.Py();

  double qout = -1.; 

  if(Pt != 0.){
    qout = abs( (Px*dPx + Py*dPy)/Pt ); 
  } 
    
  return qout;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double AliAnalysisTaskNanoFemtoProtonPion::GetQSideLCMS(const TLorentzVector Particle1, const TLorentzVector Particle2)
{
  const double Px = Particle1.Px() + Particle2.Px();
  const double Py = Particle1.Py() + Particle2.Py();
  const double Pt = sqrt(Px*Px + Py*Py);

  const double dPx = Particle1.Px() - Particle2.Px();
  const double dPy = Particle1.Py() - Particle2.Py();

  double qside = -1.; 

  if(Pt != 0.){
    qside = abs( (Px*dPy - Py*dPx)/Pt ); 
  } 
    
  return qside;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

double AliAnalysisTaskNanoFemtoProtonPion::GetQLongLCMS(const TLorentzVector Particle1, const TLorentzVector Particle2)
{
  const double E = Particle1.E() + Particle2.E();
  const double Pz = Particle1.Pz() + Particle2.Pz();
  const double Mt = sqrt(E*E - Pz*Pz);

  const double dE = Particle1.E() - Particle2.E();
  const double dPz = Particle1.Pz() - Particle2.Pz();

  double qlong = -1.; 

  if(Mt != 0.){
    qlong = abs( (E*dPz - Pz*dE)/Mt ); 
  } 
    
  return qlong;
}

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

bool AliAnalysisTaskNanoFemtoProtonPion::CommonAncestors(AliFemtoDreamBasePart& part1, AliFemtoDreamBasePart& part2) {
    bool IsCommon = false;
    if(part1.GetMotherID() == part2.GetMotherID()){
      IsCommon = true;
    }else if(part1.GetMotherID() != part2.GetMotherID()){
      IsCommon = false;
    }
    return IsCommon;
}

//==================================================================================================================================================

bool AliAnalysisTaskNanoFemtoProtonPion::CommonMotherResonance(AliFemtoDreamBasePart& part1, AliFemtoDreamBasePart& part2) {

  if(part1.GetMotherID() != part2.GetMotherID()) {
    AliFatal("AliAnalysisTaskNanoFemtoProtonPion::CommonMotherResonance: The two particle should have a common mother"); 
  }

  if(part1.GetMotherPDG() != part2.GetMotherPDG()) { //the ID is the same, but the PDG different -> Two tracks from same hard scattering but different resonances.
    return false; 
  }

  bool HasCommonMotherResonance = true;
  HasCommonMotherResonance = IsResonance(part1.GetMotherPDG()); //the resonance is of the type that should be removed 
  return HasCommonMotherResonance; 
}

//==================================================================================================================================================

bool AliAnalysisTaskNanoFemtoProtonPion::IsResonance(int PDG) {

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

//==================================================================================================================================================

void AliAnalysisTaskNanoFemtoProtonPion::InitializeArrays()
{
  std::istringstream issCombination(fCombinationInput);
  int Counter = 0;
  int PassedCombinations = 0;

  for(int i=0; i<10; i++)
  {
    fCombinations[i][0] = -1; 
    fCombinations[i][1] = -1;   
  }

  do {
    if(Counter>9){
     AliFatal("AliAnalysisTaskNanoFemtoProtonPion::InitializeArrays: Max. 10 entries of combinations allowed"); 
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
    PassedCombinations++; 
   } while (issCombination);


  //-----------------------------------------------------------------------
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
      AliFatal("AliAnalysisTaskNanoFemtoProtonPion::InitializeArrays: Passed NameTag has to have exaclty 4 entries!");
   }

  //-----------------------------------------------------------------------
  for(int i=0; i<10; i++)
  {
    fClosePairRejection[i] = false;   
  }

  std::istringstream issClosePair(fClosePairRejectionInput);
  Counter = 0;

  while (issClosePair) {
    if(Counter>9){
     AliFatal("AliAnalysisTaskNanoFemtoProtonPion::InitializeArrays: Max. 10 entries of combinations allowed"); 
     break; 
    }

    std::string subs;
    issClosePair >> subs;

    if(subs.size() == 0)
       continue;

    if(subs == "true"){
       fClosePairRejection[Counter] = true;
    } else if (subs == "false"){
       fClosePairRejection[Counter] = false;
    } else {
       AliFatal(Form("AliAnalysisTaskNanoFemtoProtonPion::InitializeArrays: Cannot read bool %s", subs.data()));
    }
  
    Counter++; 
   } 

   if(Counter!=PassedCombinations){
      AliFatal("AliAnalysisTaskNanoFemtoProtonPion::InitializeArrays: Passed Combinations do not match passed ClosePairRejectiosn!");
   }

} // void AliAnalysisTaskNanoFemtoProtonPion::InitializeArrays()

//==================================================================================================================================================
