//
// Creation date: 2015/10/01
// Author: Ionut-Cristian Arsene, iarsene@cern.ch, i.c.arsene@fys.uio.no

#include "AliReducedAnalysisTest.h"

#include <iostream>
using std::cout;
using std::endl;

#include <TClonesArray.h>

#include "AliReducedVarManager.h"
#include "AliReducedEventInfo.h"
#include "AliReducedBaseEvent.h"
#include "AliReducedTrackInfo.h"
#include "AliReducedCaloClusterInfo.h"
#include "AliReducedPairInfo.h"
#include "AliHistogramManager.h"

ClassImp(AliReducedAnalysisTest);


//___________________________________________________________________________
AliReducedAnalysisTest::AliReducedAnalysisTest() :
  AliReducedAnalysisTaskSE(),
  fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fEventCuts(),
  fTrackCuts(),
  fPairCuts(),
  fTrackFilterBitNames(""),
  fMCBitsNames("")
{
  //
  // default constructor
  //
   
}


//___________________________________________________________________________
AliReducedAnalysisTest::AliReducedAnalysisTest(const Char_t* name, const Char_t* title) :
  AliReducedAnalysisTaskSE(name,title),
  fHistosManager(new AliHistogramManager("Histogram Manager", AliReducedVarManager::kNVars)),
  fEventCuts(),
  fTrackCuts(),
  fPairCuts(),
  fTrackFilterBitNames(""),
  fMCBitsNames("")
{
  //
  // named constructor
  //
   fEventCuts.SetOwner(kTRUE);
   fTrackCuts.SetOwner(kTRUE);
   fPairCuts.SetOwner(kTRUE);
}


//___________________________________________________________________________
AliReducedAnalysisTest::~AliReducedAnalysisTest() 
{
  //
  // destructor
  //
   fEventCuts.Clear("C"); fTrackCuts.Clear("C"); fPairCuts.Clear("C");
   if(fHistosManager) delete fHistosManager;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisTest::IsEventSelected(AliReducedBaseEvent* event) {
  //
  // apply event cuts
  //
  if(fEventCuts.GetEntries()==0) return kTRUE;
  // loop over all the cuts and make a logical and between all cuts in the list
  // TODO: Cut masks or more complicated cut configurations can be implemented here
  for(Int_t i=0; i<fEventCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fEventCuts.At(i);
    if(!cut->IsSelected(event)) return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisTest::IsTrackSelected(AliReducedBaseTrack* track) {
  //
  // apply event cuts
  //
  if(fTrackCuts.GetEntries()==0) return kTRUE;
  // loop over all the cuts and make a logical and between all cuts in the list
  // TODO: Cut masks or more complicated cut configurations can be implemented here
  for(Int_t i=0; i<fTrackCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fTrackCuts.At(i);
    if(!cut->IsSelected(track)) return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________________
Bool_t AliReducedAnalysisTest::IsPairSelected(AliReducedBaseTrack* pair) {
  //
  // apply event cuts
  //
  if(fPairCuts.GetEntries()==0) return kTRUE;
  // loop over all the cuts and make a logical and between all cuts in the list
  // TODO: Cut masks or more complicated cut configurations can be implemented here
  for(Int_t i=0; i<fPairCuts.GetEntries(); ++i) {
    AliReducedInfoCut* cut = (AliReducedInfoCut*)fPairCuts.At(i);
    if(!cut->IsSelected(pair)) return kFALSE;
  }
  return kTRUE;
}

//___________________________________________________________________________
void AliReducedAnalysisTest::Init() {
  //
  // initialize stuff
  //
   AliReducedVarManager::SetDefaultVarNames();
   fHistosManager->SetUseDefaultVariableNames(kTRUE);
   fHistosManager->SetDefaultVarNames(AliReducedVarManager::fgVariableNames,AliReducedVarManager::fgVariableUnits);
}


//___________________________________________________________________________
void AliReducedAnalysisTest::Process() {
  //
  // process the current event
  //  
  if(!fEvent) return;
  
  AliReducedVarManager::SetEvent(fEvent);
  
  AliReducedVarManager::FillEventInfo(fEvent, fValues);
  fHistosManager->FillHistClass("Event_NoCuts", fValues);
  if(fEvent->IsA()==AliReducedEventInfo::Class()) {
    for(UShort_t ibit=0; ibit<64; ++ibit) {
      AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
      fHistosManager->FillHistClass("OnlineTriggers_NoCuts", fValues);
      for(UShort_t ibit2=0; ibit2<64; ++ibit2) {
         AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues, ibit2);
         fHistosManager->FillHistClass("TriggerCorrelation_NoCuts", fValues);
      }
    }
  }
  if(!IsEventSelected(fEvent)) return;
  
  AliReducedEventInfo* eventInfo = NULL;
  if(fEvent->IsA()==AliReducedEventInfo::Class()) eventInfo = (AliReducedEventInfo*)fEvent;
  
  if(eventInfo) {
    for(UShort_t ibit=0; ibit<64; ++ibit) {
      AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues);
      fHistosManager->FillHistClass("OnlineTriggers_AfterCuts", fValues);
      for(UShort_t ibit2=0; ibit2<64; ++ibit2) {
         AliReducedVarManager::FillEventOnlineTrigger(ibit, fValues, ibit2);
         fHistosManager->FillHistClass("TriggerCorrelation_AfterCuts", fValues);
      }
      for(UShort_t i=0; i<32; ++i) {
        AliReducedVarManager::FillL0TriggerInputs(eventInfo, i, fValues);
        fHistosManager->FillHistClass("OnlineTriggers_vs_L0TrigInputs", fValues);
      }
      for(UShort_t i=0; i<32; ++i) {
        AliReducedVarManager::FillL1TriggerInputs(eventInfo, i, fValues);
        fHistosManager->FillHistClass("OnlineTriggers_vs_L1TrigInputs", fValues);
      }
      for(UShort_t i=0; i<16; ++i) {
        AliReducedVarManager::FillL2TriggerInputs(eventInfo, i, fValues);
        fHistosManager->FillHistClass("OnlineTriggers_vs_L2TrigInputs", fValues);
      }
    }
  }
    
  for(UShort_t ibit=0; ibit<64; ++ibit) {
    AliReducedVarManager::FillEventTagInput(fEvent, ibit, fValues);
    fHistosManager->FillHistClass("EvtTags", fValues);
  }
  
  if(eventInfo) {
    for(UShort_t ibit=0; ibit<32; ++ibit) {
      AliReducedVarManager::FillL0TriggerInputs(eventInfo, ibit, fValues);
      fHistosManager->FillHistClass("L0TriggerInput", fValues);
      for(UShort_t ibit2=0; ibit2<32; ++ibit2) {
         AliReducedVarManager::FillL0TriggerInputs(eventInfo, ibit, fValues, ibit2);
         fHistosManager->FillHistClass("L0InputCorrelation", fValues);
      }
    }
    for(UShort_t ibit=0; ibit<32; ++ibit) {
      AliReducedVarManager::FillL1TriggerInputs(eventInfo, ibit, fValues);
      fHistosManager->FillHistClass("L1TriggerInput", fValues);
      for(UShort_t ibit2=0; ibit2<32; ++ibit2) {
         AliReducedVarManager::FillL1TriggerInputs(eventInfo, ibit, fValues, ibit2);
         fHistosManager->FillHistClass("L1InputCorrelation", fValues);
      }
    }
    for(UShort_t ibit=0; ibit<16; ++ibit) {
      AliReducedVarManager::FillL2TriggerInputs(eventInfo, ibit, fValues);
      fHistosManager->FillHistClass("L2TriggerInput", fValues);
      for(UShort_t ibit2=0; ibit2<16; ++ibit2) {
         AliReducedVarManager::FillL2TriggerInputs(eventInfo, ibit, fValues, ibit2);
         fHistosManager->FillHistClass("L2InputCorrelation", fValues);
      }
    }
    
    for(Int_t icl=0; icl<eventInfo->GetNCaloClusters(); ++icl) {
      AliReducedVarManager::FillCaloClusterInfo(eventInfo->GetCaloCluster(icl), fValues);
      fHistosManager->FillHistClass("CaloClusters", fValues);
    }
  }
     
  // fill track histograms for the first track array
  FillTrackHistograms(fEvent->GetTracks());
  // fill track histograms for the second track array
  if(fEvent->NTracks2()>0)
     FillTrackHistograms(fEvent->GetTracks2());
  
    
  TClonesArray* pairList = fEvent->GetPairs();
  TIter nextPair(pairList);
  AliReducedPairInfo* pair = 0x0;
  fValues[AliReducedVarManager::kNpairsSelected] = 0.;
  if(pairList) {
    for(Int_t ip=0; ip<pairList->GetEntries(); ++ip) {
      pair = (AliReducedPairInfo*)nextPair();
    
      if(!IsPairSelected(pair)) continue;
    
      TString pairTypeStr = "";
      if(pair->PairType()==0) pairTypeStr = "Offline";
      if(pair->PairType()==1) pairTypeStr = "OnTheFly";
      
      for(UShort_t iflag=0; iflag<32; ++iflag) {
        AliReducedVarManager::FillPairQualityFlag(pair, iflag, fValues);
        fHistosManager->FillHistClass(Form("PairQualityFlags_%s",pairTypeStr.Data()), fValues);
        for(UShort_t iflag2=0; iflag2<32; ++iflag2) {
           AliReducedVarManager::FillPairQualityFlag(pair, iflag, fValues, iflag2);
           fHistosManager->FillHistClass(Form("CorrelationQualityFlagsPairs_%s",pairTypeStr.Data()), fValues);
        }
      }
      TString type[3]={"PP","PM","MM"};
      
      AliReducedVarManager::FillPairInfo(pair, fValues);
      switch (pair->CandidateId()) {
        case AliReducedPairInfo::kGammaConv :
          fHistosManager->FillHistClass(Form("PairQA_%sGamma",pairTypeStr.Data()),fValues);
	  if(pair->IsPureV0Gamma()) fHistosManager->FillHistClass(Form("PairQA_%sPureGamma",pairTypeStr.Data()),fValues);
          break;
        case AliReducedPairInfo::kK0sToPiPi :
	  fHistosManager->FillHistClass(Form("PairQA_%sK0s",pairTypeStr.Data()),fValues);
	  if(pair->IsPureV0K0s()) fHistosManager->FillHistClass(Form("PairQA_%sPureK0s",pairTypeStr.Data()),fValues);
	  break;
        case AliReducedPairInfo::kLambda0ToPPi :
	  fHistosManager->FillHistClass(Form("PairQA_%sLambda", pairTypeStr.Data()),fValues);
	  if(pair->IsPureV0Lambda()) fHistosManager->FillHistClass(Form("PairQA_%sPureLambda",pairTypeStr.Data()),fValues);
	  break;
        case AliReducedPairInfo::kALambda0ToPPi :
	  fHistosManager->FillHistClass(Form("PairQA_%sALambda", pairTypeStr.Data()),fValues);
	  if(pair->IsPureV0ALambda()) fHistosManager->FillHistClass(Form("PairQA_%sPureALambda",pairTypeStr.Data()),fValues);
	  break;
        case AliReducedPairInfo::kJpsiToEE :
           fHistosManager->FillHistClass(Form("PairQA_Jpsi2EE_%s", type[Int_t(pair->PairType())].Data()),fValues);
           break;  
        case AliReducedPairInfo::kADzeroToKplusPiminus :
           fHistosManager->FillHistClass(Form("PairQA_ADzeroToKplusPiminus_%s", type[Int_t(pair->PairType())].Data()),fValues);
           break;     
      };
      fValues[AliReducedVarManager::kNpairsSelected] += 1.0;
    }  // end loop over pairs
  }  // end if(pairList)
    
  fHistosManager->FillHistClass("Event_AfterCuts", fValues);
}

//___________________________________________________________________________
void AliReducedAnalysisTest::FillTrackHistograms(TClonesArray* trackList) {
   //
   // fill tracks histograms
   //
   AliReducedBaseTrack* track = 0x0;
   TIter nextTrack(trackList);
   if(trackList) {
      for(Int_t it=0; it<trackList->GetEntries(); ++it) {
         track = (AliReducedBaseTrack*)nextTrack();
         if(!track) continue;
         // reset track variables
         for(Int_t i=AliReducedVarManager::kNEventVars; i<AliReducedVarManager::kEMCALmatchedEOverP; ++i) fValues[i]=-9999.;
         
         if(fProcessMCInfo && track->GetMCFlags()) {
            AliReducedVarManager::FillTrackInfo(track,fValues);
            TObjArray* namesArr = fMCBitsNames.Tokenize(";");
            for(Int_t iflag = 0; iflag<namesArr->GetEntries(); ++iflag)
               if(track->TestMCFlag(iflag)) fHistosManager->FillHistClass(Form("PureMCqa_%s", namesArr->At(iflag)->GetName()), fValues);
            for(UShort_t iflag=0; iflag<32; ++iflag) {
               AliReducedVarManager::FillTrackMCFlag(track, iflag, fValues);
               fHistosManager->FillHistClass("PureMCflags", fValues);
               for(UShort_t iflag2=0; iflag2<32; ++iflag2) {
                  AliReducedVarManager::FillTrackMCFlag(track, iflag, fValues, iflag2);
                  fHistosManager->FillHistClass("CorrelationMCflags", fValues);
               }
            }
            continue;
         }
         
         if(!IsTrackSelected(track)) continue;
         AliReducedVarManager::FillTrackInfo(track,fValues);
         fHistosManager->FillHistClass("TrackQA_AllTracks", fValues);
         TObjArray* namesBitArr = fTrackFilterBitNames.Tokenize(";");
         for(Int_t iflag = 0; iflag<namesBitArr->GetEntries(); ++iflag)
            if(track->TestQualityFlag(iflag+32)) fHistosManager->FillHistClass(Form("TrackQA_%s", namesBitArr->At(iflag)->GetName()), fValues);
         
         
         AliReducedTrackInfo* trackInfo = NULL;
         if(track->IsA()==AliReducedTrackInfo::Class()) trackInfo = (AliReducedTrackInfo*)track;
         
         if(trackInfo) {
            for(UInt_t iflag=0; iflag<AliReducedVarManager::kNTrackingFlags; ++iflag) {
               AliReducedVarManager::FillTrackingFlag(trackInfo, iflag, fValues);
               fHistosManager->FillHistClass("TrackingFlags", fValues);
            }
         }
         for(UShort_t iflag=0; iflag<64; ++iflag) {
            AliReducedVarManager::FillTrackQualityFlag(track, iflag, fValues);
            fHistosManager->FillHistClass("TrackQualityFlags", fValues);
            for(UShort_t iflag2=0; iflag2<64; ++iflag2) {
               AliReducedVarManager::FillTrackQualityFlag(track, iflag, fValues,iflag2);
               fHistosManager->FillHistClass("CorrelationQualityFlagsTracks", fValues);
            }
         }
         if(trackInfo) {
            for(Int_t iLayer=0; iLayer<6; ++iLayer) {
               AliReducedVarManager::FillITSlayerFlag(trackInfo, iLayer, fValues);
               fHistosManager->FillHistClass("ITSclusterMap", fValues);
            }
            for(Int_t iLayer=0; iLayer<8; ++iLayer) {
               AliReducedVarManager::FillTPCclusterBitFlag(trackInfo, iLayer, fValues);
               fHistosManager->FillHistClass("TPCclusterMap", fValues);
            }
         }
         if(track->IsGammaLeg()) {
            fHistosManager->FillHistClass("TrackQA_GammaLeg", fValues);
            for(UShort_t iflag=0; iflag<64; ++iflag) {
               AliReducedVarManager::FillTrackQualityFlag(track, iflag, fValues);
               fHistosManager->FillHistClass("TrackQualityFlags_GammaLeg", fValues);
            }
         }
         if(track->IsPureGammaLeg()) fHistosManager->FillHistClass("TrackQA_PureGammaLeg", fValues);
         if(track->IsK0sLeg()) {
            fHistosManager->FillHistClass("TrackQA_K0sLeg", fValues);
            for(UShort_t iflag=0; iflag<64; ++iflag) {
               AliReducedVarManager::FillTrackQualityFlag(track, iflag, fValues);
               fHistosManager->FillHistClass("TrackQualityFlags_K0sLeg", fValues);
            }
         }
         if(track->IsPureK0sLeg()) fHistosManager->FillHistClass("TrackQA_PureK0sLeg", fValues);
         if(track->IsLambdaLeg()) {
            if(track->Charge()>0) fHistosManager->FillHistClass("TrackQA_LambdaPosLeg", fValues);
            else fHistosManager->FillHistClass("TrackQA_LambdaNegLeg", fValues);
         }
         if(track->IsPureLambdaLeg()) {
            if(track->Charge()>0) fHistosManager->FillHistClass("TrackQA_PureLambdaPosLeg", fValues);
            else fHistosManager->FillHistClass("TrackQA_PureLambdaNegLeg", fValues);
         }
         if(track->IsALambdaLeg()) {
            if(track->Charge()>0) fHistosManager->FillHistClass("TrackQA_ALambdaPosLeg", fValues);
            else fHistosManager->FillHistClass("TrackQA_ALambdaNegLeg", fValues);
         }
         if(track->IsPureALambdaLeg()) {
            if(track->Charge()>0) fHistosManager->FillHistClass("TrackQA_PureALambdaPosLeg", fValues);
            else fHistosManager->FillHistClass("TrackQA_PureALambdaNegLeg", fValues);
         }  
      }  // end loop over tracks
   }  // end if(trackList)  
}

//___________________________________________________________________________
void AliReducedAnalysisTest::Finish() {
  //
  // run stuff after the event loop
  //
}
