/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
*                                                                        *
* Author: The ALICE Off-line Project.                                    *
* Contributors are mentioned in the code where appropriate.              *
*                                                                        *
* Permission to use, copy, modify and distribute this software and its   *
* documentation strictly for non-commercial purposes is hereby granted   *
* without fee, provided that the above copyright notice appears in all   *
* copies and that both the copyright notice and this permission notice   *
* appear in the supporting documentation. The authors make no claims     *
* about the suitability of this software for any purpose. It is          *
* provided "as is" without express or implied warranty.                  *
**************************************************************************/

///////////////////////////////////////////////////////////////////////////
//                                                                       //
//     Analysis Task for Random Rejection, i.e. Prefilter Efficiency     //
//       based on AliAnalysisTaskMultiDielectron (04.2015)               //
///////////////////////////////////////////////////////////////////////////

#include <TChain.h>
#include <TH1D.h>
#include <TRandom3.h> // added

#include <AliCFContainer.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliVEvent.h>
#include <AliTriggerAnalysis.h>
#include <AliPIDResponse.h>
#include <AliTPCPIDResponse.h>

#include "AliDielectron.h"
#include "AliDielectronHistos.h"
#include "AliDielectronCF.h"
#include "AliDielectronMC.h"
#include "AliDielectronMixingHandler.h"
#include "AliAnalysisTaskRandomRejection.h"

ClassImp(AliAnalysisTaskRandomRejection)

//_________________________________________________________________________________
AliAnalysisTaskRandomRejection::AliAnalysisTaskRandomRejection() :
  AliAnalysisTaskMultiDielectron(),
  fPtFunc(0x0),
  fPtExpr("exp(-x/3.)"),
  fRndmPtMin(0.2),
  fRndmPtMax(10.),
  fRndmEtaMax(0.9),
  fNTestpartPerEle(200),
  fTestparticles(0x0)//,
{
  //
  // Constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskRandomRejection::AliAnalysisTaskRandomRejection(const char *name) :
  AliAnalysisTaskMultiDielectron(name),
  fPtFunc(0x0),
  fPtExpr("exp(-x/3.)"),
  fRndmPtMin(0.2),
  fRndmPtMax(10.),
  fRndmEtaMax(0.9),
  fNTestpartPerEle(200),
  fTestparticles(0x0)//,
{
  //
  // Constructor
  //
  DefineInput(0,TChain::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TH1D::Class());
  fListHistos.SetName("Dielectron_Histos_Multi");
  fListCF.SetName("Dielectron_CF_Multi");
  fListDielectron.SetOwner();
  fListHistos.SetOwner();
  fListCF.SetOwner();
}

//_________________________________________________________________________________
AliAnalysisTaskRandomRejection::~AliAnalysisTaskRandomRejection()
{
  //
  // Destructor
  //

  //histograms and CF are owned by the dielectron framework.
  //however they are streamed to file, so in the first place the
  //lists need to be owner...
  fListHistos.SetOwner(kFALSE);
  fListCF.SetOwner(kFALSE);

  //  if(fPairArray)       { delete fPairArray;       fPairArray=0; }
  // try to reduce memory issues
  if(fEventStat)       { delete fEventStat;       fEventStat=0; }
  if(fTriggerAnalysis) { delete fTriggerAnalysis; fTriggerAnalysis=0; }
  ///
  /// _____ extension compared to AliAnalysisTaskMultiDielectron _____
  if(fPtFunc)          { delete fPtFunc;          fPtFunc=0; }
  //like done in 'AliDielectron::ClearArrays()':
  fFinalTracks[0].Clear();
  fFinalTracks[1].Clear(); // or Delete()? ("Remove all objects from the array AND delete all heap based objects.")
  if(fTestparticles) fTestparticles->Delete();
  /// _____
}

//_________________________________________________________________________________
void AliAnalysisTaskRandomRejection::UserCreateOutputObjects()
{
  //
  // Add all histogram manager histogram lists to the output TList
  //

  if (!fListHistos.IsEmpty()||!fListCF.IsEmpty()) return; //already initialised

//   AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
//   Bool_t isESD=man->GetInputEventHandler()->IsA()==AliESDInputHandler::Class();
//   Bool_t isAOD=man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  
  TIter nextDie(&fListDielectron);
  AliDielectron *die=0;
  while ( (die=static_cast<AliDielectron*>(nextDie())) ){
    die->Init();
    if (die->GetHistogramList())    fListHistos.Add(const_cast<THashList*>(die->GetHistogramList()));
    if (die->GetHistogramArray())   fListHistos.Add(const_cast<TObjArray*>(die->GetHistogramArray()));
    if (die->GetQAHistArray())      fListHistos.Add(const_cast<TObjArray*>(die->GetQAHistArray()));
    if (die->GetCFManagerPair())    fListCF.Add(const_cast<AliCFContainer*>(die->GetCFManagerPair()->GetContainer()));
  }

  Int_t cuts=fListDielectron.GetEntries();
  Int_t nbins=kNbinsEvent+2*cuts;
  if (!fEventStat){
    fEventStat=new TH1D("hEventStat_RandomRejection","Event statistics",nbins,0,nbins);
    fEventStat->GetXaxis()->SetBinLabel(1,"Before Phys. Sel.");
    fEventStat->GetXaxis()->SetBinLabel(2,"After Phys. Sel.");

    //default names
    fEventStat->GetXaxis()->SetBinLabel(3,"Bin3 not used");
    fEventStat->GetXaxis()->SetBinLabel(4,"Bin4 not used");
    fEventStat->GetXaxis()->SetBinLabel(5,"Bin5 not used");
    
    if(fTriggerOnV0AND) fEventStat->GetXaxis()->SetBinLabel(3,"V0and triggers");
    if (fEventFilter) fEventStat->GetXaxis()->SetBinLabel(4,"After Event Filter");
    if (fRejectPileup) fEventStat->GetXaxis()->SetBinLabel(5,"After Pileup rejection");
    
    for (Int_t i=0; i<cuts; ++i){
      /// _____ modification/extension compared to AliAnalysisTaskMultiDielectron _____
      fEventStat->GetXaxis()->SetBinLabel((kNbinsEvent+1)+2*i,Form("#splitline{w/ prefilter ele}{%s}",fListDielectron.At(i)->GetName()));
      fEventStat->GetXaxis()->SetBinLabel((kNbinsEvent+2)+2*i,Form("#splitline{+ w/ final ele}{%s}",fListDielectron.At(i)->GetName()));
      /// _____
    }
  }

  if (!fTriggerAnalysis) fTriggerAnalysis=new AliTriggerAnalysis;
  fTriggerAnalysis->EnableHistograms();
  fTriggerAnalysis->SetAnalyzeMC(AliDielectronMC::Instance()->HasMC());
  
  PostData(1, &fListHistos);
  PostData(2, &fListCF);
  PostData(3, fEventStat);
}

//_________________________________________________________________________________
void AliAnalysisTaskRandomRejection::UserExec(Option_t *)
{
  //
  // Main loop. Called for every event
  //

  if (fListHistos.IsEmpty()&&fListCF.IsEmpty()) return;

  AliAnalysisManager *man=AliAnalysisManager::GetAnalysisManager();
  Bool_t isESD=man->GetInputEventHandler()->IsA()==AliESDInputHandler::Class();
  Bool_t isAOD=man->GetInputEventHandler()->IsA()==AliAODInputHandler::Class();
  
  AliInputEventHandler* inputHandler = (AliInputEventHandler*) (man->GetInputEventHandler());
  if (!inputHandler) return;
  
//   AliPIDResponse *pidRes=inputHandler->GetPIDResponse();
  if ( inputHandler->GetPIDResponse() ){
    // for the 2.76 pass2 MC private train. Together with a sigma shift of -0.169
//    pidRes->GetTPCResponse().SetSigma(4.637e-3,2.41332105409873257e+04);
    AliDielectronVarManager::SetPIDResponse( inputHandler->GetPIDResponse() );
  } else {
    AliFatal("This task needs the PID response attached to the input event handler!");
  }
  
  // Was event selected ?
  ULong64_t isSelected = AliVEvent::kAny;
  Bool_t isRejected = kFALSE;
  if( fSelectPhysics && inputHandler){
    if((isESD && inputHandler->GetEventSelection()) || isAOD){
      isSelected = inputHandler->IsEventSelected();
      if (fExcludeTriggerMask && (isSelected&fExcludeTriggerMask)) isRejected=kTRUE;
      if (fTriggerLogic==kAny) isSelected&=fTriggerMask;
      else if (fTriggerLogic==kExact) isSelected=((isSelected&fTriggerMask)==fTriggerMask);
   
      TString firedTriggerClasses=InputEvent()->GetFiredTriggerClasses();
      if(!fFiredTrigger.IsNull()) isSelected=(firedTriggerClasses.Contains(fFiredTrigger))^fFiredExclude;
    }
   }
 
 
  //Before physics selection
  fEventStat->Fill(kAllEvents);
  if (isSelected==0||isRejected) {
    PostData(3,fEventStat);
    return;
  }
  //after physics selection
  fEventStat->Fill(kSelectedEvents);

  //V0and
  if(fTriggerOnV0AND){
  if(isESD){if (!fTriggerAnalysis->IsOfflineTriggerFired(static_cast<AliESDEvent*>(InputEvent()), AliTriggerAnalysis::kV0AND))
            return;}
  if(isAOD){if(!((static_cast<AliAODEvent*>(InputEvent()))->GetVZEROData()->GetV0ADecision() == AliVVZERO::kV0BB &&
            (static_cast<AliAODEvent*>(InputEvent()))->GetVZEROData()->GetV0CDecision() == AliVVZERO::kV0BB) )
            return;}
   }
  

  fEventStat->Fill(kV0andEvents);

  //Fill Event histograms before the event filter
  TIter nextDie(&fListDielectron);
  AliDielectron *die=0;
  Bool_t hasMC=AliDielectronMC::Instance()->HasMC();
  while ( (die=static_cast<AliDielectron*>(nextDie())) ){
    AliDielectronHistos *h=die->GetHistoManager();
    if (h){
      AliDielectronVarManager::SetFillMap(h->GetUsedVars());
      if (hasMC && AliDielectronMC::Instance()->ConnectMCEvent() && h->GetHistogramList()->FindObject("MCEvent_noCuts")) {
	AliDielectronVarManager::SetEvent(AliDielectronMC::Instance()->GetMCEvent());
        h->FillClass("MCEvent_noCuts",AliDielectronVarManager::kNMaxValues,AliDielectronVarManager::GetData());
      }
      if (h->GetHistogramList()->FindObject("Event_noCuts")) {
	AliDielectronVarManager::SetEvent(InputEvent());
        h->FillClass("Event_noCuts",AliDielectronVarManager::kNMaxValues,AliDielectronVarManager::GetData());
      }
    }
  }
  nextDie.Reset();
  
  //event filter
  if (fEventFilter) {
    if (!fEventFilter->IsSelected(InputEvent())) return;
  }
  fEventStat->Fill(kFilteredEvents);
  
  //pileup
  if (fRejectPileup){
    if (InputEvent()->IsPileupFromSPD(3,0.8,3.,2.,5.)) return;
  }
  fEventStat->Fill(kPileupEvents);
  
  //bz for AliKF
  Double_t bz = InputEvent()->GetMagneticField();
  AliKFParticle::SetField( bz );

  AliDielectronPID::SetCorrVal((Double_t)InputEvent()->GetRunNumber());
  AliDielectronPair::SetBeamEnergy(InputEvent(), fBeamEnergy);
  AliDielectronPair::SetRandomizeDaughters(fRandomizeDaughters);
  
  //Process event in all AliDielectron instances
  //   TIter nextDie(&fListDielectron);
  //   AliDielectron *die=0;
  Bool_t sel=kFALSE;
  Int_t idie=0;
  while ( (die=static_cast<AliDielectron*>(nextDie())) ){
    ///
    /// _____ modification/extension compared to AliAnalysisTaskMultiDielectron _____
    ///
    //AliInfo(Form("process idie = %i.", idie));
    
    if (die->GetPairPreFilter().GetCuts()->GetEntries()<1) {
      AliInfo(Form(" die #%d : no pair prefilter found, skip any calculation. (task should be deactivated on a higher level...)", idie));
      continue;
    }
    
    /// We need the (track) arrays later.
    die->SetDontClearArrays(kTRUE); // AliDielectron will take care of cleaning the arrays when the next event is executed.
    
    /// Pairing and rejection need to be switched off, otherwise the track arrays are already filtered!
    die->SetNoPairing(kTRUE);
    die->SetPreFilterAllSigns(kFALSE);
    die->SetPreFilterUnlikeOnly(kFALSE);
    /// Deactivate event mixing
    delete die->GetMixingHandler();
    die->SetMixingHandler(0x0);
    
    // Calling die->Process() returns a boolean if event was selected.
    // It also fills the track arrays after the fTrackFilter cuts, which in case of active prefiltering contain the prefilter electrons.
    // Pay attention if switching 'DoEventProcess' off: in AliDielectron::Init(): if(!fEventProcess) {... // move all track cuts (if any) into pair leg cuts // add pair leg cuts to pair filter }
    // die->SetEventProcess(kFALSE); // only meaningful with the internal train, see 'AliAnalysisTaskMultiDielectron.cxx'
    
    if (die->DoEventProcess()) {
      
      sel = die->Process(InputEvent());
      if (!sel) continue;
      /// Skip further calculation if track arrays are empty, i.e. no electron passed the prefilter cuts.
      if ((die->GetTrackArray(0)->GetEntriesFast()<1) && (die->GetTrackArray(1)->GetEntriesFast()<1)) continue;  //  0: Event1, positive particles //  1: Event1, negative particles
      fEventStat->Fill((kNbinsEvent)+2*idie); // events w/ prefilter ele
      
      FillFinalTrackArrays(InputEvent(), die);
      /// Skip further calculation if event does not contain any electron passing the final analysis cuts.
      if ((fFinalTracks[0].GetEntriesFast()<1) && (fFinalTracks[1].GetEntriesFast()<1)) {
        //printf(" Info: event has prefilter electrons, but no global ones! skip event. \n");
        continue;
      }
      fEventStat->Fill((kNbinsEvent+1)+2*idie); // events also w/ final ele
      
      CalcRandomPairs(die);
      
    }
    
    ++idie;
  }
  
  PostData(1, &fListHistos);
  PostData(2, &fListCF);
  PostData(3,fEventStat);
}


//_________________________________________________________________________________
void AliAnalysisTaskRandomRejection::CalcRandomPairs(AliDielectron* die)
{
  ///
  /// Random pairing for prefilter efficiency.
  /// much inspired by void AliDielectron::PairPreFilter(...)
  ///
  //AliInfo("begin.");
  
  UInt_t selectedMask=(1<<die->GetPairPreFilter().GetCuts()->GetEntries())-1;
  AliDielectronPair candidate;
  candidate.SetKFUsage(kFALSE);
  
  // construct or extend array of testparticles if needed
  Int_t nFinalAnaElePosi = fFinalTracks[0].GetEntriesFast() + fFinalTracks[1].GetEntriesFast();
  Int_t nNeededTestPart  = fNTestpartPerEle * nFinalAnaElePosi;
  InitTestparticles(nNeededTestPart);
  
  TObjArray* arrTracks1 = fTestparticles;
  Int_t      ntrack1    = nNeededTestPart; // may be less than the size of fTestparticles!!!
  
  /// remove a random final analysis electron from one of the prefilter electron arrays. (it may stay in the final analysis electron array)
  /// this avoids a slight multiplicity bias, because the total number of electrons is increased by 1 due to the testparticle.
  TRandom3 rnd;
  rnd.SetSeed(0);
  Int_t whichEleToRemove = Int_t(rnd.Rndm()*nFinalAnaElePosi);
  Int_t fromWhichArray   = (whichEleToRemove<fFinalTracks[0].GetEntriesFast())?0:1; // remove from first array, if it is within the arrays size, otherwise from second array.
  if (fromWhichArray==1) whichEleToRemove -= fFinalTracks[0].GetEntriesFast(); // this is needed to not go out of bounds if it is in second array.
  TObject *eleToRemove=fFinalTracks[fromWhichArray].At(whichEleToRemove);
  TObjArray* arrToReduce = const_cast<TObjArray*>(die->GetTrackArray(fromWhichArray));
  if (arrToReduce->FindObject(eleToRemove)) {
    arrToReduce->AddAt(0x0, arrToReduce->IndexOf(eleToRemove)); // remove the track from the array. (don't delete the track object itself, because it's in the input data rootfile!)
    arrToReduce->Compress(); // compress the array
  } else {
    AliWarning(Form("WARNING: Did not find track to delete from array die->GetTrackArray(%i). Not deleting this track... (this should never happen, but rare cases can be ignored. nFinalAnaElePosi=%i, whichEleToRemove=%i, array sizes: %i, %i.)", fromWhichArray, nFinalAnaElePosi, whichEleToRemove, fFinalTracks[0].GetEntriesFast(), fFinalTracks[1].GetEntriesFast()));
    // return;
  }
  // 'eleToRemove' and 'arrToReduce' are just pointers to existing objects, so they don't need to be deleted.
  
  /// The pairing needs to be done between testparticles and prefilter electrons. Get them from track arrays:
  TObjArray* arrTracks2 = const_cast<TObjArray*>(die->GetTrackArray(0));  //  0: Event1, positive particles
  Int_t      ntrack2    = (*arrTracks2).GetEntriesFast();
  TObjArray* arrTracks3 = const_cast<TObjArray*>(die->GetTrackArray(1));  //  1: Event1, negative particles
  Int_t      ntrack3    = (*arrTracks3).GetEntriesFast();
  //printf(" ntrack1 = %d \t ntrack2 = %d \n", ntrack1, ntrack2);
  
  // flag array for rejected testparticles
  Bool_t* bTracks1 = new Bool_t[ntrack1];
  for (Int_t itrack1=0; itrack1<ntrack1; ++itrack1) bTracks1[itrack1]=kFALSE;
  // flag array for rejected prefilter electrons (+ , just for curiosity)
  Bool_t* bTracks2 = new Bool_t[ntrack2];
  for (Int_t itrack2=0; itrack2<ntrack2; ++itrack2) bTracks2[itrack2]=kFALSE;
  // flag array for rejected prefilter electrons (- , just for curiosity)
  Bool_t* bTracks3 = new Bool_t[ntrack3];
  for (Int_t itrack3=0; itrack3<ntrack3; ++itrack3) bTracks3[itrack3]=kFALSE;
  
  
  Int_t nRejPasses = 2;
  for (Int_t iRP=0; iRP < nRejPasses; ++iRP) {
    TObjArray *arrTracks1RP=arrTracks1;
    TObjArray *arrTracks2RP=arrTracks2;
    Bool_t *bTracks1RP = bTracks1;
    Bool_t *bTracks2RP = bTracks2;
    Int_t fPdgLeg1 = ((static_cast<AliVTrack*>((*arrTracks1RP).UncheckedAt(0)))->Charge() > 0.)?-11:+11;
    Int_t fPdgLeg2 = -11;  //  0: Event1, positive particles
    switch (iRP) {
      case 1:
				arrTracks1RP=arrTracks1; // stays the same
				arrTracks2RP=arrTracks3;
				bTracks1RP = bTracks1; // stays the same
				bTracks2RP = bTracks3;
        fPdgLeg2 = 11;  //  1: Event1, negative particles
				break;
      default: ;//nothing to do
    }
    
    Int_t ntrack1RP=(*arrTracks1RP).GetEntriesFast();
    Int_t ntrack2RP=(*arrTracks2RP).GetEntriesFast();
    
    if (ntrack1RP<nNeededTestPart) {
      AliWarning(Form("WARNING: Size of testparticle array is smaller than needed (%i < %i). Using all available testparticles... (this should never happen, it will give less weight to high multiplicity events.)", ntrack1RP, nNeededTestPart));
      nNeededTestPart = ntrack1RP;
      //return;
    }
    //AliInfo(Form("RP %i: nNeededTestPart = %i", iRP, nNeededTestPart));
    for (Int_t itrack1=0; itrack1<nNeededTestPart; ++itrack1){
      
      for (Int_t itrack2=0; itrack2<ntrack2RP; ++itrack2){
        
        TObject *track1=(*arrTracks1RP).UncheckedAt(itrack1);
        TObject *track2=(*arrTracks2RP).UncheckedAt(itrack2);
        if (!track1 || !track2) continue;
        //create the pair
        candidate.SetTracks(static_cast<AliVTrack*>(track1), fPdgLeg1,
                            static_cast<AliVTrack*>(track2), fPdgLeg2);
        
        candidate.SetType(iRP); // could think of nicer/unique names for these pairs (instead of "ev1+_ev1+", "ev1+_ev1-")
        //candidate.SetLabel(AliDielectronMC::Instance()->GetLabelMotherWithPdg(&candidate,fPdgMother));
        
        //pair cuts
        UInt_t cutMask=die->GetPairPreFilter().IsSelected(&candidate);
        
        Bool_t wasRejected=kFALSE;
        wasRejected = (cutMask==selectedMask);
        // rejection occurs if the prefilter cuts are fulfilled!
        // comment in AliDielectron::PairPreFilter(): // remove all tracks from the Single track arrays that pass the cuts in this filter
        
        if (wasRejected) {
          //set flags for rejected tracks
          bTracks1RP[itrack1]=kTRUE;
          bTracks2RP[itrack2]=kTRUE;
        }
        
        //fill histograms
        // currently we assume same behaviour for e- and e+, so all random pairs are filled into the same histograms.
        FillHistogramsRandomPairs(die, &candidate, wasRejected);
        
      } // pairing nested loop
      
    } // pairing main loop
    
  } // rejection passes
  
  
  //fill track histograms
  FillHistogramsTestpart(die, arrTracks1, bTracks1, nNeededTestPart);
  FillHistogramsDataEle(die, arrTracks2, bTracks2);
  FillHistogramsDataEle(die, arrTracks3, bTracks3);
  
  // clean up
  delete [] bTracks1;
  delete [] bTracks2;
  delete [] bTracks3;
  
  //AliInfo("end.");
}


//_________________________________________________________________________________
void AliAnalysisTaskRandomRejection::FillHistogramsRandomPairs(AliDielectron* die, AliDielectronPair* pair, Bool_t wasRejected)
{
  ///
  /// Fill histograms from random pairing, mainly for curiosity and to double-check cuts.
  /// The histogram classes "Rand_Pair" and "Rand_RejPair" must be available.
  ///
  //printf("FillHistogramsRandomPairs( wasRejected = %s ) \n", wasRejected?"kTRUE":"kFALSE");
  
  AliDielectronHistos *h=die->GetHistoManager();
  if (h) {
    if (h->GetHistogramList()->FindObject("Rand_Pair")) {
      Double_t values[AliDielectronVarManager::kNMaxValues]={0};
      AliDielectronVarManager::SetFillMap(h->GetUsedVars()); // needs to be done before every filling!
      AliDielectronVarManager::Fill(pair, values);
      
      h->FillClass("Rand_Pair",AliDielectronVarManager::kNMaxValues,values);
      if (wasRejected) {
        h->FillClass("Rand_RejPair",AliDielectronVarManager::kNMaxValues,values);
      }
    }
  }
}


//_________________________________________________________________________________
void AliAnalysisTaskRandomRejection::FillHistogramsTestpart(AliDielectron* die, TObjArray* arrTracks1, Bool_t* bTracks1, Int_t nNeededTestPart)
{
  ///
  /// Fill histograms of testparticles.
  /// The histogram classes "Random_Testpart" and "Random_RejTestpart" must be available.
  /// The ratio RejTestpart/Testpart (e.g. in dimensions eta, phi, pt) represents the random rejection probability
  /// of final analysis electrons, which has to be applied track-by-track to the efficiency correction.
  ///
  AliDielectronHistos *h=die->GetHistoManager();
  if (h) {
    if (h->GetHistogramList()->FindObject("Random_Testpart")) {
      Double_t values[AliDielectronVarManager::kNMaxValues]={0};
      AliDielectronVarManager::SetFillMap(h->GetUsedVars());
      
      for (Int_t itrack1=0; itrack1<nNeededTestPart; ++itrack1){ // loop until 'nNeededTestPart', which is the amount of testparticles used for the current event.
        TObject *track1=(*arrTracks1).UncheckedAt(itrack1);
        AliDielectronVarManager::Fill(static_cast<AliVParticle*>(track1), values);
        
        h->FillClass("Random_Testpart",AliDielectronVarManager::kNMaxValues,values);
        if (bTracks1[itrack1]) { // was rejected // 'bTracks1[]' only has the size 'nNeededTestPart'.
          h->FillClass("Random_RejTestpart",AliDielectronVarManager::kNMaxValues,values);
        }
      }
    }
  }
}


//_________________________________________________________________________________
void AliAnalysisTaskRandomRejection::FillHistogramsDataEle(AliDielectron* die, TObjArray* arrTracks1, Bool_t* bTracks1)
{
  ///
  /// Fill histograms of prefilter electron sample, mainly for curiosity and to double-check cuts.
  /// The histogram classes "Random_DataEle" and "Random_RejDataEle" must be available.
  ///
  AliDielectronHistos *h=die->GetHistoManager();
  if (h) {
    if (h->GetHistogramList()->FindObject("Random_DataEle")) {
      Double_t values[AliDielectronVarManager::kNMaxValues]={0};
      AliDielectronVarManager::SetFillMap(h->GetUsedVars());
      
      for (Int_t itrack1=0; itrack1<(*arrTracks1).GetEntriesFast(); ++itrack1){
        TObject *track1=(*arrTracks1).UncheckedAt(itrack1);
        AliDielectronVarManager::Fill(static_cast<AliVParticle*>(track1), values);
        
        h->FillClass("Random_DataEle",AliDielectronVarManager::kNMaxValues,values);
        if (bTracks1[itrack1]) { // was rejected
          h->FillClass("Random_RejDataEle",AliDielectronVarManager::kNMaxValues,values);
        }
      }
    }
  }
}


//_________________________________________________________________________________
void AliAnalysisTaskRandomRejection::InitTestparticles(Int_t nNeededTestPart)
{
  ///
  /// Initialize array of Testparticles for random pairing.
  /// Since they probe the rejection of final analysis electrons, they should cover the corresponding kinematic region.
  /// pt and eta ranges could be extracted from trackcuts, but they may vary between the attached Dielectron objects, so it's a bit tricky.
  /// Instead, please use SetPtRange(Double_t min, Double_t max), SetPtExpr(const char *expr/*="exp(-x/3.)"*/) and SetEtaMax(Double_t val) 
  /// in your AddTask if the defaults don't fit you.
  /// The number of testparticles used per event is proportional to the number of final analysis electrons in this event.
  /// This accounts for multiplicity-dependent rejection <-> correct weighting of events according to their multiplicity.
  /// (nNeededTestPart = 'fNTestpartPerEle*[N final analysis electrons]'. May be modified with SetNTestpartPerEle(Int_t ntest).)
  ///
  //AliInfo("begin.");
  
  Int_t N_newTestpart;
  if (!fTestparticles) { // allows to call this function multiple times to increase the amount of testparticles if more are needed.
    // initially, produce enough testparticles for events with at least 4 final analysis electrons.
    N_newTestpart = TMath::Max(nNeededTestPart, fNTestpartPerEle*4);
    fTestparticles = new TObjArray();
    fTestparticles->SetOwner(kTRUE);
  } else {
    if (nNeededTestPart <= fTestparticles->GetEntriesFast()) return;
    // subsequently, extend the testparticle array if the current event has more analysis electrons than the previous maximum.
    N_newTestpart = nNeededTestPart - fTestparticles->GetEntriesFast();
    if (N_newTestpart<1) {
      AliFatal(Form("need more testparticles, but N_newTestpart<1 (=%i). this should never happen.", N_newTestpart));
      return;
    }
  }
  AliInfo(Form("creating %i new testparticles. (nNeededTestPart = %i, fTestparticles->GetEntriesFast() = %i)", N_newTestpart, nNeededTestPart, fTestparticles->GetEntriesFast()));
  
  //  Double_t ptmin=1.0, ptmax=2.0;
  Double_t pt=-999., eta=-999., phi=-999., theta=-999.;
  TRandom3 rnd;
  rnd.SetSeed(0);
  //gRandom = &rnd; // doesnt work somehow [ produces crash in unrelated place. seems to mess up 'rnd' as well. ]
  gRandom->SetSeed(0); // gRandom is used by fPtFunc->GetRandom();
  
  if (!fPtFunc) fPtFunc = new TF1("fPtFunc", fPtExpr.Data(), fRndmPtMin, fRndmPtMax);
  if (!fPtFunc) {
    AliFatal(Form("could not create function for random pt-distribution with expression: \"%s\", range: %f - %f GeV/c. Using default...", fPtExpr.Data(), fRndmPtMin, fRndmPtMax));
    fPtFunc = new TF1("fPtFunc", "exp(-x/3.)", 0.2, 10.);
  }
  AliInfo(Form("function used for random pt-distribution:  %s, %s, %f - %f GeV/c", fPtFunc->GetName(), fPtFunc->GetTitle(), fPtFunc->GetXmin(), fPtFunc->GetXmax()));
  
  for (int itest=0; itest<N_newTestpart; ++itest) {
    pt = fPtFunc->GetRandom();
    eta = fRndmEtaMax*(rnd.Rndm()*2.-1.);
    phi = TMath::TwoPi()*rnd.Rndm();
    //    pt = ptmin + (ptmax-ptmin)*(itest+0.5)/N_newTestpart; //uniform for debugging
    //    do { pt = rnd.Exp(3.); } while (pt<ptmin); //tau~3 reasonable
    //    eta = fRndmEtaMax*(itest+0.5)/N_newTestpart; //uniform for debugging (only from 0 to fRndmEtaMax!)
    //    phi = TMath::TwoPi()*(itest+0.5)/N_newTestpart; //uniform for debugging
    
    theta = 2.*TMath::ATan( TMath::Exp(-eta) ); // theta is the intrinsic AODTrack variable
    
    /// @TODO: create AODTrack or ESDtrack in consistence with the analysed data type!
    AliAODTrack* testpart = new AliAODTrack();
    testpart->SetPt(pt);        //void SetPt(Double_t pt) { fMomentum[0] = pt; };
    testpart->SetPhi(phi);      //void SetPhi(Double_t phi) { fMomentum[1] = phi; }
    testpart->SetTheta(theta);  //void SetTheta(Double_t theta) { fMomentum[2] = theta; }
    testpart->SetCharge(1);
    testpart->SetPIDForTracking(AliAODTrack::kElectron);
    
    fTestparticles->Add(testpart);
  }
  //printf(" fTestparticles->UncheckedAt(1)->M() = %f \t ->Charge() = %i \n", (static_cast<AliVTrack*>(fTestparticles->UncheckedAt(0)))->M(), (static_cast<AliVTrack*>(fTestparticles->UncheckedAt(0)))->Charge());
  //AliInfo("end.");
}


//_________________________________________________________________________________
void AliAnalysisTaskRandomRejection::FillFinalTrackArrays(AliVEvent * const ev, AliDielectron* die)
{
  ///
  /// Fill track arrays of final analysis electrons.
  /// Their presence is needed to determine if random rejection needs to be tested in this event.
  /// (taken from void AliDielectron::FillTrackArrays(AliVEvent * const ev, Int_t eventNr))
  ///
  //AliInfo("begin.");
  
  fFinalTracks[0].Clear();
  fFinalTracks[1].Clear();
  
  Int_t ntracks=ev->GetNumberOfTracks();
  
  UInt_t selectedMask=(1<<die->GetPairPreFilterLegs().GetCuts()->GetEntries())-1;
  for (Int_t itrack=0; itrack<ntracks; ++itrack){
    //get particle
    AliVParticle *particle=ev->GetTrack(itrack);
    
    //apply track cuts
    UInt_t cutmask=die->GetPairPreFilterLegs().IsSelected(particle);
    if (cutmask!=selectedMask) continue;
    
    //fill selected particle into the corresponding track arrays
    Short_t charge=particle->Charge();
    if (charge>0)      fFinalTracks[0].Add(particle);
    else if (charge<0) fFinalTracks[1].Add(particle);
  }
  //AliInfo("end.");
}


//_________________________________________________________________________________
void AliAnalysisTaskRandomRejection::FinishTaskOutput()
{
  //
  // Write debug tree
  //
  //@TODO: in AliAnalysisTaskMultiDielectron: if event mixing is off, the internal train will not work, or will it?
  //@TODO: in AliAnalysisTaskRandomRejection: this function may not be needed at all in this task?
  
  TIter nextDie(&fListDielectron);
  Int_t ic=0;
  AliDielectron *die=0;
  AliDielectron *die2=0;
  fPairArray=0x0;

  // main loop
  while ( (die=static_cast<AliDielectron*>(nextDie())) ){
    ic++;

    // debug tree
    die->SaveDebugTree();

    // skip internal train tasks in main loop
    if(!die->DoEventProcess()) continue;

    // mix remaining
    AliDielectronMixingHandler *mix=die->GetMixingHandler();
    if (!mix || !mix->GetMixUncomplete()) continue;

    // loop over all pools
    for (Int_t ipool=0; ipool<mix->GetNumberOfBins(); ++ipool){
      //      printf("mix remaining %04d/%04d \n",ipool,mix->GetNumberOfBins());
      if(! mix->MixRemaining(die, ipool) ) { fPairArray=0x0;  continue; }

      fPairArray = (*(die->GetPairArraysPointer()));
      if(!fPairArray) continue;

      // loop over internal train task candidates
      for(Int_t i=ic; i<fListDielectron.GetEntries(); i++) {
	die2 = static_cast<AliDielectron*>(fListDielectron.At(i));
	// abort if tasks following are not internal wagons
	if(die2->DoEventProcess()) break;
	// fill internal train output
	die2->SetPairArraysPointer(fPairArray);
	//	printf(" --> fill internal train output %s \n",die2->GetName());
	die2->FillHistogramsFromPairArray(kTRUE);
      }
      // printf("\n\n\n===============\ncall mix in Terminate: %p (%p)\n=================\n\n",mix,die);

    }

  }

  PostData(1, &fListHistos);
  PostData(2, &fListCF);
}

