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
  fRndmEtaMax(0.9),
  fNRndmPt(8),
  fNRndmEta(8),
  fNRndmPhi(8),
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
  fRndmEtaMax(0.9),
  fNRndmPt(8),
  fNRndmEta(8),
  fNRndmPhi(8),
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
  if(fTestparticles)   { delete fTestparticles;   fTestparticles=0; }
  if(fFinalTracks[1])  { fFinalTracks[1]->Clear(); } // or ->Delete()? ("Remove all objects from the array AND delete all heap based objects.")
  if(fFinalTracks[0])  { fFinalTracks[0]->Clear(); }
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
  
  //Process event in all AliDielectron instances
  //   TIter nextDie(&fListDielectron);
  //   AliDielectron *die=0;
  Bool_t sel=kFALSE;
  Int_t idie=0;
  while ( (die=static_cast<AliDielectron*>(nextDie())) ){
    ///
    /// _____ modification/extension compared to AliAnalysisTaskMultiDielectron _____
    ///
    AliInfo(Form(" **************** start of modified code! die #%d **************************", idie));
    
    if (die->GetPairPreFilter().GetCuts()->GetEntries()<1) {
      AliInfo(Form(" die #%d : no pair prefilter found, skip any calculation. (task should be deactivated on a higher level...)", idie));
      continue;
    }
    
    /// We need the (track) arrays later.
    die->SetDontClearArrays(kTRUE);
    
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
      /// Skip further calculation if event does not contain any electron passing the analysis cuts.
      if ((fFinalTracks[0]->GetEntriesFast()<1) && (fFinalTracks[1]->GetEntriesFast()<1)) {
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
  /// Random pairing for prefilter efficiency
  /// much inspired by void AliDielectron::PairPreFilter(...)
  ///
  printf("CalcRandomPairs() \n");
  
  UInt_t selectedMask=(1<<die->GetPairPreFilter().GetCuts()->GetEntries())-1;
  AliDielectronPair candidate;
  candidate.SetKFUsage(kFALSE);
  
  // construct array of testparticles
  if (!fTestparticles) InitTestparticles();
  TObjArray* arrTracks1 = fTestparticles;
  Int_t      ntrack1    = (*arrTracks1).GetEntriesFast();
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
				arrTracks1RP=arrTracks1;
				arrTracks2RP=arrTracks3;
				bTracks1RP = bTracks1;
				bTracks2RP = bTracks3;
        fPdgLeg2 = 11;  //  1: Event1, negative particles
				break;
      default: ;//nothing to do
    }
    
    Int_t ntrack1RP=(*arrTracks1RP).GetEntriesFast();
    Int_t ntrack2RP=(*arrTracks2RP).GetEntriesFast();
    
    for (Int_t itrack1=0; itrack1<ntrack1RP; ++itrack1){
      
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
  FillHistogramsTestpart(die, arrTracks1, bTracks1);
  FillHistogramsDataEle(die, arrTracks2, bTracks2);
  FillHistogramsDataEle(die, arrTracks3, bTracks3);
  
  // clean up
  delete [] bTracks1;
  delete [] bTracks2;
  delete [] bTracks3;
  
}


//_________________________________________________________________________________
void AliAnalysisTaskRandomRejection::FillHistogramsRandomPairs(AliDielectron* die, AliDielectronPair* pair, Bool_t wasRejected)
{
  ///
  /// Fill histograms from random pairing
  ///
  //printf("FillHistogramsRandomPairs( wasRejected = %s ) \n", wasRejected?"kTRUE":"kFALSE");
  
  AliVParticle *d1=pair->GetFirstDaughterP();
  AliVParticle *d2=pair->GetSecondDaughterP();
  if (!d1 || !d2) {
    printf(" Error: FillHistogramsRandomPairs(): one daughter is not available! return. \n");
    return;
  }
//  Bool_t d1_IsDataEle = ((AliAODTrack*)d1)->GetDetPid()?kTRUE:kFALSE;
//  Bool_t d2_IsDataEle = ((AliAODTrack*)d2)->GetDetPid()?kTRUE:kFALSE;
  
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
void AliAnalysisTaskRandomRejection::FillHistogramsTestpart(AliDielectron* die, TObjArray* arrTracks1, Bool_t* bTracks1)
{
  AliDielectronHistos *h=die->GetHistoManager();
  if (h) {
    if (h->GetHistogramList()->FindObject("Random_Testpart")) {
      Double_t values[AliDielectronVarManager::kNMaxValues]={0};
      AliDielectronVarManager::SetFillMap(h->GetUsedVars());
      
      for (Int_t itrack1=0; itrack1<(*arrTracks1).GetEntriesFast(); ++itrack1){
        TObject *track1=(*arrTracks1).UncheckedAt(itrack1);
        AliDielectronVarManager::Fill(static_cast<AliVParticle*>(track1), values);
        
        h->FillClass("Random_Testpart",AliDielectronVarManager::kNMaxValues,values);
        if (bTracks1[itrack1]) { // was rejected
          h->FillClass("Random_RejTestpart",AliDielectronVarManager::kNMaxValues,values);
        }
      }
    }
  }
}


//_________________________________________________________________________________
void AliAnalysisTaskRandomRejection::FillHistogramsDataEle(AliDielectron* die, TObjArray* arrTracks1, Bool_t* bTracks1)
{
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
void AliAnalysisTaskRandomRejection::InitTestparticles()
{
  ///
  /// Initialize array of Testparticles for random pairing.
  /// Since they probe the rejection of final analysis electrons, they should cover the corresponding kinematic region.
  /// (Could extract pt and eta ranges from trackcuts, but they may vary between the attached Dielectron objects, so it's a bit tricky.
  ///  Instead, use SetPtFunc(TF1* func) and SetEtaMax(Double_t val) in your AddTask. )
  /// The number of testparticles per event is 'fNRndmPt*fNRndmEta*fNRndmPhi'. Modify with SetNPtEtaPhi(UInt_t npt, UInt_t neta, UInt_t nphi).
  ///
  //printf("InitTestparticles() \n");
  
  fTestparticles = new TObjArray();
  fTestparticles->SetOwner(kTRUE);
  
  //  UInt_t fNRndmPt  = 8;
  //  UInt_t fNRndmEta = 8;
  //  UInt_t fNRndmPhi = 8;
  //  Double_t ptmin=1.0, ptmax=2.0;
  Double_t pt=-999., eta=-999., phi=-999., theta=-999.;
  TRandom3 rnd;
  rnd.SetSeed(0);
  //gRandom = &rnd; // doesnt work somehow [ produces crash in unrelated place. seems to mess up 'rnd' as well. ]
  gRandom->SetSeed(0); // gRandom is used by fPtFunc->GetRandom();
  
  if (!fPtFunc) fPtFunc = new TF1("fPtFunc", "exp(-x/3.)", 0.2, 10.);
  
  for (int ipt=0; ipt<fNRndmPt; ++ipt) {
    //pt = fPtFunc->GetRandom(); // better to sample pt also more often.
    //    pt = ptmin + (ptmax-ptmin)*(ipt+0.5)/fNRndmPt; //uniform for debugging
    //    do { pt = rnd.Exp(3.); } while (pt<ptmin); //tau~3 reasonable
    
    for (int ieta=0; ieta<fNRndmEta; ++ieta) {
      pt = fPtFunc->GetRandom(); // better to sample pt also more often.
      eta = fRndmEtaMax*(rnd.Rndm()*2.-1.);
      //eta = fRndmEtaMax*(ieta+0.5)/fNRndmEta; //uniform for debugging (from 0 to fRndmEtaMax)
      
      for (int iphi=0; iphi<fNRndmPhi; ++iphi) {
        phi = TMath::TwoPi()*rnd.Rndm();
        //phi = TMath::TwoPi()*(iphi+0.5)/fNRndmPhi; //uniform for debugging
        
        theta = 2.*TMath::ATan( TMath::Exp(-eta) ); // theta is the intrinsic AODTrack variable
        
        AliAODTrack* testpart = new AliAODTrack();
        testpart->SetPt(pt);        //void SetPt(Double_t pt) { fMomentum[0] = pt; };
        testpart->SetPhi(phi);      //void SetPhi(Double_t phi) { fMomentum[1] = phi; }
        testpart->SetTheta(theta);  //void SetTheta(Double_t theta) { fMomentum[2] = theta; }
        testpart->SetCharge(1);
        testpart->SetPIDForTracking(AliAODTrack::kElectron);
        
        fTestparticles->Add(testpart);
      }
    }
  }
  //printf(" fTestparticles->UncheckedAt(1)->M() = %f \t ->Charge() = %i \n", (static_cast<AliVTrack*>(fTestparticles->UncheckedAt(0)))->M(), (static_cast<AliVTrack*>(fTestparticles->UncheckedAt(0)))->Charge());
}


//_________________________________________________________________________________
void AliAnalysisTaskRandomRejection::FillFinalTrackArrays(AliVEvent * const ev, AliDielectron* die)
{
  ///
  /// Presence of final analysis electrons is needed to determine if random rejection needs to be tested in this event.
  /// (taken from void AliDielectron::FillTrackArrays(AliVEvent * const ev, Int_t eventNr))
  ///
  //printf("FillFinalTrackArrays() \n");
  
  if (!fFinalTracks[0]) {
    fFinalTracks[0] = new TObjArray();
    fFinalTracks[1] = new TObjArray();
  } else {
    fFinalTracks[0]->Clear();
    fFinalTracks[1]->Clear();
  }
  
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
    if (charge>0)      fFinalTracks[0]->Add(particle);
    else if (charge<0) fFinalTracks[1]->Add(particle);
  }
}


//_________________________________________________________________________________
void AliAnalysisTaskRandomRejection::FinishTaskOutput()
{
  //
  // Write debug tree
  //
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

