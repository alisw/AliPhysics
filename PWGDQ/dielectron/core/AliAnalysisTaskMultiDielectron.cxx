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
//                        Basic Analysis Task                            //
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include <TChain.h>
#include <TH1D.h>

#include <AliCFContainer.h>
#include <AliInputEventHandler.h>
#include <AliESDInputHandler.h>
#include <AliAODInputHandler.h>
#include <AliAnalysisManager.h>
#include <AliVEvent.h>
#include <AliTriggerAnalysis.h>
#include <AliPIDResponse.h>
#include <AliTPCPIDResponse.h>
#include <AliAnalysisUtils.h>

#include "AliDielectron.h"
#include "AliDielectronHistos.h"
#include "AliDielectronCF.h"
#include "AliDielectronMC.h"
#include "AliDielectronMixingHandler.h"
#include "AliAnalysisTaskMultiDielectron.h"

ClassImp(AliAnalysisTaskMultiDielectron)

//_________________________________________________________________________________
AliAnalysisTaskMultiDielectron::AliAnalysisTaskMultiDielectron() :
  AliAnalysisTaskSE(),
	fPairArray(0x0),
  fListDielectron(),
  fListHistos(),
  fListCF(),
  fSelectPhysics(kFALSE),
  fTriggerMask(AliVEvent::kMB),
  fExcludeTriggerMask(0),
  fTriggerOnV0AND(kFALSE),
  fFiredTrigger(""),
  fFiredExclude(kFALSE),
  fRejectPileup(kFALSE),
  fPileUpRejTool(AliDielectronEventCuts::kSPD),
  fBeamEnergy(-1.),
  fRandomizeDaughters(kFALSE),
  fTriggerLogic(kAny),
  fTriggerAnalysis(0x0),
  fRequireTRDtrigger(kFALSE),
  fRequireMatchedTrack(kFALSE),
  fEvtVsTrkHistExists(kFALSE),
  fTRDTriggerClass(AliDielectronEventCuts::kSEorQU),
  fEventFilter(0x0),
  fEventStat(0x0),
  fEventStatTRDTrigger(0x0)
{
  //
  // Constructor
  //
}

//_________________________________________________________________________________
AliAnalysisTaskMultiDielectron::AliAnalysisTaskMultiDielectron(const char *name) :
  AliAnalysisTaskSE(name),
  fPairArray(0x0),
  fListDielectron(),
  fListHistos(),
  fListCF(),
  fSelectPhysics(kFALSE),
  fTriggerMask(AliVEvent::kMB),
  fExcludeTriggerMask(0),
  fTriggerOnV0AND(kFALSE),
  fFiredTrigger(""),
  fFiredExclude(kFALSE),
  fRejectPileup(kFALSE),
  fPileUpRejTool(AliDielectronEventCuts::kSPD),
  fBeamEnergy(-1.),
  fRandomizeDaughters(kFALSE),
  fTriggerLogic(kAny),
  fTriggerAnalysis(0x0),
  fRequireTRDtrigger(kFALSE),
  fRequireMatchedTrack(kFALSE),
  fEvtVsTrkHistExists(kFALSE),
  fTRDTriggerClass(AliDielectronEventCuts::kSEorQU),
  fEventFilter(0x0),
  fEventStat(0x0),
  fEventStatTRDTrigger(0x0)
{
  //
  // Constructor
  //
  DefineInput(0,TChain::Class());
  if((AliAnalysisDataContainer*)(AliAnalysisManager::GetAnalysisManager())->GetContainers()->FindObject("ZDCEPExchangeContainer"))
     DefineInput(1,AliFlowEventSimple::Class());
  DefineOutput(1, TList::Class());
  DefineOutput(2, TList::Class());
  DefineOutput(3, TH1D::Class());
  DefineOutput(4, TH1D::Class());
  fListHistos.SetName("Dielectron_Histos_Multi");
  fListCF.SetName("Dielectron_CF_Multi");
  fListDielectron.SetOwner();
  fListHistos.SetOwner();
  fListCF.SetOwner();
}

//_________________________________________________________________________________
AliAnalysisTaskMultiDielectron::~AliAnalysisTaskMultiDielectron()
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
  if(fEventStatTRDTrigger){ delete fEventStatTRDTrigger;fEventStatTRDTrigger=0; }
  if(fTriggerAnalysis) { delete fTriggerAnalysis; fTriggerAnalysis=0; }
}
//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectron::UserCreateOutputObjects()
{
  //
  // Add all histogram manager histogram lists to the output TList
  //

  // ---| set owner ship of lists |---------------------------------------------
  fListDielectron.SetOwner();
  fListHistos.SetOwner();
  fListCF.SetOwner();

  // ---| skip list initialisation if already done |----------------------------
  if (!fListHistos.IsEmpty()||!fListCF.IsEmpty()) return;

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
  TString name = "hEventStat";
  if (!fEventStat){
    if(fRequireTRDtrigger){
      switch(fTRDTriggerClass) {
        case AliDielectronEventCuts::kSE:
          name= "hEventStatSE";
          break;
        case AliDielectronEventCuts::kQU:
          name= "hEventStatQU";
          break;
        case AliDielectronEventCuts::kSEorQU:
          name= "hEventStatSEorQU";
          break;
        case AliDielectronEventCuts::kSEandQU:
          name= "hEventStatSEandQU";
          break;
      }
    }
    fEventStat=new TH1D(name.Data(),"Event statistics",nbins,0,nbins);
    fEventStat->GetXaxis()->SetBinLabel(1,"Before Phys. Sel.");
    fEventStat->GetXaxis()->SetBinLabel(2,"After Phys. Sel.");

    //default names
    fEventStat->GetXaxis()->SetBinLabel(3,"Bin3 not used");
    fEventStat->GetXaxis()->SetBinLabel(4,"Bin4 not used");
    fEventStat->GetXaxis()->SetBinLabel(5,"Bin5 not used");
    fEventStat->GetXaxis()->SetBinLabel(6,"Bin6 not used");
    fEventStat->GetXaxis()->SetBinLabel(7,"Bin7 not used");

    if(fTriggerOnV0AND) fEventStat->GetXaxis()->SetBinLabel(3,"V0and triggers");
    if(fRequireTRDtrigger) fEventStat->GetXaxis()->SetBinLabel(4,"TRD triggered");
    if(fRequireMatchedTrack) fEventStat->GetXaxis()->SetBinLabel(5,"track matched");
    if (fEventFilter) fEventStat->GetXaxis()->SetBinLabel(6,"After Event Filter");
    if (fRejectPileup) fEventStat->GetXaxis()->SetBinLabel(7,"After Pileup rejection");

    for (Int_t i=0; i<cuts; ++i){
      fEventStat->GetXaxis()->SetBinLabel((kNbinsEvent+1)+2*i,Form("#splitline{1 candidate}{%s}",fListDielectron.At(i)->GetName()));
      fEventStat->GetXaxis()->SetBinLabel((kNbinsEvent+2)+2*i,Form("#splitline{With >1 candidate}{%s}",fListDielectron.At(i)->GetName()));
    }
  }
   if (fRequireTRDtrigger && !fEventStatTRDTrigger){
    fEventStatTRDTrigger=new TH1D("hEventStatTRDTrigger","TRD trigger Statisitcs",10,0.,10.);
    fEventStatTRDTrigger->GetXaxis()->SetBinLabel(1,"No trigger");
    fEventStatTRDTrigger->GetXaxis()->SetBinLabel(2,"only SE fired, not matched");
    fEventStatTRDTrigger->GetXaxis()->SetBinLabel(3,"only SE fired, matched");
    fEventStatTRDTrigger->GetXaxis()->SetBinLabel(4,"only QU fired, not matched");
    fEventStatTRDTrigger->GetXaxis()->SetBinLabel(5,"only QU fired, matched");
    fEventStatTRDTrigger->GetXaxis()->SetBinLabel(6,"both fired, not matched");
    fEventStatTRDTrigger->GetXaxis()->SetBinLabel(7,"both fired, matched");
    fEventStatTRDTrigger->GetXaxis()->SetBinLabel(8,"both fired, SE not matched, QU matched");
    fEventStatTRDTrigger->GetXaxis()->SetBinLabel(9,"both fired, SE matched, QU not matched");
    fEventStatTRDTrigger->GetXaxis()->SetBinLabel(10,"error");
  }


  if (!fTriggerAnalysis) fTriggerAnalysis=new AliTriggerAnalysis;
  fTriggerAnalysis->EnableHistograms();
  fTriggerAnalysis->SetAnalyzeMC(AliDielectronMC::Instance()->HasMC());

  PostData(1, &fListHistos);
  PostData(2, &fListCF);
  PostData(3, fEventStat);
  if(fRequireTRDtrigger){
    PostData(4, fEventStatTRDTrigger);
  }
}

//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectron::UserExec(Option_t *)
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
  if( inputHandler){
    if(fSelectPhysics){
      if((isESD && inputHandler->GetEventSelection()) || isAOD){
        isSelected = inputHandler->IsEventSelected();
        if (fExcludeTriggerMask && (isSelected&fExcludeTriggerMask)) isRejected=kTRUE;
        if (fTriggerLogic==kAny) isSelected&=fTriggerMask;
        else if (fTriggerLogic==kExact) isSelected=((isSelected&fTriggerMask)==fTriggerMask);
      }
    }
    if( isSelected && !fFiredTrigger.IsNull() ){
      TString firedTriggerClasses=InputEvent()->GetFiredTriggerClasses();
      isSelected=(firedTriggerClasses.Contains(fFiredTrigger))^fFiredExclude;
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

  // TRD trigger
  if(fRequireTRDtrigger){
    Bool_t trackMatched;
    Int_t bin;
    Bool_t trdTriggered = AliDielectronEventCuts::IsTRDTriggerFired( InputEvent(), fTRDTriggerClass, trackMatched, bin );
    fEventStatTRDTrigger->Fill(bin);
    if(!trdTriggered){
      return;
    }
    fEventStat->Fill( kTrdTriggeredEvents);
    if(fRequireMatchedTrack && !trackMatched){
        return;
    }
    fEventStat->Fill(kTrdTriggeredEventsMatched);
  }

  //event filter
  if (fEventFilter) {
    if (!fEventFilter->IsSelected(InputEvent())) return;
  }
  fEventStat->Fill(kFilteredEvents);

  //pileup
  if(fRejectPileup){
    switch(fPileUpRejTool) {
      case AliDielectronEventCuts::kSPD:
        if (InputEvent()->IsPileupFromSPD(3,0.8,3.,2.,5.)) return;
        break;
      case AliDielectronEventCuts::kSPDInMultBins:
        if ( (InputEvent()->GetMultiplicity()) && (InputEvent()->IsPileupFromSPDInMultBins() == kTRUE) ) return;
        break;
      case AliDielectronEventCuts::kMultiVertexer:
        AliAnalysisUtils utils;
        if (utils.IsPileUpMV(InputEvent())) return;
        break;
    }
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
    if(die->DoEventProcess()) {
      sel= die->Process(InputEvent());
      // input for internal train
      if(die->DontClearArrays()) {
        fPairArray = (*(die->GetPairArraysPointer())); // the pair arrays from the current 'die' object are stored so they can be used by the next one(s). saves computing time from pairing.
      }
    }
    else {
      // internal train
      if(sel) die->Process(fPairArray);
    }

    if (die->HasCandidates()){
      Int_t ncandidates=die->GetPairArray(1)->GetEntriesFast();
      if (ncandidates==1) fEventStat->Fill((kNbinsEvent)+2*idie);
      else if (ncandidates>1) fEventStat->Fill((kNbinsEvent+1)+2*idie);
    }
    ++idie;
  }

  PostData(1, &fListHistos);
  PostData(2, &fListCF);
  PostData(3,fEventStat);
  if(fRequireTRDtrigger){
    PostData(4,fEventStatTRDTrigger);
  }

}

//_________________________________________________________________________________
void AliAnalysisTaskMultiDielectron::FinishTaskOutput()
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

    // Event vs Track Histo Class
    if(fEvtVsTrkHistExists) die->FinishEvtVsTrkHistoClass();

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
