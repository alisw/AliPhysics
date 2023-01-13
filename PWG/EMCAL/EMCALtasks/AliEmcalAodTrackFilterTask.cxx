//
// Class to filter Aod tracks
//
// Author: C.Loizides
// Updated for ROOT6 by C.Pliatskas in 2022

#include "AliEmcalAodTrackFilterTask.h"
#include <TClonesArray.h>
#include <TRandom3.h>
#include <AliAODEvent.h>
#include <AliAODTrack.h>
#include "AliParticleContainer.h"
#include <AliAnalysisManager.h>
#include <AliEMCALRecoUtils.h>
#include <AliLog.h>
#include <AliVEventHandler.h>
#include "AliEmcalContainerUtils.h"
ClassImp(AliEmcalAodTrackFilterTask)

//________________________________________________________________________
AliEmcalAodTrackFilterTask::AliEmcalAodTrackFilterTask() : 
//  AliAnalysisTaskSE("AliEmcalAodTrackFilterTask"),
 AliAnalysisTaskEmcal("AliEmcalAodTrackFilterTask"),
  fTracksOutName("PicoTracks"),
  fTracksInName("tracks"),
  fIncludeNoITS(kTRUE),
  fCutMaxFrShTPCClus(0),
  fUseNegativeLabels(kTRUE),
  fIsMC(kFALSE),
  fDoPropagation(kFALSE),
  fAttemptProp(kFALSE),
  fAttemptPropMatch(kFALSE),
  fKeepInvMassTag(kFALSE),
  fDist(440),
  fTrackEfficiency(0),
  fTracksIn(0),
  fTracksOut(0)
{
  // Constructor.

  fAODfilterBits[0] = -1;
  fAODfilterBits[1] = -1;
}

//________________________________________________________________________
AliEmcalAodTrackFilterTask::AliEmcalAodTrackFilterTask(const char *name) : 
//  AliAnalysisTaskSE(name),
 AliAnalysisTaskEmcal(name),
  fTracksOutName("PicoTracks"),
  fTracksInName("tracks"),
  fIncludeNoITS(kTRUE),
  fCutMaxFrShTPCClus(0),
  fUseNegativeLabels(kTRUE),
  fIsMC(kFALSE),
  fDoPropagation(kFALSE),
  fAttemptProp(kFALSE),
  fAttemptPropMatch(kFALSE),
  fKeepInvMassTag(kFALSE),
  fDist(440),
  fTrackEfficiency(0),
  fTracksIn(0),
  fTracksOut(0)
{
  // Constructor.

  fAODfilterBits[0] = -1;
  fAODfilterBits[1] = -1;
  fBranchNames = "AOD:tracks";
}

//________________________________________________________________________
AliEmcalAodTrackFilterTask::~AliEmcalAodTrackFilterTask()
{
  // Destructor.
}

//________________________________________________________________________
void AliEmcalAodTrackFilterTask::UserCreateOutputObjects()
{
  // Create my user objects.

  fTracksOut = new TClonesArray("AliAODTrack");
  fTracksOut->SetName(fTracksOutName);
}

//________________________________________________________________________
AliEmcalAodTrackFilterTask* AliEmcalAodTrackFilterTask::AddTaskEmcalAodTrackFilter(
                                    const char *name,
                                    const char *inname,
                                    const char *runperiod, 
                                    const char *taskName,
                                    const char *suffix){
    // Get the pointer to the existing analysis manager via the static access method.
    //==============================================================================
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr)
    {
      ::Error("AddTaskAodTrackFilter", "No analysis manager to connect to.");
      return NULL;
    }
  
    // Check the analysis type using the event handlers connected to the analysis manager.
    //==============================================================================
    if (!mgr->GetInputEventHandler())
    {
      ::Error("AddTaskAodTrackFilter", "This task requires an input event handler");
      return NULL;
    }
  
    TString inputDataType = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
    if (inputDataType != "AOD") {
      ::Info("AddTaskEmcalAodTpcTrack", "This task is only needed for AOD analysis. No task added.");
      return NULL;
    }

    //-------------------------------------------------------
    // Init the task and do settings
    //-------------------------------------------------------
  
    AliEmcalAodTrackFilterTask *aodTask = new AliEmcalAodTrackFilterTask(taskName);
    aodTask->SetTracksOutName(name);
    aodTask->SetTracksInName(inname);
    aodTask->SetMC(kFALSE);
    TString trackname(inname);
   if(!trackname.IsNull())aodTask->AddParticleContainer(trackname);

    TString OutTrackName(name);
   if(!OutTrackName.IsNull())aodTask->AddParticleContainer(OutTrackName);

    Bool_t includeNoITS  = kFALSE;
    Bool_t doProp        = kFALSE; //force propagation of all tracks to EMCal
    Bool_t doAttemptProp = kTRUE;  //only propagate the tracks which were not propagated during AOD filtering
    Bool_t isMC          = kFALSE;
  
    TString runPeriod(runperiod);
    runPeriod.ToLower();
    if (runPeriod == "lhc10b" || runPeriod == "lhc10c" || runPeriod == "lhc10d" ||
        runPeriod == "lhc10e" || runPeriod == "lhc10h" ||
        runPeriod == "lhc11h" || runPeriod == "lhc12a" || runPeriod == "lhc12b" ||
        runPeriod == "lhc12c" || runPeriod == "lhc12d" || runPeriod == "lhc12e" ||
        runPeriod == "lhc12f" || runPeriod == "lhc12g" || runPeriod == "lhc12h" ||
        runPeriod == "lhc12i" || runPeriod == "lhc13b" || runPeriod == "lhc13c" ||
        runPeriod == "lhc13d" || runPeriod == "lhc13e" || runPeriod == "lhc13f" ||
        runPeriod == "lhc13g"
        ) {
      aodTask->SetAODfilterBits(256,512); // hybrid tracks
      if (runPeriod == "lhc10b" || runPeriod == "lhc10c" || runPeriod == "lhc10d" || runPeriod == "lhc10e" || runPeriod == "lhc10h") {
        includeNoITS = kTRUE;
      }
    } else if (runPeriod == "lhc10f7a"    || runPeriod == "lhc12a15e"   || runPeriod.Contains("lhc12a17") || runPeriod == "lhc13b4" ||
               runPeriod == "lhc13b4_fix" || runPeriod == "lhc13b4_plus"    || runPeriod.Contains("lhc14a1") || runPeriod.Contains("lhc13b2_efix")
               ) {
      aodTask->SetAODfilterBits(256,512); // hybrid tracks
      isMC = kTRUE;
      if (runPeriod == "lhc10f7a") {
        includeNoITS = kTRUE;
      }
    } else if (runPeriod == "lhc11a" || runPeriod == "lhc10hold") {
      aodTask->SetAODfilterBits(256,16); // hybrid tracks
      includeNoITS = kTRUE;
    } else if(runPeriod == "lhc11d") {
      aodTask->SetAODfilterBits(256,16); // hybrid tracks (MV: not 100% sure)
      includeNoITS = kFALSE;
    } else if (runPeriod.Contains("lhc12a15a") || runPeriod == "lhc12a15f" || runPeriod == "lhc12a15g") {
      aodTask->SetAODfilterBits(256,16); // hybrid tracks
      isMC = kTRUE;
      includeNoITS = kTRUE;
    } else if (runPeriod.Contains("lhc11a1")){
      aodTask->SetAODfilterBits(256, 16);
      isMC = kTRUE;
      includeNoITS=kTRUE;
    }  else if (runPeriod.Contains(":")) {
      TString runPeriodToken(runperiod);
      TObjArray *arr = runPeriodToken.Tokenize(":");
      TString arg1(arr->At(0)->GetName());

      TString arg2("-1");
      if (arr->GetEntries()>1)
        arg2 = arr->At(1)->GetName();
      if (arr->GetEntries()>2) {
        TString arg3 = arr->At(2)->GetName();
        if (arg3.Contains("includeNoITS=kTRUE"))
          includeNoITS=kTRUE;
        if (arg3.Contains("doProp=kTRUE"))
          doProp=kTRUE;
      if (arg3.Contains("doAttemptProp=kFALSE"))
          doAttemptProp=kFALSE;
        if (arg3.Contains("isMC=kTRUE"))
          isMC = kTRUE;
      }
      aodTask->SetAODfilterBits(arg1.Atoi(),arg2.Atoi());
      delete arr;
    } else {
  
    if (!runPeriod.IsNull())
        ::Warning("AddTaskEmcalAodTrackFilter", Form("Run period %s not known. It will use IsHybridGlobalConstrainedGlobal.", runPeriod.Data()));
    aodTask->SetIncludeNoITS(includeNoITS);
    aodTask->SetDoPropagation(doProp);
    aodTask->SetAttemptProp(doAttemptProp);
    aodTask->SetMC(isMC);
  
    //-------------------------------------------------------
    // Final settings, pass to manager and set the containers
    //-------------------------------------------------------
    mgr->AddTask(aodTask);
  
    // Create containers for input/output
    AliAnalysisDataContainer *cinput1 = mgr->GetCommonInputContainer();
    mgr->ConnectInput(aodTask, 0,  cinput1 );
  
    return aodTask;
}}
void AliEmcalAodTrackFilterTask::ExecOnce(){

  AliParticleContainer* OutputCont= static_cast<AliParticleContainer*>(fParticleCollArray.Last());
    // add tracks to event if not yet there
    fTracksOut->Delete();
    if(OutputCont->GetIsEmbedding()){
          AliVEvent *event = AliEmcalContainerUtils::GetEvent(InputEvent(),kTRUE);
          if(!(event->FindListObject(fTracksOutName))){
          event->AddObject(fTracksOut);
          OutputCont->SetArray(event);}
          }
    else if (!(InputEvent()->FindListObject(fTracksOutName))) {
      InputEvent()->AddObject(fTracksOut);
          OutputCont->SetArray(InputEvent());
    }
   AliAnalysisTaskEmcal::ExecOnce();
}

Bool_t AliEmcalAodTrackFilterTask::Run(){

    // Main loop, called for each event.
    AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
    if (!am) {
      AliError("Manager zero, returning");
      return kFALSE;
    }

AliParticleContainer* partCont = 0;
      TIter next(&fParticleCollArray);
      while ((partCont = static_cast<AliParticleContainer*>(next()))) {
      Int_t nacc = 0;
      Int_t count =0;
          if(TString(partCont->GetTitle()).Contains(fTracksOutName))continue;
    for(auto part : partCont->accepted()) {
          if (!part) continue;
      AliAODTrack *track = dynamic_cast<AliAODTrack*>(part);
      if (!track) continue;
      Int_t type = -1;
      if (fAODfilterBits[0] < 0) {
        if (track->IsHybridGlobalConstrainedGlobal())
          type = 3;
        else //*not a good track*
          continue;
      } else {
        if (track->TestFilterBit(fAODfilterBits[0])) {
          type = 0;
        } else if (fAODfilterBits[1]>-1 && track->TestFilterBit(fAODfilterBits[1])) {
          if ((track->GetStatus()&AliVTrack::kITSrefit)==0) {
            if (fIncludeNoITS)
              type = 2;
            else
              continue;
          } else {
            type = 1;
          }
        }
        else { //*not a good track*
          continue;
        }
      }
if (fCutMaxFrShTPCClus > 0) {
        Double_t frac = Double_t(track->GetTPCnclsS()) / Double_t(track->GetTPCncls());
        if (frac > fCutMaxFrShTPCClus) {
          continue;
        }
      }
  
      if (fTrackEfficiency) {
        Double_t r = gRandom->Rndm();
        if (fTrackEfficiency->Eval(track->Pt()) < r)
          continue;
      }
  
      AliAODTrack *newt = new ((*fTracksOut)[nacc]) AliAODTrack(*track);
      newt->SetUniqueID(0);
      newt->ResetBit(TObject::kHasUUID);
      newt->ResetBit(TObject::kIsReferenced);
  
      Bool_t propthistrack = kFALSE;
      if (fDoPropagation)
        propthistrack = kTRUE;
      else if (!newt->IsExtrapolatedToEMCAL()) {
        if (fAttemptProp)
          propthistrack = kTRUE;
        else if (fAttemptPropMatch && newt->IsEMCAL())
          propthistrack = kTRUE;
      }
      if (propthistrack)
        AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(newt,fDist);
  
      Int_t label = 0;
      if (fIsMC) {
        if (fUseNegativeLabels)
          label = track->GetLabel();
        else
          label = TMath::Abs(track->GetLabel());
        if (label == 0)
          AliDebug(2,Form("Track %d with label==0", count));
      }
      if(fKeepInvMassTag && !fIsMC && (track->GetLabel() == 1011000 ||
          track->GetLabel() == 1012000 ||
          track->GetLabel() == 1021000 ||
          track->GetLabel() == 1022000 ||
          track->GetLabel() == 1031000 ||
          track->GetLabel() == 1032000)){
        newt->SetLabel(track->GetLabel());
      }
      else
        newt->SetLabel(label);
      if (type==0) {
        newt->SetBit(BIT(22),0);
        newt->SetBit(BIT(23),0);
      } else if (type==1) {
        newt->SetBit(BIT(22),1);
        newt->SetBit(BIT(23),0);
      } else if (type==2) {
        newt->SetBit(BIT(22),0);
        newt->SetBit(BIT(23),1);
      } else if (type==3) {
        newt->SetBit(BIT(22),1);
        newt->SetBit(BIT(23),1);
      }
      ++nacc;++count;
    }
  }
return kTRUE;
}

/*
//________________________________________________________________________
void AliEmcalAodTrackFilterTask::UserExec(Option_t *) 
{
  // Main loop, called for each event.
  AliAnalysisManager *am = AliAnalysisManager::GetAnalysisManager();
  if (!am) {
    AliError("Manager zero, returning");
    return;
  }

  // retrieve tracks from input.
  if (!fTracksIn) { 
    fTracksIn = dynamic_cast<TClonesArray*>(InputEvent()->FindListObject(fTracksInName));
    if (!fTracksIn) {
      AliError(Form("Could not retrieve tracks %s!", fTracksInName.Data())); 
      return;
    }
    if (!fTracksIn->GetClass()->GetBaseClass("AliVParticle")) {
      AliError(Form("%s: Collection %s does not contain AliVParticle objects!", GetName(), fTracksInName.Data())); 
      return;
    }
  }

AliParticleContainer* OutputCont= static_cast<AliParticleContainer*>(fParticleCollArray.Last());
  // add tracks to event if not yet there
  fTracksOut->Delete();
  if(OutputCont->GetIsEmbedding()){
        AliVEvent *event = AliEmcalContainerUtils::GetEvent(InputEvent(),kTRUE);
        if(!(event->FindListObject(fTracksOutName))){
        event->AddObject(fTracksOut);
        OutputCont->SetArray(event);}
        }
  else if (!(InputEvent()->FindListObject(fTracksOutName))) {
    InputEvent()->AddObject(fTracksOut);
        OutputCont->SetArray(InputEvent());
  }
 AliAnalysisTaskEmcal::ExecOnce();
  // loop over tracks
AliParticleContainer* partCont = 0;
    TIter next(&fParticleCollArray);
    while ((partCont = static_cast<AliParticleContainer*>(next()))) {
    Int_t nacc = 0;
    Int_t count =0;
        if(TString(partCont->GetTitle()).Contains(fTracksOutName))continue;
  for(auto part : partCont->accepted()) {
        if (!part) continue;
    AliAODTrack *track = dynamic_cast<AliAODTrack*>(part);
    if (!track) continue;
    Int_t type = -1;
    if (fAODfilterBits[0] < 0) {
      if (track->IsHybridGlobalConstrainedGlobal())
        type = 3;
      else /*not a good track
        continue;
    } else {
      if (track->TestFilterBit(fAODfilterBits[0])) {
        type = 0;
      } else if (fAODfilterBits[1]>-1 && track->TestFilterBit(fAODfilterBits[1])) {
        if ((track->GetStatus()&AliVTrack::kITSrefit)==0) {
          if (fIncludeNoITS)
            type = 2;
          else
            continue;
        } else {
          type = 1;
        }
      }
      else {/*not a good track
        continue;
      }
    }

    if (fCutMaxFrShTPCClus > 0) {
      Double_t frac = Double_t(track->GetTPCnclsS()) / Double_t(track->GetTPCncls());
      if (frac > fCutMaxFrShTPCClus) {
        continue;
      }
    }

    if (fTrackEfficiency) {
      Double_t r = gRandom->Rndm();
      if (fTrackEfficiency->Eval(track->Pt()) < r)
        continue;
    }

    AliAODTrack *newt = new ((*fTracksOut)[nacc]) AliAODTrack(*track);
    newt->SetUniqueID(0);
    newt->ResetBit(TObject::kHasUUID);
    newt->ResetBit(TObject::kIsReferenced);

    Bool_t propthistrack = kFALSE;
    if (fDoPropagation)
      propthistrack = kTRUE;
    else if (!newt->IsExtrapolatedToEMCAL()) {
      if (fAttemptProp)
        propthistrack = kTRUE;
      else if (fAttemptPropMatch && newt->IsEMCAL())
        propthistrack = kTRUE;
    }
    if (propthistrack)
      AliEMCALRecoUtils::ExtrapolateTrackToEMCalSurface(newt,fDist);

    Int_t label = 0;
    if (fIsMC) {
      if (fUseNegativeLabels)
        label = track->GetLabel();
      else 
        label = TMath::Abs(track->GetLabel());
      if (label == 0) 
        AliDebug(2,Form("Track %d with label==0", count));
    }
    if(fKeepInvMassTag && !fIsMC && (track->GetLabel() == 1011000 ||
        track->GetLabel() == 1012000 ||
        track->GetLabel() == 1021000 ||
        track->GetLabel() == 1022000 ||
        track->GetLabel() == 1031000 ||
        track->GetLabel() == 1032000)){
      newt->SetLabel(track->GetLabel());
    }
    else
      newt->SetLabel(label);
    if (type==0) {
      newt->SetBit(BIT(22),0);
      newt->SetBit(BIT(23),0);
    } else if (type==1) {
      newt->SetBit(BIT(22),1);
      newt->SetBit(BIT(23),0);
    } else if (type==2) {
      newt->SetBit(BIT(22),0);
      newt->SetBit(BIT(23),1);
    } else if (type==3) {
      newt->SetBit(BIT(22),1);
      newt->SetBit(BIT(23),1);
    }
    ++nacc;++count;
  }
}

}
*/
