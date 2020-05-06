/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 	*
 *																			*
 * Author: Daniel Mühlheim													*
 * Version 1.0																*
 *																			*
 * Permission to use, copy, modify and distribute this software and its	 	*
 * documentation strictly for non-commercial purposes is hereby granted	 	*
 * without fee, provided that the above copyright notice appears in all	 	*
 * copies and that both the copyright notice and this permission notice	 	*
 * appear in the supporting documentation. The authors make no claims		*
 * about the suitability of this software for any purpose. It is			*
 * provided "as is" without express or implied warranty.					*
 **************************************************************************/

////////////////////////////////////////////////
//---------------------------------------------
// Basic Track Matching Class
//---------------------------------------------
////////////////////////////////////////////////


#include "AliAnalysisManager.h"
#include "AliAODEvent.h"
#include "AliCaloTriggerMimickHelper.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "TChain.h"

class iostream;

using namespace std;


ClassImp(AliCaloTriggerMimickHelper)

//________________________________________________________________________
AliCaloTriggerMimickHelper::AliCaloTriggerMimickHelper(const char *name, Int_t clusterType, Bool_t isMC) : AliAnalysisTaskSE(name),
  fClusterType(clusterType),
  fRunNumber(-1),
  fPHOSTrigUtils(0x0),
  fPHOSTrigger(kPHOSAny),
  fForceRun(kFALSE),
  fIsMC(isMC),
  fEventChosenByTrigger(kFALSE),
  fDoLightOutput(kFALSE)
{
    cout<<"Debug Output; AliCaloTriggerMimickHelper.C, Line: "<<__LINE__<<endl;
    // Default constructor
    DefineInput(0, TChain::Class());
}

//________________________________________________________________________
AliCaloTriggerMimickHelper::AliCaloTriggerMimickHelper(const char *name, Int_t clusterType, Int_t isMC) : AliAnalysisTaskSE(name),
  fClusterType(clusterType),
  fRunNumber(-1),
  fPHOSTrigUtils(0x0),
  fPHOSTrigger(kPHOSAny),
  fForceRun(kFALSE),
  fIsMC(isMC),
  fEventChosenByTrigger(kFALSE),
  fDoLightOutput(kFALSE)
{
    cout<<"Debug Output; AliCaloTriggerMimickHelper.C, Line: "<<__LINE__<<endl;
    // Default constructor
    DefineInput(0, TChain::Class());
}

//________________________________________________________________________
AliCaloTriggerMimickHelper::~AliCaloTriggerMimickHelper(){
    cout<<"Debug Output; AliCaloTriggerMimickHelper.C, Line: "<<__LINE__<<endl;
    // default deconstructor
}

//________________________________________________________________________
void AliCaloTriggerMimickHelper::Terminate(Option_t *){

}

//________________________________________________________________________
void AliCaloTriggerMimickHelper::UserCreateOutputObjects(){
    cout<<"Debug Output; AliCaloTriggerMimickHelper.C, Line: "<<__LINE__<<endl;
  //Prepare PHOS trigger utils if necessary
  fPHOSTrigUtils = new AliPHOSTriggerUtils("PHOSTrig") ;
  if(fForceRun){
    printf("Force run %d \n", fRunNumber) ;
    fPHOSTrigUtils->ForseUsingRun(fRunNumber) ;
  }
  if(fDoLightOutput) return;
  return;
}


//________________________________________________________________________
void AliCaloTriggerMimickHelper::UserExec(Option_t *){
    cout<<"Debug Output; AliCaloTriggerMimickHelper.C, Line: "<<__LINE__<<endl;
  // main method of AliCaloTriggerMimickHelper, first initialize and then process event
  if(!fForceRun)
    fRunNumber=fInputEvent->GetRunNumber() ;
  // do processing only for PHOS (2) clusters; for EMCal (1), DCal (3), EMCal with DCal (4) or  otherwise do nothing
  if(fClusterType == 2){
      SetEventChosenByTrigger(kFALSE);
      fPHOSTrigUtils->SetEvent(fInputEvent) ;
      Int_t nclus = 0;
      nclus = fInputEvent->GetNumberOfCaloClusters();
      // return if no Clusters in the event
      if(nclus == 0)  return;
      for(Int_t i = 0; i < nclus; i++){
          if (GetEventChosenByTrigger()){
              break;
          }
          AliVCluster* clus = NULL;
          clus = fInputEvent->GetCaloCluster(i);
          if (!clus) {
            continue;
          }
          SetTriggerDataOrMC(clus, fIsMC);
      }
  }
}


void AliCaloTriggerMimickHelper::SetTriggerDataOrMC(AliVCluster * clu, Bool_t isMCPhoton){
    cout<<"Debug Output; AliCaloTriggerMimickHelper.C, Line: "<<__LINE__<<endl;
    //Mark photons fired trigger
    if(isMCPhoton){
      SetEventChosenByTrigger(fPHOSTrigUtils->IsFiredTriggerMC(clu)&(1<<(fPHOSTrigger))) ;
    } else {
      SetEventChosenByTrigger(fPHOSTrigUtils->IsFiredTrigger(clu)) ;
    }
}

