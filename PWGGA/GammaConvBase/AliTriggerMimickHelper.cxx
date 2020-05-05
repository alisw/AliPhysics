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
#include "AliTriggerMimickHelper.h"
#include "AliEMCALRecoUtils.h"
#include "AliESDEvent.h"
#include "AliESDtrack.h"
#include "AliESDtrackCuts.h"
#include "AliInputEventHandler.h"
#include "AliTrackerBase.h"
#include "AliV0ReaderV1.h"
#include "AliPHOSTriggerUtils.h"

#include "TAxis.h"
#include "TChain.h"
#include "TH1F.h"
#include "TF1.h"

#include <vector>
#include <map>
#include <utility>

class iostream;

using namespace std;


ClassImp(AliTriggerMimickHelper)

//________________________________________________________________________
AliTriggerMimickHelper::AliTriggerMimickHelper(const char *name, Int_t clusterType, Bool_t isMC) : AliAnalysisTaskSE(name),
  fClusterType(clusterType),
  fRunNumber(-1),
  fPHOSTrigUtils(0x0),
  fPHOSTrigger(kPHOSAny),
  fForceRun(kFALSE),
  fIsMC(isMC),
  fEventChosenByTrigger(kFALSE),
  fDoLightOutput(kFALSE)
{
    // Default constructor
}

//________________________________________________________________________
AliTriggerMimickHelper::~AliTriggerMimickHelper(){
    // default deconstructor
}

//________________________________________________________________________
void AliTriggerMimickHelper::Terminate(Option_t *){

}

//________________________________________________________________________
void AliTriggerMimickHelper::UserCreateOutputObjects(){
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
void AliTriggerMimickHelper::UserExec(Option_t *){
  // main method of AliTriggerMimickHelper, first initialize and then process event
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


void AliTriggerMimickHelper::SetTriggerDataOrMC(AliVCluster * clu, Bool_t isMCPhoton){
    //Mark photons fired trigger
    if(isMCPhoton){
      SetEventChosenByTrigger(fPHOSTrigUtils->IsFiredTriggerMC(clu)&(1<<(fPHOSTrigger))) ;
    } else {
      SetEventChosenByTrigger(fPHOSTrigUtils->IsFiredTrigger(clu)) ;
    }
}

