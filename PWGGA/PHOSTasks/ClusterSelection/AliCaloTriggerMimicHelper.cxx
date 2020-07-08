/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. 	*
 *																			*
 * Author: Jens Lühder  													*
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
#include "AliCaloTriggerMimicHelper.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "TChain.h"

class iostream;

using namespace std;


ClassImp(AliCaloTriggerMimicHelper)

//________________________________________________________________________
AliCaloTriggerMimicHelper::AliCaloTriggerMimicHelper(const char *name, Int_t clusterType, Int_t isMC) : AliAnalysisTaskSE(name),
    fNameOfClassObject(name),
    fClusterType(clusterType),
    fRunNumber(-1),
    nModules(4),
    fNMaxPHOSModules(0),
    nMaxCellsPHOS(0),
    fPHOSTrigger(kPHOSAny),
    fPHOSTrigUtils(0x0),
    fGeomPHOS(NULL),
    fDoLightOutput(kFALSE),
    fForceRun(kFALSE),
    fIsMC(isMC),
    fEventChosenByTrigger(kFALSE),
    fDoDebugOutput(0),
    fOutputList(NULL),
    fOutputList_Debug(NULL),
    fHist_Event_Accepted(NULL),
    fHist_Cluster_Accepted(NULL),
    fHist_nModues(NULL),
    fHist_cellID_All(NULL),
    fHist_cellID_isAccepted(NULL),
    fHist_relID0_All(NULL),
    fHist_relID0_cellIDwasAccepted(NULL),
    fHist_relID0_isAccepted(NULL)
{
    // Default constructor
    DefineInput(0, TChain::Class());
}

//________________________________________________________________________
AliCaloTriggerMimicHelper::~AliCaloTriggerMimicHelper(){
    if (fDoDebugOutput>=3){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, AliCaloTriggerMimicHelper Line: "<<__LINE__<<endl;}
    // default deconstructor
}

//________________________________________________________________________
void AliCaloTriggerMimicHelper::Terminate(Option_t *){

}

//________________________________________________________________________
void AliCaloTriggerMimicHelper::UserCreateOutputObjects(){
  //SetDebugOutput(2);
  if (fDoDebugOutput>=3){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserCreateOutputObjects Line: "<<__LINE__<<endl;}
  fNMaxPHOSModules=4;
  nMaxCellsPHOS = (fNMaxPHOSModules*3584); //56*64=3584
  //Prepare PHOS trigger utils if necessary
  fPHOSTrigUtils = new AliPHOSTriggerUtils("PHOSTrig") ;
  if(fForceRun){
    printf("Force run %d \n", fRunNumber) ;
    fPHOSTrigUtils->ForseUsingRun(fRunNumber) ;
  }
  if (fDoDebugOutput>=3){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserCreateOutputObjects Line: "<<__LINE__<<endl;}
  if(fOutputList != NULL){
    delete fOutputList;
    fOutputList = NULL;
  }
  if(fOutputList == NULL){
    fOutputList = new TList();
    fOutputList->SetOwner(kTRUE);
    fOutputList->SetName(Form("%s_Output", fNameOfClassObject.Data()));
  }

  fHist_Event_Accepted              = new TH1I("fHist_Event_Accepted","fHist_Event_Accepted",4,0.5,4.5);
  fOutputList->Add(fHist_Event_Accepted);
  fHist_Event_Accepted->GetXaxis()->SetBinLabel(1,"All Events");
  fHist_Event_Accepted->GetXaxis()->SetBinLabel(2,"Accepted Events");
  fHist_Event_Accepted->GetXaxis()->SetBinLabel(3,"noCluster");
  fHist_Event_Accepted->GetXaxis()->SetBinLabel(4,"Not triggered");

  if(fDoLightOutput>=1) return;

  fHist_Cluster_Accepted              = new TH1I("fHist_Cluster_Accepted","fHist_Cluster_Accepted",6,0.5,6.5);
  fOutputList->Add(fHist_Cluster_Accepted);
  fHist_Cluster_Accepted->GetXaxis()->SetBinLabel(1,"All Clusters");
  fHist_Cluster_Accepted->GetXaxis()->SetBinLabel(2,"All Clusters Checked");
  fHist_Cluster_Accepted->GetXaxis()->SetBinLabel(3,"Cluster not good");
  fHist_Cluster_Accepted->GetXaxis()->SetBinLabel(4,"Cluster good");
  fHist_Cluster_Accepted->GetXaxis()->SetBinLabel(5,"Triggered clusters");
  fHist_Cluster_Accepted->GetXaxis()->SetBinLabel(6,"Not triggered clusters");

  fHist_cellID_All                  = new TH1I("fHist_cellID","fHist_cellID",20000,-0.5,19999.5);
  fOutputList->Add(fHist_cellID_All);
  fHist_cellID_isAccepted           = new TH1I("fHist_cellID_isOK","fHist_cellID_isOK",20000,-0.5,19999.5);
  fOutputList->Add(fHist_cellID_isAccepted);
  fHist_relID0_All                  = new TH1I("fHist_relID0","fHist_relID0",10,-0.5,9.5);
  fOutputList->Add(fHist_relID0_All);
  fHist_relID0_cellIDwasAccepted    = new TH1I("fHist_relID0_cellIDisOK","fHist_relID0_cellIDisOK",10,-0.5,9.5);
  fOutputList->Add(fHist_relID0_cellIDwasAccepted);
  fHist_relID0_isAccepted           = new TH1I("fHist_relID0_isOK","fHist_relID0_isOK",10,-0.5,9.5);
  fOutputList->Add(fHist_relID0_isAccepted);
  return;
}


//________________________________________________________________________
void AliCaloTriggerMimicHelper::UserExec(Option_t *){
    if (fDoDebugOutput>=3){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec Start, Line: "<<__LINE__<<"; GetPHOSTrigger: "<<GetPHOSTrigger()<<"; fRunNumber: "<<fRunNumber<<endl;}
  // main method of AliCaloTriggerMimicHelper, first initialize and then process event
  Double_t minEnergy_Debug=4.0;
  Int_t minEnergy_Reached_Debug=0;
  if(!fForceRun)
    fRunNumber=fInputEvent->GetRunNumber() ;
  // do processing only for PHOS (2) clusters; for EMCal (1), DCal (3), EMCal with DCal (4) or  otherwise do nothing
  if(fClusterType == 2){
      fHist_Event_Accepted->Fill(1);
      SetEventChosenByTrigger(kFALSE);
      Int_t  relid[4];
      Bool_t isClusterGood;
      Int_t cellAbsId;
      Int_t nCellsPHOS;
      fGeomPHOS = AliPHOSGeometry::GetInstance();
      nModules = fGeomPHOS->GetNModules();
      nCellsPHOS=((nModules-1)*3584); //56*64=3584
      fPHOSTrigUtils->SetEvent(fInputEvent) ;
      Int_t nclus = 0;
      nclus = fInputEvent->GetNumberOfCaloClusters();
      // return if no Clusters in the event
      if(nclus == 0)  {
          fHist_Event_Accepted->Fill(3);
          return;
      }

      if (fDoLightOutput<1){
        fHist_Cluster_Accepted->Fill(1, nclus);
      }
      for(Int_t i = 0; i < nclus; i++){
          if (fDoLightOutput<1){
            fHist_Cluster_Accepted->Fill(2);
          }
          if (fDoDebugOutput>=5){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec, Line: "<<__LINE__<<"; ClusterLoop i="<<i<<"; (nclus=="<<nclus<<")"<<endl;}
          if (GetEventChosenByTrigger()){
              break;
          }
          AliVCluster* clus = NULL;
          clus = fInputEvent->GetCaloCluster(i);
          if (!clus) {
            continue;
          }
          if (fDoDebugOutput>=5){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec, Line: "<<__LINE__<<"; cluster E:"<<clus->E()<<endl;}

          cellAbsId = 0;
          isClusterGood =kTRUE;
          for (Int_t iDig=0; iDig< clus->GetNCells(); iDig++){
                cellAbsId = clus->GetCellAbsId(iDig);
                fGeomPHOS->AbsToRelNumbering(cellAbsId,relid);
                if (fDoLightOutput<1){
                    fHist_cellID_All->Fill(cellAbsId);
                    fHist_relID0_All->Fill(relid[0]);
                }
                if ((cellAbsId<0)||(cellAbsId>=nCellsPHOS)){
                    isClusterGood=kFALSE;
                    break;
                }
                if (fDoLightOutput<1){
                    fHist_cellID_isAccepted->Fill(cellAbsId);
                    fHist_relID0_cellIDwasAccepted->Fill(relid[0]);
                }
                if ((relid[0]>=nModules)||(relid[0]<0)){
                    isClusterGood=kFALSE;
                    break;
                }
                if (fDoLightOutput<1){
                    fHist_relID0_isAccepted->Fill(relid[0]);
                }
          }
          if (isClusterGood){
              if (fDoLightOutput<1){
                  fHist_Cluster_Accepted->Fill(4);
              }
              if (fDoDebugOutput>=1) {if (clus->E()>=minEnergy_Debug){minEnergy_Reached_Debug=1;}}
              if ((fDoDebugOutput>=2)&&(minEnergy_Reached_Debug>=1)){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec, Line: "<<__LINE__<<"; cluster E:"<<clus->E()<<"; ClusterLoop i="<<i<<"; (nclus=="<<nclus<<")"<<endl;}
              SetTriggerDataOrMC(clus, fIsMC);
              if (fDoLightOutput<1){
                if (GetEventChosenByTrigger()){
                    fHist_Cluster_Accepted->Fill(5);
                } else {
                    fHist_Cluster_Accepted->Fill(6);
                }
              }
          } else {
              if (fDoLightOutput<1){
                fHist_Cluster_Accepted->Fill(3);
              }
              if ((fDoDebugOutput>=3)&&(clus->E()>=minEnergy_Debug)){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec, Line: "<<__LINE__<<"; !isClusterGood"<<endl;}
          }
      }   
      if (GetEventChosenByTrigger()){
          fHist_Event_Accepted->Fill(2);
      } else {
          fHist_Event_Accepted->Fill(4);
      }
  }
  if ((fDoDebugOutput>=1)&&(minEnergy_Reached_Debug>=1)){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec End, Line: "<<__LINE__<<"; fIsMC: "<<fIsMC<<"; GetEventChosenByTrigger(): "<<GetEventChosenByTrigger()<<endl;}
}


void AliCaloTriggerMimicHelper::SetTriggerDataOrMC(AliVCluster * clu, Bool_t isMCPhoton){
    if (fDoDebugOutput>=3){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, SetTriggerDataOrMC Start, Line: "<<__LINE__<<"; isMCPhoton: "<<isMCPhoton<<endl;}
    //Mark photons fired trigger
    if(isMCPhoton){
      SetEventChosenByTrigger(fPHOSTrigUtils->IsFiredTriggerMC(clu)&(1<<(fPHOSTrigger))) ;
      if (fDoDebugOutput>=4){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, Line: "<<__LINE__<<"; GetEventChosenByTrigger(): "<<GetEventChosenByTrigger()<<endl;}
    } else {
      SetEventChosenByTrigger(fPHOSTrigUtils->IsFiredTrigger(clu)) ;
      if (fDoDebugOutput>=4){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, Line: "<<__LINE__<<"; GetEventChosenByTrigger(): "<<GetEventChosenByTrigger()<<endl;}
    }
}

