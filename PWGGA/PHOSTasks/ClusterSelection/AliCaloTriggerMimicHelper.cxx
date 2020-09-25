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
#include "AliInputEventHandler.h"
#include "AliVEvent.h"
#include "AliCaloTriggerMimicHelper.h"
#include "AliESDEvent.h"
#include "AliInputEventHandler.h"
#include "AliVCluster.h"
#include "AliVCaloCells.h"
#include "TChain.h"

class iostream;

using namespace std;


ClassImp(AliCaloTriggerMimicHelper)

//================================================================================================================================================================
AliCaloTriggerMimicHelper::AliCaloTriggerMimicHelper(const char *name, Int_t clusterType, Int_t isMC) : AliAnalysisTaskSE(name),
    fNameOfClassObject(name),
    fClusterType(clusterType),
    fRunNumber(-1),
    nModules(4),
    fNMaxPHOSModules(0),
    nMaxCellsPHOS(0),
    maxRows(64),
    maxColumns(56),
    maxCellsModule(0),
    fPHOSTrigger(kPHOSAny),
    fPHOSTrigUtils(0x0),
    fGeomPHOS(NULL),
    fDoLightOutput(kFALSE),
    fForceRun(kFALSE),
    fIsMC(isMC),
    fTriggerHelperRunMode(0),
    fEventChosenByTrigger(kFALSE),
    fEventChosenByTriggerTrigUtils(kFALSE),
    fCurrentClusterTriggered(0),
    fCurrentClusterTriggeredTrigUtils(0),
    fCurrentClusterTriggerBadMapResult(0),
    fCurrentTriggeredClusterInBadDDL(0),
    fDoDebugOutput(0),
    fMapClusterToTriggered(),
    fMapClusterToTriggerMap(),
    fMapTriggeredClusterInBadDDL(),
    fOutputList(NULL),
    fdo_fHist_Event_Accepted(0),
    fHist_Event_Accepted(NULL),
    fdo_fHist_Triggered_wEventFlag(0),
    fHist_Triggered_wEventFlag(NULL),
    fdo_fHist_Cluster_Accepted(0),
    fHist_Cluster_Accepted(NULL),
    fdo_fHist_cellID(0),
    fHist_cellID_All(NULL),
    fHist_cellID_isAccepted(NULL),
    fdo_fHist_relID(0),
    fHist_relID0_All(NULL),
    fHist_relID0_cellIDwasAccepted(NULL),
    fHist_relID0_isAccepted(NULL),
    fdo_fHist_GammaClusE(0),
    fHist_GammaClusE_Trig(NULL),
    fHist_GammaClusE_notTrig(NULL),
    fdo_TriggeredClusters_ColumnVsRow_overThresh(0),
    fHist_TriggeredClusters_ColumnVsRow_overThresh(NULL),
    fdo_TriggeredClusters_ColumnVsRow_underThresh(0),
    fHist_TriggeredClusters_ColumnVsRow_underThresh(NULL),
    fEnergyThreshold_ColumnVsRow(1.)
{
    // Default constructor
    DefineInput(0, TChain::Class());
}

//================================================================================================================================================================
AliCaloTriggerMimicHelper::~AliCaloTriggerMimicHelper(){
    if (fDoDebugOutput>=3){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, AliCaloTriggerMimicHelper Line: "<<__LINE__<<endl;}
    // default deconstructor
    fMapClusterToTriggered.clear();
    fMapClusterToTriggerMap.clear();
    fMapTriggeredClusterInBadDDL.clear();
}

//================================================================================================================================================================
void AliCaloTriggerMimicHelper::Terminate(Option_t *){
    fMapClusterToTriggered.clear();
    fMapClusterToTriggerMap.clear();
    fMapTriggeredClusterInBadDDL.clear();

}

//================================================================================================================================================================
void AliCaloTriggerMimicHelper::UserCreateOutputObjects(){
    //SetDebugOutput(6);
    if (fDoDebugOutput>=3){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserCreateOutputObjects Line: "<<__LINE__<<endl;}
    fNMaxPHOSModules=4;
    maxCellsModule = maxColumns*maxRows; //56*64=3584
    nMaxCellsPHOS = (fNMaxPHOSModules*maxCellsModule); //56*64=3584
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
        fOutputList->SetName(Form("%s", fNameOfClassObject.Data()));
    }


    fdo_fHist_Event_Accepted             = 1;
    fdo_fHist_Triggered_wEventFlag       = 1;
    fdo_fHist_GammaClusE                 = 1;
    fdo_TriggeredClusters_ColumnVsRow_overThresh = 1;
    fdo_TriggeredClusters_ColumnVsRow_underThresh = 1;
    if ( fDoLightOutput == 0 ){
        fdo_fHist_Cluster_Accepted       = 1;
        fdo_fHist_cellID                 = 1;
        fdo_fHist_relID                  = 1;
    }

    if (fdo_fHist_Event_Accepted){
        fHist_Event_Accepted              = new TH1I("fHist_Event_Accepted","fHist_Event_Accepted",5,0.5,5.5);
        fOutputList->Add(fHist_Event_Accepted);
        fHist_Event_Accepted->GetXaxis()->SetBinLabel(1,"All Events");
        fHist_Event_Accepted->GetXaxis()->SetBinLabel(2,"Accepted Events");
        fHist_Event_Accepted->GetXaxis()->SetBinLabel(3,"noCluster");
        fHist_Event_Accepted->GetXaxis()->SetBinLabel(4,"Not triggered");
        fHist_Event_Accepted->GetXaxis()->SetBinLabel(5,"No L0");
    }

    if (fdo_fHist_Triggered_wEventFlag){
        fHist_Triggered_wEventFlag        = new TH1I("fHist_Triggered_wEventFlag","fHist_Triggered_wEventFlag",9,0.5,9.5);
        fOutputList->Add(fHist_Triggered_wEventFlag);
        fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(1,"mimickedTrigger");
        fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(2,"mimickedTrigger w L0");
        fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(3,"mimickedTrigger wo L0");
        fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(4,"L0 wo mimickedTrigger");
        fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(5,"All L0");
        fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(6,"All Events");
        fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(7,"Trig. good DDL");
        fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(8,"Trig. bad DDL");
        fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(9,"Trig. maybe bad DDL");
    }


    if (fdo_fHist_GammaClusE){
        fHist_GammaClusE_Trig             = new TH1D("fHist_GammaClusE_Trig","fHist_GammaClusE_Trig",100, 0.0, 50.);
        fOutputList->Add(fHist_GammaClusE_Trig);

        fHist_GammaClusE_notTrig          = new TH1D("fHist_GammaClusE_notTrig","fHist_GammaClusE_notTrig",100, 0.0, 50.);
        fOutputList->Add(fHist_GammaClusE_notTrig);
    }

    if (fdo_fHist_Cluster_Accepted){
        fHist_Cluster_Accepted              = new TH1I("fHist_Cluster_Accepted","fHist_Cluster_Accepted",6,0.5,6.5);
        fOutputList->Add(fHist_Cluster_Accepted);
        fHist_Cluster_Accepted->GetXaxis()->SetBinLabel(1,"All Clusters");
        fHist_Cluster_Accepted->GetXaxis()->SetBinLabel(2,"All Clusters Checked");
        fHist_Cluster_Accepted->GetXaxis()->SetBinLabel(3,"Cluster not good");
        fHist_Cluster_Accepted->GetXaxis()->SetBinLabel(4,"Cluster good");
        fHist_Cluster_Accepted->GetXaxis()->SetBinLabel(5,"Triggered clusters");
        fHist_Cluster_Accepted->GetXaxis()->SetBinLabel(6,"Not triggered clusters");
    }

    if (fdo_fHist_cellID){
        fHist_cellID_All                  = new TH1I("fHist_cellID","fHist_cellID",20000,-0.5,19999.5);
        fOutputList->Add(fHist_cellID_All);
        fHist_cellID_isAccepted           = new TH1I("fHist_cellID_isOK","fHist_cellID_isOK",20000,-0.5,19999.5);
        fOutputList->Add(fHist_cellID_isAccepted);
    }
    if (fdo_fHist_relID){
        fHist_relID0_All                  = new TH1I("fHist_relID0","fHist_relID0",10,-0.5,9.5);
        fOutputList->Add(fHist_relID0_All);
        fHist_relID0_cellIDwasAccepted    = new TH1I("fHist_relID0_cellIDisOK","fHist_relID0_cellIDisOK",10,-0.5,9.5);
        fOutputList->Add(fHist_relID0_cellIDwasAccepted);
        fHist_relID0_isAccepted           = new TH1I("fHist_relID0_isOK","fHist_relID0_isOK",10,-0.5,9.5);
        fOutputList->Add(fHist_relID0_isAccepted);
    }
    if (fdo_TriggeredClusters_ColumnVsRow_overThresh){
        fHist_TriggeredClusters_ColumnVsRow_overThresh = new TH2I*[fNMaxPHOSModules];
        for (Int_t iNofModules=0; iNofModules<fNMaxPHOSModules; iNofModules++){
            fHist_TriggeredClusters_ColumnVsRow_overThresh[iNofModules] = new TH2I(Form("TrigClusters_ColRow_ovThr_Mod%d", iNofModules ), Form("TrigClusters_ColRow_ovThr_Mod%d", iNofModules ), maxRows, 0.5, maxRows+0.5, maxColumns, 0.5, maxColumns+0.5);
            fOutputList->Add(fHist_TriggeredClusters_ColumnVsRow_overThresh[iNofModules]);
        }
    }
    if (fdo_TriggeredClusters_ColumnVsRow_underThresh){
        fHist_TriggeredClusters_ColumnVsRow_underThresh = new TH2I*[fNMaxPHOSModules];
        for (Int_t iNofModules=0; iNofModules<fNMaxPHOSModules; iNofModules++){
            fHist_TriggeredClusters_ColumnVsRow_underThresh[iNofModules] = new TH2I(Form("TrigClusters_ColRow_unThr_Mod%d", iNofModules ),Form("TrigClusters_ColRow_unThr_Mod%d", iNofModules ), maxRows, 0.5, maxRows+0.5, maxColumns, 0.5, maxColumns+0.5);
            fOutputList->Add(fHist_TriggeredClusters_ColumnVsRow_underThresh[iNofModules]);
        }
    }
    return;
}


//================================================================================================================================================================
void AliCaloTriggerMimicHelper::UserExec(Option_t *){
    if (fDoDebugOutput>=3){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec Start, Line: "<<__LINE__<<"; GetPHOSTrigger: "<<GetPHOSTrigger()<<"; fRunNumber: "<<fRunNumber<<endl;}
    // main method of AliCaloTriggerMimicHelper, first initialize and then process event
    fMapClusterToTriggered.clear();
    fMapClusterToTriggerMap.clear();
    fMapTriggeredClusterInBadDDL.clear();
    Double_t minEnergy_Debug=4.0;
    Int_t minEnergy_Reached_Debug=0;
    if(!fForceRun)
        fRunNumber=fInputEvent->GetRunNumber() ;
    AliInputEventHandler *fInputHandler=(AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());
    Bool_t isL0TriggerFlag=(fInputHandler->IsEventSelected() & AliVEvent::kPHI7);
    if (fdo_fHist_Triggered_wEventFlag){fHist_Triggered_wEventFlag->Fill(6);} //All Events
    if (isL0TriggerFlag) {
        if (fdo_fHist_Triggered_wEventFlag){fHist_Triggered_wEventFlag->Fill(5);} //All L0
    }
    // do processing only for PHOS (2) clusters; for EMCal (1), DCal (3), EMCal with DCal (4) or  otherwise do nothing
    if(fClusterType == 2){
        if (fdo_fHist_Event_Accepted){fHist_Event_Accepted->Fill(1);} //All Events
        SetEventChosenByTrigger(kFALSE);
        SetEventChosenByTriggerTrigUtils(kFALSE);
        if ((!isL0TriggerFlag)&&(fTriggerHelperRunMode == 0)) {
            if (fdo_fHist_Event_Accepted){fHist_Event_Accepted->Fill(5);} //No L0
            return;
        }
        Int_t  relid[4];
        Int_t maxId=-1;
        Double_t eMax = -111;
        Bool_t isClusterGood;
        Int_t cellAbsId;
        Int_t nCellsPHOS;
        Double_t eCell;
        Int_t mod;
        Int_t ix; //Rows: 64
        Int_t iz; //Columns: 56
        Int_t CurrentClusterID;
        Int_t CurrentDDL=0;
        fGeomPHOS = AliPHOSGeometry::GetInstance();
        nModules = fGeomPHOS->GetNModules();
        nCellsPHOS=((nModules-1)*maxCellsModule); //56*64=3584
        fPHOSTrigUtils->SetEvent(fInputEvent) ;
        Int_t nclus = 0;
        nclus = fInputEvent->GetNumberOfCaloClusters();
        // return if no Clusters in the event
        if(nclus == 0)  {
            if (fdo_fHist_Event_Accepted){fHist_Event_Accepted->Fill(3);} //noCluster
            return;
        }
        if (fdo_fHist_Cluster_Accepted){fHist_Cluster_Accepted->Fill(1, nclus);} //All Clusters
        AliVCaloCells * phsCells=fInputEvent->GetPHOSCells() ;
        fEventChosenByTrigger=0;
        fEventChosenByTriggerTrigUtils=0;
        //----------------------------------------------------------------------------------------------------
        for(Int_t i = 0; i < nclus; i++){
            if (fdo_fHist_Cluster_Accepted){fHist_Cluster_Accepted->Fill(2);} // All Clusters Checked
            if (fDoDebugOutput>=5){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec, Line: "<<__LINE__<<"; ClusterLoop i="<<i<<"; (nclus=="<<nclus<<")"<<endl;}
            //if (GetEventChosenByTriggerTrigUtils()){break;}
            AliVCluster* clus = NULL;
            clus = fInputEvent->GetCaloCluster(i);
            if (!clus) {
                continue;
            }
            if (fDoDebugOutput>=5){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec, Line: "<<__LINE__<<"; cluster E:"<<clus->E()<<endl;}
            cellAbsId = 0;
            isClusterGood =kTRUE;
            maxId=-1;
            eMax = -111;
            eCell=0;
            CurrentClusterID=i;
            fCurrentClusterTriggered=0;
            fCurrentClusterTriggeredTrigUtils=0;
            fCurrentClusterTriggerBadMapResult=0;
            fCurrentTriggeredClusterInBadDDL=0;
            //--------------------------------------------------
            for (Int_t iDig=0; iDig< clus->GetNCells(); iDig++){
                cellAbsId = clus->GetCellAbsId(iDig);
                fGeomPHOS->AbsToRelNumbering(cellAbsId,relid);
                if (fdo_fHist_cellID){fHist_cellID_All->Fill(cellAbsId);}
                if (fdo_fHist_relID){fHist_relID0_All->Fill(relid[0]);}
                if ((cellAbsId<0)||(cellAbsId>=nCellsPHOS)){
                    isClusterGood=kFALSE;
                    break;
                }
                if (fdo_fHist_cellID){fHist_cellID_isAccepted->Fill(cellAbsId);}
                if (fdo_fHist_relID){fHist_relID0_cellIDwasAccepted->Fill(relid[0]);}
                if ((relid[0]>=nModules)||(relid[0]<0)){
                    isClusterGood=kFALSE;
                    break;
                }
                eCell = phsCells->GetCellAmplitude(cellAbsId)*clus->GetCellAmplitudeFraction(iDig);
                if(eCell>eMax)
                {
                    eMax = eCell;
                    maxId = cellAbsId;
                }
                if (fdo_fHist_relID){fHist_relID0_isAccepted->Fill(relid[0]);}
            } //Loop over 0<=iDig<clus->GetNCells() ends
            //--------------------------------------------------
            //check BadMap
            fGeomPHOS->AbsToRelNumbering(maxId, relid);
            mod= relid[0]; //Module Number
            ix = relid[2]; //Row Number: 64
            iz = relid[3]; //Column Number: 56
            if (fDoDebugOutput>=5){cout<<"mod: "<<mod<<", ix: "<<ix<<"; iz: "<<iz<<endl;}
            fCurrentClusterTriggerBadMapResult=(Int_t)fPHOSTrigUtils->TestBadMap(mod,ix,iz);
            if (fCurrentClusterTriggerBadMapResult == 0){
                isClusterGood=kFALSE;
            } else {
                fMapClusterToTriggered[CurrentClusterID]=fCurrentClusterTriggerBadMapResult;
            }
            if (isClusterGood){
                if (fdo_fHist_Cluster_Accepted){fHist_Cluster_Accepted->Fill(4);} //Cluster good
                if (fDoDebugOutput>=1) {if (clus->E()>=minEnergy_Debug){minEnergy_Reached_Debug=1;}}
                if ((fDoDebugOutput>=2)&&(minEnergy_Reached_Debug>=1)){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec, Line: "<<__LINE__<<"; cluster E:"<<clus->E()<<"; ClusterLoop i="<<i<<"; (nclus=="<<nclus<<")"<<endl;}
                SetTriggerDataOrMC(clus, fIsMC);
                if ( (isL0TriggerFlag)&&(fCurrentClusterTriggeredTrigUtils>0) ){
                    if (fDoDebugOutput>=6){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec, Line: "<<__LINE__<<endl;}
                    fCurrentClusterTriggered=fCurrentClusterTriggeredTrigUtils;
                    fMapClusterToTriggerMap[CurrentClusterID]=fCurrentClusterTriggered;
                    CurrentDDL=WhichDDL(mod,ix);
                    fCurrentTriggeredClusterInBadDDL=IsDDLBad(CurrentDDL, fRunNumber);
                    if (fCurrentTriggeredClusterInBadDDL == 2 ){ //bad DDLs
                        fMapTriggeredClusterInBadDDL[CurrentClusterID]=fCurrentTriggeredClusterInBadDDL;
                        if (fdo_fHist_Triggered_wEventFlag){fHist_Triggered_wEventFlag->Fill(8);} //Trig. bad DDL
                    } else if (fCurrentTriggeredClusterInBadDDL == 1 ){ // maybe bad DDLs
                        fMapTriggeredClusterInBadDDL[CurrentClusterID]=fCurrentTriggeredClusterInBadDDL;
                        if (fdo_fHist_Triggered_wEventFlag){fHist_Triggered_wEventFlag->Fill(9);} //Trig. maybe bad DDL
                    } else {
                        if (fdo_fHist_Triggered_wEventFlag){fHist_Triggered_wEventFlag->Fill(7);} //Trig. good DDL
                    }
                    if (fdo_fHist_GammaClusE){fHist_GammaClusE_Trig->Fill(clus->E());}
                    if (fdo_TriggeredClusters_ColumnVsRow_overThresh){
                        if (clus->E()>=fEnergyThreshold_ColumnVsRow){fHist_TriggeredClusters_ColumnVsRow_overThresh[mod-1]->Fill(ix, iz, 1.);}
                    }
                    if (fdo_TriggeredClusters_ColumnVsRow_underThresh){
                        if (clus->E()<fEnergyThreshold_ColumnVsRow){fHist_TriggeredClusters_ColumnVsRow_underThresh[mod-1]->Fill(ix, iz, 1.);}
                    }
                } else {
                    if (fDoDebugOutput>=6){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec, Line: "<<__LINE__<<endl;}
                    if (fdo_fHist_GammaClusE){fHist_GammaClusE_notTrig->Fill(clus->E());}
                }
                if (fDoDebugOutput>=6){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec, Line: "<<__LINE__<<endl;}
                if (fCurrentClusterTriggeredTrigUtils){
                    if (fdo_fHist_Cluster_Accepted){fHist_Cluster_Accepted->Fill(5);} //Triggered clusters
                } else {
                    if (fDoDebugOutput>=6){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec, Line: "<<__LINE__<<endl;}
                    if (fdo_fHist_Cluster_Accepted){fHist_Cluster_Accepted->Fill(6);} //Not triggered clusters
                }
            } else {
                if (fDoDebugOutput>=6){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec, Line: "<<__LINE__<<endl;}
                if (fdo_fHist_Cluster_Accepted){fHist_Cluster_Accepted->Fill(3);} //Cluster not good
                if ((fDoDebugOutput>=3)&&(clus->E()>=minEnergy_Debug)){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec, Line: "<<__LINE__<<"; !isClusterGood"<<endl;}
            }
        }   // Loop over 0<=i<nclus ends
        //----------------------------------------------------------------------------------------------------
        if (fEventChosenByTriggerTrigUtils){
            if (fdo_fHist_Event_Accepted){fHist_Event_Accepted->Fill(2);} //Accepted Events
            if (fdo_fHist_Triggered_wEventFlag){fHist_Triggered_wEventFlag->Fill(1);} //mimickedTrigger
            if (isL0TriggerFlag) {
                fEventChosenByTrigger=fEventChosenByTriggerTrigUtils;
                if (fdo_fHist_Triggered_wEventFlag){fHist_Triggered_wEventFlag->Fill(2);} //mimickedTrigger w L0
            }
            else {
                if (fdo_fHist_Triggered_wEventFlag){fHist_Triggered_wEventFlag->Fill(3);} //mimickedTrigger wo L0
                if (fdo_fHist_Event_Accepted){fHist_Event_Accepted->Fill(5);} //No L0
            } //mimickedTrigger wo L0
        } else {
            if (fdo_fHist_Event_Accepted){fHist_Event_Accepted->Fill(4);} //Not triggered
            if (isL0TriggerFlag) {if (fdo_fHist_Triggered_wEventFlag){fHist_Triggered_wEventFlag->Fill(4);}} //L0 wo mimickedTrigger
        }
    }
    if ((fDoDebugOutput>=1)&&(minEnergy_Reached_Debug>=1)){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec End, Line: "<<__LINE__<<"; fIsMC: "<<fIsMC<<"; GetEventChosenByTrigger(): "<<GetEventChosenByTrigger()<<endl;}
}
//================================================================================================================================================================
void AliCaloTriggerMimicHelper::SetTriggerDataOrMC(AliVCluster * clu, Bool_t isMCPhoton){
    if (fDoDebugOutput>=3){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, SetTriggerDataOrMC Start, Line: "<<__LINE__<<"; isMCPhoton: "<<isMCPhoton<<endl;}
    //Mark photons fired trigger
    if(isMCPhoton){
        if(fPHOSTrigger==kPHOSAny || fPHOSTrigger==kPHOSL0) {
            fCurrentClusterTriggeredTrigUtils=fPHOSTrigUtils->IsFiredTriggerMC(clu)&1;
        } else {
            fCurrentClusterTriggeredTrigUtils=fPHOSTrigUtils->IsFiredTriggerMC(clu)&(1<<(fPHOSTrigger-1));
        }
        if ( (fEventChosenByTriggerTrigUtils==0)&&(fCurrentClusterTriggeredTrigUtils>=1) ){fEventChosenByTriggerTrigUtils=fCurrentClusterTriggeredTrigUtils;}
        if (fDoDebugOutput>=4){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, Line: "<<__LINE__<<"; GetEventChosenByTrigger(): "<<GetEventChosenByTrigger()<<endl;}
    } else {
        fCurrentClusterTriggeredTrigUtils=fPHOSTrigUtils->IsFiredTrigger(clu);
        if ( (fEventChosenByTriggerTrigUtils==0)&&(fCurrentClusterTriggeredTrigUtils>=1) ){fEventChosenByTriggerTrigUtils=fCurrentClusterTriggeredTrigUtils;}
        if (fDoDebugOutput>=4){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, Line: "<<__LINE__<<"; GetEventChosenByTrigger(): "<<GetEventChosenByTrigger()<<endl;}
    }
}
//================================================================================================================================================================
Int_t AliCaloTriggerMimicHelper::WhichDDL(Int_t module, Int_t cellx)
{
  const Int_t Nmod=5;//totally, 5 PHOS modules are designed.
  Int_t ddl = -1;

  if(cellx<1 || 64<cellx) return -1;

  if(module<1 || 4<module){
    return -1;
  }
  else{
    ddl = (Nmod-module) * 4 + (cellx-1)/16;//convert offline module numbering to online.
    return ddl;
  }
}
//================================================================================================================================================================
Int_t  AliCaloTriggerMimicHelper::IsDDLBad(Int_t iDDL, Int_t iRun){ //returns 0 for good DDLs, 1 for maybe bad DDLs and 2 for bad DDLs
    if(iDDL<6 || iDDL>19) return 2.;
    if((iRun>=252603 && iRun<=264347)){// LHC16_AllPeriods_pp_NomB
      switch(iDDL){
          case 11 : return 1;
          case 14 : return 1;
          case 16 : return 1;
          case 18 : return 2;
         default : return 0;
      }
    } else if((iRun>=270531 && iRun<=282704)) {// LHC17_AllPeriods_pp_NomB
      switch(iDDL){
          case 6 : return 2;
          case 7 : return 1;
          case 18 : return 2;
         default : return 0;
      }
    } else if((iRun>=284706 && iRun<=295232)) {// LHC18_AllPeriods_pp_NomB
      switch(iDDL){
          case 6 : return 2;
          case 18 : return 2;
         default : return 0;
      }
    }
    return 0; //all other Runranges are set to be good
  }
