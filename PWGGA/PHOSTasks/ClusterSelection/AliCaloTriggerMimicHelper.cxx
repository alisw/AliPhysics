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
    startDDLNumber(0),
    endDDLNumber(0),
    maxNumberOfDDLs(0),
    startTRU_Number(0),
    endTRU_Number(0),
    maxNumberOfTRUs(0),
    fPHOSTrigger(kPHOSAny),
    fPHOSTrigUtils(0x0),
    fGeomPHOS(NULL),
    fDoLightOutput(kFALSE),
    fForceRun(kFALSE),
    fIsMC(isMC),
    fTriggerHelperRunMode(0),
    fEventFlagPassed(0),
    fEventChosenByTrigger(kFALSE),
    fEventChosenByTriggerTrigUtils(kFALSE),
    fCurrentClusterTriggered(0),
    fCurrentClusterTriggeredTrigUtils(0),
    fCurrentClusterTriggerBadMapResult(0),
    fCurrentTriggeredClusterInBadDDL(0),
    fDoDebugOutput(0),
    minEnergyToTrigger(0),
    minCellsToTrigger(2),
    fEfficiencyChoiceOption_TriggerHelper(0),
    fMapClusterIDToHaveTriggered(),
    fMapClusterIDToIsInTriggerMap(),
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
    fdo_Any_4x4_Distance(0),
    fdo_4x4_Distance_All(0),
    fHist_4x4_Distance_All(NULL),
    fdo_Tr4x4_Distance_Triggered(0),
    fHist_Tr4x4_Distance_Triggered(NULL),
    fdo_Tr4x4_Distance_notTriggered(0),
    fHist_Tr4x4_Distance_notTriggered(NULL),
    fdo_Any_TRU(0),
    fdo_ClusEVsTiming_TRU(0),
    fHist_ClusEVsTiming_TRU(NULL),
    fdo_ClusEVsTiming_TRU_Trig(0),
    fHist_ClusEVsTiming_TRU_Trig(NULL),
    fdo_ClusEVsTiming_TRU_notTrig(0),
    fHist_ClusEVsTiming_TRU_notTrig(NULL),
    fdo_TRU_Numbers(0),
    fHist_TRU_Numbers(NULL),
    fdo_TRU_Channels(0),
    fHist_TRU_Channels(NULL),
    fdo_TRU_ChannelsXZ(0),
    fHist_TRU_ChannelsXZ(NULL),
    fEnergyThreshold_ColumnVsRow(1.)
{
    // Default constructor
    DefineInput(0, TChain::Class());
}

//================================================================================================================================================================
AliCaloTriggerMimicHelper::~AliCaloTriggerMimicHelper(){
    if (fDoDebugOutput>=3){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, AliCaloTriggerMimicHelper Line: "<<__LINE__<<endl;}
    // default deconstructor
    fMapClusterIDToHaveTriggered.clear();
    fMapClusterIDToIsInTriggerMap.clear();
    fMapTriggeredClusterInBadDDL.clear();
}

//================================================================================================================================================================
void AliCaloTriggerMimicHelper::Terminate(Option_t *){
    fMapClusterIDToHaveTriggered.clear();
    fMapClusterIDToIsInTriggerMap.clear();
    fMapTriggeredClusterInBadDDL.clear();

}

//================================================================================================================================================================
void AliCaloTriggerMimicHelper::UserCreateOutputObjects(){
    //SetDebugOutput(6);
    if (fDoDebugOutput>=3){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserCreateOutputObjects Line: "<<__LINE__<<endl;}
    fNMaxPHOSModules=5;
    maxCellsModule = maxColumns*maxRows; //56*64=3584
    nMaxCellsPHOS = (fNMaxPHOSModules*maxCellsModule); //56*64=3584
    startDDLNumber = 6;
    endDDLNumber = 19;
    maxNumberOfDDLs = endDDLNumber-startDDLNumber+1;
    startTRU_Number = 1;
    endTRU_Number   = 8;
    maxNumberOfTRUs = endTRU_Number-startTRU_Number+1;
    //Prepare PHOS trigger utils if necessary
    fPHOSTrigUtils = new AliPHOSTriggerUtils("PHOSTrig") ;
    if(fForceRun){
        printf("Force run %d \n", fRunNumber) ;
        fPHOSTrigUtils->ForseUsingRun(fRunNumber) ;
    }
    minEnergyToTrigger=0.1;
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
        if (fTriggerHelperRunMode==1){
            fHist_Event_Accepted->GetXaxis()->SetBinLabel(5,"L0");
        } else if (fTriggerHelperRunMode==2){
            fHist_Event_Accepted->GetXaxis()->SetBinLabel(5,"No (L0&INT7)");
        } else if (fTriggerHelperRunMode==3){
            fHist_Event_Accepted->GetXaxis()->SetBinLabel(5,"No INT7");
        }
    }

    if (fdo_fHist_Triggered_wEventFlag){
        Int_t iCurrentDDL_BinIndexBad;
        Int_t iCurrentDDL_BinIndexMaybeBad;
        Int_t iCurrentDDL_BinLabel;
        fHist_Triggered_wEventFlag        = new TH1I("fHist_Triggered_wEventFlag","fHist_Triggered_wEventFlag",9+(maxNumberOfDDLs*2),0.5,9.5+(maxNumberOfDDLs*2));
        fOutputList->Add(fHist_Triggered_wEventFlag);
        fHist_Triggered_wEventFlag->GetXaxis()->LabelsOption("v");
        fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(1,"mim.Trig.");
        fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(2,"mim.Trig.wL0");
        fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(3,"mim.Trig.woL0");
        fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(4,"L0womim.Trig.");
        fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(5,"All L0");
        fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(6,"All Events");
        fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(7,"good DDL");
        fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(8,"bad DDL");
        fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(9,"may.bad DDL");
        for (Int_t iCurrentDDL_BinIndex_Loop=0; iCurrentDDL_BinIndex_Loop<maxNumberOfDDLs; iCurrentDDL_BinIndex_Loop++){
            iCurrentDDL_BinIndexBad=10+iCurrentDDL_BinIndex_Loop;
            iCurrentDDL_BinIndexMaybeBad=(10+iCurrentDDL_BinIndex_Loop)+maxNumberOfDDLs;
            iCurrentDDL_BinLabel=iCurrentDDL_BinIndex_Loop+startDDLNumber;
            fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(iCurrentDDL_BinIndexBad,Form("bad DDL%d", iCurrentDDL_BinLabel));
            fHist_Triggered_wEventFlag->GetXaxis()->SetBinLabel(iCurrentDDL_BinIndexMaybeBad,Form("may.bad DDL%d", iCurrentDDL_BinLabel));
        }
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

    //Tr4x4 Information
    if (fdo_4x4_Distance_All){
        fHist_4x4_Distance_All = new TH2I("Tr4x4_Dist_All", "Tr4x4_Dist_All", 31, -15.5, 15.5, 31, -15.5, 15.5);
        fHist_4x4_Distance_All->GetXaxis()->SetTitle("Distance in x");
        fHist_4x4_Distance_All->GetYaxis()->SetTitle("Distance in z");
        fOutputList->Add(fHist_4x4_Distance_All);
    }
    if (fdo_Tr4x4_Distance_Triggered){
        fHist_Tr4x4_Distance_Triggered = new TH2I("Tr4x4_Dist_Trig", "Tr4x4_Dist_Trig", 31, -15.5, 15.5, 31, -15.5, 15.5);
        fHist_Tr4x4_Distance_Triggered->GetXaxis()->SetTitle("Distance in x");
        fHist_Tr4x4_Distance_Triggered->GetYaxis()->SetTitle("Distance in z");
        fOutputList->Add(fHist_Tr4x4_Distance_Triggered);
    }
    if (fdo_Tr4x4_Distance_notTriggered){
        fHist_Tr4x4_Distance_notTriggered = new TH2I("Tr4x4_Dist_notTrig", "Tr4x4_Dist_notTrig", 31, -15.5, 15.5, 31, -15.5, 15.5);
        fHist_Tr4x4_Distance_notTriggered->GetXaxis()->SetTitle("Distance in x");
        fHist_Tr4x4_Distance_notTriggered->GetYaxis()->SetTitle("Distance in z");
        fOutputList->Add(fHist_Tr4x4_Distance_notTriggered);
    }
    if ((fdo_4x4_Distance_All)||(fdo_Tr4x4_Distance_Triggered)||(fdo_Tr4x4_Distance_notTriggered)){
        fdo_Any_4x4_Distance            = 1;
    }

    //TRU Information
    if (fdo_ClusEVsTiming_TRU){
        fHist_ClusEVsTiming_TRU = new TH2D**[fNMaxPHOSModules];
        for (Int_t iNofModules=0; iNofModules<fNMaxPHOSModules; iNofModules++){
            fHist_ClusEVsTiming_TRU[iNofModules] = new TH2D*[maxNumberOfTRUs];
            for (Int_t iCurrentTRU_BinIndex_Loop=0; iCurrentTRU_BinIndex_Loop<maxNumberOfTRUs; iCurrentTRU_BinIndex_Loop++){
                fHist_ClusEVsTiming_TRU[iNofModules][iCurrentTRU_BinIndex_Loop] = new TH2D( Form("fHist_ClusEVsTiming_Mod%d_TRU%d", iNofModules, iCurrentTRU_BinIndex_Loop+startTRU_Number), Form("fHist_ClusEVsTiming_Mod%d_TRU%d", iNofModules, iCurrentTRU_BinIndex_Loop+startTRU_Number), 100, 0, 50, 100, 0.0e-9, 1000e-9);
                fHist_ClusEVsTiming_TRU[iNofModules][iCurrentTRU_BinIndex_Loop]->GetXaxis()->SetTitle("Energy in GeV");
                fHist_ClusEVsTiming_TRU[iNofModules][iCurrentTRU_BinIndex_Loop]->GetYaxis()->SetTitle("Timing in s");
                fOutputList->Add(fHist_ClusEVsTiming_TRU[iNofModules][iCurrentTRU_BinIndex_Loop]);
            }
        }
    }
    if (fdo_ClusEVsTiming_TRU_Trig){
        fHist_ClusEVsTiming_TRU_Trig = new TH2D**[fNMaxPHOSModules];
        for (Int_t iNofModules=0; iNofModules<fNMaxPHOSModules; iNofModules++){
            fHist_ClusEVsTiming_TRU_Trig[iNofModules] = new TH2D*[maxNumberOfTRUs];
            for (Int_t iCurrentTRU_BinIndex_Loop=0; iCurrentTRU_BinIndex_Loop<maxNumberOfTRUs; iCurrentTRU_BinIndex_Loop++){
                fHist_ClusEVsTiming_TRU_Trig[iNofModules][iCurrentTRU_BinIndex_Loop] = new TH2D( Form("fHist_ClusEVsTiming_Trig_Mod%d_TRU%d", iNofModules, iCurrentTRU_BinIndex_Loop+startTRU_Number), Form("fHist_ClusEVsTiming_Trig_Mod%d_TRU%d", iNofModules, iCurrentTRU_BinIndex_Loop+startTRU_Number), 100, 0, 50, 100, 0.0e-9, 1000e-9);
                fHist_ClusEVsTiming_TRU_Trig[iNofModules][iCurrentTRU_BinIndex_Loop]->GetXaxis()->SetTitle("Energy in GeV");
                fHist_ClusEVsTiming_TRU_Trig[iNofModules][iCurrentTRU_BinIndex_Loop]->GetYaxis()->SetTitle("Timing in s");
                fOutputList->Add(fHist_ClusEVsTiming_TRU_Trig[iNofModules][iCurrentTRU_BinIndex_Loop]);
            }
        }
    }
    if (fdo_ClusEVsTiming_TRU_notTrig){
        fHist_ClusEVsTiming_TRU_notTrig = new TH2D**[fNMaxPHOSModules];
        for (Int_t iNofModules=0; iNofModules<fNMaxPHOSModules; iNofModules++){
            fHist_ClusEVsTiming_TRU_notTrig[iNofModules] = new TH2D*[maxNumberOfTRUs];
            for (Int_t iCurrentTRU_BinIndex_Loop=0; iCurrentTRU_BinIndex_Loop<maxNumberOfTRUs; iCurrentTRU_BinIndex_Loop++){
                fHist_ClusEVsTiming_TRU_notTrig[iNofModules][iCurrentTRU_BinIndex_Loop] = new TH2D( Form("fHist_ClusEVsTiming_notTrig_Mod%d_TRU%d", iNofModules, iCurrentTRU_BinIndex_Loop+startTRU_Number), Form("fHist_ClusEVsTiming_notTrig_Mod%d_TRU%d", iNofModules, iCurrentTRU_BinIndex_Loop+startTRU_Number), 100, 0, 50, 100, 0.0e-9, 1000e-9);
                fHist_ClusEVsTiming_TRU_notTrig[iNofModules][iCurrentTRU_BinIndex_Loop]->GetXaxis()->SetTitle("Energy in GeV");
                fHist_ClusEVsTiming_TRU_notTrig[iNofModules][iCurrentTRU_BinIndex_Loop]->GetYaxis()->SetTitle("Timing in s");
                fOutputList->Add(fHist_ClusEVsTiming_TRU_notTrig[iNofModules][iCurrentTRU_BinIndex_Loop]);
            }
        }
    }
    if (fdo_TRU_Numbers){
        fHist_TRU_Numbers = new TH1I*[fNMaxPHOSModules];
        for (Int_t iNofModules=0; iNofModules<fNMaxPHOSModules; iNofModules++){
            fHist_TRU_Numbers[iNofModules] = new TH1I(Form("fHist_TRU_Numbers_Mod%d", iNofModules), Form("fHist_TRU_Numbers_Mod%d", iNofModules), 100, 0.5, 100.5);
            fOutputList->Add(fHist_TRU_Numbers[iNofModules]);
        }
    }
    if (fdo_TRU_Channels){
        fHist_TRU_Channels = new TH1I*[fNMaxPHOSModules];
        for (Int_t iNofModules=0; iNofModules<fNMaxPHOSModules; iNofModules++){
            fHist_TRU_Channels[iNofModules] = new TH1I(Form("fHist_TRU_Channels_Mod%d", iNofModules), Form("fHist_TRU_Channels_Mod%d", iNofModules), 200, 0.5, 200.5);
            fOutputList->Add(fHist_TRU_Channels[iNofModules]);
        }
    }
    if (fdo_TRU_ChannelsXZ){
        fHist_TRU_ChannelsXZ = new TH2I*[fNMaxPHOSModules];
        for (Int_t iNofModules=0; iNofModules<fNMaxPHOSModules; iNofModules++){
            fHist_TRU_ChannelsXZ[iNofModules] = new TH2I(Form("fHist_TRU_ChannelsXZ_Mod%d", iNofModules), Form("fHist_TRU_ChannelsXZ_Mod%d", iNofModules), 21, -0.5, 20.5, 21, -0.5, 20.5);
            fOutputList->Add(fHist_TRU_ChannelsXZ[iNofModules]);
        }
    }
    if ((fdo_ClusEVsTiming_TRU)||(fdo_ClusEVsTiming_TRU_Trig)||(fdo_ClusEVsTiming_TRU_notTrig)||(fdo_TRU_Numbers)||(fdo_TRU_Channels)||(fdo_TRU_ChannelsXZ)){
        fdo_Any_TRU=1;
    }

    return;
}


//================================================================================================================================================================
void AliCaloTriggerMimicHelper::UserExec(Option_t *){
    if (fDoDebugOutput>=3){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec Start, Line: "<<__LINE__<<"; GetPHOSTrigger: "<<GetPHOSTrigger()<<"; fRunNumber: "<<fRunNumber<<endl;}
    // main method of AliCaloTriggerMimicHelper, first initialize and then process event
    fMapClusterIDToHaveTriggered.clear();
    fMapClusterIDToIsInTriggerMap.clear();
    fMapTriggeredClusterInBadDDL.clear();
    fEventChosenByTrigger=0;
    fEventChosenByTriggerTrigUtils=0;
    fCurrentClusterTriggerBadMapResult=0;
    fCurrentTriggeredClusterInBadDDL=0;
    Double_t minEnergy_Debug=4.0;
    Int_t minEnergy_Reached_Debug=0;
    if(!fForceRun)
        fRunNumber=fInputEvent->GetRunNumber() ;
    AliInputEventHandler *fInputHandler=(AliInputEventHandler*)(AliAnalysisManager::GetAnalysisManager()->GetInputEventHandler());

    Bool_t isL0TriggerFlag;
    Bool_t isINT7TriggerFlag;
    if (fIsMC){
        isL0TriggerFlag=(fInputHandler->IsEventSelected() & AliVEvent::kAny);
        isINT7TriggerFlag=(fInputHandler->IsEventSelected() & AliVEvent::kAny);
    } else {
        isL0TriggerFlag=(fInputHandler->IsEventSelected() & AliVEvent::kPHI7);
        isINT7TriggerFlag=(fInputHandler->IsEventSelected() & AliVEvent::kINT7);
    }
    fEventFlagPassed=kFALSE;
    if (fdo_fHist_Triggered_wEventFlag){fHist_Triggered_wEventFlag->Fill(6);} //All Events
    if (isL0TriggerFlag) {
        if (fdo_fHist_Triggered_wEventFlag){fHist_Triggered_wEventFlag->Fill(5);} //All L0
    }
    // do processing only for PHOS (2) clusters; for EMCal (1), DCal (3), EMCal with DCal (4) or  otherwise do nothing
    if(fClusterType == 2){
        if (fdo_fHist_Event_Accepted){fHist_Event_Accepted->Fill(1);} //All Events
        if ((!isL0TriggerFlag)&&(fTriggerHelperRunMode == 0)) {//Triggered events need L0 flag; Only events with following event flags pass: L0
            if (fdo_fHist_Event_Accepted){fHist_Event_Accepted->Fill(5);} //No L0
            return;
        } else if ((isL0TriggerFlag)&&(fTriggerHelperRunMode == 1)){//MB Event Option; Only events with following event flags pass: No L0
            if (fdo_fHist_Event_Accepted){fHist_Event_Accepted->Fill(5);} //L0
            return;
        } else if ((!(isL0TriggerFlag&&isINT7TriggerFlag))&&(fTriggerHelperRunMode == 2)){//MB Event Option; Only events with following event flags pass: L0 and INT7
            if (fdo_fHist_Event_Accepted){fHist_Event_Accepted->Fill(5);} //No (L0&INT7)
            return;
        } else if ((!isINT7TriggerFlag)&&(fTriggerHelperRunMode == 3)){//MB Event Option; Only events with following event flags pass: INT7
            if (fdo_fHist_Event_Accepted){fHist_Event_Accepted->Fill(5);} //No INT7
            return;
        }
    
        fEventFlagPassed=kTRUE;
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
        Int_t CurrentClusterID=-1;
        Int_t CurrentDDL=-1;
        Int_t CurrentTRU=-1;
        Int_t CurrentTRUChannel=-1;
        Int_t CurrentTRUChannelX=-1;
        Int_t CurrentTRUChannelZ=-1;
        fGeomPHOS = AliPHOSGeometry::GetInstance();
        nModules = fGeomPHOS->GetNModules();
        nCellsPHOS=((nModules-1)*maxCellsModule); //56*64=3584
        fPHOSTrigUtils->SetEvent(fInputEvent) ;
        fPHOSTrigUtils->SetEfficiencyChoiceOption(fEfficiencyChoiceOption_TriggerHelper);
        Int_t nclus = 0;
        nclus = fInputEvent->GetNumberOfCaloClusters();
        // return if no Clusters in the event
        if(nclus == 0)  {
            if (fdo_fHist_Event_Accepted){fHist_Event_Accepted->Fill(3);} //noCluster
            return;
        }
        if (fdo_fHist_Cluster_Accepted){fHist_Cluster_Accepted->Fill(1, nclus);} //All Clusters
        AliVCaloCells * phsCells=fInputEvent->GetPHOSCells() ;
        //----------------------------------------------------------------------------------------------------
        for(Int_t i = 0; i < nclus; i++){
            fCurrentClusterTriggered=0;
            fCurrentClusterTriggeredTrigUtils=0;
            fCurrentClusterTriggerBadMapResult=0;
            fCurrentTriggeredClusterInBadDDL=0;
            cellAbsId = 0;
            isClusterGood =kTRUE;
            maxId=-1;
            eMax = -111;
            eCell=0;
            CurrentClusterID=i;
            if (fdo_fHist_Cluster_Accepted){fHist_Cluster_Accepted->Fill(2);} // All Clusters Checked
            if (fDoDebugOutput>=5){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec, Line: "<<__LINE__<<"; ClusterLoop i="<<i<<"; (nclus=="<<nclus<<")"<<endl;}
            //if (GetEventChosenByTriggerTrigUtils()){break;}
            AliVCluster* clus = NULL;
            clus = fInputEvent->GetCaloCluster(i);
            if (!clus) {
                continue;
            }
            if (fDoDebugOutput>=5){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec, Line: "<<__LINE__<<"; cluster E:"<<clus->E()<<endl;}
            if (fDoDebugOutput>=2){
                Int_t CurrentClusterID_ByCluster = clus->GetID();
                cout<<"Cluster Index by Loop: "<<i<<"; Cluster Index by GetID(): "<<CurrentClusterID_ByCluster<<endl;
            }
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
            if (fCurrentClusterTriggerBadMapResult == 0){ //Clusters which are marked bad by the TriggerBadMap have their fMapClusterIDToIsInTriggerMap value set to 0; empty map entries == 0
                isClusterGood=kFALSE;
            } else { //Clusters which are marked good by the TriggerBadMap have their fMapClusterIDToIsInTriggerMap value set to 1 and are able to trigger
                fMapClusterIDToIsInTriggerMap[CurrentClusterID]=fCurrentClusterTriggerBadMapResult;
            }
            if (clus->E()<minEnergyToTrigger){ //Do not use clusters, flagged by BadMap; Sets bad clusters energy to 0
                isClusterGood=kFALSE;
            }
            if (clus->GetNCells() < minCellsToTrigger){ //minimum amount of Cells required to trigger
                isClusterGood=kFALSE;
            }
            if (isClusterGood){
                if (fdo_fHist_Cluster_Accepted){fHist_Cluster_Accepted->Fill(4);} //Cluster good
                if (fDoDebugOutput>=1) {if (clus->E()>=minEnergy_Debug){minEnergy_Reached_Debug=1;}}
                if ((fDoDebugOutput>=2)&&(minEnergy_Reached_Debug>=1)){cout<<"Debug Output; AliCaloTriggerMimicHelper.C, UserExec, Line: "<<__LINE__<<"; cluster E:"<<clus->E()<<"; ClusterLoop i="<<i<<"; (nclus=="<<nclus<<")"<<endl;}
                SetTriggerDataOrMC(clus, fIsMC);
                if (fdo_Any_TRU){
                    CurrentTRU=WhichTRU(ix,iz);
                    if ((fdo_TRU_Channels)||(fdo_TRU_ChannelsXZ)){
                        CurrentTRUChannel=WhichTRUChannel(ix,iz, CurrentTRUChannelX, CurrentTRUChannelZ);
                        if (fdo_TRU_Channels){fHist_TRU_Channels[mod]->Fill(CurrentTRUChannel);}
                        if (fdo_TRU_ChannelsXZ){fHist_TRU_ChannelsXZ[mod]->Fill(CurrentTRUChannelX, CurrentTRUChannelZ);}
                    }
                    if (fdo_TRU_Numbers){fHist_TRU_Numbers[mod]->Fill(CurrentTRU);}
                    if (fdo_ClusEVsTiming_TRU){
                        if ((CurrentTRU>=startTRU_Number)&&(CurrentTRU<=endTRU_Number)&&(CurrentTRU!=-1)){
                            fHist_ClusEVsTiming_TRU[mod][CurrentTRU-startTRU_Number]->Fill(clus->E(), clus->GetTOF());
                        }
                    }
                }
                if (fdo_Any_4x4_Distance){//Shall only be active if you really want to see these for debugging
                    Int_t CurrentNtrg4x4        = (Int_t)fPHOSTrigUtils->GetNtrg4x4();
                    Int_t CurrentTrMod4x4       = 0;
                    Int_t CurrentTrX4x4         = 0;
                    Int_t CurrentTrZ4x4         = 0;
                    Int_t CurrentDistanceX      = 0;
                    Int_t CurrentDistanceZ      = 0;
                    for (Int_t itr=0; itr < CurrentNtrg4x4; itr++ ){
                        CurrentTrMod4x4         = (Int_t)fPHOSTrigUtils->GetTrMod4x4(itr);
                        if(CurrentTrMod4x4 != mod ) continue; //trigger fired in other module
                        CurrentTrX4x4           = (Int_t)fPHOSTrigUtils->GetTrX4x4(itr);
                        CurrentTrZ4x4           = (Int_t)fPHOSTrigUtils->GetTrZ4x4(itr);
                        CurrentDistanceX        = ix - CurrentTrX4x4;
                        CurrentDistanceZ        = iz - CurrentTrZ4x4;
                        if (fdo_4x4_Distance_All){
                            fHist_4x4_Distance_All->Fill(CurrentDistanceX, CurrentDistanceZ);
                        }
                        if (fCurrentClusterTriggeredTrigUtils>0) {
                            if (fdo_Tr4x4_Distance_Triggered){
                                fHist_Tr4x4_Distance_Triggered->Fill(CurrentDistanceX, CurrentDistanceZ);
                            }
                        } else {
                            if (fdo_Tr4x4_Distance_notTriggered){
                                fHist_Tr4x4_Distance_notTriggered->Fill(CurrentDistanceX, CurrentDistanceZ);
                            }
                        }
                    }
                }
                if (fCurrentClusterTriggeredTrigUtils>0){
                    fCurrentClusterTriggered=fCurrentClusterTriggeredTrigUtils;
                    fMapClusterIDToHaveTriggered[CurrentClusterID]=fCurrentClusterTriggered;
                    CurrentDDL=WhichDDL(mod,ix);
                    fCurrentTriggeredClusterInBadDDL=IsDDLBad(CurrentDDL, fRunNumber);
                    if (fCurrentTriggeredClusterInBadDDL >= 1 ){ //1==maybe bad; 2==bad DDLs
                        fMapTriggeredClusterInBadDDL[CurrentClusterID]=fCurrentTriggeredClusterInBadDDL;
                        if (fdo_fHist_Triggered_wEventFlag){fHist_Triggered_wEventFlag->Fill(9);} //Trig. maybe bad DDL
                        if (fdo_fHist_Triggered_wEventFlag){fHist_Triggered_wEventFlag->Fill(10+maxNumberOfDDLs+CurrentDDL-startDDLNumber);} //Trig. maybe bad DDL, specific
                        if (fCurrentTriggeredClusterInBadDDL == 2 ){ //2==bad DDLs
                            if (fdo_fHist_Triggered_wEventFlag){fHist_Triggered_wEventFlag->Fill(8);} //Trig. bad DDL
                            if (fdo_fHist_Triggered_wEventFlag){fHist_Triggered_wEventFlag->Fill(10+CurrentDDL-startDDLNumber);} //Trig. bad DDL, specific
                        }
                    } else {
                        if (fdo_fHist_Triggered_wEventFlag){fHist_Triggered_wEventFlag->Fill(7);} //Trig. good DDL
                    }
                    if (fdo_fHist_GammaClusE){fHist_GammaClusE_Trig->Fill(clus->E());}
                    if (fdo_TriggeredClusters_ColumnVsRow_overThresh){
                        if (clus->E()>=fEnergyThreshold_ColumnVsRow){fHist_TriggeredClusters_ColumnVsRow_overThresh[mod]->Fill(ix, iz, 1.);}
                    }
                    if (fdo_TriggeredClusters_ColumnVsRow_underThresh){
                        if (clus->E()<fEnergyThreshold_ColumnVsRow){fHist_TriggeredClusters_ColumnVsRow_underThresh[mod]->Fill(ix, iz, 1.);}
                    }
                    if (fdo_ClusEVsTiming_TRU_Trig){
                        if ((CurrentTRU>=startTRU_Number)&&(CurrentTRU<=endTRU_Number)&&(CurrentTRU!=-1)){
                            fHist_ClusEVsTiming_TRU_Trig[mod][CurrentTRU-startTRU_Number]->Fill(clus->E(), clus->GetTOF());
                        }
                    }
                } else {
                    if (fdo_fHist_GammaClusE){fHist_GammaClusE_notTrig->Fill(clus->E());}
                    if (fdo_ClusEVsTiming_TRU_notTrig){
                        if ((CurrentTRU>=startTRU_Number)&&(CurrentTRU<=endTRU_Number)&&(CurrentTRU!=-1)){
                            fHist_ClusEVsTiming_TRU_notTrig[mod][CurrentTRU-startTRU_Number]->Fill(clus->E(), clus->GetTOF());
                        }
                    }
                }
                if (fCurrentClusterTriggeredTrigUtils){
                    if (fdo_fHist_Cluster_Accepted){fHist_Cluster_Accepted->Fill(5);} //Triggered clusters
                } else {
                    if (fdo_fHist_Cluster_Accepted){fHist_Cluster_Accepted->Fill(6);} //Not triggered clusters
                }
            } else {
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
        }
        //Not Triggered Event
        else {
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
Int_t AliCaloTriggerMimicHelper::WhichTRU(Int_t cellx, Int_t cellz) //TRU go from x to x
{
    Int_t tru = -1;
    if(cellx<1 || 64<cellx){
      AliError("cellx is wrong! tru=-1 will return.");
      return -1;
    }
    if(cellz<1 || 56<cellz){
      AliError("cellz is wrong! tru=-1 will return.");
      return -1;
    }

    Int_t XID = (cellx -1) / 16 + 1;
    Int_t ZID = (cellz -1) / 28;

    tru = 2*XID - ZID;
    //cout << "cellx = " << cellx << " , cellz = " << cellz << " , tru = " << tru << endl;

    return tru;
}
//================================================================================================================================================================
Int_t AliCaloTriggerMimicHelper::WhichTRUChannel(Int_t cellx, Int_t cellz, Int_t &chX, Int_t &chZ)
{
  //this will return TRU channel 0-111.
  Int_t ch = -1;
  if(cellx<1 || 64<cellx){
    AliError("cellx is wrong! tru=-1 will return.");
    return -1;
  }
  if(cellz<1 || 56<cellz){
    AliError("cellz is wrong! tru=-1 will return.");
    return -1;
  }

  chX = ((cellx -1)/2) %  8;
  chZ = ((cellz -1)/2) % 14;

  ch = 8*chZ + chX;
  //printf("cellx = %d , cellz = %d , chX = %d , chZ = %d , ch = %d.\n",cellx,cellz,chX,chZ,ch);

  return ch;
}
//================================================================================================================================================================
Int_t  AliCaloTriggerMimicHelper::IsDDLBad(Int_t iDDL, Int_t iRun){ //returns 0 for good DDLs, 1 for maybe bad DDLs and 2 for bad DDLs
    if(iDDL<startDDLNumber || iDDL>endDDLNumber) return 2.;
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
