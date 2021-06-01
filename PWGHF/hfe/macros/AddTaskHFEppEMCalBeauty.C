
///////////////////////////////////////////////////////////////////
//                                                               //            
// AddTaskHFEppEMCalBeauty                                       //
// Author: Vivek Singh                                           //
//                                                               //
///////////////////////////////////////////////////////////////////

class AliAnalysisDataContainer;

AliAnalysisHFEppEMCalBeauty* AddTaskHFEppEMCalBeauty(

TString name = "", 
Bool_t isMC=kFALSE,
AliVEvent::EOfflineTriggerTypes trigger=AliVEvent::kINT7,
Bool_t isEG1=kFALSE,
Bool_t useTender = kTRUE,
Bool_t ClsTypeEMC = kTRUE,
Bool_t ClsTypeDCAL = kTRUE,
Bool_t fSwitchRIP=kTRUE,
Double_t Etarange= 0.7, 
Int_t TPCNCrRows=70,
Double_t RatioCrossedRowOverFindable=0.8,
Int_t ITSNclus= 3,
Int_t TPCNclusPID= 60,
Bool_t SPDBoth= kFALSE,
Bool_t SPDAny= kTRUE,
Bool_t SPDFirst= kFALSE,
Double_t DCAxyCut= 1.0,
Double_t DCAzCut=2.0,
Double_t TPCnsigmin= -1.0,
Double_t TPCnsigmax= 3.0,
Double_t EopEMin= 0.9,     
Double_t EopEMax= 1.2,     
Double_t  M02Min= 0.02,      
Double_t M02Max1= 0.9,
Double_t M02Max2= 0.7,
Double_t M02Max3= 0.5,
Double_t InvmassCut= 0.14,   
Int_t AssoTPCCluster= 60,
Bool_t AssoITSRefit= kTRUE ,
Double_t AssopTMin= 0.1  ,
Double_t AssoEtarange= 0.9 ,
Double_t AssoTPCnsig=  3.0

)

{

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) { ::Error("AliAnalysisHFEppEMCalBeauty", "No analysis manager to connect to.");
    return NULL;
    }

    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":MyTask";      // create a subfolder in the file

    TString finDirname         = "_";
    TString outBasicname = "MyOutputContainer";
    TString profname       = "coutputProf";     
    
    finDirname        += name.Data();
    outBasicname      += finDirname.Data();

   //+++++++++++++++++++++++++++++++++++++++++++++++++++++
    TString taskname="ElecAnalysis";
    AliAnalysisHFEppEMCalBeauty *HFeTask = new AliAnalysisHFEppEMCalBeauty(taskname.Data());
    HFeTask->SetDebugLevel(2);

    HFeTask->SelectCollisionCandidates(trigger);
    HFeTask->SetTrigger(trigger);
    HFeTask->SetMCAnalysis(isMC);
    HFeTask->SetTenderSwitch(useTender);
    HFeTask->SetClusterTypeEMC(ClsTypeEMC);
    HFeTask->SetClusterTypeDCAL(ClsTypeDCAL);
    HFeTask->SwitchRecalImpPar(fSwitchRIP);

    HFeTask->SetEtaRange(Etarange);            
    HFeTask->SetMinTPCCluster(TPCNCrRows);     
    HFeTask->SetMinRatioCrossedRowOverFindable(RatioCrossedRowOverFindable); 
    HFeTask->SetMinITSCluster(ITSNclus);
    HFeTask->SetMinTPCClusterPID(TPCNclusPID);
    HFeTask->SetHitsOnSPDLayers(SPDBoth,SPDAny,SPDFirst);
    HFeTask->SetDCACut(DCAxyCut,DCAzCut);
    
    HFeTask->SetTPCnsigma(TPCnsigmin,TPCnsigmax);
    HFeTask->SetEopE(EopEMin,EopEMax);
    HFeTask->SetShowerShapeEM02(M02Min,M02Max1,M02Max2,M02Max3);
    
    HFeTask->SetInvMassCut(InvmassCut);
    HFeTask->SetAssoTPCclus(AssoTPCCluster);
    HFeTask->SetAssoITSrefit(AssoITSRefit);
    HFeTask->SetAssopTMin(AssopTMin);
    HFeTask->SetAssoEtarange(AssoEtarange);
    HFeTask->SetAssoTPCnsig(AssoTPCnsig);


  if(trigger==AliVEvent::kINT7){

    isEG1=kFALSE;
    HFeTask->SetEMCalTriggerEG1(kFALSE);
    HFeTask->SetEMCalTriggerDG1(kFALSE);
    HFeTask->SetEMCalTriggerEG2(kFALSE);
    HFeTask->SetEMCalTriggerDG2(kFALSE);
    }


  if(trigger==AliVEvent::kEMCEGA &&  isEG1==kTRUE){

    HFeTask->SetEMCalTriggerEG2(kFALSE);
    HFeTask->SetEMCalTriggerDG2(kFALSE);

    if(ClsTypeEMC && ClsTypeDCAL){
        HFeTask->SetEMCalTriggerEG1(kTRUE);
        HFeTask->SetEMCalTriggerDG1(kTRUE);
    }
    if(ClsTypeEMC && !ClsTypeDCAL){
        HFeTask->SetEMCalTriggerEG1(kTRUE);
        HFeTask->SetEMCalTriggerDG1(kFALSE);
    }
    if(!ClsTypeEMC && ClsTypeDCAL){
        HFeTask->SetEMCalTriggerEG1(kFALSE);
        HFeTask->SetEMCalTriggerDG1(kTRUE);
    }
        
  }
 
  if(trigger==AliVEvent::kEMCEGA && isEG1==kFALSE){

     HFeTask->SetEMCalTriggerEG1(kFALSE);
     HFeTask->SetEMCalTriggerDG1(kFALSE);
               
     if(ClsTypeEMC && ClsTypeDCAL){
        HFeTask->SetEMCalTriggerEG2(kTRUE);
        HFeTask->SetEMCalTriggerDG2(kTRUE);
     }
     if(ClsTypeEMC && !ClsTypeDCAL){
        HFeTask->SetEMCalTriggerEG2(kTRUE);
        HFeTask->SetEMCalTriggerDG2(kFALSE);
     }
     if(!ClsTypeEMC && ClsTypeDCAL){
        HFeTask->SetEMCalTriggerEG2(kFALSE);
        HFeTask->SetEMCalTriggerDG2(kTRUE);
     }

    }

    mgr->AddTask(HFeTask);    

    //_________Structure of Task O/P
    AliAnalysisDataContainer *cinput = mgr->GetCommonInputContainer();
    AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(outBasicname,TList::Class(),AliAnalysisManager::kOutputContainer, fileName.Data());

    // your HFeTask needs input: here we connect the manager to your HFeTask
    mgr->ConnectInput(HFeTask,0,mgr->GetCommonInputContainer());
    mgr->ConnectOutput(HFeTask,1,coutput1);
  //  mgr->ConnectOutput(HFetaskINT7,2,mgr->CreateContainer(profname, TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));    


  return HFeTask;
}
