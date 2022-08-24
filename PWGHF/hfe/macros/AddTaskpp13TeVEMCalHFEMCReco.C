///////////////////////////////////////////////////////////////////
//                                                               //            
// AddTaskpp13TeVEMCalHFEMCReco                                       //
// Author: Vivek Singh                                           //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysispp13TeVEMCalHFEMCReco* AddTaskpp13TeVEMCalHFEMCReco(

TString name = "", 
Bool_t PhysSelINT7 =  kTRUE,
Bool_t isEG1=kFALSE,
Bool_t isMC=kFALSE,
Bool_t SwitchPi0EtaWeight=kTRUE,
Bool_t SwitchNHFEeffi = kTRUE,
Bool_t SwitchEleRecoEffi = kTRUE,
Bool_t SwitchMCTempWeight= kTRUE,
Bool_t SwitchFillMCTemp = kTRUE,
Bool_t useTender = kTRUE,
Bool_t ClsTypeEMC = kTRUE,
Bool_t ClsTypeDCAL = kTRUE,
Bool_t fSwitchRIP=kTRUE,
Double_t Etarange= 0.6, 
Int_t TPCNCrRows=70,
Double_t RatioCrossedRowOverFindable=0.8,
Int_t ITSNclus= 3,
Int_t TPCNclusPID= 60,
Bool_t SPDBoth= kFALSE,
Bool_t SPDAny= kTRUE,
Bool_t SPDFirst= kFALSE,
Double_t DCAxyCut= 1,
Double_t DCAzCut=2,
Double_t TPCnsigmin= -1,
Double_t TPCnsigmax= 3,
Double_t EopEMin= 0.9,     
Double_t EopEMax= 1.2,     
Double_t  M02Min= 0.02,      
Double_t M02Max1= 0.9,
Double_t M02Max2= 0.7,
Double_t M02Max3= 0.5,
Double_t InvmassCut= 0.14,   
Int_t AssoTPCCluster= 60,
Bool_t AssoITSRefit= kTRUE,
Double_t AssopTMin= 0.1,
Double_t AssoEtarange= 0.9,
Double_t AssoTPCnsig=  3.0

)

{
    AliVEvent::EOfflineTriggerTypes trigger;
    if(PhysSelINT7)  trigger=AliVEvent::kINT7;
    if(!PhysSelINT7) trigger=AliVEvent::kEMCEGA;

  cout<<"  PhysSelINT7  ===    "<< PhysSelINT7   <<"  trigger  ===    "<< trigger <<endl;


    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) { ::Error("AliAnalysispp13TeVEMCalHFEMCReco", "No analysis manager to connect to.");
    return 0x0;
    }

    if (!mgr->GetInputEventHandler()) {
        return 0x0;
    }

    // by default, a file is open for writing. here, we get the filename
    TString fileName = AliAnalysisManager::GetCommonFileName();
    fileName += ":MyTask";      // create a subfolder in the file

    TString finDirname         = "_";
    TString outBasicname = "MyOutputContainer";   
    TString profname       = "coutputProf";     
    
    finDirname        += name.Data();
    outBasicname      += finDirname.Data();
    
    TString outBasicname1 = "Weights";
    outBasicname1      += finDirname.Data();
   //+++++++++++++++++++++++++++++++++++++++++++++++++++++   

    TString taskname="ElecAnalysis";
    AliAnalysispp13TeVEMCalHFEMCReco *HFeTask = new AliAnalysispp13TeVEMCalHFEMCReco(name.Data());
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

if(trigger==AliVEvent::kINT7)
{

  cout<<"  AliVEvent::kINT7  ===    "<< trigger <<" isEG1  "<<isEG1<<endl;


    isEG1=kFALSE;
    HFeTask->SetEMCalTriggerEG1(kFALSE);
    HFeTask->SetEMCalTriggerDG1(kFALSE);
    HFeTask->SetEMCalTriggerEG2(kFALSE);
    HFeTask->SetEMCalTriggerDG2(kFALSE);

    HFeTask->SwitchPi0EtaWeightCalc(SwitchPi0EtaWeight);
    HFeTask->SetNonHFEEffi(SwitchNHFEeffi);
    HFeTask->SetElecRecoEffi(SwitchEleRecoEffi);
    HFeTask->SwitchMCTemplateWeightCalc(SwitchMCTempWeight);
    HFeTask->SwitchFillMCTemplate(SwitchFillMCTemp);

    if(SwitchFillMCTemp)
    {
        TString DMesonWeightMaps, BMesonWeightMaps, CharmpTWeightMaps;
        
        BMesonWeightMaps = "alien:///alice/cern.ch/user/v/vksingh/DandBmesonpTweightCorrectionFiles/BMesonpTWeight.root";
        DMesonWeightMaps = "alien:///alice/cern.ch/user/v/vksingh/DandBmesonpTweightCorrectionFiles/DMesonpTWeight.root";
    
        CharmpTWeightMaps = "alien:///alice/cern.ch/user/v/vksingh/DandBmesonpTweightCorrectionFiles/CharmpTWeight.root";

        printf("\n### reading file %s ...\n",DMesonWeightMaps.Data());
        printf("\n### reading file %s ...\n",BMesonWeightMaps.Data());
        printf("\n### reading file %s ...\n",CharmpTWeightMaps.Data());
        
        TFile* f2 = TFile::Open(DMesonWeightMaps.Data());
        if(f2)
        {
            TH1 *D1 = (TH1*)f2->Get("RatD0");
            TH1 *D2 = (TH1*)f2->Get("RatD0Up");
            TH1 *D3 = (TH1*)f2->Get("RatD0Down");
            HFeTask->SetDmesonWeightHist(D1,D2,D3);
        }
            //  f2->Close();
        TFile* f3 = TFile::Open(BMesonWeightMaps.Data());
        if(f3)
        {
            TH1 *B1 = (TH1*)f3->Get("RatBMes");
            TH1 *B2 = (TH1*)f3->Get("RatBMesMin");
            TH1 *B3 = (TH1*)f3->Get("RatBMesMax");
            HFeTask->SetBmesonWeightHist(B1,B2,B3);
        }
            //  f3->Close();
        TFile* f4 = TFile::Open(CharmpTWeightMaps.Data());
        if(f4){
            TH1 *D0 = (TH1*)f4->Get("WeightD0");
            TH1 *DPlus = (TH1*)f4->Get("WeightDPlus");
            TH1 *Ds = (TH1*)f4->Get("WeightDs");
            TH1 *Lc = (TH1*)f4->Get("WeightLc");
                
            HFeTask->SetDmesonWeightHistPbPb(D0,DPlus,Ds,Lc);
        }

    }

}


  if(trigger==AliVEvent::kEMCEGA &&  isEG1==kTRUE){

  cout<<"  AliVEvent::kEMCEGA and EG1  ===    "<< trigger <<" isEG1  "<<isEG1<<endl;


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

  cout<<"  AliVEvent::kEMCEGA and EG1  ===    "<< trigger <<" isEG1  "<<isEG1<<endl;


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
