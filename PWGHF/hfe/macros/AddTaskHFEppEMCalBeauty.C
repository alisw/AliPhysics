///////////////////////////////////////////////////////////////////
//                                                               //            
// AddTaskHFEppEMCalBeauty                                       //
// Author: Vivek Singh                                           //
//                                                               //
///////////////////////////////////////////////////////////////////
class AliAnalysisDataContainer;

AliAnalysisHFEppEMCalBeauty* AddTaskHFEppEMCalBeauty(

TString name = "", 
Bool_t PhysSelINT7 =  kTRUE,
Bool_t isEG1=kFALSE,
Bool_t isMC=kFALSE,
Bool_t SwitchPi0EtaWeight=kFALSE,
Bool_t SwitchNHFEeffi = kFALSE,
Bool_t SwitchEleRecoEffi = kFALSE,
Bool_t SwitchMCTempWeight= kFALSE,
Bool_t SwitchFillTemp = kFALSE,
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
Double_t DeltaEta =0.01,
Double_t DeltaPhi =0.01,
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
    printf("=======================================\n");
    printf("==========ADD TASK PARAMETERS==========\n");
    printf("=======================================\n");
    printf("Task name: %s \n", name.Data());
 
    printf("PhysSelINT7 Flag: %i \n",PhysSelINT7);
    printf("isEG1 Flag: %i \n",isEG1);
    printf("isMC Flag: %i \n",isMC);
  
    printf("SwitchPi0EtaWeight Flag: %i \n",SwitchPi0EtaWeight);
    printf("SwitchNHFEeffi Flag: %i \n",SwitchNHFEeffi);
    printf("SwitchEleRecoEffi Flag: %i \n",SwitchEleRecoEffi);
    printf("SwitchMCTempWeight Flag: %i \n",SwitchMCTempWeight);
    printf("SwitchFillTemp Flag: %i \n",SwitchFillTemp);
    printf("useTender Flag: %i \n",useTender);
    printf("ClsTypeEMC Flag: %i \n",ClsTypeEMC);
    printf("ClsTypeDCAL Flag: %i \n",ClsTypeDCAL);
 
    printf("fSwitchRIP Flag: %i \n",fSwitchRIP);
 
    printf("Min. Etarange: %f \n", Etarange);
    printf("Min. TPCNCrRows: %i \n", TPCNCrRows);
    printf("Min. RatioCrossedRowOverFindable: %f \n", RatioCrossedRowOverFindable);
    printf("Min. ITSNclus: %i \n", ITSNclus);
    printf("Max. TPCNclusPID: %i \n", TPCNclusPID);
    printf("ITS SPDBoth: %i \n",SPDBoth);
    printf("ITS SPDAny: %i \n",SPDAny);
    printf("ITS SPDFirst: %i \n",SPDFirst);
    printf("Min. DCAxyCut: %f \n",DCAxyCut);
    printf("Min. DCAzCut: %f \n",DCAzCut);
    printf("Min. DeltaEta: %f \n",DeltaEta);
    printf("Min. DeltaPhi: %f \n",DeltaPhi);
    printf("Min. TPC nsigma: %f \n", TPCnsigmin);
    printf("Max. TPC nsigma: %f \n", TPCnsigmax);
    printf("Min. EopEMin: %f \n", EopEMin);
    printf("Max. EopEMax: %f \n", EopEMax);
    printf("Min. M02Min: %f \n", M02Min);
    printf("Max1. M02Max1: %f \n", M02Max1);
    printf("Min2. M02Max2: %f \n", M02Max2);
    printf("Max3. M02Max3: %f \n", M02Max3);
 
    printf("InvmassCut: %f \n", InvmassCut);
    printf("AssoTPCCluster: %i \n", AssoTPCCluster);
    printf("AssoITSRefit: %i \n", AssoITSRefit);
    printf("AssopTMin: %f \n", AssopTMin);
    printf("AssoEtarange: %f \n", AssoEtarange);
    printf("AssoTPCnsig: %f \n", AssoTPCnsig);
    printf("===================================\n\n");

    AliVEvent::EOfflineTriggerTypes trigger;
    if(PhysSelINT7)  trigger=AliVEvent::kINT7;
    if(!PhysSelINT7) trigger=AliVEvent::kEMCEGA;

    cout<<"  PhysSelINT7  ===    "<< PhysSelINT7   <<"  trigger  ===    "<< trigger <<endl;

    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) { ::Error("AliAnalysisHFEppEMCalBeauty", "No analysis manager to connect to.");
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
    AliAnalysisHFEppEMCalBeauty *HFeTask = new AliAnalysisHFEppEMCalBeauty(taskname.Data());
    //HFeTask->SetDebugLevel(2);

    HFeTask->SelectCollisionCandidates(trigger);
    HFeTask->SetTrigger(trigger);
    HFeTask->SetMCAnalysis(isMC);

    HFeTask->SetPi0EtaWeightCalc(SwitchPi0EtaWeight);
    HFeTask->SetNonHFEEffi(SwitchNHFEeffi);
    HFeTask->SetElecRecoEffi(SwitchEleRecoEffi);
    HFeTask->SwitchMCTemplateWeightCalc(SwitchMCTempWeight);
    HFeTask->SwitchFillMCTemplate(SwitchFillTemp);


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
    HFeTask->SetEMCalMatching(DeltaEta,DeltaPhi);
    
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
     
      if(isMC && SwitchFillTemp){    
      //------------------- Jonghan's Method -------------------------------------------------------------------------
      // B meson pt correction
      TF1 *BmesonCentLow = new TF1("BmesonCentLow","(0.604848+0.323547*x-0.0530093*x*x)/1.48462",0.1,3.5);
      TF1 *BmesonCentHigh = new TF1("BmesonCentHigh","(1.39587-0.119037*x+0.0102669*x*x-(3.65110e-04)*x*x*x+(6.44063e-06)*x*x*x*x-(5.51122e-08)*x*x*x*x*x+(1.81335e-10)*x*x*x*x*x*x)/1.48462",3.5,100);
      TF1 *BmesonMinLow = new TF1("BmesonMinLow","(0.568465+0.305052*x-0.0472978*x*x)/1.5673",0.1,3.5);
      TF1 *BmesonMinHigh = new TF1("BmesonMinHigh","(1.29080-0.0910133*x+0.00815947*x*x-(2.81802e-04)*x*x*x+(4.76347e-06)*x*x*x*x-(3.89458e-08)*x*x*x*x*x+(1.22428e-10)*x*x*x*x*x*x)/1.5673",3.5,100);
      TF1 *BmesonMaxLow = new TF1("BmesonMaxLow","(0.634004+0.338291*x-0.0575622*x*x)/1.42542",0.1,3.5);
      TF1 *BmesonMaxHigh = new TF1("BmesonMaxHigh","(1.47864-0.140949*x+0.0118876*x*x-(4.28193e-04)*x*x*x+(7.69184e-06)*x*x*x*x-(6.69947e-08)*x*x*x*x*x+(2.23972e-10)*x*x*x*x*x*x)/1.42542",3.5,100);
      HFeTask->SetBcorrCentLow(BmesonCentLow);
      HFeTask->SetBcorrCentHigh(BmesonCentHigh);
      HFeTask->SetBcorrMinLow(BmesonMinLow);
      HFeTask->SetBcorrMinHigh(BmesonMinHigh);
      HFeTask->SetBcorrMaxLow(BmesonMaxLow);
      HFeTask->SetBcorrMaxHigh(BmesonMaxHigh);
      cout<<"---------------------------------------------------------------------------------------"<<endl;
      cout<<"------------------------------B meson pt correction------------------------------------"<<endl;
      cout<<"---------------------------------------------------------------------------------------"<<endl<<endl;

      // D meson pt correction
      TF1 *DmesonCorr = new TF1("DmesonCorr","(1/1.73303)*(0.552699*TMath::Power(1+(1.18107-1)*x/1.4446,1.18107/(1-1.18107))+0.158337*TMath::Exp(1.8436-x/1.28754)) / ((3.33487*((3.54275-1.)*(3.54275-2.))/(3.54275*0.501107*(3.54275*0.501107+9.24732*(3.54275-2.)))*pow(1.+(sqrt(9.24732*9.24732 + x*x)-9.24732)/(3.54275*0.501107),-3.54275)))",0.1,100);
      TF1 *DmesonCorrVar1 = new TF1("DmesonCorrVar1","(1/1.27055)*(0.552699*TMath::Power(1+(1.18107-1)*x/1.5446,1.18107/(1-1.18107))+0.088337*TMath::Exp(1.8436-x/1.28754)) / ((3.33487*((3.54275-1.)*(3.54275-2.))/(3.54275*0.501107*(3.54275*0.501107+9.24732*(3.54275-2.)))*pow(1.+(sqrt(9.24732*9.24732 + x*x)-9.24732)/(3.54275*0.501107),-3.54275)))",0.1,100);
      TF1 *DmesonCorrVar2 = new TF1("DmesonCorrVar2","(1/2.193)*(0.552699*TMath::Power(1+(1.18107-1)*x/1.3446,1.18107/(1-1.18107))+0.228337*TMath::Exp(1.8436-x/1.28754)) / ((3.33487*((3.54275-1.)*(3.54275-2.))/(3.54275*0.501107*(3.54275*0.501107+9.24732*(3.54275-2.)))*pow(1.+(sqrt(9.24732*9.24732 + x*x)-9.24732)/(3.54275*0.501107),-3.54275)))",0.1,100);
      HFeTask->SetDcorrFtn(DmesonCorr);
      HFeTask->SetDcorrFtnVar1(DmesonCorrVar1);
      HFeTask->SetDcorrFtnVar2(DmesonCorrVar2);
      cout<<"---------------------------------------------------------------------------------------"<<endl;
      cout<<"------------------------------D meson pt correction------------------------------------"<<endl;
      cout<<"---------------------------------------------------------------------------------------"<<endl<<endl;
      
      // Lc meson pt correction
      //TF1 *LcCorr = new TF1("LcCorr","12.6935*TMath::Exp(-0.568250*x)+0.00681196",1., 30.);
      TF1 *LcCorr = new TF1("LcCorr","1", .5, 30.);
      HFeTask->SetLccorrFtn(LcCorr);
      cout<<"---------------------------------------------------------------------------------------"<<endl;
      cout<<"------------------------------LambdaC pt correction for most-central analysis----------"<<endl;
      cout<<"---------------------------------------------------------------------------------------"<<endl<<endl;

      //------------------- Erin's Method --------------------------------------------------------------------------
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
