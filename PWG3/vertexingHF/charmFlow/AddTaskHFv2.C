AliAnalysisTaskSEHFv2 *AddTaskHFv2(TString filename="DplustoKpipiCuts.root", AliAnalysisTaskSEHFv2::DecChannel decCh=AliAnalysisTaskSEHFv2::kD0toKpi,Bool_t readMC=kFALSE,TString name="",Int_t flagep=1 /*0=tracks,1=V0*/)
{
  //
  // Test macro for the AliAnalysisTaskSE for  D 
  // mesons v2 analysis with event plane method
  // Authors: Chiara Bianchin, cbianchi@pd.infn.it, 
  //          Robert Grajcarek, grajcarek@physi.uni-heidelberg.de
  //          Giacomo Ortona, ortona@to.infn.it,
  //          Carlos Perez Lara, carlos.eugenio.perez.lara@cern.ch
  //          Francesco Prino, prino@to.infn.it
  // Get the pointer to the existing analysis manager via the static access method.
  //============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskHFv2", "No analysis manager to connect to.");
    return NULL;
  }
  Bool_t stdcuts=kFALSE;
  TFile* filecuts=new TFile(filename.Data());
  if(!filecuts->IsOpen()){
    cout<<"Input file not found:  using std cut object"<<endl;
    stdcuts=kTRUE;
  }
  
  AliRDHFCuts *analysiscuts=0x0;
  TString suffix="";

  TString cutsobjname="loosercuts";
  //Analysis cuts
  switch (decCh){
  case 0:
    cutsobjname="AnalysisCuts"; 
    if(stdcuts){
      analysiscuts = new AliRDHFCutsDplustoKpipi();
      analysiscuts->SetStandardCutsPbPb2010();
    }
    else analysiscuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get(cutsobjname);
    suffix="Dplus";
    break;
  case 1:
    cutsobjname="D0toKpiCuts";
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsD0toKpi();
      analysiscuts->SetStandardCutsPbPb2010();
    }
    else analysiscuts = (AliRDHFCutsD0toKpi*)filecuts->Get(cutsobjname);
    suffix="D0";
    break;
  case 2:
    cutsobjname="DStartoKpipiCuts";
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDStartoKpipi();
      analysiscuts->SetStandardCutsPbPb2010();
    }
    else analysiscuts = (AliRDHFCutsDStartoKpipi*)filecuts->Get(cutsobjname);
    suffix="Dstar";
    break;
  default:
    cout<<"Not available"<<endl;
    break;
  }

  if(!analysiscuts){
    cout<<"Specific AliRDHFCuts not found"<<endl;
    return;
  }

  suffix+=name;
  const Int_t nphibins=4;
  const Int_t nphibinlimits=nphibins+1;
  Float_t pi=TMath::Pi();
  Float_t philimits[nphibinlimits]={0., pi/4.,pi/2., 3./4.*pi, pi};

  // Analysis task    
  AliAnalysisTaskSEHFv2 *v2Task = new AliAnalysisTaskSEHFv2("HFv2Analysis",analysiscuts,decCh,nphibinlimits,philimits);
  v2Task->SetReadMC(readMC);

  v2Task->SetEtaGapFeatureForEventplaneFromTracks(kTRUE);
  
  v2Task->SetDebugLevel(0);
  
  if(flagep){
    //histogram for V0
    TFile *fpar = TFile::Open("VZEROParHist.root");
    TH2D *hh[6];
    for(Int_t i=0;i<6;i++){
      TString hhname;hhname.Form("parhist%d_%d",(i+2)*10,(i+3)*10);
      hh[i]=(TH2D*)fpar->Get(hhname.Data());
    } 
    v2Task->SetVZEROParHist(hh);
  }
  mgr->AddTask(v2Task);

  // Create containers for input/output

  TString contname=Form("cinputv2%s",suffix.Data());
  AliAnalysisDataContainer *cinputv2 = mgr->CreateContainer(contname.Data(),TChain::Class(),AliAnalysisManager::kInputContainer);

  TString outputfile = AliAnalysisManager::GetCommonFileName();
  TString outputhistos = outputfile += ":PWG3_D2H_HFv2"; 

  contname=Form("hEventsInfo%s",suffix.Data());
  AliAnalysisDataContainer *coutputstat = mgr->CreateContainer(contname.Data(),TH1F::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

  contname=Form("coutputv2%s",suffix.Data());
  AliAnalysisDataContainer *coutputv2 = mgr->CreateContainer(contname.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

  contname=Form("cutobj%s",suffix.Data());
  AliAnalysisDataContainer *cutobj = mgr->CreateContainer(contname.Data(),AliRDHFCuts::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());

  mgr->ConnectInput(v2Task,0,mgr->GetCommonInputContainer());
  
  mgr->ConnectOutput(v2Task,1,coutputstat);
  
  mgr->ConnectOutput(v2Task,2,coutputv2);

  mgr->ConnectOutput(v2Task,3,cutobj);
 
  if(flagep){
    contname=Form("coutputVZEROpar%s",suffix.Data());
    AliAnalysisDataContainer *coutputpar = mgr->CreateContainer(contname.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    mgr->ConnectOutput(v2Task,4,coutputpar);
  }

  return v2Task;
}
