AliAnalysisTaskSESignificance *AddTaskSignificance(TString filename="cuts4SignifMaxim.root",Int_t decCh=1,Bool_t readMC=kFALSE)
{
  //
  // Test macro for the AliAnalysisTaskSE for D meson candidates
  // Invariant mass histogram and
  // association with MC truth (using MC info in AOD)
  //  R. Bala, bala@to.infn.it
  // C. Bianchin, cbianchi@pd.infn.it
  // Get the pointer to the existing analysis manager via the static access method.
  //============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskSignificance", "No analysis manager to connect to.");
    return NULL;
  }

  TFile* filecuts=new TFile(filename.Data());
  if(!filecuts->IsOpen()){
    cout<<"Input file not found: exit"<<endl;
    return;
  }
  
  AliRDHFCuts *analysiscuts=0x0;
  TString suffix="";

  TString cutsobjname="loosercuts";
  //Analysis cuts
  switch (decCh){
  case 0:
    analysiscuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get(cutsobjname);
    suffix="Dplus";
    break;
  case 1:
    analysiscuts = (AliRDHFCutsD0toKpi*)filecuts->Get(cutsobjname);
    suffix="D0";
    break;
  case 2:
    analysiscuts = (AliRDHFCutsDstartoKpipi*)filecuts->Get(cutsobjname);
    suffix="Dstar";
    break;
  case 3:
    analysiscuts = (AliRDHFCutsDstoKKpi*)filecuts->Get(cutsobjname);
    suffix="Ds";
    break;
  case 4:
    analysiscuts = (AliRDHFCutsD0toKpipipi*)filecuts->Get(cutsobjname);
    suffix="D04";
    break;
  case 5:
    analysiscuts = (AliRDHFCutsLctopKpi*)filecuts->Get(cutsobjname);
    suffix="Lc";
    break;
  }
  //ptbins

  if(!analysiscuts){
    cout<<"Specific AliRDHFCuts not found"<<endl;
    return;
  }
  
  const Int_t nptbins=analysiscuts->GetNPtBins();
  Float_t* ptbins=analysiscuts->GetPtBinLimits();

  // Analysis task    
  const Int_t npars=analysiscuts->GetNVarsForOpt();//numbers of var for opt
  Bool_t* varsforopt=analysiscuts->GetVarsForOpt();
  cout<<"pt bins= "<<nptbins<<"  varsforopt= "<<npars<<"  nvars= "<<analysiscuts->GetNVars()<<endl;
  Int_t* nofcells=new Int_t[npars];//={4,4,4,4};
  Float_t** looses;
  looses=new Float_t*[npars];
  Float_t** tights;
  tights=new Float_t*[npars];
  TString axisTitle[npars];
  TString parname="";

  for (Int_t ivop=0;ivop<npars;ivop++){
    looses[ivop]=new Float_t[nptbins];
    tights[ivop]=new Float_t[nptbins];
    nofcells[ivop]=4;
  }

  Int_t count=0;
  
  for (Int_t ip=0;ip<nptbins;ip++){
    for(Int_t iv=0;iv<analysiscuts->GetNVars();iv++){
      if(varsforopt[iv]){
	looses[count][ip]=analysiscuts->GetCutValue(iv,ip);
	axisTitle[count]=(analysiscuts->GetVarNames())[iv];
	parname=Form("par%dptbin%d",count,ip);
	TParameter<float>* par=(TParameter<float>*)filecuts->Get(parname.Data());
	tights[count][ip]=par->GetVal();
	count++;
      }
    }
    count=0;
  }
  
  //creation TList of AliMultiDimVector
  TList *listMDV=new TList();
  listMDV->SetOwner();
  listMDV->SetName("listMDV");

  for (Int_t ip=0;ip<nptbins;ip++){
    Float_t *loosescut=new Float_t[npars];
    Float_t *tightscut=new Float_t[npars];
    for(Int_t i=0;i<npars;i++){
      loosescut[i]=looses[i][ip];
      tightscut[i]=tights[i][ip];
    }
    Float_t ptbincut[2]={ptbins[ip],ptbins[ip+1]};
    TString mdvname=Form("multiDimVectorPtBin%d",ip);
    AliMultiDimVector *mv=new AliMultiDimVector(mdvname.Data(),"MultiDimVector",1,ptbincut,npars,nofcells,loosescut,tightscut,axisTitle);
    listMDV->Add(mv);
  }
  
  AliAnalysisTaskSESignificance *sigTask = new AliAnalysisTaskSESignificance("SignificanceAnalysis",listMDV,analysiscuts,decCh,AliRDHFCuts::kAll);//AliRDHFCuts::kCandidate
  sigTask->SetReadMC(readMC);
  //sigTask->SetDoLikeSign(kTRUE);
  sigTask->SetDebugLevel(3);
  mgr->AddTask(sigTask);

  TString contname=Form("cinputSig%s",suffix.Data());
  // Create containers for input/output
  AliAnalysisDataContainer *cinputSig = mgr->CreateContainer(contname.Data(),TChain::Class(),AliAnalysisManager::kInputContainer);
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  TString outputhistos = outputfile += ":PWG3_D2H_Significance"; 
  contname=Form("coutputSig%s",suffix.Data());
  AliAnalysisDataContainer *coutputSig = mgr->CreateContainer(contname.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  contname=Form("coutputmv%s",suffix.Data());
  AliAnalysisDataContainer *coutputmv = mgr->CreateContainer(contname.Data(),TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());


  mgr->ConnectInput(sigTask,0,mgr->GetCommonInputContainer());
  
  mgr->ConnectOutput(sigTask,1,coutputSig);
  
  mgr->ConnectOutput(sigTask,2,coutputmv);
 
  return sigTask;
}

