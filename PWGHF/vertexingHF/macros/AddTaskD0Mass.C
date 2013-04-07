AliAnalysisTaskSED0Mass *AddTaskD0Mass(Int_t flag=0/*0 = D0,1 = LS*/,Bool_t readMC=kFALSE,
				       Bool_t filldistr=kFALSE,Bool_t cutOnDistr=kFALSE,
				       Int_t system=0/*0=pp,1=PbPb*/,Int_t flagD0D0bar=0,
				       Float_t minC=0, Float_t maxC=0,
				       TString finDirname="Loose",
				       TString finname="",TString finObjname="D0toKpiCuts", Bool_t flagAOD049=kFALSE,
				       Bool_t FillMassPt=false, Bool_t FillImpPar=false,
					   Bool_t DrawDetSignal=false, Bool_t PIDCheck=false)
{
  //
  // AddTask for the AliAnalysisTaskSE for D0 candidates
  // invariant mass histogram and association with MC truth 
  // (using MC info in AOD) and cut variables distributions
  // C.Bianchin  chiara.bianchin@pd.infn.it
  //


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskD0Distr", "No analysis manager to connect to.");
    return NULL;
  }   

  TString filename="",out1name="",out2name="",out3name="",out4name="",out5name="",out6name="",out7name="",out8name="", inname="";
  filename = AliAnalysisManager::GetCommonFileName();
  filename += ":PWG3_D2H_";
  if(flag==0){
    filename+="D0InvMass";
    if(cutOnDistr) filename+="C"; 
    if(flagD0D0bar==1)filename+="D0";
    if(flagD0D0bar==2)filename+="D0bar";
    //list mass
    out1name="coutputmassD0Mass";
    if(cutOnDistr) out1name+="C"; 
    if(flagD0D0bar==1)out1name+="D0";
    if(flagD0D0bar==2)out1name+="D0bar";
    //list distr
    out2name="coutputmassD0distr";
    if(cutOnDistr) out2name+="C"; 
    if(flagD0D0bar==1)out2name+="D0";
    if(flagD0D0bar==2)out2name+="D0bar";
   //hist entries
    out3name="nEntriesD0";
    if(cutOnDistr) out3name+="C"; 
    if(flagD0D0bar==1)out3name+="D0";
    if(flagD0D0bar==2)out3name+="D0bar";
   //cuts object
    out4name="cutsD0";
    if(cutOnDistr) out4name+="C"; 
    if(flagD0D0bar==1)out4name+="D0";
    if(flagD0D0bar==2)out4name+="D0bar";

    //AliNormalizationCounter
    out5name="normalizationCounter";
    if(cutOnDistr) out5name+="C"; 
    if(flagD0D0bar==1)out5name+="D0";
    if(flagD0D0bar==2)out5name+="D0bar";

    // mass, pt, imp param distr
    out6name="coutputmassD0MassPt";
    if(cutOnDistr) out6name+="C"; 
    if(flagD0D0bar==1)out6name+="D0";
    if(flagD0D0bar==2)out6name+="D0bar";

    out7name ="coutputVarTree";

   //detector signal hists

    out8name="detectorSignals";
    if(cutOnDistr) out8name+="C"; 
    if(flagD0D0bar==1)out8name+="D0";
    if(flagD0D0bar==2)out8name+="D0bar";

    


    inname="cinputmassD0_0";
    if(cutOnDistr) inname+="C"; 
    if(flagD0D0bar==1)inname+="D0";
    if(flagD0D0bar==2)inname+="D0bar";

  } else {
    filename+="D0InvMassLikeSign";
    if(cutOnDistr) filename+="C"; 
    if(flagD0D0bar==1)filename+="D0";
    if(flagD0D0bar==2)filename+="D0bar";
    //list mass
    out1name="coutputmassLSMass";
    if(cutOnDistr) out1name+="C"; 
    if(flagD0D0bar==1)out1name+="D0";
    if(flagD0D0bar==2)out1name+="D0bar";
    //list distr
    out2name="coutputmassLSdistr";
    if(cutOnDistr) out2name+="C"; 
    if(flagD0D0bar==1)out2name+="D0";
    if(flagD0D0bar==2)out2name+="D0bar";
   //hist entries
    out3name="nEntriesLS";
    if(cutOnDistr) out3name+="C"; 
    if(flagD0D0bar==1)out3name+="D0";
    if(flagD0D0bar==2)out3name+="D0bar";
   //cuts object
    out4name="cutsLS";
    if(cutOnDistr) out4name+="C"; 
    if(flagD0D0bar==1)out4name+="D0";
    if(flagD0D0bar==2)out4name+="D0bar";

    inname="cinputmassD0_1";
    if(cutOnDistr) inname+="C"; 
    if(flagD0D0bar==1)inname+="D0";
    if(flagD0D0bar==2)inname+="D0bar";
  }
  filename += finDirname.Data();
  out1name += finDirname.Data();
  out2name += finDirname.Data();
  out3name += finDirname.Data(); 
  out4name += finDirname.Data(); 
  out5name += finDirname.Data(); 
  out6name += finDirname.Data();
  out7name += finDirname.Data();
  out8name += finDirname.Data();
  inname += finDirname.Data();

   //setting my cut values

    //cuts order
    //       printf("    |M-MD0| [GeV]    < %f\n",fD0toKpiCuts[0]);
    //     printf("    dca    [cm]  < %f\n",fD0toKpiCuts[1]);
    //     printf("    cosThetaStar     < %f\n",fD0toKpiCuts[2]);
    //     printf("    pTK     [GeV/c]    > %f\n",fD0toKpiCuts[3]);
    //     printf("    pTpi    [GeV/c]    > %f\n",fD0toKpiCuts[4]);
    //     printf("    |d0K|  [cm]  < %f\n",fD0toKpiCuts[5]);
    //     printf("    |d0pi| [cm]  < %f\n",fD0toKpiCuts[6]);
    //     printf("    d0d0  [cm^2] < %f\n",fD0toKpiCuts[7]);
    //     printf("    cosThetaPoint    > %f\n",fD0toKpiCuts[8]);

  Bool_t stdcuts=kFALSE;
  TFile* filecuts;
  if( finname.EqualTo("") ) {
    stdcuts=kTRUE; 
  } else {
      filecuts=TFile::Open(finname.Data());
      if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
	AliFatal("Input file not found : check your cut object");
      }
  }

  AliRDHFCutsD0toKpi* RDHFD0toKpi=new AliRDHFCutsD0toKpi();
  if(stdcuts) {
    if(system==0) RDHFD0toKpi->SetStandardCutsPP2010();
    else {
      RDHFD0toKpi->SetStandardCutsPbPb2011();
      if(minC!=0 && maxC!=0) { //if centrality 0 and 0 leave the values in the cut object
	RDHFD0toKpi->SetMinCentrality(minC);
	RDHFD0toKpi->SetMaxCentrality(maxC);
      }
      if(flagAOD049)RDHFD0toKpi->SetUseAOD049(kTRUE);
      RDHFD0toKpi->SetUseCentrality(AliRDHFCuts::kCentV0M);
    }
  }
  else   {
    RDHFD0toKpi = (AliRDHFCutsD0toKpi*)filecuts->Get(finObjname.Data());
    if(!RDHFD0toKpi){
      AliFatal("Specific AliRDHFCuts not found");
      return;
    }
    if(flagAOD049)RDHFD0toKpi->SetUseAOD049(kTRUE);
    if(minC!=0 && maxC!=0) { //if centrality 0 and 0 leave the values in the cut object
      RDHFD0toKpi->SetMinCentrality(minC);
      RDHFD0toKpi->SetMaxCentrality(maxC);
    } 
  }
  //  RDHFD0toKpi->SetName(Form("D0toKpiCuts%d",flag));

  TString centr="";
  if(minC!=0 && maxC!=0) centr = Form("%.0f%.0f",minC,maxC);
  else centr = Form("%.0f%.0f",RDHFD0toKpi->GetMinCentrality(),RDHFD0toKpi->GetMaxCentrality());
  out1name+=centr;
  out2name+=centr;
  out3name+=centr;
  out4name+=centr;
  out5name+=centr;
  out6name+=centr;
  out7name+=centr;
  out8name+=centr;
  inname+=centr;

  // Aanalysis task    
  TString taskname="MassAndDistrAnalysis";
  if (flag==0)taskname.Prepend("D0");
  else taskname.Prepend("LS");
  AliAnalysisTaskSED0Mass *massD0Task = new AliAnalysisTaskSED0Mass(taskname.Data(),RDHFD0toKpi);
  massD0Task->SetDebugLevel(0);
  massD0Task->SetArray(flag);
  massD0Task->SetReadMC(readMC);
  massD0Task->SetCutOnDistr(cutOnDistr);
  massD0Task->SetUsePid4Distr(kFALSE);
  massD0Task->SetFillOnlyD0D0bar(flagD0D0bar);
  massD0Task->SetSystem(system); //0=pp, 1=PbPb
  massD0Task->SetFillVarHists(filldistr); // default is FALSE if System=PbPb

  massD0Task->SetFillPtHistos(FillMassPt);
  massD0Task->SetFillImpactParameterHistos(FillImpPar);
  massD0Task->SetDrawDetSignal(DrawDetSignal);
  massD0Task->SetPIDCheck(PIDCheck);
  //  massD0Task->SetRejectSDDClusters(kTRUE);

  //   massD0Task->SetWriteVariableTree(kTRUE);

  mgr->AddTask(massD0Task);
  
  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinputmassD0 = mgr->CreateContainer(inname,TChain::Class(), 
							  AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *coutputmassD01 = mgr->CreateContainer(out1name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //mass
  AliAnalysisDataContainer *coutputmassD02 = mgr->CreateContainer(out2name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //distr
  AliAnalysisDataContainer *coutputmassD03 = mgr->CreateContainer(out3name,TH1F::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //nev
  AliAnalysisDataContainer *coutputmassD04 = mgr->CreateContainer(out4name,AliRDHFCutsD0toKpi::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //cuts
  AliAnalysisDataContainer *coutputmassD05 = mgr->CreateContainer(out5name,AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //counter
  AliAnalysisDataContainer *coutputmassD06 = mgr->CreateContainer(out6name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //mass vs pt vs impt par
  AliAnalysisDataContainer *coutputmassD07 = mgr->CreateContainer(out7name,TTree::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //mass vs pt vs impt par
  AliAnalysisDataContainer *coutputmassD08 = mgr->CreateContainer(out8name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //dedx
  

  
  mgr->ConnectInput(massD0Task,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(massD0Task,1,coutputmassD01);
  mgr->ConnectOutput(massD0Task,2,coutputmassD02);
  mgr->ConnectOutput(massD0Task,3,coutputmassD03);
  mgr->ConnectOutput(massD0Task,4,coutputmassD04);
  mgr->ConnectOutput(massD0Task,5,coutputmassD05);
  mgr->ConnectOutput(massD0Task,6,coutputmassD06);
  mgr->ConnectOutput(massD0Task,7,coutputmassD07);
  mgr->ConnectOutput(massD0Task,8,coutputmassD08);


  return massD0Task;
}
