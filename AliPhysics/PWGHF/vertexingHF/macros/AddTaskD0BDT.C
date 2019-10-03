AliAnalysisTaskSED0BDT *AddTaskD0BDT(Bool_t readMC=kFALSE, Int_t system=0/*0=pp,1=PbPb*/,
								     Float_t minC=0, Float_t maxC=0,
								     TString finDirname="Loose", TString finname="",TString finObjname="D0toKpiCuts",
								     TString BDTfilename="", TString BDTobjnamepre="BDT",
								     Float_t BDTRespCut = -1., Bool_t DoSidebndSample=kFALSE, Bool_t GetRespTree = kTRUE, Float_t SBndSampleFrac = 0.1,
								     Float_t LeftSBndCut = 1.792, Float_t RightSBndCut = 1.942)
{
  //
  // AddTask for the AliAnalysisTaskSE for D0 candidates
  // invariant mass histogram and association with MC truth 
  // (using MC info in AOD) and cut variables distributions
  // C.Bianchin  chiara.bianchin@pd.infn.it
  //


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  Int_t flag=0;
  Int_t flagD0D0bar=0;
  Bool_t cutOnDistr = kFALSE;
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskD0Distr", "No analysis manager to connect to.");
    return NULL;
  }   

  TString filename="",out1name="",out2name="",out3name="",out4name="",out5name="",out6name="",out7name="",out8name="",out9name="", inname="";
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

    // mass, y distr
    out9name="coutputmassD0MassY";
    if(cutOnDistr) out9name+="C"; 
    if(flagD0D0bar==1)out9name+="D0";
    if(flagD0D0bar==2)out9name+="D0bar";


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

    //AliNormalizationCounter
    out5name="normalizationCounterLS";
    if(cutOnDistr) out5name+="C"; 
    if(flagD0D0bar==1)out5name+="D0";
    if(flagD0D0bar==2)out5name+="D0bar";

    // mass, pt, imp param distr
    out6name="coutputmassD0MassPtLS";
    if(cutOnDistr) out6name+="C"; 
    if(flagD0D0bar==1)out6name+="D0";
    if(flagD0D0bar==2)out6name+="D0bar";

    out7name ="coutputVarTreeLS";

    //detector signal hists
    out8name="detectorSignalsLS";
    if(cutOnDistr) out8name+="C"; 
    if(flagD0D0bar==1)out8name+="D0";
    if(flagD0D0bar==2)out8name+="D0bar";

    // mass, y distr
    out9name="coutputmassD0MassYLS";
    if(cutOnDistr) out9name+="C"; 
    if(flagD0D0bar==1)out9name+="D0";
    if(flagD0D0bar==2)out9name+="D0bar";

    
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
  out9name += finDirname.Data();
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
	::Fatal("AddTaskD0BDT", "Input file not found : check your cut object");
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
      RDHFD0toKpi->SetUseCentrality(AliRDHFCuts::kCentV0M);
    }
  }
  else   {
    RDHFD0toKpi = (AliRDHFCutsD0toKpi*)filecuts->Get(finObjname.Data());
    if(!RDHFD0toKpi){
      ::Fatal("AddTaskD0Mass", "Specific AliRDHFCuts not found");
      return NULL;
    }
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
  out9name+=centr;
  inname+=centr;
  
  Int_t Nptbins = RDHFD0toKpi->GetNPtBins();
  Float_t *ptbin = RDHFD0toKpi->GetPtBinLimits();
  
  TFile *fileBDT = TFile::Open(BDTfilename);
  if(!fileBDT ||(fileBDT&& !fileBDT->IsOpen())) ::Fatal("AddTaskD0BDT", "BDT file not found : check your BDT object");

  // Aanalysis task    
  TString taskname="BDTAnalysis";
  if (flag==0)taskname.Prepend("D0");
  else taskname.Prepend("LS");
  AliAnalysisTaskSED0BDT *massD0Task = new AliAnalysisTaskSED0BDT(taskname.Data(),RDHFD0toKpi);
  massD0Task->SetDebugLevel(0);
  massD0Task->SetArray(flag);
  massD0Task->SetReadMC(readMC);
  massD0Task->SetCutOnDistr(cutOnDistr);
  massD0Task->SetUsePid4Distr(kFALSE);
  massD0Task->SetFillOnlyD0D0bar(flagD0D0bar);
  massD0Task->SetSystem(system); //0=pp, 1=PbPb
  massD0Task->SetFillVarHists(kFALSE); // default is FALSE if System=PbPb

  massD0Task->SetAODMismatchProtection(1);
  massD0Task->SetFillPtHistos(kFALSE);
  massD0Task->SetFillImpactParameterHistos(kFALSE);
  massD0Task->SetFillYHistos(kFALSE);
  massD0Task->SetDrawDetSignal(false);
  massD0Task->SetPIDCheck(kFALSE);
  massD0Task->SetDoMCAcceptanceHistos(false);
  massD0Task->SetRejectSDDClusters(kFALSE);
  massD0Task->SetWriteVariableTree(kFALSE);
  
  TList *bdtlist = new TList();
  for(Int_t i=0;i<Nptbins;i++){
	  TString BDTobjname = BDTobjnamepre;
	  BDTobjname += Form("1_%.0f_%.0f",ptbin[i],ptbin[i+1]);
	  AliRDHFBDT *thisbdt = (AliRDHFBDT*)(fileBDT->Get(BDTobjname)->Clone(Form("_%s",BDTobjname.Data())));
	  if(!thisbdt) ::Fatal("AddTaskD0BDT", Form("Failed to find BDT named %s",BDTobjname.Data()));
	  //~ std::cout<<thisbdt->GetDesc()<<endl;
	  bdtlist->Add(thisbdt);
	  if(!DoSidebndSample){
		  TString BDT2objname1 = BDTobjnamepre; TString BDT2objname2 = BDTobjnamepre; TString BDT2objname3 = BDTobjnamepre;
		  TString BDT2objname4 = BDTobjnamepre; TString BDT2objname5 = BDTobjnamepre; TString BDT2objname6 = BDTobjnamepre;
		  BDT2objname1 += Form("2_%.0f_%.0f_0",ptbin[i],ptbin[i+1]);
		  BDT2objname2 += Form("2_%.0f_%.0f_1",ptbin[i],ptbin[i+1]);
		  BDT2objname3 += Form("2_%.0f_%.0f_2",ptbin[i],ptbin[i+1]);
		  BDT2objname4 += Form("2_%.0f_%.0f_3",ptbin[i],ptbin[i+1]);
		  BDT2objname5 += Form("2_%.0f_%.0f_4",ptbin[i],ptbin[i+1]);
		  BDT2objname6 += Form("2_%.0f_%.0f_5",ptbin[i],ptbin[i+1]);
		  AliRDHFBDT *thisbdt2_0 = (AliRDHFBDT*)(fileBDT->Get(BDT2objname1)->Clone(Form("_%s",BDT2objname1.Data())));
		  AliRDHFBDT *thisbdt2_1 = (AliRDHFBDT*)(fileBDT->Get(BDT2objname2)->Clone(Form("_%s",BDT2objname2.Data())));
		  AliRDHFBDT *thisbdt2_2 = (AliRDHFBDT*)(fileBDT->Get(BDT2objname3)->Clone(Form("_%s",BDT2objname3.Data())));
		  AliRDHFBDT *thisbdt2_3 = (AliRDHFBDT*)(fileBDT->Get(BDT2objname4)->Clone(Form("_%s",BDT2objname4.Data())));
		  AliRDHFBDT *thisbdt2_4 = (AliRDHFBDT*)(fileBDT->Get(BDT2objname5)->Clone(Form("_%s",BDT2objname5.Data())));
		  AliRDHFBDT *thisbdt2_5 = (AliRDHFBDT*)(fileBDT->Get(BDT2objname6)->Clone(Form("_%s",BDT2objname6.Data())));
		  if(!thisbdt2_0) ::Fatal("AddTaskD0BDT", Form("Failed to find BDT named %s",BDT2objname1.Data()));
		  if(!thisbdt2_1) ::Fatal("AddTaskD0BDT", Form("Failed to find BDT named %s",BDT2objname2.Data()));
		  if(!thisbdt2_2) ::Fatal("AddTaskD0BDT", Form("Failed to find BDT named %s",BDT2objname3.Data()));
		  if(!thisbdt2_3) ::Fatal("AddTaskD0BDT", Form("Failed to find BDT named %s",BDT2objname4.Data()));
		  if(!thisbdt2_4) ::Fatal("AddTaskD0BDT", Form("Failed to find BDT named %s",BDT2objname5.Data()));
		  if(!thisbdt2_5) ::Fatal("AddTaskD0BDT", Form("Failed to find BDT named %s",BDT2objname6.Data()));
		  bdtlist->Add(thisbdt2_0);
		  bdtlist->Add(thisbdt2_1);
		  bdtlist->Add(thisbdt2_2);
		  bdtlist->Add(thisbdt2_3);
		  bdtlist->Add(thisbdt2_4);
		  bdtlist->Add(thisbdt2_5);
	  }
  }
  fileBDT->Close();
  massD0Task->SetBDTGetRespTree(GetRespTree);
  massD0Task->SetBDTRespCut(BDTRespCut);
  massD0Task->SetBDTSidebandCut(LeftSBndCut,RightSBndCut);
  massD0Task->SetBDTSampleSideband(DoSidebndSample);
  massD0Task->SetBDTSidebandSamplingFraction(SBndSampleFrac);
  massD0Task->SetBDTList(bdtlist);
  
  mgr->AddTask(massD0Task);
  
  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinputmassD0 = mgr->CreateContainer(inname,TChain::Class(), 
							  AliAnalysisManager::kInputContainer);
							  
  out6name = "D0Ntuples" + finDirname;
  out7name = "D0BDTResponses" + finDirname;

  AliAnalysisDataContainer *coutputmassD01 = mgr->CreateContainer(out1name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //mass
  AliAnalysisDataContainer *coutputmassD02 = mgr->CreateContainer(out2name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //distr
  AliAnalysisDataContainer *coutputmassD03 = mgr->CreateContainer(out3name,TH1F::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //nev
  AliAnalysisDataContainer *coutputmassD04 = mgr->CreateContainer(out4name,AliRDHFCutsD0toKpi::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //cuts
  AliAnalysisDataContainer *coutputmassD05 = mgr->CreateContainer(out5name,AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //counter
  AliAnalysisDataContainer *coutputmassD06 = mgr->CreateContainer(out6name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data());
  AliAnalysisDataContainer *coutputmassD07 = mgr->CreateContainer(out7name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data());

    
  mgr->ConnectInput(massD0Task,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(massD0Task,1,coutputmassD01);
  mgr->ConnectOutput(massD0Task,2,coutputmassD02);
  mgr->ConnectOutput(massD0Task,3,coutputmassD03);
  mgr->ConnectOutput(massD0Task,4,coutputmassD04);
  mgr->ConnectOutput(massD0Task,5,coutputmassD05);
  mgr->ConnectOutput(massD0Task,6,coutputmassD06);
  mgr->ConnectOutput(massD0Task,7,coutputmassD07);

  return massD0Task;
}
