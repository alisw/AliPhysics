AliAnalysisTaskSED0BDT *AddTaskD0BDT(Bool_t readMC=kFALSE, Int_t system=0/*0=pp,1=PbPb*/,
								     Float_t minC=0, Float_t maxC=0,
								     TString finDirname="Loose", TString finname="",TString finObjname="D0toKpiCuts_pp",
								     TString BDTfilename="", Bool_t DoSidebndSample=kFALSE, Float_t SBndSampleFrac = 0.1)
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
  
  if(!readMC&&!DoSidebndSample){
	  TFile *fileBDT = TFile::Open(BDTfilename);
	  if(!fileBDT ||(fileBDT&& !fileBDT->IsOpen())) ::Fatal("AddTaskD0BDT", "BDT file not found : check your BDT object");
	  AliRDHFCutsD0toKpi* cut4bdt = (AliRDHFCutsD0toKpi*)fileBDT->Get("Cut4BDTptbin")->Clone();	// An simple cut file for trained BDT pT binning
	  //~ cut4bdt->SetDirectory(0);
	  Int_t Nptbins = cut4bdt->GetNPtBins();
	  TDirectory *initdir = (TDirectory*)fileBDT->Get("pT_0");
	  TList *BDTNamelist = (TList*)initdir->GetListOfKeys()->Clone("BDTNamelist");		// TKey list, only fname used
	  TList *bdtlist = new TList();														// to be saved BDT list

	  for(Int_t i=0;i<Nptbins;i++){
		  TDirectory *thisdir = (TDirectory*)fileBDT->Get(Form("pT_%d",i));
		  for(Int_t j=0;j<BDTNamelist->GetEntries();j++){
			  TString BDTobjname = BDTNamelist->At(j)->GetName();
			  AliRDHFBDT *thisbdt = (AliRDHFBDT*)(thisdir->Get(BDTobjname)->Clone(Form("pT_%d_%s",i,BDTobjname.Data())));
			  if(!thisbdt) ::Fatal("AddTaskD0BDT", Form("Failed to find BDT named %s",BDTobjname.Data()));
			  bdtlist->Add(thisbdt);
		  }
	  }
	  massD0Task->SetBDTNamesList(BDTNamelist);
	  massD0Task->SetBDTPtbins(cut4bdt);
	  massD0Task->SetBDTList(bdtlist);
	  fileBDT->Close();
  }
  if(DoSidebndSample){
	  massD0Task->SetBDTSampleSideband(DoSidebndSample);
	  massD0Task->SetBDTSidebandSamplingFraction(SBndSampleFrac);
  }
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
