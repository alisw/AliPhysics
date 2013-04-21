AliAnalysisTaskSED0Correlations *AddTaskD0Correlations(Bool_t readMC=kFALSE, Bool_t mixing=kFALSE, Bool_t recoTrMC=kFALSE, Bool_t recoD0MC = kFALSE, Bool_t effOn=kTRUE, Double_t etacorr=1.5, Int_t system=0/*0=pp,1=PbPb*/, Int_t flagD0D0bar=0, Float_t minC=0, Float_t maxC=0, TString finDirname="Output", TString cutsfilename="D0toKpiCuts.root", TString cutsfilename2="AssocPartCuts_fBit0_woITS.root", TString 
cutsD0name="D0toKpiCuts", TString cutsTrkname="AssociatedTrkCuts", TString effName = "3D_eff_wo_ITScls2_f0_p8eta.root", Bool_t flagAOD049=kFALSE, Int_t standardbins=1, Bool_t stdcuts=kFALSE)
{
  //
  // AddTask for the AliAnalysisTaskSE for D0 candidates
  // invariant mass histogram and association with MC truth 
  // (using MC info in AOD) and cut variables distributions
  // C.Bianchin  chiara.bianchin@pd.infn.it
  // F.Colamaria fabio.colamaria@ba.infn.it


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskD0Distr", "No analysis manager to connect to.");
    return NULL;
  }   

  TString filename="",out1name="",out2name="",out3name="",out4name="",out5name="",out6name="",out7name="",inname="";
  filename = AliAnalysisManager::GetCommonFileName();
  filename += ":PWG3_D2H_";
  filename+="D0InvMass";
  if(flagD0D0bar==1)filename+="D0";
  if(flagD0D0bar==2)filename+="D0bar";
  //list mass
  out1name="coutputmassD0Mass";
  if(flagD0D0bar==1)out1name+="D0";
  if(flagD0D0bar==2)out1name+="D0bar";
  //hist entries
  out2name="nEntriesD0";
  if(flagD0D0bar==1)out2name+="D0";
  if(flagD0D0bar==2)out2name+="D0bar";
  //cuts object
  out3name="cutsD0"; 
  if(flagD0D0bar==1)out3name+="D0";
  if(flagD0D0bar==2)out3name+="D0bar";
  //AliNormalizationCounter
  out4name="normalizationCounter";
  if(flagD0D0bar==1)out4name+="D0";
  if(flagD0D0bar==2)out4name+="D0bar";
  //correlation outputs
  out5name ="correlations";
  if(flagD0D0bar==1)out5name+="D0";
  if(flagD0D0bar==2)out5name+="D0bar";
  //correlation further studies
  out6name ="MCStudyPlots";
  if(flagD0D0bar==1)out6name+="D0";
  if(flagD0D0bar==2)out6name+="D0bar";
  //correlated trk cuts
  out7name ="cutsTracks";
  if(flagD0D0bar==1)out7name+="D0";
  if(flagD0D0bar==2)out7name+="D0bar";
  inname="cinputmassD0_0";
  if(flagD0D0bar==1)inname+="D0";
  if(flagD0D0bar==2)inname+="D0bar";

  filename += finDirname.Data();
  out1name += finDirname.Data();
  out2name += finDirname.Data(); 
  out3name += finDirname.Data(); 
  out4name += finDirname.Data(); 
  out5name += finDirname.Data();
  out6name += finDirname.Data();
  out7name += finDirname.Data();
  inname += finDirname.Data();

    //setting my cut values

    //cuts order
    //     printf("    |M-MD0| [GeV]    < %f\n",fD0CorrCuts[0]);
    //     printf("    dca    [cm]  < %f\n",fD0CorrCuts[1]);
    //     printf("    cosThetaStar     < %f\n",fD0CorrCuts[2]);
    //     printf("    pTK     [GeV/c]    > %f\n",fD0CorrCuts[3]);
    //     printf("    pTpi    [GeV/c]    > %f\n",fD0CorrCuts[4]);
    //     printf("    |d0K|  [cm]  < %f\n",fD0CorrCuts[5]);
    //     printf("    |d0pi| [cm]  < %f\n",fD0CorrCuts[6]);
    //     printf("    d0d0  [cm^2] < %f\n",fD0CorrCuts[7]);
    //     printf("    cosThetaPoint    > %f\n",fD0CorrCuts[8]);

  TFile* filecuts=new TFile(cutsfilename.Data());
  if(!filecuts->IsOpen()){
    cout<<"Input file not found for D0 cuts: using std cut object"<<endl;
    stdcuts=kTRUE;
  }
  TFile* filecuts2=new TFile(cutsfilename2.Data());
  if(!filecuts2->IsOpen()){
    cout<<"Input file not found for tracks cuts!"<<endl;
    return;
  }

  //Cuts for D0
  AliRDHFCutsD0toKpi* RDHFD0Corrs=new AliRDHFCutsD0toKpi();
//RDHFD0Corrs->SetUsePhysicsSelection(kFALSE);
  RDHFD0Corrs->SetLowPt(kFALSE); //low-pt special PID disabled
  if(stdcuts) {
    if(system==0) RDHFD0Corrs->SetStandardCutsPP2010();
    else {
      RDHFD0Corrs->SetStandardCutsPbPb2010();
      if(minC!=0 && maxC!=0) { //if centrality 0 and 0 leave the values in the cut object
	RDHFD0Corrs->SetMinCentrality(minC);
	RDHFD0Corrs->SetMaxCentrality(maxC);
      }
      if(flagAOD049)RDHFD0Corrs->SetUseAOD049(kTRUE);
      RDHFD0Corrs->SetUseCentrality(AliRDHFCuts::kCentV0M);
    }
  }
  else {
    RDHFD0Corrs = (AliRDHFCutsD0toKpi*)filecuts->Get(cutsD0name.Data());
    if(!RDHFD0Corrs){
      cout<<"Specific AliRDHFCuts not found"<<endl;
      return;
    }
    if(flagAOD049)RDHFD0Corrs->SetUseAOD049(kTRUE);
    if(minC!=0 && maxC!=0) { //if centrality 0 and 0 leave the values in the cut object
      RDHFD0Corrs->SetMinCentrality(minC);
      RDHFD0Corrs->SetMaxCentrality(maxC);
    } 
  }

  if(effOn) {
    //Load efficiency map
    TFile* fileeff=new TFile(effName.Data());
    if(!fileeff->IsOpen()){
      cout<<"Input file not found for efficiency! Exiting..."<<endl;
      return;
    }  
    TCanvas *c = (TCanvas*)fileeff->Get("c");
    TH3D *h3D = (TH3D*)c->FindObject("heff_rebin");
  }

  //Cuts for correlated tracks/K0
  AliHFAssociatedTrackCuts* corrCuts=new AliHFAssociatedTrackCuts();
  corrCuts = (AliHFAssociatedTrackCuts*)filecuts2->Get(cutsTrkname.Data());
  if(!corrCuts){
      cout<<"Specific AliHFAssociatedTrackCuts not found"<<endl;
      return;
  }
  corrCuts->SetTrackCutsNames();
  corrCuts->SetvZeroCutsNames();
  corrCuts->SetEfficiencyWeightMap(h3D);
  corrCuts->PrintAll();

  TString centr="";
  if(minC!=0 && maxC!=0) centr = Form("%.0f%.0f",minC,maxC);
  else centr = Form("%.0f%.0f",RDHFD0Corrs->GetMinCentrality(),RDHFD0Corrs->GetMaxCentrality());
  out1name+=centr;
  out2name+=centr;
  out3name+=centr;
  out4name+=centr;
  out5name+=centr;
  out6name+=centr;
  out7name+=centr;
  inname+=centr;

  // Aanalysis task    
  TString taskname="MassAndDistrAnalysis"; taskname.Prepend("D0");
  AliAnalysisTaskSED0Correlations *massD0Task = new AliAnalysisTaskSED0Correlations(taskname.Data(),RDHFD0Corrs,corrCuts);
  massD0Task->SetDebugLevel(2);
  massD0Task->SetReadMC(readMC);
  massD0Task->SetMCReconstructedTracks(recoTrMC);
  massD0Task->SetMCReconstructedD0(recoD0MC);
  massD0Task->SetEvMixing(mixing);
  massD0Task->SetFillOnlyD0D0bar(flagD0D0bar);
  massD0Task->SetSystem(system); //0=pp, 1=PbPb
  massD0Task->SetEtaForCorrel(etacorr);

//*********************
//correlation settings
//*********************

  if(standardbins==1) {
    printf("Standard bins (from D0Mass cuts object)\n");
    massD0Task->SetNPtBinsCorr(RDHFD0Corrs->GetNPtBins()); 
    massD0Task->SetPtBinsLimsCorr(RDHFD0Corrs->GetPtBinLimits());
  } else {
    Double_t ptlimits[15] = {0,0.5,1,2,3,4,5,6,7,8,12,16,20,24,9999};
    massD0Task->SetNPtBinsCorr(15);
    massD0Task->SetPtBinsLimsCorr(ptlimits);
  }
  Double_t pttreshlow[15] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t pttreshup[15] = {999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.};
  massD0Task->SetPtTreshLow(pttreshlow);
  massD0Task->SetPtTreshUp(pttreshup);
  massD0Task->PrintBinsAndLimits();

  //  massD0Task->SetRejectSDDClusters(kTRUE);
  //  massD0Task->SetWriteVariableTree(kTRUE);

  mgr->AddTask(massD0Task);
  
  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinputmassD0 = mgr->CreateContainer(inname,TChain::Class(),AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *coutputmassD01 = mgr->CreateContainer(out1name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //mass
  AliAnalysisDataContainer *coutputmassD02 = mgr->CreateContainer(out2name,TH1F::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //nev
  AliAnalysisDataContainer *coutputmassD03 = mgr->CreateContainer(out3name,AliRDHFCutsD0toKpi::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //cuts
  AliAnalysisDataContainer *coutputmassD04 = mgr->CreateContainer(out4name,AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //counter
  AliAnalysisDataContainer *coutputmassD05 = mgr->CreateContainer(out5name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //correlations
  AliAnalysisDataContainer *coutputmassD06 = mgr->CreateContainer(out6name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //MC study plots (for corrs)
  AliAnalysisDataContainer *coutputmassD07 = mgr->CreateContainer(out7name,AliHFAssociatedTrackCuts::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //cuts for tracks/K0

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
