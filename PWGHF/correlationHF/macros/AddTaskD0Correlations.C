AliAnalysisTaskSED0Correlations *AddTaskD0Correlations(Bool_t readMC=kFALSE, Bool_t mixing=kFALSE, Bool_t recoTrMC=kFALSE, Bool_t recoD0MC = kFALSE,  Bool_t flagsoftpicut = kTRUE, Bool_t MEthresh = kFALSE, Bool_t pporpPb_lims=kFALSE /*0=pp,1=pPb limits*/, TString cutsfilename="D0toKpiCuts.root", TString cutsfilename2="AssocPartCuts_Std_NewPools.root", TString effD0namec="D0Eff_From_c_wLimAcc_2D.root", TString effD0nameb="D0Eff_From_b_wLimAcc_2D.root", TString effName = "3D_eff_Std.root", TString cutsD0name="D0toKpiCuts", TString cutsTrkname="AssociatedTrkCuts", Double_t etacorr=1.5, Int_t system=0/*0=useMultipl(pp),1=useCentral(PbPb,pA depends)-*/, Int_t flagD0D0bar=0, Float_t minC=0, Float_t maxC=0, TString finDirname="Output", Bool_t flagAOD049=kFALSE, Int_t standardbins=1, Bool_t stdcuts=kFALSE, Bool_t analyszeKaon=kFALSE, Bool_t speed=kTRUE, Bool_t mergepools=kFALSE, Bool_t useDeff=kTRUE, Bool_t useTrackeff=kTRUE, Bool_t useCutFileSBRanges=kFALSE, Double_t ptAssocLim=1., Int_t fillTrees=AliAnalysisTaskSED0Correlations::kNoTrees, Double_t fractAccME=100., Double_t minDPt=2., Int_t AODprot=1, Bool_t puritystudies=kFALSE)
{
  //
  // AddTask for the AliAnalysisTaskSE for D0 candidates
  // invariant mass histogram and association with MC truth 
  // (using MC info in AOD) and cut variables distributions
  // F.Colamaria fabio.colamaria@ba.infn.it


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskD0Distr", "No analysis manager to connect to.");
    return NULL;
  }   

  TString filename="",out1name="",out2name="",out3name="",out4name="",out5name="",out6name="",out7name="",out8name="",out9name="",inname="";
  filename = AliAnalysisManager::GetCommonFileName();
  filename += ":";
  filename+="D0hCorrel";
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
  out6name ="debugPlots";
  if(flagD0D0bar==1)out6name+="D0";
  if(flagD0D0bar==2)out6name+="D0bar";
  //correlated trk cuts
  out7name ="cutsTracks";
  if(flagD0D0bar==1)out7name+="D0";
  if(flagD0D0bar==2)out7name+="D0bar";
  //D0 tree
  out8name ="TreeD0";
  if(flagD0D0bar==1)out8name+="D0";
  if(flagD0D0bar==2)out8name+="D0bar";
  //Assoc track tree
  out9name ="TreeTracks";
  if(flagD0D0bar==1)out9name+="D0";
  if(flagD0D0bar==2)out9name+="D0bar";

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
  out8name += finDirname.Data();
  out9name += finDirname.Data();
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

  TFile* filecuts=TFile::Open(cutsfilename.Data());
  if(!filecuts->IsOpen()){
    cout<<"Input file not found for D0 cuts: using std cut object"<<endl;
    stdcuts=kTRUE;
  }
  TFile* filecuts2=TFile::Open(cutsfilename2.Data());
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
      if(minC!=0 || maxC!=0) { //if centrality 0 and 0 leave the values in the cut object
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
    if(minC!=0 || maxC!=0) { //if centrality 0 and 0 leave the values in the cut object
      RDHFD0Corrs->SetMinCentrality(minC);
      RDHFD0Corrs->SetMaxCentrality(maxC);
    } 
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

  //Load efficiency map
  if(!effName.EqualTo("")) {
    TFile* fileeff=TFile::Open(effName.Data());
    if(!fileeff->IsOpen()){
      cout<<"Input file not found for efficiency! Exiting..."<<endl;
      return;
    }  
    TCanvas *c = (TCanvas*)fileeff->Get("c");
    TH3D *h3D = (TH3D*)c->FindObject("heff_rebin");
    if(recoD0MC) corrCuts->SetEfficiencyWeightMap(h3D); //data and MC Reco
  } else cout<<"*** WARNING! No tracking efficiency map set! ***"<<endl;

  //Load D0 efficiency map
  if(!effD0namec.EqualTo("")) {
    TFile* fileeffD0c=TFile::Open(effD0namec.Data());
    if(!fileeffD0c->IsOpen()){
      cout<<"Input file not found for efficiency! Exiting..."<<endl;
      return;
    }
    TH2D *hEffD0c = (TH2D*)fileeffD0c->Get("h_Eff");
    if(recoD0MC) corrCuts->SetTriggerEffWeightMap(hEffD0c); //data and MC Reco
  } else cout<<"*** WARNING! No prompt trigger efficiency map set! ***"<<endl;

  //Load D0 efficiency map from b
  if(readMC) {
    if(!effD0nameb.EqualTo("")) {
      TFile* fileeffD0b=TFile::Open(effD0nameb.Data());
      if(!fileeffD0b->IsOpen()){
        cout<<"Input file not found for efficiency! Exiting..."<<endl;
        return;
      }
      TH2D *hEffD0b = (TH2D*)fileeffD0b->Get("h_Eff");
      if(recoD0MC && readMC) corrCuts->SetTriggerEffWeightMapB(hEffD0b); //MC Reco
    } else cout<<"*** WARNING! No feed-down trigger efficiency map set! ***"<<endl;
  }

  corrCuts->PrintAll();

  TString centr="";
  if(minC!=0 || maxC!=0) centr = Form("_%.0f%.0f",minC,maxC);
  else centr = Form("_%.0f%.0f",RDHFD0Corrs->GetMinCentrality(),RDHFD0Corrs->GetMaxCentrality());
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
  TString taskname="MassAndDistrAnalysis"; taskname.Prepend("D0");
  AliAnalysisTaskSED0Correlations *massD0Task = new AliAnalysisTaskSED0Correlations(taskname.Data(),RDHFD0Corrs,corrCuts);
  massD0Task->SetDebugLevel(2);
  massD0Task->SetReadMC(readMC);
  massD0Task->SetMCReconstructedTracks(recoTrMC);
  massD0Task->SetMCReconstructedD0(recoD0MC);
  massD0Task->SetEvMixing(mixing);
  massD0Task->SetFillOnlyD0D0bar(flagD0D0bar);
  massD0Task->SetSystem(system); //0=use multiplicity (pp), 1=use centrality (PbPb). For pA you can choose how to behave
  massD0Task->SetEtaForCorrel(etacorr);
  massD0Task->SetSoftPiFlag(flagsoftpicut);
  massD0Task->SetMEAxisThresh(MEthresh);
  massD0Task->SetUseDeff(useDeff); 
  massD0Task->SetUseTrackeff(useTrackeff); 
  massD0Task->SetMinDPt(minDPt);
  massD0Task->SetFillTrees(fillTrees,fractAccME);
  massD0Task->SetAODMismatchProtection(AODprot);
  massD0Task->SetPurityStudies(puritystudies);
  if(analyszeKaon) massD0Task->SetKaonCorrelations(kTRUE);

//*********************
//correlation settings
//*********************

  if(standardbins==1) {
    printf("Standard bins (from D0Mass cuts object)\n");
    massD0Task->SetNPtBinsCorr(RDHFD0Corrs->GetNPtBins()); 
    massD0Task->SetPtBinsLimsCorr(RDHFD0Corrs->GetPtBinLimits());
  } else {
    Double_t ptlimits[15] = {0,0.5,1,2,3,4,5,6,7,8,12,16,20,24,9999};
    massD0Task->SetNPtBinsCorr(14);
    massD0Task->SetPtBinsLimsCorr(ptlimits);
  }
  Double_t pttreshlow[15] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
  Double_t pttreshup[15] = {999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.};
  massD0Task->SetPtTreshLow(pttreshlow);
  massD0Task->SetPtTreshUp(pttreshup);

  Double_t LeftSignReg_LowPt = 1.7968;
  Double_t RightSignReg_LowPt = 1.9528;
  Double_t LeftSignReg_HighPt = 1.7488;
  Double_t RightSignReg_HighPt = 2.0008;

  massD0Task->SetLeftSignReg_LowPt(LeftSignReg_LowPt);
  massD0Task->SetRightSignReg_LowPt(RightSignReg_LowPt);
  massD0Task->SetLeftSignReg_HighPt(LeftSignReg_HighPt);
  massD0Task->SetRightSignReg_HighPt(RightSignReg_HighPt); 
  
  if(useCutFileSBRanges) { //use SB ranges from cut file

    TVectorD *LSBLow = (TVectorD*)filecuts2->Get("vLSBLow");
    TVectorD *LSBUpp = (TVectorD*)filecuts2->Get("vLSBUpp");
    TVectorD *RSBLow = (TVectorD*)filecuts2->Get("vRSBLow");
    TVectorD *RSBUpp = (TVectorD*)filecuts2->Get("vRSBUpp");      

    if(!LSBLow||!LSBUpp||!RSBLow||!RSBUpp) {printf("Error! No SB ranges found in the Associated track cut file, but useCutFileSBRanges==kTRUE! Exiting...\n"); return;}

    massD0Task->SetLSBLowLim(LSBLow->GetMatrixArray());
    massD0Task->SetLSBHighLim(LSBUpp->GetMatrixArray());
    massD0Task->SetRSBLowLim(RSBLow->GetMatrixArray());
    massD0Task->SetRSBHighLim(RSBUpp->GetMatrixArray());

  } else { //use SB ranges from AddTask (following settings)
    if(!pporpPb_lims) { //pp limits
				//      1-2    2-3    3-4    4-5    5-6    6-7    7-8   8-12   12-16  16-20  20-24   24+
      Double_t LSBLowLim[14] = {0.,0.,1.7688,1.7688,1.7488,1.7368,1.7088,1.7168,1.7168,1.7008,1.7088,1.7088,1.7088,1.7088}; //to be filled looking at results from invariant mass fits!
      Double_t LSBUppLim[14] = {0.,0.,1.8168,1.8168,1.8088,1.8008,1.7888,1.7928,1.7928,1.7528,1.7648,1.7648,1.7648,1.7648};
      Double_t RSBLowLim[14] = {0.,0.,1.9168,1.9168,1.9248,1.9288,1.9448,1.9448,1.9488,1.9728,1.9768,1.9768,1.9768,1.9768};
      Double_t RSBUppLim[14] = {0.,0.,1.9688,1.9688,1.9848,1.9928,2.0248,2.0208,2.0248,2.0848,2.0808,2.0808,2.0808,2.0808};
 
      massD0Task->SetLSBLowLim(LSBLowLim);
      massD0Task->SetLSBHighLim(LSBUppLim);
      massD0Task->SetRSBLowLim(RSBLowLim);
      massD0Task->SetRSBHighLim(RSBUppLim);

    } else { //pPb limits
				//      1-2    2-3    3-4    4-5    5-6    6-7    7-8   8-12   12-16  16-20  20-24   24+
      Double_t LSBLowLim[14] = {0.,0.,1.7928,1.7928,1.7768,1.7728,1.7648,1.7488,1.7448,1.7728,1.7048,1.7048,1.7048,1.7048}; //to be filled looking at results from invariant mass fits!
      Double_t LSBUppLim[14] = {0.,0.,1.8288,1.8288,1.8208,1.8208,1.8128,1.8088,1.8048,1.8048,1.7568,1.7568,1.7568,1.7568};
      Double_t RSBLowLim[14] = {0.,0.,1.9008,1.9008,1.9088,1.9128,1.9168,1.9288,1.9288,1.9288,1.9728,1.9728,1.9728,1.9728};
      Double_t RSBUppLim[14] = {0.,0.,1.9408,1.9408,1.9528,1.9608,1.9688,1.9848,1.9888,1.9928,2.0768,2.0768,2.0768,2.0768}; 
 
      massD0Task->SetLSBLowLim(LSBLowLim);
      massD0Task->SetLSBHighLim(LSBUppLim);
      massD0Task->SetRSBLowLim(RSBLowLim);
      massD0Task->SetRSBHighLim(RSBUppLim);
    }
  }

  //  massD0Task->SetRejectSDDClusters(kTRUE);
  //  massD0Task->SetWriteVariableTree(kTRUE);

  massD0Task->SetPtAssocLim(ptAssocLim);
  massD0Task->SetSpeed(speed);
  massD0Task->SetMergePools(mergepools);
  massD0Task->PrintBinsAndLimits();

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
  AliAnalysisDataContainer *coutputmassD08 = mgr->CreateContainer(out8name,TTree::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //TTree D0
  AliAnalysisDataContainer *coutputmassD09 = mgr->CreateContainer(out9name,TTree::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //TTree Tracks

  mgr->ConnectInput(massD0Task,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(massD0Task,1,coutputmassD01);
  mgr->ConnectOutput(massD0Task,2,coutputmassD02);
  mgr->ConnectOutput(massD0Task,3,coutputmassD03);
  mgr->ConnectOutput(massD0Task,4,coutputmassD04);
  mgr->ConnectOutput(massD0Task,5,coutputmassD05);
  mgr->ConnectOutput(massD0Task,6,coutputmassD06);
  mgr->ConnectOutput(massD0Task,7,coutputmassD07);
  mgr->ConnectOutput(massD0Task,8,coutputmassD08);
  mgr->ConnectOutput(massD0Task,9,coutputmassD09);

  return massD0Task;
}
