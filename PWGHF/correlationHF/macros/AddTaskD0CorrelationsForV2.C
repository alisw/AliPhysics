AliAnalysisTaskSED0Correlations *AddTaskD0CorrelationsForV2(Bool_t readMC=kFALSE, Bool_t mixing=kFALSE, Bool_t recoTrMC=kFALSE, Bool_t recoD0MC = kFALSE,  Bool_t flagsoftpicut = kTRUE, Bool_t MEthresh = kFALSE, Bool_t pporpPb_lims=kFALSE /*0=pp,1=pPb limits*/, TString cutsfilename="D0toKpiCuts.root", TString cutsfilename2="AssocPartCuts_Std_NewPools.root", TString effD0namec="D0Eff_From_c_wLimAcc_2D.root", TString effD0nameb="D0Eff_From_b_wLimAcc_2D.root", TString effName = "3D_eff_Std.root", TString cutsD0name="D0toKpiCuts", TString cutsTrkname="AssociatedTrkCuts", Int_t system=0/*0=useMultipl(pp),1=useCentral(PbPb,pA depends)-*/, TString finDirname="Output", Int_t speed=AliAnalysisTaskSED0Correlations::kOneBinSB, Bool_t mergepools=kFALSE, Bool_t useDeff=kTRUE, Bool_t useTrackeff=kTRUE, Double_t ptAssocLim=4., Double_t minDPt=2., TString multSelEstimator="V0M", Float_t multmin=0., Float_t multmax=0., Int_t fillTrees=AliAnalysisTaskSED0Correlations::kNoTrees, Double_t fractAccME=100., Bool_t puritystudies=kFALSE, Bool_t equalizeTracklets = kFALSE, Double_t refmult = 12.524, TString sample = "pp13TeV", TString fileprofiles2016="", TString fileprofiles2017="", TString fileprofiles2018="", Bool_t reweighMC=kFALSE, TString filenameWeights="")
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
  //list mass
  out1name="coutputmassD0Mass";
  //hist entries
  out2name="nEntriesD0";
  //cuts object
  out3name="cutsD0"; 
  //AliNormalizationCounter
  out4name="normalizationCounter";
  //correlation outputs
  out5name ="correlations";
  //correlation further studies
  out6name ="debugPlots";
  //correlated trk cuts
  out7name ="cutsTracks";
  //D0 tree
  out8name ="TreeD0";
  //Assoc track tree
  out9name ="TreeTracks";

  inname="cinputmassD0_0";

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
    cout<<"Input file not found for D0 cuts!"<<endl;
    return;
  }
  TFile* filecuts2=TFile::Open(cutsfilename2.Data());
  if(!filecuts2->IsOpen()){
    cout<<"Input file not found for tracks cuts!"<<endl;
    return;
  }

  //Cuts for D0
  AliRDHFCutsD0toKpi* RDHFD0Corrs=new AliRDHFCutsD0toKpi();
  RDHFD0Corrs = (AliRDHFCutsD0toKpi*)filecuts->Get(cutsD0name.Data());
  if(!RDHFD0Corrs){
    cout<<"Specific AliRDHFCuts not found"<<endl;
    return;
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
  if(multmin!=0 || multmax!=0) centr = Form("_%.0f%.0f",multmin,multmax);
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
  massD0Task->SetFillOnlyD0D0bar(0);
  massD0Task->SetSystem(system); //0=use multiplicity (pp), 1=use centrality (PbPb). For pA you can choose how to behave
  massD0Task->SetEtaForCorrel(1.5);
  massD0Task->SetSoftPiFlag(flagsoftpicut);
  massD0Task->SetMEAxisThresh(MEthresh);
  massD0Task->SetUseDeff(useDeff); 
  massD0Task->SetUseTrackeff(useTrackeff); 
  massD0Task->SetMinDPt(minDPt);
  massD0Task->SetFillTrees(fillTrees,fractAccME);
  massD0Task->SetAODMismatchProtection(kTRUE);
  massD0Task->SetPurityStudies(puritystudies);
  if(multSelEstimator=="SPD") massD0Task->SetTrackletRange(multmin,multmax);
  if(multSelEstimator=="V0M") massD0Task->SetCentralityV0(multmin,multmax);

  massD0Task->SetAnalysisVsMult(kTRUE);
  
  // multiplicity weights
  if(reweighMC) {
      printf(" Loading Ntracklet weights...\n");
      TFile* filewgt = TFile::Open(filenameWeights);
      TH1D* hWeight = (TH1D*)filewgt->Get("hNtrUnCorrEvWithCandWeight");
      if(!hWeight) {
          printf("Error! Weights for Ntracklets not correctly loaded!");
          return;
      }
      ->SetUseNtrklWeight(kTRUE);
      massD0Task->SetHistNtrklWeight(hWeight);
  }

  //tracklet equalisation, load profiles
  if(equalizeTracklets) {
    massD0Task->SetEqualizeTracklets(kTRUE);
    massD0Task->SetReferenceMultiplicity(refmult);
    TFile* fprof1 = TFile::Open(fileprofiles2016);
    TFile* fprof2 = TFile::Open(fileprofiles2017);
    TFile* fprof3 = TFile::Open(fileprofiles2018);
    if(!fileprofiles2016 || !fileprofiles2017 || !fileprofiles2018) {
      printf("Error! Tracklet profile file not correctly loaded!");
      return;
    }

    if(sample=="pp13TeV") {
      const Char_t* periodNames[33]={"LHC16d","LHC16e","LHC16g","LHC16h","LHC16j","LHC16k","LHC16l","LHC16o","LHC16p",
                                     "LHC17c","LHC17e","LHC17f","LHC17h","LHC17i","LHC17j","LHC17k","LHC17l","LHC17m","LHC17o","LHC17r",
                                     "LHC18b","LHC18d","LHC18e","LHC18f","LHC18g","LHC18h","LHC18i","LHC18k","LHC18l","LHC18m","LHC18n","LHC18o","LHC18p"};
      TProfile *multEstimatorAvg[33];
      for(Int_t ip=0;ip<9; ip++){
        multEstimatorAvg[ip] = (TProfile*)(fprof1->Get(Form("SPDmult10_%s",periodNames[ip]))->Clone(Form("SPDmult10_%s_clone",periodNames[ip])));
        if (!multEstimatorAvg[ip]) {
          printf("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]);
          return;
        }
        massD0Task->SetMultiplVsZProfile(multEstimatorAvg[ip],ip);
      }
      for(Int_t ip=9;ip<20; ip++){
        multEstimatorAvg[ip] = (TProfile*)(fprof2->Get(Form("SPDmult10_%s",periodNames[ip]))->Clone(Form("SPDmult10_%s_clone",periodNames[ip])));
        if (!multEstimatorAvg[ip]) {
          printf("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]);
          return;
        }
        massD0Task->SetMultiplVsZProfile(multEstimatorAvg[ip],ip);        
      }
      for(Int_t ip=20;ip<33; ip++){
        multEstimatorAvg[ip] = (TProfile*)(fprof3->Get(Form("SPDmult10_%s",periodNames[ip]))->Clone(Form("SPDmult10_%s_clone",periodNames[ip])));
        if (!multEstimatorAvg[ip]) {
          printf("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]);
          return;
        }
        massD0Task->SetMultiplVsZProfile(multEstimatorAvg[ip],ip);     
      }        
      for(Int_t ip=0;ip<33;ip++) printf("Estimator for %d = %s\n",ip,multEstimatorAvg[ip]->GetName());//just a cross-check
    } else { //not pp@13 TeV
      printf("Year not configured! Exiting...\n");
      return;
    }
  } //end equalize tracklets

//*********************
//correlation settings
//*********************

  printf("Bins pT(D) from D0Mass cuts object)\n");
  Int_t nbins=RDHFD0Corrs->GetNPtBins();
  massD0Task->SetNPtBinsCorr(nbins); 
  massD0Task->SetPtBinsLimsCorr(RDHFD0Corrs->GetPtBinLimits());
  Double_t pttreshlow[50];
  Double_t pttreshup[50];
  for(int k=0;k<nbins;k++) {pttreshlow[k]=0.; pttreshup[k]=999.;}
  massD0Task->SetPtTreshLow(pttreshlow);
  massD0Task->SetPtTreshUp(pttreshup);

  //needed for fSpeed==0 and ==1 (i.e. when inv mass axis of THnSparse has the normal signal bins)
  Double_t LeftSignReg_LowPt = 1.7968;
  Double_t RightSignReg_LowPt = 1.9528;
  Double_t LeftSignReg_HighPt = 1.7488;
  Double_t RightSignReg_HighPt = 2.0008;

  massD0Task->SetLeftSignReg_LowPt(LeftSignReg_LowPt);
  massD0Task->SetRightSignReg_LowPt(RightSignReg_LowPt);
  massD0Task->SetLeftSignReg_HighPt(LeftSignReg_HighPt);
  massD0Task->SetRightSignReg_HighPt(RightSignReg_HighPt); 
 
  //needed for fSpeed==2

  //use SB ranges from cut file
  TVectorD *LSBLow = (TVectorD*)filecuts2->Get("vLSBLow");
  TVectorD *LSBUpp = (TVectorD*)filecuts2->Get("vLSBUpp");
  TVectorD *RSBLow = (TVectorD*)filecuts2->Get("vRSBLow");
  TVectorD *RSBUpp = (TVectorD*)filecuts2->Get("vRSBUpp");      

  if(!LSBLow||!LSBUpp||!RSBLow||!RSBUpp) {printf("Error! No SB ranges found in the Associated track cut file, but useCutFileMassRanges==kTRUE! Exiting...\n"); return;}

  massD0Task->SetLSBLowLim(LSBLow->GetMatrixArray());
  massD0Task->SetLSBHighLim(LSBUpp->GetMatrixArray());
  massD0Task->SetRSBLowLim(RSBLow->GetMatrixArray());
  massD0Task->SetRSBHighLim(RSBUpp->GetMatrixArray());

  if(speed==2) {
    TVectorD *SignLow = (TVectorD*)filecuts2->Get("vSignLow");
    TVectorD *SignUpp = (TVectorD*)filecuts2->Get("vSignUpp");

    if(!SignLow||!SignUpp) {printf("Error! No Signal ranges found in the Associated track cut file, but useCutFileMassRanges==kTRUE, and fSpeed==2! Exiting...\n"); return;}

    massD0Task->SetSignLowLim(SignLow->GetMatrixArray());
    massD0Task->SetSignHighLim(SignUpp->GetMatrixArray());
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
