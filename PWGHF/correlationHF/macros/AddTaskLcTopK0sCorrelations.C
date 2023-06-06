AliAnalysisTaskSELcTopK0sCorrelations *AddTaskLcTopK0sCorrelations(Bool_t readMC=kFALSE, Bool_t mixing=kFALSE, Bool_t recoTrMC=kFALSE, Bool_t recoLambdacMC = kFALSE, Bool_t MEthresh = kFALSE, Bool_t pporpPb_lims=kFALSE /*0=pp,1=pPb limits*/, TString cutsfilename="cuts_pidOFF.root", TString cutsfilename2="AssocPartCuts_Std_NewPools.root", TString effLambdacnamec="LambdacEff_From_c_wLimAcc_2D.root", TString effLambdacnameb="LambdacEff_From_b_wLimAcc_2D.root", TString effName = "3D_eff_Std.root", TString cutsLambdacname="LctopK0sAnalysisCuts", TString cutsTrkname="AssociatedTrkCuts", Double_t etacorr=1.5, Int_t system=0/*0=useMultipl(pp),1=useCentral(PbPb,pA depends)-*/, Int_t flagLambdacLambdacbar=0, Float_t minC=0, Float_t maxC=0, TString finDirname="Output", Bool_t flagAOD049=kFALSE, Int_t standardbins=1, Bool_t stdcuts=kFALSE, Bool_t analyszeKaon=kFALSE, Int_t speed=AliAnalysisTaskSELcTopK0sCorrelations::kOneBinSB, Bool_t mergepools=kFALSE, Bool_t UseLceff=kTRUE, Bool_t useTrackeff=kTRUE, Bool_t useCutFileMassRanges=kFALSE, Double_t ptAssocLim=1., Int_t fillTrees=AliAnalysisTaskSELcTopK0sCorrelations::kNoTrees, Double_t fractAccME=100., Double_t minDPt=2., Int_t AODprot=0, Bool_t puritystudies=kFALSE, Bool_t reweighMC=kFALSE, TString filenameWeights="",Bool_t flag=kFALSE, Bool_t isMultiClass = kFALSE, TString path = "", Bool_t flagsoftpicut = kFALSE)
{
  //
  // AddTask for the AliAnalysisTaskSE for Lambdac candidates
  // invariant mass histogram and association with MC truth 
  // (using MC info in AOD) and cut variables distributions
  // Antonio Palasciano  antonio.palasciano@ba.infn.it 


  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskLambdacDistr", "No analysis manager to connect to.");
    return NULL;
  }   

  TString filename="",out1name="",out2name="",out3name="",out4name="",out5name="",out6name="",out7name="",out8name="",out9name="",inname="";
  filename = AliAnalysisManager::GetCommonFileName();
  filename += ":";
  filename+="LambdachCorrel";
  if(flagLambdacLambdacbar==1)filename+="Lambdac";
  if(flagLambdacLambdacbar==2)filename+="Lambdacbar";
  //list mass
  out1name="coutputmassLambdacMass";
  if(flagLambdacLambdacbar==1)out1name+="Lambdac";
  if(flagLambdacLambdacbar==2)out1name+="Lambdacbar";
  //hist entries
  out2name="nEntriesLambdac";
  if(flagLambdacLambdacbar==1)out2name+="Lambdac";
  if(flagLambdacLambdacbar==2)out2name+="Lambdacbar";
  //cuts object
  out3name="cutsLambdac"; 
  if(flagLambdacLambdacbar==1)out3name+="Lambdac";
  if(flagLambdacLambdacbar==2)out3name+="Lambdacbar";
  //AliNormalizationCounter
  out4name="normalizationCounter";
  if(flagLambdacLambdacbar==1)out4name+="Lambdac";
  if(flagLambdacLambdacbar==2)out4name+="Lambdacbar";
  //correlation outputs
  out5name ="correlations";
  if(flagLambdacLambdacbar==1)out5name+="Lambdac";
  if(flagLambdacLambdacbar==2)out5name+="Lambdacbar";
  //correlation further studies
  out6name ="debugPlots";
  if(flagLambdacLambdacbar==1)out6name+="Lambdac";
  if(flagLambdacLambdacbar==2)out6name+="Lambdacbar";
  //correlated trk cuts
  out7name ="cutsTracks";
  if(flagLambdacLambdacbar==1)out7name+="Lambdac";
  if(flagLambdacLambdacbar==2)out7name+="Lambdacbar";
  //Lambdac tree
  out8name ="TreeLambdac";
  if(flagLambdacLambdacbar==1)out8name+="Lambdac";
  if(flagLambdacLambdacbar==2)out8name+="Lambdacbar";
  //Assoc track tree
  out9name ="TreeTracks";
  if(flagLambdacLambdacbar==1)out9name+="Lambdac";
  if(flagLambdacLambdacbar==2)out9name+="Lambdacbar";

  inname="cinputmassLambdac_0";
  if(flagLambdacLambdacbar==1)inname+="Lambdac";
  if(flagLambdacLambdacbar==2)inname+="Lambdacbar";

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
    //     printf("    |M-MLambdac| [GeV]    < %f\n",fLambdacCorrCuts[0]);
    //     printf("    dca    [cm]  < %f\n",fLambdacCorrCuts[1]);
    //     printf("    cosThetaStar     < %f\n",fLambdacCorrCuts[2]);
    //     printf("    pTK     [GeV/c]    > %f\n",fLambdacCorrCuts[3]);
    //     printf("    pTpi    [GeV/c]    > %f\n",fLambdacCorrCuts[4]);
    //     printf("    |LambdacK|  [cm]  < %f\n",fLambdacCorrCuts[5]);
    //     printf("    |Lambdacpi| [cm]  < %f\n",fLambdacCorrCuts[6]);
    //     printf("    LambdacLambdac  [cm^2] < %f\n",fLambdacCorrCuts[7]);
    //     printf("    cosThetaPoint    > %f\n",fLambdacCorrCuts[8]);

  TFile* filecuts=TFile::Open(cutsfilename.Data());
  if(!filecuts->IsOpen()){
    std::cout<<"Input file not found for Lambdac cuts: using std cut object"<<std::endl;
    stdcuts=kTRUE;
  }
  TFile* filecuts2=TFile::Open(cutsfilename2.Data());
  if(!filecuts2->IsOpen()){
    std::cout<<"Input file not found for tracks cuts!"<<std::endl;
    return 0;
  }

  //Cuts for Lambdac
  AliRDHFCutsLctoV0* RDHFLambdacCorrs=new AliRDHFCutsLctoV0();
//RDHFLambdacCorrs->SetUsePhysicsSelection(kFALSE);
 // RDHFLambdacCorrs->SetLowPt(kFALSE); //low-pt special PID disabled
  if(stdcuts) {
    if(system==0) RDHFLambdacCorrs->SetStandardCutsPP2010();
    else {
      RDHFLambdacCorrs->SetStandardCutsPbPb2010();
      if(minC!=0 || maxC!=0) { //if centrality 0 and 0 leave the values in the cut object
	      RDHFLambdacCorrs->SetMinCentrality(minC);
	      RDHFLambdacCorrs->SetMaxCentrality(maxC);
      }
      if(flagAOD049)RDHFLambdacCorrs->SetUseAOD049(kTRUE);
    }
  }
  else {
    RDHFLambdacCorrs = (AliRDHFCutsLctoV0*)filecuts->Get(cutsLambdacname.Data());
    if(!RDHFLambdacCorrs){
      std::cout<<"Specific AliRDHFCuts not found"<<std::endl;
      return 0;
    }
    if(flagAOD049)RDHFLambdacCorrs->SetUseAOD049(kTRUE);
    if(minC!=0 || maxC!=0) { //if centrality 0 and 0 leave the values in the cut object
      RDHFLambdacCorrs->SetMinCentrality(minC);
      RDHFLambdacCorrs->SetMaxCentrality(maxC); //****you need to define the estimator in the cut file, though!!!****
    } 
  }

  //Cuts for correlated tracks/K0
  AliHFAssociatedTrackCuts* corrCuts=new AliHFAssociatedTrackCuts();
  corrCuts = (AliHFAssociatedTrackCuts*)filecuts2->Get(cutsTrkname.Data());
  if(!corrCuts){
      std::cout<<"Specific AliHFAssociatedTrackCuts not found"<<std::endl;
      return 0;
  }
  corrCuts->SetTrackCutsNames();
  corrCuts->SetvZeroCutsNames();

  //Load efficiency map
  if(!effName.EqualTo("")) {
    TFile* fileeff=TFile::Open(effName.Data());
    if(!fileeff->IsOpen()){
      std::cout<<"Input file not found for efficiency! Exiting..."<<std::endl;
      return 0;
    }  
    TCanvas *c = (TCanvas*)fileeff->Get("c");
    TH3D *h3D = (TH3D*)c->FindObject("heff_rebin");
    if(recoLambdacMC) corrCuts->SetEfficiencyWeightMap(h3D); //data and MC Reco
  } else std::cout<<"*** WARNING! No tracking efficiency map set! ***"<<std::endl;

  //Load Lambdac efficiency map
  if(!effLambdacnamec.EqualTo("")) {
    TFile* fileeffLambdacc=TFile::Open(effLambdacnamec.Data());
    if(!fileeffLambdacc->IsOpen()){
      std::cout<<"Input file not found for efficiency! Exiting..."<<std::endl;
      return 0;
    }
    TH2D *hEffLambdacc = (TH2D*)fileeffLambdacc->Get("h_Eff");
    if(recoLambdacMC) corrCuts->SetTriggerEffWeightMap(hEffLambdacc); //data and MC Reco
  } else std::cout<<"*** WARNING! No prompt trigger efficiency map set! ***"<<std::endl;

  //Load Lambdac efficiency map from b
  if(readMC) {
    if(!effLambdacnameb.EqualTo("")) {
      TFile* fileeffLambdacb=TFile::Open(effLambdacnameb.Data());
      if(!fileeffLambdacb->IsOpen()){
        std::cout<<"Input file not found for efficiency! Exiting..."<<std::endl;
        return 0;
      }
      TH2D *hEffLambdacb = (TH2D*)fileeffLambdacb->Get("h_Eff");
      if(recoLambdacMC && readMC) corrCuts->SetTriggerEffWeightMapB(hEffLambdacb); //MC Reco
    } else std::cout<<"*** WARNING! No feed-down trigger efficiency map set! ***"<<std::endl;
  }

  corrCuts->PrintAll();

  TString centr="";
  if(minC!=0 || maxC!=0) centr = Form("_%.0f%.0f",minC,maxC);
  else centr = Form("_%.0f%.0f",RDHFLambdacCorrs->GetMinCentrality(),RDHFLambdacCorrs->GetMaxCentrality());
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
  TString taskname="MassAndDistrAnalysis"; taskname.Prepend("Lambdac");
  AliAnalysisTaskSELcTopK0sCorrelations *massLambdacTask = new AliAnalysisTaskSELcTopK0sCorrelations(taskname.Data(),RDHFLambdacCorrs,corrCuts);
  massLambdacTask->SetDebugLevel(2);
  massLambdacTask->SetReadMC(readMC);
  massLambdacTask->SetMCReconstructedTracks(recoTrMC);
  massLambdacTask->SetMCReconstructedLambdac(recoLambdacMC);
  massLambdacTask->SetEvMixing(mixing);
  massLambdacTask->SetFillOnlyLambdacLambdacbar(flagLambdacLambdacbar);
  massLambdacTask->SetSystem(system); //0=use multiplicity (pp), 1=use centrality (PbPb). For pA you can choose how to behave
  massLambdacTask->SetEtaForCorrel(etacorr);
  massLambdacTask->SetSoftPiFlag(flagsoftpicut);
  massLambdacTask->SetMEAxisThresh(MEthresh);
  massLambdacTask->SetUseLceff(UseLceff); 
  massLambdacTask->SetUseTrackeff(useTrackeff); 
  massLambdacTask->SetMinDPt(minDPt);
  massLambdacTask->SetFillTrees(fillTrees,fractAccME);
  massLambdacTask->SetAODMismatchProtection(AODprot);
  massLambdacTask->SetPurityStudies(puritystudies);
  if(analyszeKaon) massLambdacTask->SetKaonCorrelations(kTRUE);

  // multiplicity weights
  if(reweighMC) {
      printf(" Loading Ntracklet weights...\n");
      TFile* filewgt = TFile::Open(filenameWeights);
      TH1D* hWeight = (TH1D*)filewgt->Get("hNtrUnCorrEvWithCandWeight");
      if(!hWeight) {
          printf("Error! Weights for Ntracklets not correctly loaded!");
          return 0;
      }
      massLambdacTask->SetUseNtrklWeight(kTRUE);
      massLambdacTask->SetHistNtrklWeight(hWeight);
  }


//*********************
//correlation settings
//*********************

  if(standardbins==1) {
    printf("Standard bins (from LambdacMass cuts object)\n");
    Int_t nbins=RDHFLambdacCorrs->GetNPtBins();
    massLambdacTask->SetNPtBinsCorr(nbins); 
    massLambdacTask->SetPtBinsLimsCorr(RDHFLambdacCorrs->GetPtBinLimits());
    Double_t pttreshlow[50];
    Double_t pttreshup[50];
    for(int k=0;k<nbins;k++) {pttreshlow[k]=0.; pttreshup[k]=999.;}
    massLambdacTask->SetPtTreshLow(pttreshlow);
    massLambdacTask->SetPtTreshUp(pttreshup);
  } else {
    Double_t ptlimits[15] = {0,0.5,1,2,3,4,5,6,7,8,12,16,20,24,9999};
    massLambdacTask->SetNPtBinsCorr(14);
    massLambdacTask->SetPtBinsLimsCorr(ptlimits);
    Double_t pttreshlow[15] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    Double_t pttreshup[15] = {999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.,999.};
    massLambdacTask->SetPtTreshLow(pttreshlow);
    massLambdacTask->SetPtTreshUp(pttreshup);
  }


  //needed for fSpeed==0 and ==1 (i.e. when inv mass axis of THnSparse has the normal signal bins)
  Double_t LeftSignReg_LowPt = 2.2224;
  Double_t RightSignReg_LowPt = 2.3944;
  Double_t LeftSignReg_HighPt = 2.2224;
  Double_t RightSignReg_HighPt = 2.3944;

  massLambdacTask->SetLeftSignReg_LowPt(LeftSignReg_LowPt);
  massLambdacTask->SetRightSignReg_LowPt(RightSignReg_LowPt);
  massLambdacTask->SetLeftSignReg_HighPt(LeftSignReg_HighPt);
  massLambdacTask->SetRightSignReg_HighPt(RightSignReg_HighPt); 
 
  //needed for fSpeed==2

  if(useCutFileMassRanges) { //use SB ranges from cut file

    TVectorD *LSBLow = (TVectorD*)filecuts2->Get("vLSBLow");
    TVectorD *LSBUpp = (TVectorD*)filecuts2->Get("vLSBUpp");
    TVectorD *RSBLow = (TVectorD*)filecuts2->Get("vRSBLow");
    TVectorD *RSBUpp = (TVectorD*)filecuts2->Get("vRSBUpp");      

    if(!LSBLow||!LSBUpp||!RSBLow||!RSBUpp) {printf("Error! No SB ranges found in the Associated track cut file, but useCutFileMassRanges==kTRUE! Exiting...\n"); return 0;}

    massLambdacTask->SetLSBLowLim(LSBLow->GetMatrixArray());
    massLambdacTask->SetLSBHighLim(LSBUpp->GetMatrixArray());
    massLambdacTask->SetRSBLowLim(RSBLow->GetMatrixArray());
    massLambdacTask->SetRSBHighLim(RSBUpp->GetMatrixArray());

    if(speed==2) {

      TVectorD *SignLow = (TVectorD*)filecuts2->Get("vSignLow");
      TVectorD *SignUpp = (TVectorD*)filecuts2->Get("vSignUpp");

      if(!SignLow||!SignUpp) {printf("Error! No Signal ranges found in the Associated track cut file, but useCutFileMassRanges==kTRUE, and fSpeed==2! Exiting...\n"); return 0;}

      massLambdacTask->SetSignLowLim(SignLow->GetMatrixArray());
      massLambdacTask->SetSignHighLim(SignUpp->GetMatrixArray());

    }

  } else { //use SB ranges from AddTask (following settings)

    if(speed==2) {printf("Error! fSpeed 2 set with useCutFileMassRanges==kTRUE, this is not allowed!\nYou have to pass the signal ranges via cut file with fSpeed=2! Exiting...\n"); return 0;}

    if(!pporpPb_lims) { //pp limits
	              //         3-4    4-5    5-6    6-7    7-8   8-12   
      Double_t LSBLowLim[7] = {2.2144,2.2144,2.1904,2.2144,2.2064,2.1664,2.1664}; //to be filled looking at results from invariant mass fits!
      Double_t LSBUppLim[7] = {2.2584,2.2584,2.2464,2.2544,2.2504,2.2384,2.2384};
      Double_t RSBLowLim[7] = {2.3264,2.3264,2.3304,2.3264,2.3304,2.3424,2.3424};
      Double_t RSBUppLim[7] = {2.3704,2.3704,2.3904,2.3764,2.3744,2.4104,2.4104};
 
      massLambdacTask->SetLSBLowLim(LSBLowLim);
      massLambdacTask->SetLSBHighLim(LSBUppLim);
      massLambdacTask->SetRSBLowLim(RSBLowLim);
      massLambdacTask->SetRSBHighLim(RSBUppLim);

    } else { //pPb limits
				//      1-2    2-3    3-4    4-5    5-6    6-7    7-8   8-12   12-16  16-20  20-24   24+
      Double_t LSBLowLim[14] = {0.,0.,1.7928,1.7928,1.7768,1.7728,1.7648,1.7488,1.7448,1.7728,1.7048,1.7048,1.7048,1.7048}; //to be filled looking at results from invariant mass fits!
      Double_t LSBUppLim[14] = {0.,0.,1.8288,1.8288,1.8208,1.8208,1.8128,1.8088,1.8048,1.8048,1.7568,1.7568,1.7568,1.7568};
      Double_t RSBLowLim[14] = {0.,0.,1.9008,1.9008,1.9088,1.9128,1.9168,1.9288,1.9288,1.9288,1.9728,1.9728,1.9728,1.9728};
      Double_t RSBUppLim[14] = {0.,0.,1.9408,1.9408,1.9528,1.9608,1.9688,1.9848,1.9888,1.9928,2.0768,2.0768,2.0768,2.0768}; 
 
      massLambdacTask->SetLSBLowLim(LSBLowLim);
      massLambdacTask->SetLSBHighLim(LSBUppLim);
      massLambdacTask->SetRSBLowLim(RSBLowLim);
      massLambdacTask->SetRSBHighLim(RSBUppLim);
    }
  }

  //  massLambdacTask->SetRejectSDDClusters(kTRUE);
  //  massLambdacTask->SetWriteVariableTree(kTRUE);

  massLambdacTask->SetPtAssocLim(ptAssocLim);
  massLambdacTask->SetSpeed(speed);
  massLambdacTask->SetMergePools(mergepools);
  massLambdacTask->PrintBinsAndLimits();

  massLambdacTask->SetDoMLApplication(flag, isMultiClass);
  massLambdacTask->SetMLConfigFile(path);
  
  mgr->AddTask(massLambdacTask);
  
  //
  // Create containers for input/output
  AliAnalysisDataContainer *cinputmassLambdac = mgr->CreateContainer(inname,TChain::Class(),AliAnalysisManager::kInputContainer);

  AliAnalysisDataContainer *coutputmassLambdac1 = mgr->CreateContainer(out1name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //mass
  AliAnalysisDataContainer *coutputmassLambdac2 = mgr->CreateContainer(out2name,TH1F::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //nev
  AliAnalysisDataContainer *coutputmassLambdac3 = mgr->CreateContainer(out3name,AliRDHFCutsLctoV0::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //cuts
  AliAnalysisDataContainer *coutputmassLambdac4 = mgr->CreateContainer(out4name,AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //counter
  AliAnalysisDataContainer *coutputmassLambdac5 = mgr->CreateContainer(out5name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //correlations
  AliAnalysisDataContainer *coutputmassLambdac6 = mgr->CreateContainer(out6name,TList::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //MC study plots (for corrs)
  AliAnalysisDataContainer *coutputmassLambdac7 = mgr->CreateContainer(out7name,AliHFAssociatedTrackCuts::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //cuts for tracks/K0
  AliAnalysisDataContainer *coutputmassLambdac8 = mgr->CreateContainer(out8name,TTree::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //TTree Lambdac
  AliAnalysisDataContainer *coutputmassLambdac9 = mgr->CreateContainer(out9name,TTree::Class(),AliAnalysisManager::kOutputContainer, filename.Data()); //TTree Tracks

  mgr->ConnectInput(massLambdacTask,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(massLambdacTask,1,coutputmassLambdac1);
  mgr->ConnectOutput(massLambdacTask,2,coutputmassLambdac2);
  mgr->ConnectOutput(massLambdacTask,3,coutputmassLambdac3);
  mgr->ConnectOutput(massLambdacTask,4,coutputmassLambdac4);
  mgr->ConnectOutput(massLambdacTask,5,coutputmassLambdac5);
  mgr->ConnectOutput(massLambdacTask,6,coutputmassLambdac6);
  mgr->ConnectOutput(massLambdacTask,7,coutputmassLambdac7);
  mgr->ConnectOutput(massLambdacTask,8,coutputmassLambdac8);
  mgr->ConnectOutput(massLambdacTask,9,coutputmassLambdac9);

  return massLambdacTask;
}
