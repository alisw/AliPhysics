AliAnalysisTaskSEDvsMultiplicity_BDT *AddTaskDvsMultiplicity_BDT(Int_t system=0,
							 Bool_t readMC=kFALSE,
							 Int_t MCOption=0,
							 Int_t pdgMeson=4122,
							 TString finDirname="Lc2pK0S",
							 TString filename="",
							 TString finAnObjname="LctoV0AnalysisCuts",
							 TString estimatorFilename="",
							 Double_t refMult=9.26,
							 Bool_t subtractDau=kFALSE,
							 Int_t NchWeight=0,
							 Int_t recoEstimator = AliAnalysisTaskSEDvsMultiplicity_BDT::kNtrk10,
							 Int_t MCEstimator = AliAnalysisTaskSEDvsMultiplicity_BDT::kEta10,
							 Bool_t isPPbData=kFALSE,
							 Int_t year = 16,
                                                         Bool_t isLcV0=kTRUE,
                                                         Bool_t applyBDT = 0,
                                                         TString BDTfilename="", TString BDTobjnamepre= "BDT",
                                                         Float_t BDTRespCut = -1., Bool_t DoSidebndSample=kFALSE,
                                                         Bool_t GetRespTree = kTRUE, Float_t SBndSampleFrac = 0.1,
                                                         Float_t LeftSBndCut = 1.792, Float_t RightSBndCut = 1.942
               )
{
  //
  // Macro for the AliAnalysisTaskSE for D candidates vs Multiplicity
  // Invariant mass histogram in pt and multiplicity bins in a 3D histogram
  //   different estimators implemented
  //==============================================================================


  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AddTaskDvsMultiplicity_BDT", "No analysis manager to connect to.");
  }
    
  Bool_t stdcuts=kFALSE;
  TFile* filecuts;
  if( filename.EqualTo("") ) {
    stdcuts=kTRUE;
  } else {
    filecuts=TFile::Open(filename.Data());
    if(!filecuts ||(filecuts&& !filecuts->IsOpen())){
      Printf("FATAL: Input file not found : check your cut object");
    }
  }
 

  //Analysis Task
  AliRDHFCuts *analysiscuts=0x0;
  Int_t Nptbins;
  Float_t *ptbin;   
 
  TString Name="";
  if(pdgMeson==411){
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDplustoKpipi();
      if (system == 0) analysiscuts->SetStandardCutsPP2010();
      else analysiscuts->SetStandardCutsPbPb2011();
    }
    else analysiscuts = (AliRDHFCutsDplustoKpipi*)filecuts->Get(finAnObjname);
    Name="Dplus";
  }else if(pdgMeson==421){
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsD0toKpi();
      if (system == 0) analysiscuts->SetStandardCutsPP2010();
      else analysiscuts->SetStandardCutsPbPb2011();
    }
    else {analysiscuts = (AliRDHFCutsD0toKpi*)filecuts->Get(finAnObjname);
       Nptbins = analysiscuts->GetNPtBins();
       ptbin = analysiscuts->GetPtBinLimits();
    }    
    Name="D0";
  }else if(pdgMeson==413){
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDStartoKpipi();
      if (system == 0) analysiscuts->SetStandardCutsPP2010();
      else analysiscuts->SetStandardCutsPbPb2011();
    }
    else analysiscuts = (AliRDHFCutsDStartoKpipi*)filecuts->Get(finAnObjname);
    Name="DStar";
  }else if(pdgMeson==431){
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsDstoKKpi();
      if (system == 0) analysiscuts->SetStandardCutsPP2010();
      else analysiscuts->SetStandardCutsPbPb2011();
    }
    else analysiscuts = (AliRDHFCutsDstoKKpi*)filecuts->Get(finAnObjname);
    Name="Ds";
  }else if(pdgMeson==4122){
    if(stdcuts) {
      analysiscuts = new AliRDHFCutsLctoV0();
      if (system == 0) analysiscuts->SetStandardCutsPP2010();
      else analysiscuts->SetStandardCutsPbPb2011();
    }
    else analysiscuts = (AliRDHFCutsLctoV0*)filecuts->Get(finAnObjname);
    Name="Lc2pK0S";
  }
    
  AliAnalysisTaskSEDvsMultiplicity_BDT *dMultTask = new AliAnalysisTaskSEDvsMultiplicity_BDT("dMultAnalysis",pdgMeson,analysiscuts,isPPbData,readMC,applyBDT);
  dMultTask->SetDebugLevel(0);
  dMultTask->SetUseBit(kTRUE);
  dMultTask->SetDoImpactParameterHistos(kFALSE); //kFALSE
  dMultTask->SetSubtractTrackletsFromDaughters(subtractDau);
  dMultTask->SetMultiplicityEstimator(recoEstimator);
  dMultTask->SetMCPrimariesEstimator(MCEstimator);
  dMultTask->SetMCOption(MCOption);
  dMultTask->SetLcToV0decay(isLcV0);
  if(isPPbData) dMultTask->SetIsPPbData();
  dMultTask->SetFillTree(kTRUE);
  
  if(applyBDT){ 
    TFile *fileBDT = TFile::Open(BDTfilename);
    if(!fileBDT ||(fileBDT&& !fileBDT->IsOpen())) ::Fatal("AddTaskD0BDT_change", "BDT file not found : check your BDT object"); 
    TList *bdtlist = new TList();
    for(Int_t i=0;i<Nptbins;i++){
          TString BDTobjname = BDTobjnamepre;
          BDTobjname += Form("1_%.0f_%.0f",ptbin[i],ptbin[i+1]);
          AliRDHFBDT *thisbdt = (AliRDHFBDT*)(fileBDT->Get(BDTobjname)->Clone(Form("_%s",BDTobjname.Data())));
          if(!thisbdt) ::Fatal("AddTaskD0BDT_change", Form("Failed to find BDT named %s",BDTobjname.Data()));
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
                  if(!thisbdt2_0) ::Fatal("AddTaskD0BDT_change", Form("Failed to find BDT named %s",BDT2objname1.Data()));
                  if(!thisbdt2_1) ::Fatal("AddTaskD0BDT_change", Form("Failed to find BDT named %s",BDT2objname2.Data()));
                  if(!thisbdt2_2) ::Fatal("AddTaskD0BDT_change", Form("Failed to find BDT named %s",BDT2objname3.Data()));
                  if(!thisbdt2_3) ::Fatal("AddTaskD0BDT_change", Form("Failed to find BDT named %s",BDT2objname4.Data()));
                  if(!thisbdt2_4) ::Fatal("AddTaskD0BDT_change", Form("Failed to find BDT named %s",BDT2objname5.Data()));
                  if(!thisbdt2_5) ::Fatal("AddTaskD0BDT_change", Form("Failed to find BDT named %s",BDT2objname6.Data()));
                  bdtlist->Add(thisbdt2_0);
                  bdtlist->Add(thisbdt2_1);
                  bdtlist->Add(thisbdt2_2);
                  bdtlist->Add(thisbdt2_3);
                  bdtlist->Add(thisbdt2_4);
                  bdtlist->Add(thisbdt2_5);
          }
  }
  fileBDT->Close();
  dMultTask->SetBDTList(bdtlist);
  dMultTask->SetBDTGetRespTree(GetRespTree); 
  dMultTask->SetBDTRespCut(BDTRespCut);
  dMultTask->SetBDTSidebandCut(LeftSBndCut,RightSBndCut);
  dMultTask->SetBDTSampleSideband(DoSidebndSample);
  dMultTask->SetBDTSidebandSamplingFraction(SBndSampleFrac);
  }
  if(pdgMeson==421) dMultTask->SetBDTFullVarString("ptD:topo1:topo2:lxy:nlxy:iscut:ispid:type:mass:d0d0:cosp:dca:ptk:ptpi:cospxy:d0k:d0pi:cosstar:ptB:pdgcode:YD0:phi:nsigmaTPCPi_0:nsigmaTPCK_0:nsigmaTOFPi_0:nsigmaTOFK_0:nsigmaTPCPi_1:nsigmaTPCK_1:nsigmaTOFPi_1:nsigmaTOFK_1:nsigmaCombPi_0:nsigmaCombK_0:nsigmaCombPi_1:nsigmaCombK_1:coutMulti");
  if(pdgMeson==431) dMultTask->SetBDTFullVarString("invMass:ptCand:imp_par_XY:d_len:multForCand:d_len_XY:norm_dl_XY:Cos_P:Cos_P_XY:dca:max_norm_d0d0exp:Sig_Vert:delta_mass_KK:cos_PiDs:cos_PiKPhi_3:imp_par_prong0:imp_par_prong1:imp_par_prong2:PtProng0:PtProng1:PtProng2:nsigTPC_Pi_0:nsigTPC_K_0:nsigTOF_Pi_0:nsigTOF_K_0:nsigTPC_Pi_1:nsigTPC_K_1:nsigTOF_Pi_1:nsigTOF_K_1:nsigTPC_Pi_2:nsigTPC_K_2:nsigTOF_Pi_2:nsigTOF_K_2:type");
   
 
  if(NchWeight){
    TH1F *hNchPrimaries = NULL;
    TH1F *hMeasNchPrimaries = NULL;
    if(NchWeight==1){
      if(isPPbData) {
	hNchPrimaries = (TH1F*)filecuts->Get("hNtrUnCorrEvWithDWeight"); // MC distribution
      }
      else hNchPrimaries = (TH1F*)filecuts->Get("hGenPrimaryParticlesInelGt0");
      if(hNchPrimaries) {
	dMultTask->UseMCNchWeight(NchWeight);
	dMultTask->SetHistoNchWeight(hNchPrimaries);
      } else {
	Printf("FATAL: Histogram for Nch multiplicity weights not found");
	return 0x0;
      }
      hMeasNchPrimaries = (TH1F*)filecuts->Get("hMeasNtrUnCorrEvWithD"); // data distribution
      if(hMeasNchPrimaries) {
	dMultTask->SetMeasuredNchHisto(hMeasNchPrimaries);
      }
    }
    else if(NchWeight==2){
      hNchPrimaries = (TH1F*)filecuts->Get("hNtrUnCorrEvWithDWeight"); // MC distribution
      hMeasNchPrimaries = (TH1F*)filecuts->Get("hMeasNtrUnCorrEvWithD"); // data distribution
      if(hNchPrimaries && hMeasNchPrimaries) {
	dMultTask->UseMCNchWeight(NchWeight);
	dMultTask->SetHistoNchWeight(hNchPrimaries);
	dMultTask->SetMeasuredNchHisto(hMeasNchPrimaries);
      } else {
	Printf("FATAL: Histogram for Ntrk multiplicity weights not found");
	return 0x0;
      }
    }
  }
    
    
  if(pdgMeson==421) {
    //dMultTask->SetMassLimits(1.5648,2.1648);
    dMultTask->SetMassLimits(0.,2.5);
    dMultTask->SetNMassBins(200);
  } else if(pdgMeson==4122) {
    dMultTask->SetMassLimits(pdgMeson,0.250);
    dMultTask->SetNMassBins(1000);
  }else if(pdgMeson==411)dMultTask->SetMassLimits(pdgMeson,0.2);
    
  if(estimatorFilename.EqualTo("") ) {
    printf("Estimator file not provided, multiplcity corrected histograms will not be filled\n");
  } else{
        
    TFile* fileEstimator=TFile::Open(estimatorFilename.Data());
    if(!fileEstimator)  {
      Printf("FATAL: File with multiplicity estimator not found\n");
      return NULL;
    }
        
    dMultTask->SetReferenceMultiplcity(refMult);
        
    const Char_t* profilebasename="SPDmult10";
    if(recoEstimator==AliAnalysisTaskSEDvsMultiplicity_BDT::kVZEROA || recoEstimator==AliAnalysisTaskSEDvsMultiplicity_BDT::kVZEROAEq) profilebasename="VZEROAmult";
    else if(recoEstimator==AliAnalysisTaskSEDvsMultiplicity_BDT::kVZERO || recoEstimator==AliAnalysisTaskSEDvsMultiplicity_BDT::kVZEROEq) profilebasename="VZEROMmult";
    cout<<endl<<endl<<" profilebasename="<<profilebasename<<endl<<endl;
        
    if (isPPbData) {
      if(year == 16) {
	    const Char_t* periodNames[4] = {"LHC16q_265499to265525_265309to265387", "LHC16q_265435","LHC16q_265388to265427","LHC16t_267163to267166"};
	    TProfile* multEstimatorAvg[4];
	    for(Int_t ip=0; ip<4; ip++) {
	    cout<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
	    multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
	    if (!multEstimatorAvg[ip]) {
	      Printf("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]);
	      return NULL;
	    }
	  }
	  dMultTask->SetMultiplVsZProfileLHC16qt1stBunch(multEstimatorAvg[0]);
	  dMultTask->SetMultiplVsZProfileLHC16qt2ndBunch(multEstimatorAvg[1]);
	  dMultTask->SetMultiplVsZProfileLHC16qt3rdBunch(multEstimatorAvg[2]);
	  dMultTask->SetMultiplVsZProfileLHC16qt4thBunch(multEstimatorAvg[3]);
      }else {
  	//Only use two profiles if pPb 2013
  	const Char_t* periodNames[2] = {"LHC13b", "LHC13c"};
	  TProfile* multEstimatorAvg[2];
	  for(Int_t ip=0; ip<2; ip++) {
	    cout<< " Trying to get "<<Form("%s_%s",profilebasename,periodNames[ip])<<endl;
	    multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
	    if (!multEstimatorAvg[ip]) {
	      Printf("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]);
	      return NULL;
	    }
	  }
	dMultTask->SetMultiplVsZProfileLHC13b(multEstimatorAvg[0]);
	dMultTask->SetMultiplVsZProfileLHC13c(multEstimatorAvg[1]);
      }
    }else {//pp data
      if(year == 10){
      const Char_t* periodNames[4] = {"LHC10b", "LHC10c", "LHC10d", "LHC10e"};
      TProfile* multEstimatorAvg[4];
      for(Int_t ip=0; ip<4; ip++) {
	multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
	if (!multEstimatorAvg[ip]) {
	  Printf("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]);
	  return NULL;
	}
      }
      dMultTask->SetMultiplVsZProfileLHC10b(multEstimatorAvg[0]);
      dMultTask->SetMultiplVsZProfileLHC10c(multEstimatorAvg[1]);
      dMultTask->SetMultiplVsZProfileLHC10d(multEstimatorAvg[2]);
      dMultTask->SetMultiplVsZProfileLHC10e(multEstimatorAvg[3]);
    }else if(year ==16){
      const Char_t* periodNames[10]={"LHC16d","LHC16e","LHC16g","LHC16h_1", "LHC16h_2","LHC16j","LHC16k","LHC16l","LHC16o","LHC16p"};
      TProfile *multEstimatorAvg[10];
      for(Int_t ip=0;ip<10; ip++){
        multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
  if (!multEstimatorAvg[ip]) {
    Printf("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]);
    return NULL;
  }
      }
      dMultTask->SetMultiplVsZProfileLHC16d(multEstimatorAvg[0]);
      dMultTask->SetMultiplVsZProfileLHC16e(multEstimatorAvg[1]);
      dMultTask->SetMultiplVsZProfileLHC16g(multEstimatorAvg[2]);
      dMultTask->SetMultiplVsZProfileLHC16h1(multEstimatorAvg[3]);
      dMultTask->SetMultiplVsZProfileLHC16h2(multEstimatorAvg[4]);
      dMultTask->SetMultiplVsZProfileLHC16j(multEstimatorAvg[5]);
      dMultTask->SetMultiplVsZProfileLHC16k(multEstimatorAvg[6]);
      dMultTask->SetMultiplVsZProfileLHC16l(multEstimatorAvg[7]);
      dMultTask->SetMultiplVsZProfileLHC16o(multEstimatorAvg[8]);
      dMultTask->SetMultiplVsZProfileLHC16p(multEstimatorAvg[9]);

    }else if(year == 17){
     const Char_t* periodNames[10]={"LHC17e","LHC17f","LHC17h","LHC17i", "LHC17j","LHC17k","LHC17l","LHC17m","LHC17o","LHC17r"};
      TProfile *multEstimatorAvg[10];
      for(Int_t ip=0;ip<10; ip++){
        multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
        if (!multEstimatorAvg[ip]) {
        Printf("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]);
        return NULL;
        }
      }
      dMultTask->SetMultiplVsZProfileLHC17e(multEstimatorAvg[0]);
      dMultTask->SetMultiplVsZProfileLHC17f(multEstimatorAvg[1]);
      dMultTask->SetMultiplVsZProfileLHC17h(multEstimatorAvg[2]);
      dMultTask->SetMultiplVsZProfileLHC17i(multEstimatorAvg[3]);
      dMultTask->SetMultiplVsZProfileLHC17j(multEstimatorAvg[4]);
      dMultTask->SetMultiplVsZProfileLHC17k(multEstimatorAvg[5]);
      dMultTask->SetMultiplVsZProfileLHC17l(multEstimatorAvg[6]);
      dMultTask->SetMultiplVsZProfileLHC17m(multEstimatorAvg[7]);
      dMultTask->SetMultiplVsZProfileLHC17o(multEstimatorAvg[8]);
      dMultTask->SetMultiplVsZProfileLHC17r(multEstimatorAvg[9]);
    }else if(year == 18){
const Char_t* periodNames[14]={"LHC18b","LHC18d","LHC18e","LHC18f", "LHC18g","LHC18h","LHC18i","LHC18j","LHC18k","LHC18l","LHC18m","LHC18n","LHC18o","LHC18p"};
      TProfile *multEstimatorAvg[14];
      for(Int_t ip=0;ip<14; ip++){
        multEstimatorAvg[ip] = (TProfile*)(fileEstimator->Get(Form("%s_%s",profilebasename,periodNames[ip]))->Clone(Form("%s_%s_clone",profilebasename,periodNames[ip])));
        if (!multEstimatorAvg[ip]) {
        Printf("Multiplicity estimator for %s not found! Please check your estimator file",periodNames[ip]);
        return NULL;
        }
      }
      dMultTask->SetMultiplVsZProfileLHC18b(multEstimatorAvg[0]);
      dMultTask->SetMultiplVsZProfileLHC18d(multEstimatorAvg[1]);
      dMultTask->SetMultiplVsZProfileLHC18e(multEstimatorAvg[2]);
      dMultTask->SetMultiplVsZProfileLHC18f(multEstimatorAvg[3]);
      dMultTask->SetMultiplVsZProfileLHC18g(multEstimatorAvg[4]);
      dMultTask->SetMultiplVsZProfileLHC18h(multEstimatorAvg[5]);
      dMultTask->SetMultiplVsZProfileLHC18i(multEstimatorAvg[6]);
      dMultTask->SetMultiplVsZProfileLHC18j(multEstimatorAvg[7]);
      dMultTask->SetMultiplVsZProfileLHC18k(multEstimatorAvg[8]);
      dMultTask->SetMultiplVsZProfileLHC18l(multEstimatorAvg[9]);
      dMultTask->SetMultiplVsZProfileLHC18m(multEstimatorAvg[10]);
      dMultTask->SetMultiplVsZProfileLHC18n(multEstimatorAvg[11]);
      dMultTask->SetMultiplVsZProfileLHC18o(multEstimatorAvg[12]);
      dMultTask->SetMultiplVsZProfileLHC18p(multEstimatorAvg[13]);
    }//18
      }
   mgr->AddTask(dMultTask);
   }
  // Create containers for input/output
    
  TString inname = "cinput";
  TString outname = "coutput";
  TString cutsname = "coutputCuts";
  TString normname = "coutputNorm";
  TString profname = "coutputProf";
  TString treenameList = "coutputD0Tree";  
  TString D0BDTResponsesList = "coutputD0BDTResponse"; 
  inname += Name.Data();
  outname += Name.Data();
  cutsname += Name.Data();
  normname += Name.Data();
  profname += Name.Data();
  treenameList += Name.Data();
  D0BDTResponsesList += Name.Data();
  inname += finDirname.Data();
  outname += finDirname.Data();
  cutsname += finDirname.Data();
  normname += finDirname.Data();
  profname += finDirname.Data();
  treenameList += finDirname.Data();
  D0BDTResponsesList += finDirname.Data();
  
  AliAnalysisDataContainer *cinput = mgr->CreateContainer(inname,TChain::Class(),AliAnalysisManager::kInputContainer);
    
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG3_D2H_DMult_";
  outputfile += Name.Data(); 
  outputfile += finDirname.Data(); 
  AliAnalysisDataContainer *coutputCuts = mgr->CreateContainer(cutsname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutput = mgr->CreateContainer(outname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutputNorm = mgr->CreateContainer(normname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutputProf = mgr->CreateContainer(profname,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutputD0Tree = mgr->CreateContainer(treenameList,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutputD0BDTResponse = mgr->CreateContainer(D0BDTResponsesList,TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data()); 
  mgr->ConnectInput(dMultTask,0,mgr->GetCommonInputContainer());

  mgr->ConnectOutput(dMultTask,1,coutput);

  mgr->ConnectOutput(dMultTask,2,coutputCuts);

  mgr->ConnectOutput(dMultTask,3,coutputNorm);

  mgr->ConnectOutput(dMultTask,4,coutputProf);
  mgr->ConnectOutput(dMultTask,5,coutputD0Tree);
  if(!readMC && applyBDT)
  {
     mgr->ConnectOutput(dMultTask,6,coutputD0BDTResponse);
  }
  return dMultTask;
}
