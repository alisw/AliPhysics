AliAnalysisTaskBEpp13TeV* AddTaskBEpp13TeV(
	TString	taskName = "beauty",
	bool	isMC = false,

    //TString estimatorFileName = "fileEstimatorAvg.root",
    //double  multRef = 11.7,

	int		minTPCnCrossedRow = 70,
	int		minTPCnClsPID = 80,
	double	maxTPCchi2 = 4,
	double	minTPCclsRatio = 0.8,
	int		minITSnCls = 3,
	int		itsLayer = 2,
	double	tpcPIDlow = -1,
	double	tpcPIDhigh = 3,
	double	tofPID = 3)
{
  //TString taskName = "beauty";
  //TString itsLayer = citsLayer;
  printf("=======================================\n");
  printf("==========ADD TASK PARAMETERS==========\n");
  printf("=======================================\n");
  printf("Task name: %s \n", taskName.Data());
  printf("MC Flag: %d \n",isMC);
  printf("Min. TPC crossed raw: %d \n", minTPCnCrossedRow);
  printf("Min. TPC clusters for PID: %d \n", minTPCnClsPID);
  printf("Max. TPC chi2/ndf: %.1f \n", maxTPCchi2);
  printf("Min. TPC cls/clsPID: %.1f \n", minTPCclsRatio);
  printf("Min. ITS clusters: %d \n", minITSnCls);
  printf("ITS requirement: %d \n",itsLayer);
  printf("Min. TPC nsigma: %.1f \n", tpcPIDlow);
  printf("Max. TPC nsigma: %.1f \n", tpcPIDhigh);
  printf("TOF nsigma: %.1f \n", tofPID);
  printf("===================================\n\n");

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if(!mgr){
    return 0x0;
  }
  
  if(!mgr->GetInputEventHandler()){
    return 0x0;
  }
  
  // now we create an instance of your task
  AliAnalysisTaskBEpp13TeV* task = new AliAnalysisTaskBEpp13TeV(taskName.Data());   
  if(!task) return 0x0;
  if(!isMC) task->SelectCollisionCandidates(AliVEvent::kINT7);
  task->SetMinTPCnCrossedRow(minTPCnCrossedRow);
  task->SetMinTPCNclsPID(minTPCnClsPID);
  task->SetMaxTPCchi2(maxTPCchi2);
  task->SetMinTPCclsRatio(minTPCclsRatio);
  task->SetMinITSNcls(minITSnCls);
  task->SetITSlayer(itsLayer);
  task->SetPIDCuts(tpcPIDlow, tpcPIDhigh, tofPID);

  // Multiplicity
  double multRef = 11.7;
  TString estimatorFileName = "fileEstimatorAvg.root";
  if(estimatorFileName.EqualTo("")){
    printf("Estimator file not provided, multiplcity corrected histograms will not be filled\n");
  }else{
    TString alienpath = "alien:///alice/cern.ch/user/j/jonghan/be_pp13TeV/";
    TString path_file = alienpath + estimatorFileName;
    //TFile* fileEstimator = TFile::Open(estimatorFileName.Data());
    TFile* fileEstimator = TFile::Open(path_file.Data());
    if(!fileEstimator){
      printf("File with multiplicity estimator not found\n");
      return 0x0;
    }
    task->SetMultReference(multRef);

    TProfile *multEstimatorAvg = (TProfile*)fileEstimator->Get("fNtrAvg_pp_pass1_all");
    if(!multEstimatorAvg) {
      printf("Multiplicity estimator not found! Please check your estimator file\n");
      return 0x0;
    }

    task->SetEstimatorAvg(multEstimatorAvg);
    //task->SetEstimatorHistogram(period);
  }
  
  
  if(isMC){
	task->SetMCanalysis();

	// B meson pt correction
	TF1 *BmesonCentLow = new TF1("BmesonCentLow","(0.604848+0.323547*x-0.0530093*x*x)/1.48462",0.1,3.5);
	TF1 *BmesonCentHigh = new TF1("BmesonCentHigh","(1.39587-0.119037*x+0.0102669*x*x-(3.65110e-04)*x*x*x+(6.44063e-06)*x*x*x*x-(5.51122e-08)*x*x*x*x*x+(1.81335e-10)*x*x*x*x*x*x)/1.48462",3.5,100);
    TF1 *BmesonMinLow = new TF1("BmesonMinLow","(0.568465+0.305052*x-0.0472978*x*x)/1.5673",0.1,3.5);
    TF1 *BmesonMinHigh = new TF1("BmesonMinHigh","(1.29080-0.0910133*x+0.00815947*x*x-(2.81802e-04)*x*x*x+(4.76347e-06)*x*x*x*x-(3.89458e-08)*x*x*x*x*x+(1.22428e-10)*x*x*x*x*x*x)/1.5673",3.5,100);
    TF1 *BmesonMaxLow = new TF1("BmesonMaxLow","(0.634004+0.338291*x-0.0575622*x*x)/1.42542",0.1,3.5);
    TF1 *BmesonMaxHigh = new TF1("BmesonMaxHigh","(1.47864-0.140949*x+0.0118876*x*x-(4.28193e-04)*x*x*x+(7.69184e-06)*x*x*x*x-(6.69947e-08)*x*x*x*x*x+(2.23972e-10)*x*x*x*x*x*x)/1.42542",3.5,100);
    task->SetBcorrCentLow(BmesonCentLow);
    task->SetBcorrCentHigh(BmesonCentHigh);
    task->SetBcorrMinLow(BmesonMinLow);
    task->SetBcorrMinHigh(BmesonMinHigh);
    task->SetBcorrMaxLow(BmesonMaxLow);
    task->SetBcorrMaxHigh(BmesonMaxHigh);
    std::cout<<"---------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"------------------------------B meson pt correction------------------------------------"<<std::endl;
    std::cout<<"---------------------------------------------------------------------------------------"<<std::endl<<std::endl;

    // D meson pt correction
    // from old D0 spectrum (measured)
		//TF1 *DmesonCorr = new TF1("DmesonCorr","(1/1.73303)*(0.552699*TMath::Power(1+(1.18107-1)*x/1.4446,1.18107/(1-1.18107))+0.158337*TMath::Exp(1.8436-x/1.28754)) / ((3.33487*((3.54275-1.)*(3.54275-2.))/(3.54275*0.501107*(3.54275*0.501107+9.24732*(3.54275-2.)))*pow(1.+(sqrt(9.24732*9.24732 + x*x)-9.24732)/(3.54275*0.501107),-3.54275)))",0.1,100);
    //TF1 *DmesonCorrVar1 = new TF1("DmesonCorrVar1","(1/1.27055)*(0.552699*TMath::Power(1+(1.18107-1)*x/1.5446,1.18107/(1-1.18107))+0.088337*TMath::Exp(1.8436-x/1.28754)) / ((3.33487*((3.54275-1.)*(3.54275-2.))/(3.54275*0.501107*(3.54275*0.501107+9.24732*(3.54275-2.)))*pow(1.+(sqrt(9.24732*9.24732 + x*x)-9.24732)/(3.54275*0.501107),-3.54275)))",0.1,100);
    //TF1 *DmesonCorrVar2 = new TF1("DmesonCorrVar2","(1/2.193)*(0.552699*TMath::Power(1+(1.18107-1)*x/1.3446,1.18107/(1-1.18107))+0.228337*TMath::Exp(1.8436-x/1.28754)) / ((3.33487*((3.54275-1.)*(3.54275-2.))/(3.54275*0.501107*(3.54275*0.501107+9.24732*(3.54275-2.)))*pow(1.+(sqrt(9.24732*9.24732 + x*x)-9.24732)/(3.54275*0.501107),-3.54275)))",0.1,100);
    // from newest D0 spectrum 
		TF1 *DmesonCorr = new TF1("DmesonCorr","(1/1.85581)*(0.0218308*TMath::Power(1+(0.974074-1)*x/4.33734,0.974074/(1-0.974074))+0.395143*TMath::Exp(1.35595-x/1.39509)) / (3.33623*((3.57495-1.)*(3.57495-2.))/(3.57495*0.508028*(3.57495*0.508028+9.15572*(3.57495-2.)))*pow(1.+(sqrt(9.15572*9.15572 + x*x)-9.15572)/(3.57495*0.508028),-3.57495))",0.1,100);
    TF1 *DmesonCorrVar1 = new TF1("DmesonCorrVar1","(1/1.41037)*(0.0288308*TMath::Power(1+(0.974074-1)*x/4.33734,0.974074/(1-0.974074))+0.295143*TMath::Exp(1.35595-x/1.39509)) / (3.33623*((3.57495-1.)*(3.57495-2.))/(3.57495*0.508028*(3.57495*0.508028+9.15572*(3.57495-2.)))*pow(1.+(sqrt(9.15572*9.15572 + x*x)-9.15572)/(3.57495*0.508028),-3.57495))",0.1,100);
    TF1 *DmesonCorrVar2 = new TF1("DmesonCorrVar2","(1/2.53073)*(0.0148308*TMath::Power(1+(0.974074-1)*x/4.33734,0.974074/(1-0.974074))+0.545143*TMath::Exp(1.35595-x/1.39509)) / (3.33623*((3.57495-1.)*(3.57495-2.))/(3.57495*0.508028*(3.57495*0.508028+9.15572*(3.57495-2.)))*pow(1.+(sqrt(9.15572*9.15572 + x*x)-9.15572)/(3.57495*0.508028),-3.57495))",0.1,100);
    task->SetDcorrFtn(DmesonCorr);
    task->SetDcorrFtnVar1(DmesonCorrVar1);
    task->SetDcorrFtnVar2(DmesonCorrVar2);
    std::cout<<"---------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"------------------------------D meson pt correction------------------------------------"<<std::endl;
    std::cout<<"---------------------------------------------------------------------------------------"<<std::endl<<std::endl;
    
    // Lc meson pt correction
    //TF1 *LcCorr = new TF1("LcCorr","12.6935*TMath::Exp(-0.568250*x)+0.00681196",1., 30.);
    TF1 *LcCorr = new TF1("LcCorr","1", .5, 30.);
    task->SetLccorrFtn(LcCorr);
    std::cout<<"---------------------------------------------------------------------------------------"<<std::endl;
    std::cout<<"------------------------------LambdaC pt correction for most-central analysis----------"<<std::endl;
    std::cout<<"---------------------------------------------------------------------------------------"<<std::endl<<std::endl;
  }
  
  
  // add your task to the manager
  mgr->AddTask(task);
  
  // by default, a file is open for writing. here, we get the filename
  TString fileName = AliAnalysisManager::GetCommonFileName();
  fileName += ":ElectronID_";      // create a subfolder in the file
  fileName += taskName;

  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,mgr->CreateContainer(Form("ccontainer0_%s",taskName.Data()), TList::Class(), AliAnalysisManager::kOutputContainer, fileName.Data()));
  
	return task;
}
