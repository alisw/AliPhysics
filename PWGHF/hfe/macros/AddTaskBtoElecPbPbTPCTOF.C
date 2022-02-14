AliAnalysisTaskBtoElecPbPbTPCTOF* AddTaskBtoElecPbPbTPCTOF( ///-> to run locally
			TString uniqueID        = "",
			bool 	isMC 			= kFALSE,
			double centMin = 30.,
			double centMax = 50.,
			int minTPCNclusters =	100,
			int minTPCclustersPID =	80,
			double maxTPCchi2 = 4.,
			double minTPCclusterRatio =	0.6,
			int  minITSNcluster = 4,
			double maxITSclusterFrac = 0.3,
			double maxITSchi2 = 10.,
			TString itsLayer = "kBoth",
			double eta = 0.8,
			double ptMin = 0.3,
			double ptMax = 30.,
			float dcaxy = 0.1,
			float dcaz = 0.2,
			double tpcPIDlow = 0.,
			double tpcPIDhigh = 3.,
			double tofPID = 3,
			int Bcorr = 2,
			int Dcorr = 2,
			int Lccorr = 2
)           
{
    
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	
	if (!mgr) {
		::Error("AddTaskBtoElecPbPbTPCTOF", "No analysis manager to connect to.");
		return NULL;
	}
	
	if (!mgr->GetInputEventHandler()) {
		::Error("AddTaskBtoElecPbPbTPCTOF", "This task requires an input event handler");
		return NULL;
	}

	///Task config
	AliAnalysisTaskBtoElecPbPbTPCTOF *task = new AliAnalysisTaskBtoElecPbPbTPCTOF(uniqueID.Data());
	printf("task ------------------------ %p\n ", task);
	task->SetCentrality(centMin, centMax);
	task->SetMinTPCNcls(minTPCNclusters);
	task->SetMinTPCNclsPID(minTPCclustersPID);
	task->SetMaxTPCchi2(maxTPCchi2);
	task->SetMinTPCclsRatio(minTPCclusterRatio);
	task->SetMinITSNcls(minITSNcluster);
	task->SetMaxITSclsFrac(maxITSclusterFrac);
	task->SetMaxITSchi2(maxITSchi2);
	task->SetITSlayer(itsLayer);
	task->SetEtaRange(eta);
	task->SetPtRange(ptMin, ptMax);
	task->SetDCARange(dcaxy, dcaz);

	task->SetPIDCuts(tpcPIDlow, tpcPIDhigh, tofPID);
  
	//Trigger
	if(!isMC){
		task->SelectCollisionCandidates(AliVEvent::kINT7); //Selecting Minumum Bias events (selected randomlly)
	}

	//In the case of a simulation
	if(isMC)
	{
		task->SetMCanalysis();

		// B meson pt correction
		if(Bcorr==1)
		{
			TF1 *BmesonCorr = new TF1("BmesonCorr","(6.44844/(TMath::Power(TMath::Exp(0.117695*x-0.0250428*x*x+1.56425/x)+x/3.89413,2.05618))+0.910764*TMath::Exp(-0.137564*x)) * (0.5/(1.+exp((x-7.)*0.7))+0.5+(x-15.)/300)", .7, 50.);
			TF1 *BmesonCorr1 = new TF1("BmesonCorr1","(1/1.35104)*(6.44844/(TMath::Power(TMath::Exp(0.117695*x-0.0250428*x*x+1.56425/x)+x/3.89413,2.05618))+0.910764*TMath::Exp(-0.137564*x)) * (0.5/(1.+exp((x-7.)*0.7))+0.5+(x-15.)/300)", .7, 50.);
			task->SetBcorrFtn(BmesonCorr);
			task->SetBcorrFtn1(BmesonCorr1);
			std::cout<<"---------------------------------------------------------------------------------------"<<std::endl;
			std::cout<<"------------------------------B meson pt correction for most-central analysis----------"<<std::endl;
			std::cout<<"---------------------------------------------------------------------------------------"<<std::endl<<std::endl;
		}
		if(Bcorr==2)
		{
			TF1 *BmesonCorr = new TF1("BmesonCorr","(5.55733/(TMath::Power(TMath::Exp(-0.0335019*x-0.0112941*x*x+0.832544/x)+x/6.75566,2.30081))+1.79595*TMath::Exp(-0.645011*x))*(.55/(.9+exp((x-6.5)*.65))+.55+(x-15.)/300)",.1,50.);
			TF1 *BmesonCorr1 = new TF1("BmesonCorr1","(1/2.51757)*(5.55733/(TMath::Power(TMath::Exp(-0.0335019*x-0.0112941*x*x+0.832544/x)+x/6.75566,2.30081))+1.79595*TMath::Exp(-0.645011*x))*(.55/(.9+exp((x-6.5)*.65))+.55+(x-15.)/300)",.1,50.);
			task->SetBcorrFtn(BmesonCorr);
			task->SetBcorrFtn1(BmesonCorr1);
			std::cout<<"---------------------------------------------------------------------------------------"<<std::endl;
			std::cout<<"------------------------------B meson pt correction for semi-central analysis----------"<<std::endl;
			std::cout<<"---------------------------------------------------------------------------------------"<<std::endl<<std::endl;
		}

		// D meson pt correction
		if(Dcorr==1)
		{
			TF1 *DmesonCorr = new TF1("DmesonCorr","(4.01931*TMath::Power(1+(1.25363-1)*x/0.769208,1.25363/(1-1.25363))+1.61721*TMath::Exp(2.86484-x/0.809167)) / (5655.77*TMath::Gaus(x,-7.79843,-0.217064)+26.0041*TMath::Landau(x,1.77719,0.608266,0))",1.,50.);
			TF1 *DmesonCorr1 = new TF1("DmesonCorr1","(1/2.82768)*((4.01931*TMath::Power(1+(1.25363-1)*x/0.769208,1.25363/(1-1.25363))+1.61721*TMath::Exp(2.86484-x/0.809167)) / (5655.77*TMath::Gaus(x,-7.79843,-0.217064)+26.0041*TMath::Landau(x,1.77719,0.608266,0)))",1.,50.);
			TF1 *DmesonCorr2 = new TF1("DmesonCorr2","(1/2.21365)*((4.01931*TMath::Power(1+(1.25363-1)*x/0.769208,1.25363/(1-1.25363))+1.61721*TMath::Exp(2.86484-x/0.809167)) / (5655.77*TMath::Gaus(x,-7.79843,-0.217064)+26.0041*TMath::Landau(x,1.77719,0.608266,0)))",1.,50.);
			TF1 *DmesonCorr3 = new TF1("DmesonCorr3","(1/1.47851)*((4.01931*TMath::Power(1+(1.25363-1)*x/0.769208,1.25363/(1-1.25363))+1.61721*TMath::Exp(2.86484-x/0.809167)) / (5655.77*TMath::Gaus(x,-7.79843,-0.217064)+26.0041*TMath::Landau(x,1.77719,0.608266,0)))",1.,50.);
			TF1 *DmesonCorr4 = new TF1("DmesonCorr4","(1/1.07738)*((4.01931*TMath::Power(1+(1.25363-1)*x/0.769208,1.25363/(1-1.25363))+1.61721*TMath::Exp(2.86484-x/0.809167)) / (5655.77*TMath::Gaus(x,-7.79843,-0.217064)+26.0041*TMath::Landau(x,1.77719,0.608266,0)))",1.,50.);
			TF1 *DmesonCorr5 = new TF1("DmesonCorr5","(1/0.617101)*((4.01931*TMath::Power(1+(1.25363-1)*x/0.769208,1.25363/(1-1.25363))+1.61721*TMath::Exp(2.86484-x/0.809167)) / (5655.77*TMath::Gaus(x,-7.79843,-0.217064)+26.0041*TMath::Landau(x,1.77719,0.608266,0)))",1.,50.);
			TF1 *DmesonCorr6 = new TF1("DmesonCorr6","(1/0.419671)*((4.01931*TMath::Power(1+(1.25363-1)*x/0.769208,1.25363/(1-1.25363))+1.61721*TMath::Exp(2.86484-x/0.809167)) / (5655.77*TMath::Gaus(x,-7.79843,-0.217064)+26.0041*TMath::Landau(x,1.77719,0.608266,0)))",1.,50.);
			TF1 *DmesonCorr7 = new TF1("DmesonCorr7","(1/0.305515)*((4.01931*TMath::Power(1+(1.25363-1)*x/0.769208,1.25363/(1-1.25363))+1.61721*TMath::Exp(2.86484-x/0.809167)) / (5655.77*TMath::Gaus(x,-7.79843,-0.217064)+26.0041*TMath::Landau(x,1.77719,0.608266,0)))",1.,50.);
			TF1 *DmesonCorr8 = new TF1("DmesonCorr8","(1/0.171142)*((4.01931*TMath::Power(1+(1.25363-1)*x/0.769208,1.25363/(1-1.25363))+1.61721*TMath::Exp(2.86484-x/0.809167)) / (5655.77*TMath::Gaus(x,-7.79843,-0.217064)+26.0041*TMath::Landau(x,1.77719,0.608266,0)))",1.,50.);
			TF1 *DmesonCorr9 = new TF1("DmesonCorr9","(1/0.0978727)*((4.01931*TMath::Power(1+(1.25363-1)*x/0.769208,1.25363/(1-1.25363))+1.61721*TMath::Exp(2.86484-x/0.809167)) / (5655.77*TMath::Gaus(x,-7.79843,-0.217064)+26.0041*TMath::Landau(x,1.77719,0.608266,0)))",1.,50.);
			TF1 *DmesonCorr10 = new TF1("DmesonCorr10","(1/0.0585945)*((4.01931*TMath::Power(1+(1.25363-1)*x/0.769208,1.25363/(1-1.25363))+1.61721*TMath::Exp(2.86484-x/0.809167)) / (5655.77*TMath::Gaus(x,-7.79843,-0.217064)+26.0041*TMath::Landau(x,1.77719,0.608266,0)))",1.,50.);
			TF1 *DmesonCorr11 = new TF1("DmesonCorr11","(1/0.027615)*((4.01931*TMath::Power(1+(1.25363-1)*x/0.769208,1.25363/(1-1.25363))+1.61721*TMath::Exp(2.86484-x/0.809167)) / (5655.77*TMath::Gaus(x,-7.79843,-0.217064)+26.0041*TMath::Landau(x,1.77719,0.608266,0)))",1.,50.);
			TF1 *DmesonCorr12 = new TF1("DmesonCorr12","(1/0.0179361)*((4.01931*TMath::Power(1+(1.25363-1)*x/0.769208,1.25363/(1-1.25363))+1.61721*TMath::Exp(2.86484-x/0.809167)) / (5655.77*TMath::Gaus(x,-7.79843,-0.217064)+26.0041*TMath::Landau(x,1.77719,0.608266,0)))",1.,50.);
			task->SetDcorrFtn(DmesonCorr);
			task->SetDcorrFtn1(DmesonCorr1);
			task->SetDcorrFtn2(DmesonCorr2);
			task->SetDcorrFtn3(DmesonCorr3);
			task->SetDcorrFtn4(DmesonCorr4);
			task->SetDcorrFtn5(DmesonCorr5);
			task->SetDcorrFtn6(DmesonCorr6);
			task->SetDcorrFtn7(DmesonCorr7);
			task->SetDcorrFtn8(DmesonCorr8);
			task->SetDcorrFtn9(DmesonCorr9);
			task->SetDcorrFtn10(DmesonCorr10);
			task->SetDcorrFtn11(DmesonCorr11);
			task->SetDcorrFtn12(DmesonCorr12);
			std::cout<<"---------------------------------------------------------------------------------------"<<std::endl;
			std::cout<<"------------------------------D meson pt correction for most-central analysis----------"<<std::endl;
			std::cout<<"---------------------------------------------------------------------------------------"<<std::endl<<std::endl;
		}
		if(Dcorr==2)
		{
			TF1 *DmesonCorr = new TF1("DmesonCorr","((0.540505*TMath::Power(1+(1.19365-1)*x/1.09029,1.19365/(1-1.19365))+0.787186*TMath::Exp(1.09960-x/0.943388)) / (1.57403*TMath::Gaus(x,0.652698,1.77447)+1.51963*TMath::Landau(x, 3.41542, 1.46549, 0)))", 0.1, 50.);
			TF1 *DmesonCorr1 = new TF1("DmesonCorr1","(1.64317)*((0.540505*TMath::Power(1+(1.19365-1)*x/1.09029,1.19365/(1-1.19365))+0.787186*TMath::Exp(1.09960-x/0.943388)) / (1.57403*TMath::Gaus(x,0.652698,1.77447)+1.51963*TMath::Landau(x, 3.41542, 1.46549, 0)))", 1., 50.);
			TF1 *DmesonCorr2 = new TF1("DmesonCorr2","(1.8145)*((0.540505*TMath::Power(1+(1.19365-1)*x/1.09029,1.19365/(1-1.19365))+0.787186*TMath::Exp(1.09960-x/0.943388)) / (1.57403*TMath::Gaus(x,0.652698,1.77447)+1.51963*TMath::Landau(x, 3.41542, 1.46549, 0)))", 1.1, 50.);
			TF1 *DmesonCorr3 = new TF1("DmesonCorr3","(2.19212)*((0.540505*TMath::Power(1+(1.19365-1)*x/1.09029,1.19365/(1-1.19365))+0.787186*TMath::Exp(1.09960-x/0.943388)) / (1.57403*TMath::Gaus(x,0.652698,1.77447)+1.51963*TMath::Landau(x, 3.41542, 1.46549, 0)))", 1.3, 50.);
			TF1 *DmesonCorr4 = new TF1("DmesonCorr4","(2.61569)*((0.540505*TMath::Power(1+(1.19365-1)*x/1.09029,1.19365/(1-1.19365))+0.787186*TMath::Exp(1.09960-x/0.943388)) / (1.57403*TMath::Gaus(x,0.652698,1.77447)+1.51963*TMath::Landau(x, 3.41542, 1.46549, 0)))", 1.5, 50.);
			TF1 *DmesonCorr5 = new TF1("DmesonCorr5","(3.85597)*((0.540505*TMath::Power(1+(1.19365-1)*x/1.09029,1.19365/(1-1.19365))+0.787186*TMath::Exp(1.09960-x/0.943388)) / (1.57403*TMath::Gaus(x,0.652698,1.77447)+1.51963*TMath::Landau(x, 3.41542, 1.46549, 0)))", 2., 50.);
			TF1 *DmesonCorr6 = new TF1("DmesonCorr6","(5.28398)*((0.540505*TMath::Power(1+(1.19365-1)*x/1.09029,1.19365/(1-1.19365))+0.787186*TMath::Exp(1.09960-x/0.943388)) / (1.57403*TMath::Gaus(x,0.652698,1.77447)+1.51963*TMath::Landau(x, 3.41542, 1.46549, 0)))", 2.5, 50.);
			TF1 *DmesonCorr7 = new TF1("DmesonCorr7","(6.78326)*((0.540505*TMath::Power(1+(1.19365-1)*x/1.09029,1.19365/(1-1.19365))+0.787186*TMath::Exp(1.09960-x/0.943388)) / (1.57403*TMath::Gaus(x,0.652698,1.77447)+1.51963*TMath::Landau(x, 3.41542, 1.46549, 0)))", 3., 50.);
			TF1 *DmesonCorr8 = new TF1("DmesonCorr8","(9.70143)*((0.540505*TMath::Power(1+(1.19365-1)*x/1.09029,1.19365/(1-1.19365))+0.787186*TMath::Exp(1.09960-x/0.943388)) / (1.57403*TMath::Gaus(x,0.652698,1.77447)+1.51963*TMath::Landau(x, 3.41542, 1.46549, 0)))", 4., 50.);
			TF1 *DmesonCorr9 = new TF1("DmesonCorr9","(13.0435)*((0.540505*TMath::Power(1+(1.19365-1)*x/1.09029,1.19365/(1-1.19365))+0.787186*TMath::Exp(1.09960-x/0.943388)) / (1.57403*TMath::Gaus(x,0.652698,1.77447)+1.51963*TMath::Landau(x, 3.41542, 1.46549, 0)))", 5., 50.);
			TF1 *DmesonCorr10 = new TF1("DmesonCorr10","(18.5169)*((0.540505*TMath::Power(1+(1.19365-1)*x/1.09029,1.19365/(1-1.19365))+0.787186*TMath::Exp(1.09960-x/0.943388)) / (1.57403*TMath::Gaus(x,0.652698,1.77447)+1.51963*TMath::Landau(x, 3.41542, 1.46549, 0)))", 6., 50.);
			TF1 *DmesonCorr11 = new TF1("DmesonCorr11","(38.536)*((0.540505*TMath::Power(1+(1.19365-1)*x/1.09029,1.19365/(1-1.19365))+0.787186*TMath::Exp(1.09960-x/0.943388)) / (1.57403*TMath::Gaus(x,0.652698,1.77447)+1.51963*TMath::Landau(x, 3.41542, 1.46549, 0)))", 8., 50.);
			TF1 *DmesonCorr12 = new TF1("DmesonCorr12","(65.55)*((0.540505*TMath::Power(1+(1.19365-1)*x/1.09029,1.19365/(1-1.19365))+0.787186*TMath::Exp(1.09960-x/0.943388)) / (1.57403*TMath::Gaus(x,0.652698,1.77447)+1.51963*TMath::Landau(x, 3.41542, 1.46549, 0)))", 10., 50.);
			task->SetDcorrFtn(DmesonCorr);
			task->SetDcorrFtn1(DmesonCorr1);
			task->SetDcorrFtn2(DmesonCorr2);
			task->SetDcorrFtn3(DmesonCorr3);
			task->SetDcorrFtn4(DmesonCorr4);
			task->SetDcorrFtn5(DmesonCorr5);
			task->SetDcorrFtn6(DmesonCorr6);
			task->SetDcorrFtn7(DmesonCorr7);
			task->SetDcorrFtn8(DmesonCorr8);
			task->SetDcorrFtn9(DmesonCorr9);
			task->SetDcorrFtn10(DmesonCorr10);
			task->SetDcorrFtn11(DmesonCorr11);
			task->SetDcorrFtn12(DmesonCorr12);
			std::cout<<"---------------------------------------------------------------------------------------"<<std::endl;
			std::cout<<"------------------------------D meson pt correction for semi-central analysis----------"<<std::endl;
			std::cout<<"---------------------------------------------------------------------------------------"<<std::endl<<std::endl;
		}
		
		// Lc meson pt correction
		if(Lccorr==1)
		{
			TF1 *LcCorr = new TF1("LcCorr","12.6935*TMath::Exp(-0.568250*x)+0.00681196",1., 30.);
			task->SetLccorrFtn(LcCorr);
			std::cout<<"---------------------------------------------------------------------------------------"<<std::endl;
			std::cout<<"------------------------------LambdaC pt correction for most-central analysis----------"<<std::endl;
			std::cout<<"---------------------------------------------------------------------------------------"<<std::endl<<std::endl;
		}
		if(Lccorr==2)
		{
			TF1 *LcCorr = new TF1("LcCorr","7.45534*TMath::Exp(-0.551280*x)+7.21259e-03", 1., 30.);
			task->SetLccorrFtn(LcCorr);
			std::cout<<"---------------------------------------------------------------------------------------"<<std::endl;
			std::cout<<"------------------------------LambdaC pt correction for semi-central analysis----------"<<std::endl;
			std::cout<<"---------------------------------------------------------------------------------------"<<std::endl<<std::endl;
		}
	}
    
	mgr->AddTask(task);
	//added to run on train
	TString fileName = AliAnalysisManager::GetCommonFileName();
	fileName += ":ElectroID_";
	fileName += uniqueID;

	//Create containers for input/output -> to run locally
	AliAnalysisDataContainer *coutput = mgr->CreateContainer(Form("ccontainer0_%s",uniqueID.Data()),TList::Class(),AliAnalysisManager::kOutputContainer,fileName);
	//Connect input/output
	mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(task, 1, coutput);
	
	return task;
}

