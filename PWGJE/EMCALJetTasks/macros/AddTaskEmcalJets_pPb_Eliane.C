// $Id$
AliAnalysisTaskSAJF* AddTaskEmcalJets_pPb_Eliane(
   const char *outfilename    = "AnalysisOutput.root",
   const char *nRhosCh        = "rhoCh",
   const char *nRhosChEm      = "rhoChEm",
   const char *nRhosEm        = "rhoEm",
   const Double_t scale       = 1.0,
   const Double_t radius      = 0.2,
   const char* usedTracks     = "PicoTracks",
   const char* outClusName    = "caloClustersCorr",
   const Double_t minTrackPt  = 0.15,
   const Double_t minClusterPt = 0.30,
   const char *CentEst         = "V0A",
   const Int_t type            =1
)
{  
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//
	//  This task calls a jet finder (EMCal Jet) and then analyzes the jet with an analyzer ()
	//
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -


	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//
	//    Run the jet finder and rho tasks first  (common task)
    //
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	  cout<<"a------------------------------here we go 11 !!!"<<endl;
	  cout<<"a------------------------------here we go 11 !!!"<<endl;
	  cout<<"a------------------------------here we go 11 !!!"<<endl;
	  cout<<"a------------------------------here we go 11 !!!"<<endl;
	  cout<<"a------------------------------here we go 11 !!!"<<endl;
	  cout<<"a------------------------------here we go 11 !!!"<<endl;
	  cout<<"a------------------------------here we go 11 !!!"<<endl;
	  cout<<"a------------------------------here we go 11 !!!"<<endl;
	  cout<<"a------------------------------here we go 11 !!!"<<endl;
	  cout<<"a------------------------------here we go 11 !!!"<<endl;
	  cout<<"a------------------------------here we go 11 !!!"<<endl;
	  cout<<"a------------------------------here we go 11 !!!"<<endl;
	  cout<<"a------------------------------here we go 11 !!!"<<endl;

	gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");
	gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoSparse.C");

	// Some constants for the jet finders
	const Int_t cKT                 = 0;
	const Int_t cANTIKT             = 1;
	const Int_t cFULLJETS           = 0;
	const Int_t cCHARGEDJETS        = 1;
	const Int_t cNEUTRALJETS        = 2;

	char *typeTPC                = "TPC";
	char *typeEMC                = "EMC";

//probably not needed SA task uses 0.6 and the radius fed into the jet finder
	//float AreaCut = 0.6*radius*radius*TMath::Pi();

	TF1 *sfunc=new TF1("sfunc","[0]*x*x+[1]*x+[2]",-1,100);
	sfunc->SetParameter(0,0.0);
	sfunc->SetParameter(1,0.0);
	sfunc->SetParameter(2,scale);

	TString nJets("");

	TString scaledname(Form("%s_Scaled", nRhosCh));
	TString newrhoname(Form("%s_All", nRhosCh));
	//TString scaledname(Form("%s_Scaled", newrhoname));

   // Here are the jet finders
	AliEmcalJetTask* jetFinderTask = AddTaskEmcalJet(AliEmcalJetTask::kAKT | AliEmcalJetTask::kFullJet,usedTracks,outClusName,minTrackPt,minClusterPt,0.005,radius,1,"Jet",0.1);
	//const UInt_t type           = AliEmcalJetTask::kAKT | AliEmcalJetTask::kFullJet | AliEmcalJetTask::kR040Jet,
	//but the last info should match the "radius" info,
	//I don'T know if this is checked within the code to avoide confusion

	jetFinderTask->GetName();


	// Here are the rho tasks
//	AliAnalysisTaskRhoSparse *rhochtask = AddTaskRhoSparse(jetFinderTaskChBack->GetName(),jetFinderTaskChSig->GetName(),usedTracks,outClusName,nRhosCh,radius,typeTPC,0.01,0,0,sfunc,0,kTRUE,nRhosCh);
//	rhochtask->SetCentralityEstimator(CentEst);

//	AliAnalysisTaskRhoSparse *rhochalltask = AddTaskRhoSparse(jetFinderTaskChBackall->GetName(),jetFinderTaskChSig->GetName(),usedTracks,outClusName,newrhoname,radius,0,0.0,0,0,sfunc,0,kTRUE,newrhoname);
//	rhochtask->SetCentralityEstimator(CentEst);

//	AliAnalysisTaskRhoSparse *rhochemtask = AddTaskRhoSparse(jetFinderTaskChEmBack->GetName(),jetFinderTask->GetName(),usedTracks,outClusName,nRhosChEm,radius,typeEMC,0.01,0,0,0,0,kTRUE,nRhosChEm);
//	rhochemtask->SetCentralityEstimator(CentEst);

	nJets=jetFinderTask->GetName();

	/*
	--->>	no idea what that is about

		gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskDeltaPt.C");

		TString deltaname(Form("DeltaPt_%s_Scaled", nRhosCh));
		AliAnalysisTaskDeltaPt* deltapt = AddTaskDeltaPt(usedTracks,outClusName,nJets,"","","","","",scaledname,radius,AreaCut,minTrackPt,minClusterPt,typeEMC,deltaname);
		deltapt->SetCentralityEstimator(CentEst);

		TString chdeltaname(Form("DeltaPt_%s", nRhosCh));
		AliAnalysisTaskDeltaPt* deltaptch = AddTaskDeltaPt(usedTracks,"",nJetsCh,"","","","","",nRhosCh,radius,AreaCut,minTrackPt,minClusterPt,typeTPC,chdeltaname);
		deltaptch->SetCentralityEstimator(CentEst);

		TString emcdeltaname(Form("DeltaPt_%s", nRhosChEm));
		AliAnalysisTaskDeltaPt* deltaptEMC = AddTaskDeltaPt(usedTracks,outClusName,nJets,"","","","","",nRhosChEm,radius,AreaCut,minTrackPt,minClusterPt,typeEMC,emcdeltaname);
		deltaptEMC->SetCentralityEstimator(CentEst);

		gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskScale.C");

		Int_t radlabel=(Int_t)floor(radius*100+0.5);
		Int_t mincluslabel=(Int_t)floor(minClusterPt*1000+0.5);
		Int_t mintracklabel=(Int_t)floor(minTrackPt*1000+0.5);
		TString scalename(Form("Scale_R0%d", radlabel));

		AliAnalysisTaskScale* scaletask = AddTaskScale(usedTracks,outClusName,minTrackPt,minClusterPt,scalename);
		scaletask->SetCentralityEstimator(CentEst);
		scaletask->SetScaleFunction(sfunc);
	 */
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//
	//    Run the jet analyzer (personal task)
	//
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    //create the according class object for the Taks
	gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskSAJF.C");
	//AliAnalysisTaskSAJF* JetAnalysisEliane;

	AliAnalysisTaskSAJF* JetAnalysisEliane = AddTaskSAJF(usedTracks,outClusName,jetFinderTask->GetName(),"",radius,1.0,0.6,"EMCAL",0,"AliAnalysisTaskEmcalJets_pPb_Eliane");


	//That might be important!!
	//spectratask->SetCentralityEstimator(CentEst);

	//what is the scalefactor about
	//--->> I need to add the correct centrality estimator!! spectratask->SetCentralityEstimator(CentEst);
//where is fNcentBins defined?  ali analysistask emcal  0-10, 10-30, 30-50


/*
	TString name(Form("SpectraMECpA_%s", nJets.Data()));
	AliAnalysisTaskEmcalJetSpectraMECpA *spectratask = new AliAnalysisTaskEmcalJetSpectraMECpA(name);
	spectratask->SetJetsName(nJets.Data());
	spectratask->SetCentralityEstimator(CentEst);

	if(type==0){
		//spectratask->SetAnaType(typeTPC);
		spectratask->SetAnaType(0);
		spectratask->SetRhoName(nRhosCh);
	}else{
		//spectratask->SetAnaType(typeEMC);
		spectratask->SetAnaType(1);
		if(!(usedTracks=="")) spectratask->SetRhoName(scaledname);
		else spectratask->SetRhoName(nRhosEm);
	}
	spectratask->SetJetPhiLimits(minPhi,maxPhi);
	spectratask->SetJetEtaLimits(minEta,maxEta);
	spectratask->SetJetAreaCut(AreaCut);
	spectratask->SetTracksName(usedTracks);

	//-------------------------------------------------------
	// Final settings, pass to manager and set the containers
	//-------------------------------------------------------


	mgr->AddTask(spectratask);

	// Create containers for input/output
	mgr->ConnectInput (spectratask, 0, mgr->GetCommonInputContainer());
	AliAnalysisDataContainer *cospectra = mgr->CreateContainer(name,
			TList::Class(),
			AliAnalysisManager::kOutputContainer,
			outfilename);
	mgr->ConnectOutput(spectratask,1,cospectra);

*/

	return JetAnalysisEliane;
}
