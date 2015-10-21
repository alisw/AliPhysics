// $Id$
AliAnalysisTaskSAJF* AddTaskEmcalJets_pPb_Eliane(
   const char *nRhosCh        = "rhoCh",
   const char *nRhosChEm      = "rhoChEm",
   const char *nRhosEm        = "rhoEm",
   const Double_t scale       = 1.0,
   const Double_t radius      = 0.2,
   const char* usedTracks     = "PicoTracks",
   const char* outClusName    = "caloClustersCorr",
   const Double_t minTrackPt  = 0.15,
   const Double_t minClusterPt = 0.30,
   const char *CentEst         = "V0A"
)
{  
	//>>  This task calls a jet finder (EMCal Jet) and then analyzes the jet with an analyzer ()

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//
	//    Run the jet finder (general tasks)
    //
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskEmcalJet.C");

    //hard coded variables for the jet finder
	const Int_t RecScheme  =1;         //RecombinationScheme: 0=E_scheme, 1=pt_scheme
	const Double_t minJetPt    = 0.1;  // is supposed to make the jet analysis faster as the jet container don't get so big


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


	TString nJets("");

   // Here are the jet finders
	AliEmcalJetTask* jetFinderTask = AddTaskEmcalJet(AliEmcalJetTask::kAKT | AliEmcalJetTask::kFullJet,usedTracks,outClusName,minTrackPt,minClusterPt,0.005,radius,RecScheme,"Jet",minJetPt);
	//const UInt_t type           = AliEmcalJetTask::kAKT | AliEmcalJetTask::kFullJet | AliEmcalJetTask::kR040Jet,
	//but the last info should match the "radius" info,
	//I don'T know if this is checked within the code to avoide confusion
	nJets=jetFinderTask->GetName();


	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	//
	//    Run the rho tasks (general tasks)
    //
	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	gROOT->LoadMacro("$ALICE_PHYSICS/PWGJE/EMCALJetTasks/macros/AddTaskRhoSparse.C");

	TF1 *sfunc=new TF1("sfunc","[0]*x*x+[1]*x+[2]",-1,100);
	sfunc->SetParameter(0,0.0);
	sfunc->SetParameter(1,0.0);
	sfunc->SetParameter(2,scale);

	TString scaledname(Form("%s_Scaled", nRhosCh));
	TString newrhoname(Form("%s_All", nRhosCh));
	//TString scaledname(Form("%s_Scaled", newrhoname));


//	AliAnalysisTaskRhoSparse *rhochtask = AddTaskRhoSparse(jetFinderTaskChBack->GetName(),jetFinderTaskChSig->GetName(),usedTracks,outClusName,nRhosCh,radius,typeTPC,0.01,0,0,sfunc,0,kTRUE,nRhosCh);
//	rhochtask->SetCentralityEstimator(CentEst);

//	AliAnalysisTaskRhoSparse *rhochalltask = AddTaskRhoSparse(jetFinderTaskChBackall->GetName(),jetFinderTaskChSig->GetName(),usedTracks,outClusName,newrhoname,radius,0,0.0,0,0,sfunc,0,kTRUE,newrhoname);
//	rhochtask->SetCentralityEstimator(CentEst);

//	AliAnalysisTaskRhoSparse *rhochemtask = AddTaskRhoSparse(jetFinderTaskChEmBack->GetName(),jetFinderTask->GetName(),usedTracks,outClusName,nRhosChEm,radius,typeEMC,0.01,0,0,0,0,kTRUE,nRhosChEm);
//	rhochemtask->SetCentralityEstimator(CentEst);



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
	AliAnalysisTaskSAJF* JetAnalysisEliane = AddTaskSAJF(usedTracks,outClusName,jetFinderTask->GetName(),"",radius,1.0,0.6,"EMCAL",0,"EmcalJets_pPb_Eliane");

	//That might be important!!
	//spectratask->SetCentralityEstimator(CentEst);

	//what is the scalefactor about
	//--->> I need to add the correct centrality estimator!! spectratask->SetCentralityEstimator(CentEst);
    //where is fNcentBins defined?  ali analysistask emcal  0-10, 10-30, 30-50

	return JetAnalysisEliane;
}
