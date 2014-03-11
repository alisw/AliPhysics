 //DEFINITION OF A FEW CONSTANTS
//----------------------------------------------------

/* $Id$ */

AliAnalysisTaskDStarCorrelations *AddTaskDStarCorrelations(AliAnalysisTaskDStarCorrelations::CollSyst syst, 
                                                            Bool_t theMCon, Bool_t mixing, Bool_t UseReco=kTRUE,Bool_t UseHadChannelinMC,Bool_t fullmode = kFALSE,Bool_t UseEffic=kFALSE,Bool_t UseDEffic = kFALSE,
                                                        AliAnalysisTaskDStarCorrelations::DEffVariable var,
                                                            Int_t trackselect =1, Int_t usedispl =0, Int_t nbins, Float_t DStarSigma, Float_t D0Sigma, Float_t D0SBSigmaLow, Float_t D0SBSigmaHigh, Float_t eta,
                                                                     TString DStarCutsFile, TString TrackCutsFile,
                                                                     TString suffix = "")
{

 
    
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error(" AliAnalysisTaskDStarCorrelations", "No analysis manager to connect to.");
    return NULL;
  } 

	cout << "==========================================================" << endl;
    cout << "Set Inputs : " << endl;
    cout << " " << endl;
    if(syst == AliAnalysisTaskDStarCorrelations::pp) cout << "Running on pp @ 7 TeV" << endl;
    if(syst == AliAnalysisTaskDStarCorrelations::pA) cout << "Running on pPb @ 5.02 TeV" << endl;
    if(syst == AliAnalysisTaskDStarCorrelations::AA) cout << "Running on PbPb @ 2.76 TeV" << endl;
    
    if(theMCon) cout << "Analysis on MonteCarlo" << endl;
    else cout << "Analysis on Data" << endl;
    
    if(mixing) cout << "Analysis on Mixed Events" << endl;
    else cout << "Analysis on Single Events" << endl;
    
    if(UseReco) cout << "Using Reconstructed objects" << endl;
    else cout << "Using Pure MC information " << endl;

    if(fullmode) cout << "Running in full mode" << endl;
    else cout << "Running in fast mode " << endl;
    
    if(UseEffic) cout << "Using single track efficiency map" << endl;
    else cout << "Not Using single track efficiency map " << endl;
    
    if(UseDEffic) cout << "Using Dmeson efficiency map" << endl;
    else cout << "Not Using Dmeson efficiency map " << endl;
    
 
    if(var == AliAnalysisTaskDStarCorrelations::kNone) cout << "Applying D Efficiency map vs pT " << endl;
    if(var == AliAnalysisTaskDStarCorrelations::kMult) cout << "Applying D Efficiency map vs pT vs Multiplicity" << endl;
    if(var == AliAnalysisTaskDStarCorrelations::kCentr) cout << "Applying D Efficiency map vs pT vs Centrality" << endl;
    if(var == AliAnalysisTaskDStarCorrelations::kRapidity) cout << "Applying D Efficiency map vs pT vs Rapidity" << endl;
    if(var == AliAnalysisTaskDStarCorrelations::kEta) cout << "Applying D Efficiency map vs pT vs Eta" << endl;

    if(trackselect == 1) cout << "Correlating with hadrons" << endl;
    if(trackselect == 2) cout << "Correlating with kaons" << endl;
    if(trackselect == 3) cout << "Correlating with kzeros" << endl;
    
    if(usedispl == 0) cout << "Not using displacement cut" << endl;
    if(usedispl == 1) cout << "Using absolute displacement cut" << endl;
    if(usedispl == 2) cout << "Using relative displacement cut" << endl;
    
    
    cout << "Number of bins in deltaphi = " << nbins << endl;
    
    cout << "N of Sigmas in D* selection =" << DStarSigma << endl;
    cout << "N of Sigmas in D0 selection = " << D0Sigma << endl;
    cout << "D0 Sidebands taken from  = " << D0SBSigmaLow << " - " << D0SBSigmaHigh << " sigmas " << endl; endl;
    
   
    
    cout << "DStar cut object:     " << DStarCutsFile << endl;
    cout << "Tracks cut object:    " << TrackCutsFile << endl;
  
    
    cout << "==========================================================" << endl;
    //gSystem->Sleep(2000);
    
//	TString DCutObjPath = "CutObjects/DStar/";
	 	
	
// ******************************** OPENING THE D* CUTS ************************************
    cout << "Getting D meson cut object from file \n" << DStarCutsFile.Data() << "\n " << endl;
    TFile* filecuts=TFile::Open(DStarCutsFile.Data());
    if(!filecuts->IsOpen()){
    cout<<"DStar cut object file not found: exit"<<endl;
    return;
    }
    
    AliRDHFCutsDStartoKpipi* RDHFDStartoKpipi=new AliRDHFCutsDStartoKpipi();
    RDHFDStartoKpipi = (AliRDHFCutsDStartoKpipi*)filecuts->Get("DStartoKpipiCuts");
    RDHFDStartoKpipi->SetName("DStartoKpipiCuts");
    
    // mm let's see if everything is ok
	if(!RDHFDStartoKpipi){
		cout<<"Specific AliRDHFCuts not found"<<endl;
		return;
	}
    
        
       // RDHFDStartoKpipi->SetTriggerClass("");
       // RDHFDStartoKpipi->SetTriggerMask(AliVEvent::kCentral);
    
    
	
// ******************************** OPENING THE ASSOCIATED TRACK CUTS ************************************
	cout << "Getting associated track cut object from file \n" << TrackCutsFile.Data() << "\n " << endl;
	TFile* filecuts2=TFile::Open(TrackCutsFile.Data());
	  if(!filecuts2->IsOpen()){
		  cout<<"Track cut object file not found: exit"<<endl;
		  return;
    }
 	AliHFAssociatedTrackCuts* corrCuts=new AliHFAssociatedTrackCuts();
	corrCuts = (AliHFAssociatedTrackCuts*)filecuts2->Get("AssociatedCuts");
	corrCuts->SetName("AssociatedCuts");
	corrCuts->PrintAll();
	if(!corrCuts){
		cout<<"Specific associated track cuts not found"<<endl;
		return;
	}
	

// ******************************** SELECTING THE MC PROCESS  ************************************

	TString selectMCproc = "";
	
	Int_t NMCevents = corrCuts->GetNofMCEventType();
	for(Int_t k=0; k<NMCevents; k++){
		Int_t * MCEventType = corrCuts->GetMCEventType();
		selectMCproc += Form("%d",MCEventType[k]);
	}

	 cout << "Select process string = " << selectMCproc << endl;
	


// ******************************** CREATING THE TASK ************************************
 
    printf("CREATE TASK \n");
  // create the task
  AliAnalysisTaskDStarCorrelations *task = new AliAnalysisTaskDStarCorrelations("AliAnalysisTaskDStarCorrelations",RDHFDStartoKpipi,corrCuts,syst,fullmode);
	
	// Setters

	if(!theMCon) {
		printf("Analysis on Data - reconstruction only!");
		UseReco = kTRUE;
        printf("Analysis on Data - hadronic channel only!");
		UseHadChannelinMC = kFALSE;
	}
    
    

	
    task->SetNofPhiBins(nbins);
	task->SetMonteCarlo(theMCon);
	task->SetUseMixing(mixing);
	task->SetCorrelator(trackselect) ;
	task->SetUseDisplacement(usedispl);
	//task->SetCollSys(syst);
	task->SetLevelOfDebug(2);
	task->SetUseReconstruction(UseReco); // set kTRUE for Using Reconstruction, kFALSe for MC Truth
    task->SetDMesonSigmas(DStarSigma,D0Sigma,D0SBSigmaLow,D0SBSigmaHigh);
	//task->SetDMesonSigmas(sigmas);
    //task->SetDMesonSigmas(sigmas);
    task->SetUseEfficiencyCorrection(UseEffic);
    task->SetUseDmesonEfficiencyCorrection(UseDEffic);
    
    task->SetEfficiencyVariable(var);
    task->SetMaxDStarEta(eta);
    task->SetUseHadronicChannelAtKineLevel(UseHadChannelinMC);
	//task->SetBkgEstimationMethod(AliAnalysisTaskDStarCorrelations::kDStarSB);

	if(trackselect == 1) Info(" AliAnalysisTaskDStarCorrelations","Correlating D* with charged hadrons \n");
	else if(trackselect == 2) Info(" AliAnalysisTaskDStarCorrelations","Correlating D* with charged kaons \n");
	else if(trackselect == 3) Info(" AliAnalysisTaskDStarCorrelations","Correlating D* with reconstructed K0s \n");
	else Fatal(" AliAnalysisTaskDStarCorrelations","Nothing to correlate with!");
	if(mixing) Info (" AliAnalysisTaskDStarCorrelations","Event Mixing Analysis\n");
	if(!mixing) Info (" AliAnalysisTaskDStarCorrelations","Single Event Analysis \n");

  // Create and connect containers for input/output
	//TString dcavalue = " ";
    
  TString contname1 = "OutputEvent";
    TString contname2 = "OutputDmeson";
    TString contname3 = "OutputTracks";
    TString contname4 = "OutputEventMixing";
    TString contname5 = "OutputCorrelations";
    TString contname6 = "OutputMC";
    TString counter = "NormCounter";
    TString cutname1 = "Dcuts" ;
	TString cutname2 = "hadroncuts" ;
    
    //if(suffix != ""){
        contname1 += suffix;
        contname2 += suffix;
        contname3 += suffix;
        contname4 += suffix;
        contname5 += suffix;
         contname6 += suffix;
        counter += suffix;
        cutname1 += suffix;
        cutname2 += suffix;
        
        
   // }
    
   // else{
    
	//if(!theMCon) TString contname = "Data";
	//if(theMCon) TString contname = "MonteCarlo";
	//TString contname2 = "MC";
	//if(trackselect ==1) TString particle = "Hadron";
	//if(trackselect ==2) TString particle = "Kaon";
	//if(trackselect ==3) TString particle = "KZero";
	
	
	TString outputfile = AliAnalysisManager::GetCommonFileName();
	//TString outputfileMC = AliAnalysisManager::GetCommonFileName();
	//TString counter = "NormCounter";
	outputfile += ":PWGHF_HFCJ_";
	//outputfileMC += ":PWGHF_HFCJ_";
	
	if(!mixing) {
		outputfile += "SE";
		//outputfileMC += "SE";
		/*contname1 += "SE";
		contname2 += "SE";
        contname3 += "SE";
        contname4 += "SE";
        contname5 += "SE";
        counter+= "SE";
        cutname1 += "SE";
		cutname2 += "SE";*/
		
	}
	if(mixing){
		outputfile += "ME";
		//outputfileMC += "ME";
		
	}
    

    
	/*outputfile += "Dphi_DStar";
	//outputfileMC += "Dphi_DStar";
	outputfile += particle;
	//outputfileMC += particle;
	*/
    if(UseEffic){
        outputfile += "_EffY_";
		//outputfileMC += "_EffY_";
	
    }
    
    if(!UseEffic){
        outputfile += "_EffN_";
		//outputfileMC += "_EffN_";
		
    }
    
    if(UseDEffic){
        TString string = "DEffY_";
        
        if(var == AliAnalysisTaskDStarCorrelations::kNone) string += "vsPt_";
        if(var == AliAnalysisTaskDStarCorrelations::kMult) string += "vsPtMult_";
        if(var == AliAnalysisTaskDStarCorrelations::kCentr) string += "vsPCentrt_";
        if(var == AliAnalysisTaskDStarCorrelations::kRapidity) string += "vsPtY_";
        if(var == AliAnalysisTaskDStarCorrelations::kEta) string += "vsPtEta_";
        
        outputfile += string;
		//outputfileMC += string;
		
    }
    
    if(!UseDEffic){
        outputfile += "DEffN_";
		//outputfileMC += "DEffN_";
		
    }
	
    
    outputfile += Form("%d_bins_",nbins);
	//outputfileMC += Form("%d_bins_",nbins);

    
	
	outputfile += suffix;
	//outputfileMC += Form("task_%d",tasknumber);
	
	
/*	outputfile += TrackCutObjNamePrefix;
	//outputfileMC += TrackCutObjNamePrefix;
	cutname += TrackCutObjNamePrefix;
	cutname2 += TrackCutObjNamePrefix;
	contname += TrackCutObjNamePrefix;
	contname2 += TrackCutObjNamePrefix;
	counter+= TrackCutObjNamePrefix;
*/	
	outputfile += selectMCproc;
	//outputfileMC += selectMCproc;
	
	
	TString reco = "";
	
	if(UseReco) reco = "_reco";
	if(!UseReco) {
        if(UseHadChannelinMC) reco = "_MCTruthHadChan";
        if(!UseHadChannelinMC) reco = "_MCTruthAllChan";
	}
	outputfile += reco;
	//outputfileMC += reco;

	
	TString nsigma = Form("_%.0f_%.0f%.0f%.0f_sigmas",DStarSigma,D0Sigma,D0SBSigmaLow,D0SBSigmaHigh);
	
	//cout << "nsigma string = "<< nsigma << endl;
	
	outputfile += nsigma;
	//outputfileMC += nsigma;

	
	/*
	cout << contname << endl;
	cout << contname2 << endl;
	cout << cutname << endl;
	cout << cutname2 << endl;
	cout << counter << endl;
	cout << outputfile << endl;*/
	//return;
	//}// end else
  mgr->AddTask(task);
  // ------ input data ------
  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
  
  // ----- output data -----
  
  // output TH1I for event counting
	
	//TLists
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname1, TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contname2, TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(contname3, TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
     AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(contname4, TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
     AliAnalysisDataContainer *coutput5 = mgr->CreateContainer(contname5, TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    AliAnalysisDataContainer *coutput6 = mgr->CreateContainer(contname6, TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
    // Normalization Counter
    AliAnalysisDataContainer *coutput7 = mgr->CreateContainer(counter,AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());
   // Cut Objects
  AliAnalysisDataContainer *coutput8 = mgr->CreateContainer(cutname1,AliRDHFCutsDStartoKpipi::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data()); //cuts D
  AliAnalysisDataContainer *coutput9 = mgr->CreateContainer(cutname2,AliHFAssociatedTrackCuts::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data()); //cuts tracks

  
	
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);
  mgr->ConnectOutput(task,4,coutput4);
  mgr->ConnectOutput(task,5,coutput5);
     mgr->ConnectOutput(task,6,coutput6);
    mgr->ConnectOutput(task,7,coutput7);
    mgr->ConnectOutput(task,8,coutput8);
    mgr->ConnectOutput(task,9,coutput9);

  return task ;

}

