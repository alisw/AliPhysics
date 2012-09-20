 //DEFINITION OF A FEW CONSTANTS
//----------------------------------------------------

/* $Id$ */

AliAnalysisTaskDStarCorrelations *AddTaskDStarCorrelations(Bool_t runOnPbPb,Bool_t theMCon, Bool_t mixing, Bool_t UseReco = kTRUE, Int_t trackselect =1, Int_t usedispl =0, TString DCutObjNamePrefix = "_corr", TString TrackCutObjNamePrefix = "")
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error(" AliAnalysisTaskDStarCorrelations", "No analysis manager to connect to.");
    return NULL;
  } 

	
	TString DCutObjPath = "~/CorrelationAnalysis/CutObjects/DStar/";
	
	
	
	TString DCutObjName = "DStartoKpipiCuts";
	DCutObjName += DCutObjNamePrefix;
	DCutObjName += ".root";
	
	DCutObjName.Prepend(DCutObjPath.Data());
	
	cout << "D* cut object is " << DCutObjName << endl;
  TFile* filecuts=new TFile(DCutObjName.Data());
  if(!filecuts->IsOpen()){
    cout<<"DStar cut object file not found: exit"<<endl;
    return;
  }  
	
	TString TrackCutObjPath = "~/CorrelationAnalysis/CutObjects/AssocTracks/";
	
	TString TrackCutObjName = "AssocPartCuts";
	TrackCutObjName += TrackCutObjNamePrefix;
	TrackCutObjName += ".root";
	
	TrackCutObjName.Prepend(TrackCutObjPath.Data());
	
	cout << "tracks cut object is " << TrackCutObjName << endl;
	  TFile* filecuts2=new TFile(TrackCutObjName.Data());
	  if(!filecuts2->IsOpen()){
		  cout<<"Track cut object file not found: exit"<<endl;
		  return;
  }

  AliRDHFCutsDStartoKpipi* RDHFDStartoKpipi=new AliRDHFCutsDStartoKpipi();
  RDHFDStartoKpipi = (AliRDHFCutsDStartoKpipi*)filecuts->Get("DStartoKpipiCuts");
  RDHFDStartoKpipi->SetName("DStartoKpipiCuts");
	
	
	AliHFAssociatedTrackCuts* corrCuts=new AliHFAssociatedTrackCuts();
	corrCuts = (AliHFAssociatedTrackCuts*)filecuts2->Get("AssociatedCuts");
	corrCuts->SetName("AssociatedCuts");
	corrCuts->PrintAll();
	
	// mm let's see if everything is ok
	if(!RDHFDStartoKpipi){
		cout<<"Specific AliRDHFCuts not found"<<endl;
		return;
	} 
	
	if(!corrCuts){
		cout<<"Specific associated track cuts not found"<<endl;
		return;
	} 
	
	TString selectMCproc = "";
	
	Int_t NMCevents = corrCuts->GetNofMCEventType();
	for(Int_t k=0; k<NMCevents; k++){
		Int_t * MCEventType = corrCuts->GetMCEventType();
		selectMCproc += Form("%d",MCEventType[k]);
	}

	cout << "Select process string = " << selectMCproc << endl;
	


  //CREATE THE TASK
  printf("CREATE TASK \n");
  // create the task
  AliAnalysisTaskDStarCorrelations *task = new AliAnalysisTaskDStarCorrelations("AliAnalysisTaskDStarCorrelations",RDHFDStartoKpipi,corrCuts);
	
	// Setters

	if(!theMCon) {
		printf("Analysis on Data - reconstruction only!");
		UseReco = kTRUE;
	}
	
	task->SetMonteCarlo(theMCon);
	task->SetUseMixing(mixing);
	task->SetCorrelator(trackselect) ;
	task->SetUseDisplacement(usedispl);
	task->SetRunPbPb(runOnPbPb);
	task->SetLevelOfDebug(2);
	task->SetUseReconstruction(UseReco); // set kTRUE for Using Reconstruction, kFALSe for MC Truth
	

	if(trackselect == 1) Info(" AliAnalysisTaskDStarCorrelations","Correlating D* with charged hadrons \n");
	else if(trackselect == 2) Info(" AliAnalysisTaskDStarCorrelations","Correlating D* with charged kaons \n");
	else if(trackselect == 3) Info(" AliAnalysisTaskDStarCorrelations","Correlating D* with reconstructed K0s \n");
	else Fatal(" AliAnalysisTaskDStarCorrelations","Nothing to correlate with!");
	if(mixing) Info (" AliAnalysisTaskDStarCorrelations","Event Mixing Analysis\n");
	if(!mixing) Info (" AliAnalysisTaskDStarCorrelations","Single Event Analysis \n");

  // Create and connect containers for input/output
	//TString dcavalue = " ";
	if(!theMCon) TString contname = "Data";
	if(theMCon) TString contname = "MonteCarlo";
	TString contname2 = "MC";
	if(trackselect ==1) TString particle = "Hadron";
	if(trackselect ==2) TString particle = "Kaon";
	if(trackselect ==3) TString particle = "KZero";
	
	TString cutname = "cuts" ;
	TString cutname2 = "hadroncuts" ;
	TString outputfile = AliAnalysisManager::GetCommonFileName();
	TString outputfileMC = AliAnalysisManager::GetCommonFileName();
	TString counter = "NormCounter";
	outputfile += ":PWGHF_D2H_";
	outputfileMC += ":PWGHF_D2H_";
	
	if(!mixing) {
		outputfile += "SE";
		outputfileMC += "SE";
		contname += "SE";
		contname2 += "SE";
		cutname += "SE";
		cutname2 += "SE";
		counter+= "SE";
	}
	if(mixing){
		outputfile += "ME";
		outputfileMC += "ME";
		contname += "ME";
		contname2 += "ME";
		cutname += "ME";
		cutname2 += "ME";
		counter+= "ME";
	}
	outputfile += "Dphi_DStar";
	outputfileMC += "Dphi_DStar";
	outputfile += particle;
	outputfileMC += particle;
	cutname += particle;
	cutname2 += particle;
	contname += particle;
	contname2 += particle;
	counter+= particle;
	
	
	
	outputfile += DCutObjNamePrefix;
	outputfileMC += DCutObjNamePrefix;
	cutname += DCutObjNamePrefix;
	cutname2 += DCutObjNamePrefix;
	contname += DCutObjNamePrefix;
	contname2 += DCutObjNamePrefix;
	counter+= DCutObjNamePrefix;
	
	outputfile += TrackCutObjNamePrefix;
	outputfileMC += TrackCutObjNamePrefix;
	cutname += TrackCutObjNamePrefix;
	cutname2 += TrackCutObjNamePrefix;
	contname += TrackCutObjNamePrefix;
	contname2 += TrackCutObjNamePrefix;
	counter+= TrackCutObjNamePrefix;
	
	outputfile += selectMCproc;
	outputfileMC += selectMCproc;
	cutname += selectMCproc;
	cutname2 += selectMCproc;
	contname += selectMCproc;
	contname2 += selectMCproc;
	counter+= selectMCproc;
	
	TString reco = "";
	
	if(UseReco) reco = "_reco";
	if(!UseReco) reco = "_MCTruth";
	
	outputfile += reco;
	outputfileMC += reco;
	cutname += reco;
	cutname2 += reco;
	contname += reco;
	contname2 += reco;
	counter+= reco;
	
	
	cout << contname << endl;
	cout << contname2 << endl;
	cout << cutname << endl;
	cout << cutname2 << endl;
	cout << counter << endl;
	cout << outputfile << endl;
	//return;
	
  mgr->AddTask(task);
  // ------ input data ------
  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
  
  // ----- output data -----
  
  // output TH1I for event counting
	
	//TLists
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname, TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(contname2, TList::Class(),AliAnalysisManager::kOutputContainer,outputfileMC.Data());
   // Cut Objects
  AliAnalysisDataContainer *coutput3 = mgr->CreateContainer(cutname,AliRDHFCutsDStartoKpipi::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data()); //cuts D
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(cutname2,AliHFAssociatedTrackCuts::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data()); //cuts tracks
   // Normalization Counter
  AliAnalysisDataContainer *coutput5 = mgr->CreateContainer(counter,AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());
  
	
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutput3);
  mgr->ConnectOutput(task,4,coutput4);
  mgr->ConnectOutput(task,5,coutput5);

  return task ;

}

