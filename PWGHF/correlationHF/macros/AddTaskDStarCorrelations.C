 //DEFINITION OF A FEW CONSTANTS
//----------------------------------------------------

/* $Id$ */

AliAnalysisTaskDStarCorrelations *AddTaskDStarCorrelations(Int_t theMCon =5, Int_t mixing=4, Int_t trackselect =1)
{

  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    ::Error("AliAnalysisTaskDStarCorrelations", "No analysis manager to connect to.");
    return NULL;
  } 

  TFile* filecuts=new TFile("DStartoKpipiCuts_corr.root");
  if(!filecuts->IsOpen()){
    cout<<"Input file not found: exit"<<endl;
    return;
  }  
	  TFile* filecuts2=new TFile("AssocPartCuts.root");
	  if(!filecuts2->IsOpen()){
		  cout<<"Input file2 not found: exit"<<endl;
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

  //CREATE THE TASK
  printf("CREATE TASK \n");
  // create the task
  AliAnalysisTaskDStarCorrelations *task = new AliAnalysisTaskDStarCorrelations("AliAnalysisTaskDStarCorrelations",RDHFDStartoKpipi,corrCuts);
	
	// Setters

	task->SetMonteCarlo(theMCon);
	task->SetUseMixing(mixing);
	task->SetCorrelator(trackselect) ;
	//task->SetDebugLevel(0);
	

	if(trackselect == 1) Info("AliAnalysisTaskDStarCorrelations","Correlating D* with charged hadrons \n");
	else if(trackselect == 2) Info("AliAnalysisTaskDStarCorrelations","Correlating D* with charged kaons \n");
	else if(trackselect == 3) Info("AliAnalysisTaskDStarCorrelations","Correlating D* with reconstructed K0s \n");
	else Fatal("AliAnalysisTaskDStarCorrelations","Nothing to correlate with!");
	if(mixing) Info ("AliAnalysisTaskDStarCorrelations","Event Mixing Analysis\n");
	if(!mixing) Info ("AliAnalysisTaskDStarCorrelations","Single Event Analysis \n");

  // Create and connect containers for input/output
	//TString dcavalue = " ";
	if(!theMCon) TString contname = "Data";
	if(theMCon) TString contname = "MonteCarlo";
	if(trackselect ==1) TString particle = "Hadron";
	if(trackselect ==2) TString particle = "Kaon";
	if(trackselect ==3) TString particle = "KZero";
	
	TString cutname = "cuts" ;
	TString cutname2 = "hadroncuts" ;
   TString outputfile = AliAnalysisManager::GetCommonFileName();
	TString counter = "NormCounter";
   outputfile += ":PWGHF_D2H_";
	if(!mixing) {
		outputfile += "SE";
		contname += "SE";
		cutname += "SE";
		cutname2 += "SE";
		counter+= "SE";
	}
	if(mixing){
		outputfile += "ME";
		contname += "ME";
		cutname += "ME";
		cutname2 += "ME";
		counter+= "ME";
	}
	outputfile += "Dphi_DStar";
	outputfile += particle;
	cutname += particle;
	cutname2 += particle;
	contname += particle;
	counter+= particle;

	
	//cout << "Contname = " << contname << endl;
  mgr->AddTask(task);
  // ------ input data ------
  AliAnalysisDataContainer *cinput0  = mgr->GetCommonInputContainer();
  
  // ----- output data -----
  
  // output TH1I for event counting
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname, TList::Class(),AliAnalysisManager::kOutputContainer,outputfile.Data());
  AliAnalysisDataContainer *coutput2 = mgr->CreateContainer(cutname,AliRDHFCutsDStartoKpipi::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data()); //cuts
  AliAnalysisDataContainer *coutputDstarNorm = mgr->CreateContainer(counter,AliNormalizationCounter::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data());
  AliAnalysisDataContainer *coutput4 = mgr->CreateContainer(cutname2,AliHFAssociatedTrackCuts::Class(),AliAnalysisManager::kOutputContainer, outputfile.Data()); //cuts
  
  mgr->ConnectInput(task,0,mgr->GetCommonInputContainer());
  mgr->ConnectOutput(task,1,coutput1);
  mgr->ConnectOutput(task,2,coutput2);
  mgr->ConnectOutput(task,3,coutputDstarNorm);
  mgr->ConnectOutput(task,4,coutput4);

  return task ;

}

