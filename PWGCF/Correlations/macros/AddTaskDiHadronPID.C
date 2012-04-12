// AddTask Macro (v 8.00).
// Updated: Apr 12th 2012.
// Author: Misha Veldhoen (m.veldhoen@cern.ch)

AliAnalysisTaskDiHadronPID *AddTaskDiHadronPID(Int_t verbose = 1,
                                               Bool_t printbuffersize = kFALSE,
                                               Bool_t mixedevents = kTRUE,
                                               TString beamtype = "PbPb",
                                               Double_t MaxEta = 0.8,
                                               Double_t MaxPlotEta = 0.8,
                                               Double_t MaxPt = 10.,
                                               Int_t NEtaBins = 36,
                                               Int_t NPhiBins = 36,
                                               Double_t VertexZMixedEvents = 2.,
                                               Bool_t zoomed = kFALSE,
                                               Bool_t DoITSCut = kFALSE,
                                               Bool_t DoDCACut = kTRUE,
                                               Bool_t DemandNoMismatch = kTRUE,
                                               Int_t trigbuffermaxsize=2500,
                                               Double_t centralitycutmax=0.,
                                               Double_t centralitycutmin=10.)

                                               

{
    	
	// Get the current analysis manager.
    AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
    if (!mgr) {
        Error("AddTaskDiHadronPID.C", "No analysis manager found.");
        return 0;
    }
    
    // Create an instance of the task.
    AliAnalysisTaskDiHadronPID *task = new AliAnalysisTaskDiHadronPID("DiHadronPID");
    
    // Configure the task.
    task->SetVerbose(verbose);
    task->SetPrintBufferSize(printbuffersize);
    task->SetCalculateMixedEvents(mixedevents);
    task->SetBeamType(beamtype);
    task->SetMaxEta(MaxEta);
    task->SetMaxPlotEta(MaxPlotEta);
    task->SetMaxPt(MaxPt);
    task->SetNEtaBins(NEtaBins);
    task->SetNPhiBins(NPhiBins);
    task->SetVertexZMixedEvents(VertexZMixedEvents);
    task->SetZoomed(zoomed);
    task->SetDoITSCut(DoITSCut);
    task->SetDoDCACut(DoDCACut);
    task->SetDemandNoMismatch(DemandNoMismatch);
    task->SetTrigBufferMaxSize(trigbuffermaxsize);
    task->SetCentralityCut(centralitycutmax,centralitycutmin);

    // Add the task.
	mgr->AddTask(task);
    
	// Data containers.
	AliAnalysisDataContainer *cinput  = mgr->GetCommonInputContainer();
	mgr->ConnectInput(task, 0, cinput); 
	
	AliAnalysisDataContainer *coutput1 = 
    mgr->CreateContainer("DiHadronPID", TList::Class(),
                         AliAnalysisManager::kOutputContainer,"DiHadronPID.root");
	
	mgr->ConnectOutput (task,  1, coutput1);
	
	return task;
	
}


