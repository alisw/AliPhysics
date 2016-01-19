AliAnalysisTaskEmcalHFeJetCorrel* AddTaskEmcalHFCJ(
  const char *ntracks            = "Tracks",
  const char *nclusters          = "CaloClusters",
  const char *njets              = "Jets",
  const char *nrho               = "Rho",
  Int_t       nCentBins          = 1,
  Double_t    jetradius          = 0.4,
  Double_t    jetptcut           = 0,
  Double_t    jetareacut         = 0.6,
  const char *type               = "EMCAL",
  Int_t       leadhadtype        = 0,
  const char *taskname           = "AliAnalysisTaskEmcalHFeJetCorrel",
  TString cutfile				 ="HFCJCuts.root",
  UInt_t triggerMask			 =-1,/*AliVEvent::kEMC1 | AliVEvent::kEMC7 | AliVEvent::kEMC8, kMB kEMC7 (kEMC8) kEMCEJE kEMCEGA*/
  Bool_t isMC					 = kFALSE,
  bool TestContainer			 = kFALSE,
  bool EGA1                      = kFALSE,
  const char *MCtracks            = "",
  const char *MCclusters          = "",
  const char *MCjets              = ""
)
{  
  // Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr)
  {
    ::Error("AddTaskEmcalHFCJ", "No analysis manager to connect to.");
    return NULL;
  }  
  
  // Check the analysis type using the event handlers connected to the analysis manager.
  //==============================================================================
  if (!mgr->GetInputEventHandler())
  {
    ::Error("AddTaskEmcalHFCJ", "This task requires an input event handler");
    return NULL;
  }
  
  //-------------------------------------------------------
  // Init the task and do settings
  //-------------------------------------------------------

  TString name(taskname);
  if (strcmp(njets,"")) {
    name += "_";
    name += njets;
  }
  if (strcmp(nrho,"")) {
    name += "_";
    name += nrho;
  }
  if (strcmp(type,"")) {
    name += "_";
    name += type;
  }

  Printf("name: %s",name.Data());

  AliAnalysisTaskEmcalHFeJetCorrel* jetTask = new AliAnalysisTaskEmcalHFeJetCorrel(name);
  jetTask->SetCentRange(0.,100.);
  jetTask->SetNCentBins(nCentBins);
  jetTask->SetReadMC(isMC);
  jetTask->SetMCParticles(MCtracks);
  jetTask->SetQA(kTRUE);//FALSE);
  jetTask->SetCheckClusterMatching(kFALSE);
  jetTask->SetDoAnalysis(kTRUE);//FALSE);
  jetTask->SetGeneralSpectra(kTRUE);//FALSE);
  jetTask->SetDetectorsTest(kTRUE);//FALSE);
  //Shingo Electron Selection
  jetTask->SetFilterBitElectron(AliAODTrack::kTrkGlobalNoDCA);

  //Defaut parameters
  jetTask->SetMinPtElectron(0.1);
  jetTask->SetEtaTracks(-0.9,0.9);
  jetTask->SetNsigmaTPCelectron(-10,3.5);
  jetTask->SetPhotonicMassCut(0.5);
  jetTask->SetMinNumberOfClusterCells(4);
  jetTask->SetMaxM20(0.9);
  jetTask->SetMaxM02(0.9);
  jetTask->SetEoverPlimits(0.0, 3.0);
  
  AliParticleContainer *trackCont  = jetTask->AddParticleContainer(ntracks);
  trackCont->SetClassName("AliVTrack");
  AliClusterContainer *clusterCont = jetTask->AddClusterContainer(nclusters);

  TString strType(type);
  AliJetContainer *jetCont = jetTask->AddJetContainer(njets,strType,jetradius);
  if(jetCont) {
    jetCont->SetRhoName(nrho);
    jetCont->ConnectParticleContainer(trackCont);
    jetCont->ConnectClusterContainer(clusterCont);
    jetCont->SetZLeadingCut(0.98,0.98);
    jetCont->SetPercAreaCut(0.6);
    jetCont->SetJetPtCut(jetptcut);    
    jetCont->SetLeadingHadronType(leadhadtype);
  }


//=========================CUTS=========================
AliRDHFJetsCuts *cuts;

bool kFileExists=kFALSE;
//if(!gSystem->AccessPathName(cutfile.Data(),kFileExists)){
  
if(gSystem->AccessPathName(cutfile.Data(),kFileExists))
{
	::Error("\n==CutObject not Defined. Verify your .root file==\n");
	return NULL;
}
	TFile *f=TFile::Open(cutfile.Data());
	//cuts= (AliRDHFCutsD0toKpi*)f->Get("EventTrackCuts");
	cuts= (AliRDHFJetsCuts*)f->Get("HFCJCuts");

	cout<<"\n==========================================\n Cutfile used:\n"<<cutfile.Data()<<endl;
	//cuts->PrintAll();
    
    if(triggerMask>0)
    {
        cuts->SetTriggerMask(triggerMask);
        if( (triggerMask == AliVEvent::kEMCEGA)&&(EGA1) )cuts->SetTriggerClass("EG1");
    }

jetTask->SetJetCuts(cuts);
delete cuts;
//========================================================== 
  
//-------------------------------------------------------
// Final settings, pass to manager and set the containers
//-------------------------------------------------------
  
  mgr->AddTask(jetTask);

  // Create containers for input/output
  AliAnalysisDataContainer *cinput1  = mgr->GetCommonInputContainer()  ;
  TString contname(name);
  contname += "_histos";
  AliAnalysisDataContainer *coutput1 = mgr->CreateContainer(contname.Data(), 
							    TList::Class(),AliAnalysisManager::kOutputContainer,
							    Form("%s", AliAnalysisManager::GetCommonFileName()));
  mgr->ConnectInput  (jetTask, 0,  cinput1 );
  mgr->ConnectOutput (jetTask, 1, coutput1 );
  
  return jetTask;
}
