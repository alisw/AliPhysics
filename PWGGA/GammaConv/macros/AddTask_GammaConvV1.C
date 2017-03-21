
void AddTask_GammaConvV1(TString mode = "test", Bool_t isMC = kFALSE){

	gSystem->Load("libCore");
	gSystem->Load("libTree");
	gSystem->Load("libGeom");
	gSystem->Load("libVMC");
	gSystem->Load("libPhysics");
	gSystem->Load("libMinuit");
	gSystem->Load("libSTEERBase");
	gSystem->Load("libESD");
	gSystem->Load("libAOD");
	gSystem->Load("libANALYSIS");
	gSystem->Load("libANALYSISalice");  
	gSystem->Load("libPWGGAGammaConv");
	gSystem->Load("libEve");
	gSystem->Load("libCDB");
	gSystem->Load("libProof");
	gSystem->Load("libRAWDatabase");
	gSystem->Load("libSTEER");
	gSystem->Load("libSTEERBase");
	gSystem->Load("libTender");
	gSystem->Load("libTRDbase");
	gSystem->Load("libVZERObase");
	gSystem->Load("libVZEROrec");
	gSystem->Load("libTenderSupplies");

	gSystem->Load("libGui.so");
	gSystem->Load("libXMLParser.so");
	gSystem->Load("libMinuit2.so");
	gSystem->Load("libOADB.so");
	gSystem->Load("libCORRFW.so");
	gSystem->Load("libPWGflowBase.so");
	gSystem->Load("libPWGflowTasks.so");
	gSystem->Load("libPWGGAGammaConv.so");
   
   Bool_t readTR = kTRUE;


   TString Pattern;
   TString Prefix;
   if(isMC){
      Pattern="*/AliESDs.root";
      Prefix ="";
   } else {
      Pattern="pass2/*/AliESDs.root";
      Prefix ="000";
   }

   Int_t run_numbers[1] = {168464};
   gSystem->Setenv("alien_CLOSE_SE","ALICE::Grenoble::SE");

   AliAnalysisManager *mgr = new AliAnalysisManager("Analysis");
   if (!mgr) {
      Error("AddTask_GammaConvV1", "No analysis manager found.");
      return 0;
   }

   gROOT->LoadMacro("CreateAlienHandlerNonPTTest.C");
   AliAnalysisAlien *alienHandler = CreateAlienHandler(mode);
   if (!alienHandler) return;

   TString WorkingDir = "";
   TString data;
   TString period;
   if(isMC){
     data = "/alice/sim/2014/LHC14a1a/";
     period = "LHC14a1a";
   } else {
     data = "/alice/data/2011/LHC11h_2/";
     period = "LHC11h";
   }
   
   AliESDInputHandler * handler;
   handler = new AliESDInputHandler(); mgr->SetInputEventHandler((AliESDInputHandler*)handler);

   AliMCEventHandler* mcHandler;
   if(isMC){
      mcHandler = new AliMCEventHandler();
      mgr->SetMCtruthEventHandler(mcHandler);
      mcHandler->SetReadTR(readTR);
   }

   //========= Add PID Reponse to ANALYSIS manager ====
   if(!(AliPIDResponse*)mgr->GetTask("PIDResponseTask")){
      gROOT->LoadMacro("$ALICE_ROOT/ANALYSIS/macros/AddTaskPIDResponse.C");
	  if(period.CompareTo("LHC11h_2"))
		AddTaskPIDResponse(isMC,kTRUE,isMC,2,kFALSE,"TPC:$ALICE_PHYSICS/OADB/COMMON/PID/data/TPCPIDResponse_special.root;TPC-Maps:$ALICE_PHYSICS/OADB/COMMON/PID/data/TPCetaMaps_special.root");
	  else
		AddTaskPIDResponse(isMC);
   }//========= Add PID Reponse to ANALYSIS manager ====
   
   gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskPhysicsSelection.C");
   AliPhysicsSelectionTask* physSelTask = AddTaskPhysicsSelection(isMC);

   AliCentralitySelectionTask *taskCentrality;
   gROOT->LoadMacro("$ALICE_PHYSICS/OADB/macros/AddTaskCentrality.C");
   taskCentrality = AddTaskCentrality();
   if (isMC){
      taskCentrality->SetMCInput();
      Info("AddTask_tender_CentralitySelection", "This task has MC.");
   }
   
   gROOT->LoadMacro("$ALICE_PHYSICS/PWGGA/GammaConv/macros/AddTask_GammaConvV1_PbPb.C");
   
   cout<<"Connecting to Alien..."<<endl;
   TGrid::Connect("alien://");
   cout<<"==============================="<<endl;
   cout<<"Successfully connected to Alien"<<endl;
   cout<<"==============================="<<endl;

   AliAnalysisTask *taskA = AddTask_GammaConvV1_PbPb(200,  //change different set of cuts
                                                     isMC,
                                                     kTRUE, //enable QA in AliAnalysisTaskGammaConvV1
                                                     kTRUE, // enable additional QA task
                                                     "alien:///alice/cern.ch/user/f/fbock/MCSpectraInput.root", // path to file for weigting input
                                                     0,  // 1 pi0 header, 2 eta header, 3 both (only for "named" boxes)
                                                     "1000000060084000001500000",
                                                     "LHC11h_2",  //name of the period for added signals and weighting
                                                     kFALSE,  //enable Weighting
                                                     kFALSE,  //use THnSparse
                                                     0); 
   if (!mgr->InitAnalysis())return;
   mgr->StartAnalysis("grid");
}
Int_t CheckLoadLibrary(const char* library)                                                                                                        
{                                                                                                                                                  
   // checks if a library is already loaded, if not loads the library                                                                               
                                                                                                                                                   
   if (strlen(gSystem->GetLibraries(library, "", kFALSE)) > 0)                                                                       
      return 1;                                                                                                                                      
                                                                                                                                                   
   return gSystem->Load(library);                                                                                                                   
}
