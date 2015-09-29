
void AddTask_GammaConvV1(TString mode){

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
   
   Bool_t isMC    = kFALSE;
   Bool_t readTR = kTRUE;

   TString Pattern="pass2/*/AliESDs.root";
   //TString Pattern="*/AliAOD.root";
   TString Prefix ="000";

//   Int_t run_numbers[100] = {137161, 137162, 137165, 137230, 137231, 137232, 137235, 137236, 137243, 137366, 137430, 137431, 137432, 137434, 137439, 137440, 137441, 137443, 137530, 137531, 137539, 137541, 137544, 137546, 137549, 137595, 137608, 137638, 137639, 137685, 137686, 137691, 137692, 137693, 137704, 137718, 137722, 137724, 137751, 137752, 137844, 137848, 138190, 138192, 138197, 138201, 138225, 138275, 138364, 138396, 138438, 138439, 138442, 138469, 138534, 138578, 138579, 138582, 138583, 138621, 138624, 138638, 138652, 138653, 138662, 138666, 138730, 138732, 138837, 138870, 138871, 138872, 139028, 139029, 139036, 139037, 139038, 139104, 139105, 139107, 139173, 139308, 139309, 139310, 139311, 139314, 139328, 139329, 139360, 139437, 139438, 139439, 139465, 139503, 139504, 139505, 139507, 139510};
   Int_t run_numbers[1] = {168464}; //, 168310, 169094, 169411, 169415, 169417, 169584, 169965, 170208, 170270};
   //Int_t run_numbers[10] = {137161, 137162, 137165, 137230, 137231, 137232, 137235, 137236, 137243, 137366};
   //Int_t run_numbers[10] = {137430, 137431, 137432, 137434, 137439, 137440, 137441, 137443, 137530, 137531};
   //Int_t run_numbers[10] = {137539, 137541, 137544, 137546, 137549, 137595, 137608, 137638, 137639, 137685};
   //Int_t run_numbers[10] = {137686, 137691, 137692, 137693, 137704, 137718, 137722, 137724, 137751, 137752};
   //Int_t run_numbers[10] = {137844, 137848, 138190, 138192, 138197, 138201, 138225, 138275, 138364, 138396};
   //Int_t run_numbers[10] = {138438, 138439, 138442, 138469, 138534, 138578, 138579, 138582, 138583, 138621};
   //Int_t run_numbers[10] = {138624, 138638, 138652, 138653, 138662, 138666, 138730, 138732, 138837, 138870};
   //Int_t run_numbers[10] = {138871, 138872, 139028, 139029, 139036, 139037, 139038, 139104, 139105, 139107};
   //Int_t run_numbers[10] = {139173, 139308, 139309, 139310, 139311, 139314, 139328, 139329, 139360, 139437};
   //Int_t run_numbers[10] = {139438, 139439, 139465, 139503, 139504, 139505, 139507, 139510};
   gSystem->Setenv("alien_CLOSE_SE","ALICE::Grenoble::SE");

   AliAnalysisManager *mgr = new AliAnalysisManager("Analysis");
   if (!mgr) {
      Error("AddTask_GammaConvV1", "No analysis manager found.");
      return 0;
   }

   gROOT->LoadMacro("CreateAlienHandlerNonPT.C");
   AliAnalysisAlien *alienHandler = CreateAlienHandler(mode);
   if (!alienHandler) return;

   TString WorkingDir = "";
   //TString data = "/alice/sim/2014/LHC14a1a/";
   TString data = "/alice/data/2011/LHC11h_2/";
   TString period = "LHC11h_2";
   
   mgr->SetGridHandler(alienHandler);
   alienHandler->SetGridWorkingDir(WorkingDir.Data());
   alienHandler->SetGridOutputDir(Form("OUTPUT_%s",WorkingDir.Data()));

   alienHandler->SetGridDataDir(data);
   alienHandler->SetDataPattern(Pattern.Data());
   alienHandler->SetRunPrefix(Prefix.Data());
   //alienHandler->SetRunRange(137161,139510);
   for (Int_t i=0; i<1; i++) {
      if (run_numbers[i]==0) break;
      alienHandler->AddRunNumber(run_numbers[i]);
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
													 0,
													 "alien:///alice/cern.ch/user/l/lleardin/InterpValuesAndFlattening.root",
													 0); 
//    AliAnalysisTask *taskB = AddTask_GammaConvV1_PbPb(194,  //change different set of cuts
//                                                      isMC,
//                                                      kTRUE, //enable QA in AliAnalysisTaskGammaConvV1
//                                                      kTRUE, // enable additional QA task
//                                                      "alien:///alice/cern.ch/user/f/fbock/MCSpectraInput.root", // path to file for weigting input
//                                                      0,  // 1 pi0 header, 2 eta header, 3 both (only for "named" boxes)
//                                                      "1000000060084000001500000",
//                                                      "LHC11h_2",  //name of the period for added signals and weighting
//                                                      kFALSE,  //enable Weighting
//                                                      kFALSE,  //use THnSparse
// 													 0,
// 													 "alien:///alice/cern.ch/user/l/lleardin/InterpValuesAndFlattening.root",
// 													 0); 
//    AliAnalysisTask *taskC = AddTask_GammaConvV1_PbPb(182,  //change different set of cuts
//                                                      isMC,
//                                                      kTRUE, //enable QA in AliAnalysisTaskGammaConvV1
//                                                      kTRUE, // enable additional QA task
//                                                      "alien:///alice/cern.ch/user/f/fbock/MCSpectraInput.root", // path to file for weigting input
//                                                      0,  // 1 pi0 header, 2 eta header, 3 both (only for "named" boxes)
//                                                      "1000000060084000001500000",
//                                                      "LHC11h_2",  //name of the period for added signals and weighting
//                                                      kFALSE,  //enable Weighting
//                                                      kFALSE,  //use THnSparse
// 													 0,
// 													 "alien:///alice/cern.ch/user/l/lleardin/CentralityFlatFile20April.root",
// 													 1); 
//    AliAnalysisTask *taskC = AddTask_GammaConvV1_PbPb(125,  //change different set of cuts
//                                                      isMC,
//                                                      kTRUE, //enable QA in AliAnalysisTaskGammaConvV1
//                                                      kTRUE, // enable additional QA task
//                                                      "alien:///alice/cern.ch/user/f/fbock/MCSpectraInput.root", // path to file for weigting input
//                                                      0,  // 1 pi0 header, 2 eta header, 3 both (only for "named" boxes)
//                                                      "1000000060084000001500000",
//                                                      "LHC11h_2",  //name of the period for added signals and weighting
//                                                      kFALSE,  //enable Weighting
//                                                      kFALSE);  //use THnSparse
//    AliAnalysisTask *taskD = AddTask_GammaConvV1_PbPb(133,  //change different set of cuts
//                                                      isMC,
//                                                      kTRUE, //enable QA in AliAnalysisTaskGammaConvV1
//                                                      kTRUE, // enable additional QA task
//                                                      "alien:///alice/cern.ch/user/f/fbock/MCSpectraInput.root", // path to file for weigting input
//                                                      0,  // 1 pi0 header, 2 eta header, 3 both (only for "named" boxes)
//                                                      "1000000060084000001500000",
//                                                      "LHC11h_2",  //name of the period for added signals and weighting
//                                                      kFALSE,  //enable Weighting
//                                                      kFALSE);  //use THnSparse
//    AliAnalysisTask *taskD = AddTask_GammaConvV1_PbPb(134,  //change different set of cuts
//                                                      isMC,
//                                                      kTRUE, //enable QA in AliAnalysisTaskGammaConvV1
//                                                      kTRUE, // enable additional QA task
//                                                      "alien:///alice/cern.ch/user/f/fbock/MCSpectraInput.root", // path to file for weigting input
//                                                      0,  // 1 pi0 header, 2 eta header, 3 both (only for "named" boxes)
//                                                      "1000000060084000001500000",
//                                                      "LHC11h_2",  //name of the period for added signals and weighting
//                                                      kFALSE,  //enable Weighting
//                                                      kFALSE);  //use THnSparse   

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
