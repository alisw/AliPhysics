//=============================================================================
//
// *** AliAnalysisTaskParticleEff.C ***
//
// This macro initialize an AnalysisTask for Particle Efficiency
//
//=============================================================================

AliAnalysisTaskEffK0ss *AddTaskEffK0s(TString containerName="femtolist",int method=3, int filterbit=128, double nsigmaDaught=5, double cosPA=0.998, double window=0.015, double dcaDaughtersMax=0.7, double dcaPVMin =0.1, bool electronRejection = kTRUE, double nsigmaErej=3.0, double minDCApos =0.1, double minDCAneg=0.1, double maxctau=26.8,  double minV0rad=10, bool ownDCA=kTRUE,  float dcaxy=0.3,  float dcaz=1)
{
  // A. Get the pointer to the existing analysis manager via the static access method.
  //==============================================================================
  AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
  if (!mgr) {
    Error("AddTaskEffDY", "No analysis manager to connect to.");
    return NULL;
  }

  // B. Check the analysis type using the event handlers connected to the analysis
  //    manager. The availability of MC handler can also be checked here.
  //==============================================================================
  if (!mgr->GetInputEventHandler()) {
    ::Error("AddTaskEffDY", "This task requires an input event handler");
    return NULL;
  }
  TString type = mgr->GetInputEventHandler()->GetDataType(); // can be "ESD" or "AOD"
  cout << "Found " <<type << " event handler" << endl;

  // C. Create the task, add it to manager.
  //===========================================================================
//  gSystem->SetIncludePath("-I$ROOTSYS/include  -I./PWG2AOD/AOD -I./PWG2femtoscopy/FEMTOSCOPY/AliFemto -I./PWG2femtoscopyUser/FEMTOSCOPY/AliFemtoUser -I$ALICE_ROOT/include");

//   if (TProofMgr::GetListOfManagers()->GetEntries()) {
// //     if (dynamic_cast<TProofLite *> gProof) {
// //       char *macrocommand[10000];
// //       sprintf(macrocommand, ".L %s", configMacroName);
// //       gProof->Exec(macrocommand);
// //     }
// //     else
//     gProof->Load(configMacroName);
//   }
  //  gROOT->LoadMacro("ConfigFemtoAnalysis.C++");

  AliAnalysisTaskEffK0ss *taskEffK0s = new AliAnalysisTaskEffK0ss("EffTaskK0s", method, filterbit, nsigmaDaught, cosPA, window, dcaDaughtersMax, dcaPVMin, electronRejection, nsigmaErej, minDCApos, minDCAneg, maxctau, minV0rad, ownDCA, dcaxy, dcaz);
  taskEffK0s->SetPidMethod(method);
  taskEffK0s->SetFB(filterbit);
  taskEffK0s->SetNsigmaDaughters(nsigmaDaught);
  taskEffK0s->SetMinCosPointingAngle(cosPA);
  taskEffK0s->SetInvMassWindow(window);
  taskEffK0s->SetMaxDCADaughters(dcaDaughtersMax);
  taskEffK0s->SetMinDCAToPrimVtx(dcaPVMin);
  taskEffK0s->SetElectronRejection(electronRejection);
  taskEffK0s->SetNsigmaElectronRejection(nsigmaErej);
  taskEffK0s->SetMinDcaPosToPrimVertex(minDCApos);
  taskEffK0s->SetMinDcaNegToPrimVertex(minDCAneg);
  taskEffK0s->SetMaxCTauK0s(maxctau);
  taskEffK0s->SetMinV0Radius(minV0rad);
  taskEffK0s->SetUseDCAcuts(ownDCA);
  taskEffK0s->SetUseDCAcutsxy(dcaxy);
  taskEffK0s->SetUseDCAcutsz(dcaz);
  mgr->AddTask(taskEffK0s);

  // D. Configure the analysis task. Extra parameters can be used via optional
  // arguments of the AddTaskXXX() function.
  //===========================================================================

  // E. Create ONLY the output containers for the data produced by the task.
  // Get and connect other common input/output containers via the manager as below
  //==============================================================================
  TString outputfile = AliAnalysisManager::GetCommonFileName();
  outputfile += ":PWG2FEMTO";
  AliAnalysisDataContainer *cout_femto  = mgr->CreateContainer(containerName,  TList::Class(),
  							       AliAnalysisManager::kOutputContainer,outputfile);


   mgr->ConnectInput(taskEffK0s, 0, mgr->GetCommonInputContainer());
   mgr->ConnectOutput(taskEffK0s, 1, cout_femto);

   // std::ofstream ofile;
   // ofile.open("test.txt", std::ofstream::app);
   // ofile<<"AddTask:" <<method<<" "<<taskEfficiency->GetPidMethod()<<std::endl;								 
   // ofile.close();
   
   // Return task pointer at the end
   return taskEffK0s;
}
