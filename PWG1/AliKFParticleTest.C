
//---------------------------------------------------------------------------------
// The example of usage of AliKFParticle & AliKFVertex classes for V0 analysis
// .
// @author  S.Gorbunov, I.Kisel
// @version 1.0
// @since   13.05.07
// 
// The AliKFParticleTest macro contains a toy V0 finder for ESD tracks.
// At the first step, event primary vertex is reconstructed. 
// At the second step, ideal PID hypothesis are assigned to all the particles.
// At the third step, V0 candidates are constructed for each pair 
// of positive and negative particles.
// V0 candidate considered as good when it passes Chi^2 cut and 
// it is placed >= 3 Sigma away from the primary vertex.
// Invariant mass distribution for all good V0 candidates is plotted.
//
//  -= Copyright &copy ALICE HLT Group =-
//_________________________________________________________________________________


void AliKFParticleTest(const char *dirname="/d/alice09/sma/v4-05-Rev-03/charmwmi/", Int_t limit = 100) 
{
  UInt_t startsample = 1000;
  cout <<"Using events : "<<dirname<<endl;
  
  gSystem->AddIncludePath("-I\"$ALICE_ROOT/include\"");   
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISRL.so");
  gROOT->LoadMacro("$ALICE_ROOT/STEER/AliKFParticleTest.h+");
  
  //
  // Setup chain
  //
  UInt_t nFile = 0;
  UInt_t iFile = startsample-1;
  TChain *chain = new TChain("esdTree");
  
  while (nFile < limit)
    {
      iFile++;       
      TString pathFile (Form("%s%3.3d/AliESDs.root", dirname, iFile));
      TFile * file = TFile::Open (pathFile.Data());
      if (file == 0x0)continue;
      
      cout<<"File "<<pathFile.Data()<<endl;
      
      if (file->IsZombie ()){
	cout<<"File "<<pathFile.Data()<<" is zombie"<<endl;
	file->Close();
	continue;
      }
      
      file->Close(); 
      chain->Add (pathFile.Data()); 	
      nFile++;       
    }
  cout << "Number of files : "<<nFile<<endl;
  
  
  AliAnalysisManager *mgr = new AliAnalysisManager("testEvent");
  
  AliAnalysisTask *task = new AliKFParticleTest("AliKFParticleTest");
  
  // Create containers for input/output
  AliAnalysisDataContainer *cinput = mgr->CreateContainer("input0", TChain::Class(), AliAnalysisManager::kInputContainer);
  
  // Connect containers to the task input/outputs

  cout << "Adding task " << task->GetName() << endl;
  
  mgr->ConnectInput(task, 0, cinput);
  
  AliAnalysisDataContainer *coutput =
    mgr->CreateContainer("output0", TObjArray::Class(), AliAnalysisManager::kOutputContainer,Form("%s.root",task->GetName()));
  mgr->ConnectOutput(task,0,coutput);
  
  // Init analysis and start event loop
  if (mgr->InitAnalysis()) {
    mgr->PrintStatus();
    mgr->StartAnalysis("local",chain);
  }
  delete mgr;   
}
