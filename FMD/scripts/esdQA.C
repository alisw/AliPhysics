//______________________________________________________________________
void 
esdQA() 
{  
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libAnalysisCheck.so");

  // create the analysis goodies object
  AliAnalysisGoodies * master     = new AliAnalysisGoodies(); 
  AliAnalysisTask*     tasks[2]   = { 0, 0 };
  TClass*              inputs[2]  = { TChain::Class(),    0 };
  TClass*              outputs[2] = { TObjArray::Class(), 0 };
  tasks[0]                        = new AliFMDQATask("FMD");
  master->SetTasks(1, tasks, inputs, outputs) ; 

  TChain* chain = new TChain("esdTree") ;
  chain->AddFile("AliESDs.root");
  master->Process(chain) ; 
  return ;
}

