void RunReducedAnalysis(TString outputdir="", TString inputfilename="", Int_t run=-1, Int_t howMany=10000000, Int_t offset=0, TString chunk="") 
{

  TString outputfiledir = outputdir+chunk;

  gROOT->LoadMacro("$ALICE_PHYSICS/PWGDQ/reducedTree/macros/RunReducedEventAnalysis.C+");
  if(run==-1) return;
  //gROOT->LoadMacro(Form("%s/macros/AddTask_reducedep.C",outputdir.Data()));
  //AliReducedAnalysisTaskSE* task = AddTask_reducedep();
  gROOT->LoadMacro("$ALICE_PHYSICS/PWGDQ/reducedTree/macros/AddTask_iarsene_testTask.C");
  AliReducedAnalysisTaskSE* task = AddTask_iarsene_testTask(kFALSE)->GetReducedTask();
  RunReducedEventAnalysis(task,outputfiledir,inputfilename,run,howMany,offset);
}
