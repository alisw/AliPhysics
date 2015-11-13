//void RunReducedAnalysis()
void RunReducedAnalysis(TString outputdir="", TString inputfilename="", Int_t run=-1, Int_t howMany=10000000, Int_t offset=0, TString chunk="") 
{


  TString codedir=outputdir;
  TString outputfiledir = "";
  if(chunk!=""){
    codedir+="/code/";
    outputfiledir=outputdir+chunk;
  }


  gROOT->ProcessLine(Form(".x %s/macros/loadClasses.C(\"%s\")",codedir.Data(),codedir.Data()));
  gROOT->LoadMacro(Form("%s/macros/RunReducedEventAnalysis.C+",codedir.Data()));
  if(run==-1) return;
  //gROOT->LoadMacro(Form("%s/macros/AddTask_reducedep.C",codedir.Data()));
  //AliReducedAnalysisTaskSE* task = AddTask_reducedep();
  gROOT->LoadMacro(Form("%s/macros/AddTask_iarsene_testTask.C",codedir.Data()));
  AliReducedAnalysisTaskSE* task = AddTask_iarsene_testTask(kFALSE);
  RunReducedEventAnalysis(task,outputfiledir,inputfilename,run,howMany,offset);


}
