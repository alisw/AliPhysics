void AddD2HTrain() {
  // 
  // Task of the D2H subgroup of PWG3 that can run in the Official Train
  //
  // They all use AOD+AOD.VertexingHF as input. 
  // They need to read the cuts from the macro 
  // $ALICE_ROOT/PWG3/vertexingHF/ConfigVertexingHF.C (trunk version).
  // This macro is loaded in the Init(), so it should be enough to have it 
  // in the local dir when the train is launched with the alien plugin.
  // They all produce only histos in the output. 
  //
  // andrea.dainese@pd.infn.it
  //

  TString taskName="",loadMacroPath="$ALICE_ROOT/PWG3/vertexingHF/";
  
  //taskName="AddTaskCompareHF.C"; taskName.Prepend(loadMacroPath.Data());
  //gROOT->LoadMacro(taskName.Data());
  //AliAnalysisTaskSECompareHF *cmpTask = AddTaskCompareHF();
  
  taskName="AddTaskD0Mass.C"; taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  AliAnalysisTaskSED0Mass *d0massTask = AddTaskD0Mass();
  AliAnalysisTaskSED0Mass *d0massLikeSignTask = AddTaskD0Mass(1); 
  
  taskName="AddTaskDplus.C"; taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  AliAnalysisTaskSEDplus *dplusTask = AddTaskDplus();
  
  //taskName="AddTaskSelectHF.C"; taskName.Prepend(loadMacroPath.Data());
  //gROOT->LoadMacro(taskName.Data());
  //AliAnalysisTaskSESelectHF *seleTask = AddTaskSelectHF();
  
  taskName="AddTaskBkgLikeSignD0.C"; taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  AliAnalysisTaskSEBkgLikeSignD0 *lsD0Task = AddTaskBkgLikeSignD0();
  
  taskName="AddTaskBkgLikeSignJPSI.C"; taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  AliAnalysisTaskSEBkgLikeSignJPSI *lsJPSITask = AddTaskBkgLikeSignJPSI();
  
  //taskName="AddTaskBtoJPSItoEle.C"; taskName.Prepend(loadMacroPath.Data());
  //gROOT->LoadMacro(taskName.Data());
  //AliAnalysisTaskSEBtoJPSItoEle *jpsiTask = AddTaskBtoJPSItoEle();
  
  taskName="AddTaskCFMultiVarMultiStep.C"; taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  AliCFHeavyFlavourTaskMultiVarMultiStep *cfmvmsTask = AddTaskCFMultiVarMultiStep();
  
  taskName="AddTaskCharmFraction.C"; taskName.Prepend(loadMacroPath.Data());
  gROOT->LoadMacro(taskName.Data());
  Int_t switchMC[5]={1,1,1,1,1};
  AliAnalysisTaskSECharmFraction *cFractTask = AddTaskCharmFraction("d0D0.root",switchMC);
  

  return;
}
