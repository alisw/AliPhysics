Int_t AddD2HTrain(Bool_t readMC=kTRUE,
		  Int_t addD0Mass=1,
		  Int_t addD0MassLS=1,
		  Int_t addDplus=1,
		  Int_t addLSD0=1,
		  Int_t addLSJpsi=1,
		  Int_t addCFD0=1,
		  Int_t addPromptD0=1,
		  Int_t addDs=1,
		  Int_t addDStar=1,
		  Int_t addDStarJets=1,
		  Int_t addCFDStar=1) {
  // 
  // Tasks of the D2H subgroup of PWG3 that can run in the Official Train
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
  Int_t ntasks=0;

  //taskName="AddTaskCompareHF.C"; taskName.Prepend(loadMacroPath.Data());
  //gROOT->LoadMacro(taskName.Data());
  //AliAnalysisTaskSECompareHF *cmpTask = AddTaskCompareHF();
  
  if(addD0Mass || addD0MassLS) {
    taskName="AddTaskD0Mass.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    if(addD0Mass) {
      AliAnalysisTaskSED0Mass *d0massTask = AddTaskD0Mass(0,readMC);
      ntasks++;
    }
    if(addD0MassLS) {
      AliAnalysisTaskSED0Mass *d0massLikeSignTask = AddTaskD0Mass(1,readMC); 
      ntasks++;
    }
  }

  if(addDplus) {
    taskName="AddTaskDplus.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSEDplus *dplusTask = AddTaskDplus(kFALSE,readMC);
    ntasks++;
  }  

  //taskName="AddTaskSelectHF.C"; taskName.Prepend(loadMacroPath.Data());
  //gROOT->LoadMacro(taskName.Data());
  //AliAnalysisTaskSESelectHF *seleTask = AddTaskSelectHF();
  
  if(addLSD0 && readMC) {
    taskName="AddTaskBkgLikeSignD0.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSEBkgLikeSignD0 *lsD0Task = AddTaskBkgLikeSignD0();
    ntasks++;
  }

  if(addLSJpsi && readMC) {
    taskName="AddTaskBkgLikeSignJPSI.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSEBkgLikeSignJPSI *lsJPSITask = AddTaskBkgLikeSignJPSI();
    ntasks++;
  }

  //taskName="AddTaskBtoJPSItoEle.C"; taskName.Prepend(loadMacroPath.Data());
  //gROOT->LoadMacro(taskName.Data());
  //AliAnalysisTaskSEBtoJPSItoEle *jpsiTask = AddTaskBtoJPSItoEle();
 
  if(addCFD0 && readMC) {
    taskName="AddTaskCFMultiVarMultiStep.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliCFHeavyFlavourTaskMultiVarMultiStep *cfmvmsTask = AddTaskCFMultiVarMultiStep();
    ntasks++;
  }

  if(addPromptD0) {
    taskName="AddTaskCharmFraction.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    Int_t switchMC[5]={1,1,1,1,1};
    AliAnalysisTaskSECharmFraction *cFractTask = AddTaskCharmFraction("d0D0.root",switchMC,readMC);
    ntasks++;
  }

  if(addDs) {
    taskName="AddTaskDs.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSEDs *dsTask = AddTaskDs(readMC);
    ntasks++;
  }  

  if(addDStar) {
    taskName="AddTaskDStarSpectra.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSEDStarSpectra *dstarTask = AddTaskDStarSpectra(readMC);
    ntasks++;
  }  

  if(addDStarJets) {
    taskName="AddTaskDStarJets.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSEDStarJets *dstarjetsTask = AddTaskDStarJets(readMC);
    ntasks++;
  }  

  if(addCFDStar && readMC) {
    taskName="AddTaskCFDStar.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliCFTaskForDStarAnalysis *cfDstarTask = AddTaskCFDStar();
    ntasks++;
  }


  return ntasks;
}
