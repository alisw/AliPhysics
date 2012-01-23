Int_t AddD2HTrain(Bool_t readMC=kTRUE,
                  Int_t addQA=1,
		  Int_t addD0Mass=1,
		  Int_t addD0MassLS=1,
		  Int_t addDplus=1,
		  Int_t addLSD0=0,
		  Int_t addCFD0=0,
		  Int_t addPromptD0=1,
		  Int_t addDs=0,
		  Int_t addDStar=1,
		  Int_t addDStarJets=0,
		  Int_t addCFDStar=0) {
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

  TString taskName="",loadMacroPath="$ALICE_ROOT/PWG3/vertexingHF/macros/";
  Int_t ntasks=0;

  //taskName="AddTaskCompareHF.C"; taskName.Prepend(loadMacroPath.Data());
  //gROOT->LoadMacro(taskName.Data());
  //AliAnalysisTaskSECompareHF *cmpTask = AddTaskCompareHF();
  
  if(addQA) {
    taskName="AddTaskHFQA.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliAnalysisTaskSEHFQA *taskQAHF = AddTaskHFQA(AliAnalysisTaskSEHFQA::kD0toKpi,"",readMC,kTRUE);
    ntasks++;
  }


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

  if(addCFD0 && readMC) {
    taskName="AddTaskCFMultiVarMultiStep.C"; taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    AliCFHeavyFlavourTaskMultiVarMultiStep *cfmvmsTask = AddTaskCFMultiVarMultiStep();
    ntasks++;
  }

  if(addPromptD0) {
    taskName="AddTaskSECharmFraction.C"; 
    taskName.Prepend(loadMacroPath.Data());
    gROOT->LoadMacro(taskName.Data());
    Int_t switchMC[5]={0,0,0,0,0};
    Int_t ppPbPb=1;// 0 for pp, 1 for PbPb, used to siwtch on/off the removal of daughters from the primary vertex
    AliAnalysisTaskSECharmFraction *cFractTask = AddTaskSECharmFraction("standard",switchMC,readMC,kTRUE,kFALSE,"D0toKpiCharmFractCuts.root","c",ppPbPb);
    // arguments: filename,switchMC,readmc,usepid,likesign,cutfilename,containerprefix
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
