/// \ingroup Macros
/// \file    QAtrendingfitExample.C
/// \brief   Demo usage of the information from the AliExternalInfo and AliNDFunctionInterface.cxx for MVA regression of TPC QA information (QA.EVS+QA calibration)
///
/// \author  Marian Ivanov

/*
  AliDrawStyle::SetDefaults()
  AliDrawStyle::ApplyStyle("figTemplate");
  gStyle->SetOptTitle(1);
  .L $AliRoot_SRC/STAT/Macros/AliNDFunctionInterface.cxx+
  .L $AliRoot_SRC/STAT/Macros/QAtrendingFitExample.C

  0.)  loadTree();
  1.)  cacheTree();      /// needed in case friend trees or arrays - which are not supported by TMVA
  3.)  RegisterFitters();
  4.)  makeMVAFits();
  5.)  loadMVAreaders();
  tree->Draw("AliNDFunctionInterface::EvalMVA(2,interactionRate, bz0,qmaxQASum,qmaxQASumR):resolutionMIP:interactionRate","run==QA.EVS.run&&meanMIP>40&&bz0<-0.3","colz");
*/



TTree * tree = 0;
TTree * treeCache=0;
/// Load tree and defining derived  information and metadata
void loadTree() {
  AliExternalInfo info;
  tree = info.GetChain("QA.TPC", "LHC17*", "cpass1_pass1", "QA.EVS;QA.rawTPC");
  tree->SetAlias("interactionRate", "QA.EVS.interactionRate");
  ///
  tree->SetAlias("qmaxQASum", "Sum$(qmaxQA.fElements*((abs(qmaxQA.fElements-40)<20)))/Sum$((abs(qmaxQA.fElements-40)<20))");
  tree->SetAlias("qmaxQASumIn", "Sum$(qmaxQA.fElements*((Iteration$<36&&abs(qmaxQA.fElements-40)<20)))/Sum$((Iteration$<36&&abs(qmaxQA.fElements-40)<20))");
  tree->SetAlias("qmaxQASumOut", "Sum$(qmaxQA.fElements*((Iteration$>=36&&abs(qmaxQA.fElements-40)<20)))/Sum$((Iteration$>=36&&abs(qmaxQA.fElements-40)<20))");
  tree->SetAlias("qmaxQASumR", "qmaxQASumIn/qmaxQASum");
  tree->SetAlias("meanMIPeleR", "meanMIPele/meanMIP");
  tree->SetAlias("bz0", "bz+rndm*0.0001");
  tree->SetMarkerStyle(21); tree->SetMarkerSize(0.5);
}

/// TMVA can not work with friend trees with indeces, respec. with array of the measurements
///  => we have to create  input "flat tree" for MVA learning
void cacheTree(){
  AliTreePlayer::MakeCacheTree(tree,"resolutionMIP:meanMIPeleR:tpcItsMatchA:bz0:interactionRate:qmaxQASum:qmaxQASumIn:qmaxQASumOut:qmaxQASumR:run:time","TMVAInput.root","MVAInput","meanMIP>30&&run==QA.EVS.run");
}

/// Register methods used for regression
void RegisterFitters(){
  /// DNN - for the moment not used as need BLASS library  - not in default AliRoot
  //TString dnnConfigString = "!H:V::VarTransform=Norm:ErrorStrategy=CROSSETROPY:WeightInitialization=XAVIERUNIFORM:";
  TString dnnConfigString = "!H:V::VarTransform=Norm:WeightInitialization=XAVIERUNIFORM:";
  TString  dnnLayoutString="Layout=TANH|100,TANH|50,TANH|10,LINEAR:";
  TString dnnTrainingString ="LearningRate=1e-5,Momentum=0.5,Repetitions=1,ConvergenceSteps=500,BatchSize=50,TestRepetitions=7,WeightDecay=0.01,Regularization=NONE,DropConfig=0.5+0.5+0.5+0.5,DropRepetitions=2";
  TString dnnConfig=dnnConfigString+dnnLayoutString+ dnnTrainingString;
  ///
  AliNDFunctionInterface::registerMethod("BDTRF25_8","!H:!V:NTrees=25:Shrinkage=0.1:UseRandomisedTrees:nCuts=20:MaxDepth=8",TMVA::Types::kBDT);
  AliNDFunctionInterface::registerMethod( "BDTRF12_16", "!H:!V:NTrees=12:Shrinkage=0.1:UseRandomisedTrees:nCuts=20:MaxDepth=16", TMVA::Types::kBDT);
  AliNDFunctionInterface::registerMethod("KNN","nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim", TMVA::Types::kKNN);
  AliNDFunctionInterface::registerMethod("MLP", "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator",TMVA::Types::kMLP);
  AliNDFunctionInterface::registerMethod("DNN",dnnConfig.Data(), TMVA::Types::kDNN);
}

void makeMVAFits(){
  TFile *f= TFile::Open("TMVAInput.root");
  f->GetObject("MVAInput",treeCache);
  AliNDFunctionInterface::FitMVA(treeCache, "resolutionMIP","interactionRate>0", "interactionRate:bz0:qmaxQASum:qmaxQASumR","BDTRF25_8:BDTRF12_16:MLP:KNN");
  AliNDFunctionInterface::FitMVA(treeCache, "meanMIPeleR","interactionRate>0", "interactionRate:bz0:qmaxQASum:qmaxQASumR","BDTRF25_8:BDTRF12_16:MLP:KNN");
  AliNDFunctionInterface::FitMVA(treeCache, "tpcItsMatchA","interactionRate>0", "interactionRate:bz0:qmaxQASum:qmaxQASumR","BDTRF25_8:BDTRF12_16:MLP:KNN");
}

/// Emulation of the bootstrap - training repeated several time
void makeMVABootstrapMI(){
  TFile *f= TFile::Open("TMVAInput.root");
  f->GetObject("MVAInput",treeCache);
  for (Int_t iBoot=0; iBoot<10; iBoot++) {
    AliNDFunctionInterface::FitMVA(treeCache, "resolutionMIP", "interactionRate>0", "interactionRate:bz0:qmaxQASum:qmaxQASumR", "BDTRF25_8:BDTRF12_16:KNN", 0, iBoot);
  }
}

///
void loadMVAreaders(){
  AliNDFunctionInterface::LoadMVAReader(0,"TMVA_RegressionOutput.root","MLP","dir_resolutionMIP");
  AliNDFunctionInterface::LoadMVAReader(1,"TMVA_RegressionOutput.root","KNN","dir_resolutionMIP");
  AliNDFunctionInterface::LoadMVAReader(2,"TMVA_RegressionOutput.root","BDTRF25_8","dir_resolutionMIP");
  AliNDFunctionInterface::LoadMVAReader(3,"TMVA_RegressionOutput.root","BDTRF12_16","dir_resolutionMIP");
}
///
void loadMVAreadersBootstrap() {
  AliNDFunctionInterface::LoadMVAReaderArray(0,"TMVA_RegressionOutput.root","BDTRF12_16",".*resolutionMIP");
  AliNDFunctionInterface::LoadMVAReaderArray(1,"TMVA_RegressionOutput.root","BDTRF25_8",".*resolutionMIP");
  AliNDFunctionInterface::LoadMVAReaderArray(2,"TMVA_RegressionOutput.root","KNN",".*resolutionMIP");
}

void drawExamplePlots(){
  loadTree();
  cacheTree();
  RegisterFitters();
  makeMVABootstrapMI();
  loadMVAreadersBootstrap();
  /// 0.) Draw correlation plot
  tree->Draw("AliNDFunctionInterface::EvalMVAStat(0,1,interactionRate, bz0,qmaxQASum,qmaxQASumR):resolutionMIP:run","run==QA.EVS.run","colz")
  gPad->SaveAs("resolutionMIPvsfit.png");
  /// 1.) Draw ration regression/fit as function of run
  tree->Draw("AliNDFunctionInterface::EvalMVAStat(0,1,interactionRate, bz0,qmaxQASum,qmaxQASumR)/resolutionMIP:run:interactionRate","run==QA.EVS.run","colz");
  gPad->SaveAs("resolutionMIPFitRatiovsRun.png");
  /// 2.) compare  of the regression RMS for method array 0 and method array1
  tree->Draw("AliNDFunctionInterface::EvalMVAStat(0,2,interactionRate, bz0,qmaxQASum,qmaxQASumR):AliNDFunctionInterface::EvalMVAStat(1,2,interactionRate, bz0,qmaxQASum,qmaxQASumR)","run==QA.EVS.run","colz")
}