/// \ingroup Macros
/// \file    QAtrendingFitExample.C
/// \brief   Demo usage of the information from the AliExternalInfo and visualization using the TStatToolkit, AliTreePlayer and AliDrawStyle
///
/// \author  Marian Ivanov
/// \code
///  AliDrawStyle::SetDefaults()
///  AliDrawStyle::ApplyStyle("figTemplate");
///  gStyle->SetOptTitle(1);
///  .L $AliRoot_SRC/STAT/Macros/AliNDFunctionInterface.cxx+
///  .L $AliRoot_SRC/STAT/Macros/QAtrendingFitExample.C
/// \endcode
///  0.)  loadTree();
///  1.)  cacheTree();      /// needed in case friend trees or arrays - which are not supported by TMVA
///  3.)  RegisterFitters();
///  4.)  makeMVAFits();
///  5.)  loadMVAReaders();
///  tree->Draw("AliNDFunctionInterface::EvalMVA(2,interactionRate, bz0,qmaxQASum,qmaxQASumR):resolutionMIP:interactionRate","run==QA.EVS.run&&meanMIP>40&&bz0<-0.3","colz");



Bool_t useDNN=kFALSE;
TTree * tree = 0;
TTree * treeCache=0;
/// Load tree and defining derived  information (TTree aliases)  and metadata
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

/// CacheTree input variables to tree format usable by TMVA
///-------------------------------
/// TMVA can not work with friend trees with indices, resp. with array of the measurements, resp aliases and functions
/// * input "flat tree" for MVA learning to be created with variables of interest
/// * AliTreePlayer::MakeCacheTree to create flat input tree for TMVA
void cacheTree(){
  AliTreePlayer::MakeCacheTree(tree,"resolutionMIP:meanMIPeleR:tpcItsMatchA:bz0:interactionRate:qmaxQASum:qmaxQASumIn:qmaxQASumOut:qmaxQASumR:run:time","TMVAInput.root","MVAInput","meanMIP>30&&run==QA.EVS.run");
}

/// Register example methods used for regression
///----------------------------------------
/// * BDT and MLP example
/// * DNN - for the moment not used as need BLASS library  - not in default AliRoot
///   * Naive adaptation of the DNN configuration from the ROOT tutorials  TMVARegression.C -
///   * slow 100 times smaller than MLP
///   * lead to floating point exception
void RegisterFitters() {
  TString layoutString("Layout=TANH|20,LINEAR");
  TString training0("LearningRate=1e-5,Momentum=0.5,Repetitions=1,ConvergenceSteps=500,BatchSize=50,"
                    "TestRepetitions=7,WeightDecay=0.01,Regularization=L1,DropConfig=0.5+0.5+0.5+0.5,"
                    "DropRepetitions=2");
  TString training1("LearningRate=1e-5,Momentum=0.9,Repetitions=1,ConvergenceSteps=170,BatchSize=30,"
                    "TestRepetitions=7,WeightDecay=0.01,Regularization=L1,DropConfig=0.1+0.1+0.1,DropRepetitions="
                    "1");
  TString trainingStrategyString("TrainingStrategy=");
  trainingStrategyString += training0 + "|" + training1;
  TString dnnOptions("!H:V:ErrorStrategy=SUMOFSQUARES:VarTransform=G:WeightInitialization=XAVIERUNIFORM:Architecture=CPU");
  dnnOptions.Append(":");
  dnnOptions.Append(layoutString);
  dnnOptions.Append(":");
  dnnOptions.Append(trainingStrategyString);
  ///
  AliNDFunctionInterface::registerMethod("BDTRF25_8","!H:!V:NTrees=25:Shrinkage=0.1:UseRandomisedTrees:nCuts=20:MaxDepth=8",TMVA::Types::kBDT);
  AliNDFunctionInterface::registerMethod("BDTRF12_16", "!H:!V:NTrees=12:Shrinkage=0.1:UseRandomisedTrees:nCuts=20:MaxDepth=16", TMVA::Types::kBDT);
  AliNDFunctionInterface::registerMethod("KNN","nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim", TMVA::Types::kKNN);
  AliNDFunctionInterface::registerMethod("MLP", "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator",TMVA::Types::kMLP);
  AliNDFunctionInterface::registerMethod("DNN_CPU",dnnOptions.Data(), TMVA::Types::kDNN);
}

/*
/// Make MVA regression example
///----------------------------
void makeMVAFits(){
  TFile *f= TFile::Open("TMVAInput.root");
  f->GetObject("MVAInput",treeCache);
  AliNDFunctionInterface::FitMVA(treeCache, "resolutionMIP","interactionRate>0", "interactionRate:bz0:qmaxQASum:qmaxQASumR","BDTRF25_8:BDTRF12_16:MLP:KNN");
  AliNDFunctionInterface::FitMVA(treeCache, "meanMIPeleR","interactionRate>0", "interactionRate:bz0:qmaxQASum:qmaxQASumR","BDTRF25_8:BDTRF12_16:MLP:KNN");
  AliNDFunctionInterface::FitMVA(treeCache, "tpcItsMatchA","interactionRate>0", "interactionRate:bz0:qmaxQASum:qmaxQASumR","BDTRF25_8:BDTRF12_16:MLP:KNN");
}
*/


/// Emulation of the TMVA bootstrap - training repeated several time with different subset of data
/// TODO - Implement real bootstrap ( random sampling with replacement + other methods) in the TMVA (to check with ROOT)
/// TODO - DNN example - TO USE DNN (Deep neural network)  BLASS or CBLASS library has to be enabled in ROOT. This is not the case for the default aliBuild recipes
void makeMVABootstrapMI(Int_t nRegression=10){
  TFile *f= TFile::Open("TMVAInput.root");
  f->GetObject("MVAInput",treeCache);
  gSystem->Unlink("TMVA_RegressionOutput.root");
  TString output="TMVA_RegressionOutput.root#";
  for (Long_t iBoot=0; iBoot<nRegression; iBoot++) {
    AliNDFunctionInterface::FitMVARegression(output+"resolutionMIP"+iBoot,treeCache, "resolutionMIP", "interactionRate>0", "interactionRate:bz0:qmaxQASum:qmaxQASumR", "BDTRF25_8:BDTRF12_16:KNN", "");
    AliNDFunctionInterface::FitMVARegression(output+"meanMIPeleR"+iBoot,treeCache, "meanMIPeleR", "interactionRate>0", "interactionRate:bz0:qmaxQASum:qmaxQASumR", "BDTRF25_8:BDTRF12_16:KNN","");
    AliNDFunctionInterface::FitMVARegression(output+"tpcItsMatchA"+iBoot,treeCache, "tpcItsMatchA", "interactionRate>0", "interactionRate:bz0:qmaxQASum:qmaxQASumR", "BDTRF25_8:BDTRF12_16:KNN","");
  }
}

/// Load regression and register it for later usage
void loadMVAReaders(){
  AliNDFunctionInterface::LoadMVAReader(0,"TMVA_RegressionOutput.root","BDTRF25_8","resolutionMIP0");
  AliNDFunctionInterface::LoadMVAReader(1,"TMVA_RegressionOutput.root","BDTRF12_16","resolutionMIP0");
  AliNDFunctionInterface::LoadMVAReader(2,"TMVA_RegressionOutput.root","KNN","resolutionMIP0");
}


/// Load array of regression  -used later in the array regression evaluation (mean, median, rms)
///-------------------------------
void loadMVAReadersBootstrap() {
  AliNDFunctionInterface::LoadMVAReaderArray(0,"TMVA_RegressionOutput.root","BDTRF12_16",".*resolutionMIP");
  AliNDFunctionInterface::LoadMVAReaderArray(1,"TMVA_RegressionOutput.root","BDTRF25_8",".*resolutionMIP");
  AliNDFunctionInterface::LoadMVAReaderArray(2,"TMVA_RegressionOutput.root","KNN",".*resolutionMIP");
}


/// QAtrendingFitExample
/// ----------------------
void QAtrendingFitExample(){
  /// 0.) Remove regression file
  gSystem->Unlink("TMVA_RegressionOutput.root");
  /// 1.) Load Input data
  loadTree();
  /// 2.) Cache tree - TMVA expect variables  - not functions and aliases
  cacheTree();
  /// 3.) Register fitters
  RegisterFitters();
  /// 4.) Make bootstrap regression
  makeMVABootstrapMI(10);
  /// 5.) Load array of regression readers
  loadMVAReaders();            /// load one of the bootstrap method
  loadMVAReadersBootstrap();   /// load bootstrap arrays
  ///
  /// 6.) Draw example plots
  /// ===============================
  tree->Draw("AliNDFunctionInterface::EvalMVAStat(0,1,interactionRate, bz0,qmaxQASum,qmaxQASumR):AliNDFunctionInterface::EvalMVA(0,interactionRate, bz0,qmaxQASum,qmaxQASumR):resolutionMIP","run==QA.EVS.run","colz");
  /// 6.1 Correlation plot - mean regression of resolutionMIP vs observed resolutionMIP
  /// ----------------------------------
  /// \code tree->Draw("AliNDFunctionInterface::EvalMVAStat(0,1,interactionRate, bz0,qmaxQASum,qmaxQASumR):resolutionMIP:run","run==QA.EVS.run","colz"); \endcode
  tree->Draw("AliNDFunctionInterface::EvalMVAStat(0,1,interactionRate, bz0,qmaxQASum,qmaxQASumR):resolutionMIP:run","run==QA.EVS.run","colz");
  gPad->SaveAs("resolutionMIPvsfit.png");
  /// 6.2 Draw ration regression/fit as function of run
  /// ----------------------------------------
  /// \code tree->Draw("AliNDFunctionInterface::EvalMVAStat(0,1,interactionRate, bz0,qmaxQASum,qmaxQASumR):resolutionMIP:run","run==QA.EVS.run","colz"); \endcode
  tree->Draw("AliNDFunctionInterface::EvalMVAStat(0,1,interactionRate, bz0,qmaxQASum,qmaxQASumR)/resolutionMIP:run:interactionRate","run==QA.EVS.run","colz");
  gPad->SaveAs("resolutionMIPFitRatiovsRun.png");
  /// 6.3  compare  of the regression RMS for method array 0 and method array1
  /// ----------------------------------------
  /// \code tree->Draw("AliNDFunctionInterface::EvalMVAStat(0,1,interactionRate, bz0,qmaxQASum,qmaxQASumR):resolutionMIP:run","run==QA.EVS.run","colz"); \endcode
  tree->Draw("AliNDFunctionInterface::EvalMVAStat(0,2,interactionRate, bz0,qmaxQASum,qmaxQASumR):AliNDFunctionInterface::EvalMVAStat(1,2,interactionRate, bz0,qmaxQASum,qmaxQASumR)","run==QA.EVS.run","colz");
}
