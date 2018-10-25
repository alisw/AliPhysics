/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
///////////////////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////
/// \file AliNDFunctionInterface.cxx
/*
  NDimensional function interpolation using THn  and using TMVA reader interface
  .L $AliRoot_SRC/STAT/Macros/AliNDFunctionInterface.cxx+

*/

#include "TObjArray.h"
#include "TH1.h"

#include "TGraph.h"
#include "TTreeStream.h"
#include "TVirtualFFT.h"
#include "TMath.h"
#include "TVector.h"
#include "TStatToolkit.h"
#include <stdio.h>
#include "TPRegexp.h"
#include "THn.h"
#include "TKey.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/DataLoader.h"
#include "TMVA/MethodBase.h"
#include "AliNDFunctionInterface.h"

//Int_t AliNDFunctionInterface::fVerbose=0;
///
/// \brief Linear interpolation of the bin content
/// \param his       - input N-dimensional histogram to interpolate
/// \param xyz       - NDimensional point where to interpolate
/// \param verbose   - verbosity flag
/// \return
Double_t AliNDFunctionInterface::GetInterpolationLinear(THn * his, Double_t *xyz, Int_t  verbose){
  const Int_t kMaxDim=100;
  Int_t nDim=his->GetNdimensions();
  Double_t deltaXYZ[kMaxDim];
  Double_t ddXYZ[kMaxDim];
  Int_t bin0, index0[kMaxDim];
  Double_t val0,value;
  bin0=his->GetBin(xyz);
  val0=his->GetBinContent(bin0,index0);
  value=val0;
  if (verbose==0) return val0;
  for (Int_t i=0;i<nDim; i++) {
    deltaXYZ[i] = (his->GetAxis(i)->GetBinCenter(index0[i]) - xyz[i])/his->GetAxis(i)->GetBinWidth(index0[i]);
    Int_t dBin = (deltaXYZ[i] > 0) ? 1 : -1;
    if (verbose>2) dBin*=-1;
    index0[i] += dBin;
    if (index0[i]<1) {index0[i]=1; dBin=0;};
    if (index0[i]>=his->GetAxis(i)->GetNbins()) {index0[i]=his->GetAxis(i)->GetNbins(); dBin=0;};
    ddXYZ[i] = -(his->GetBinContent(index0)-val0)*dBin;
    index0[i] -= dBin;
    value += ddXYZ[i]*deltaXYZ[i];
    if (verbose>1){ printf("%d\t%f\t%f\t%f\n", i, val0, deltaXYZ[i], ddXYZ[i]); }
  }
  return value;
}

/// \brief Linear interpolation of the numerical derivative
/// dV/dx \approx (V(i+1)-V(i-1))/(2(delta)
/// \param his
/// \param xyz
/// \param dIndex
/// \param verbose
/// \return
Double_t AliNDFunctionInterface::GetDeltaInterpolationLinear(THn * his, Double_t *xyz, Int_t dIndex, Int_t  verbose){
  const Int_t kMaxDim=100;
  Int_t nDim=his->GetNdimensions();
  Double_t deltaXYZ[kMaxDim];
  Double_t ddXYZ[kMaxDim];
  Int_t bin0, index0[kMaxDim],index1[kMaxDim];
  Double_t val0,value;
  bin0=his->GetBin(xyz);                   // Get nearest bin
  his->GetBinContent(bin0,index0);         // get Nd index of closest bin
  memcpy(index1,index0,nDim*sizeof(Int_t));
  if (index1[dIndex]<=1)  index1[dIndex]=2;
  if (index1[dIndex]>=his->GetAxis(dIndex)->GetNbins()-1)  index1[dIndex]=his->GetAxis(dIndex)->GetNbins();
  index1[dIndex]+=1;
  Double_t derivative=his->GetBinContent(index1);
  Double_t dx=his->GetAxis(dIndex)->GetBinWidth(index1[dIndex]);
  index1[dIndex]-=2;
  derivative-=his->GetBinContent(index1);
  dx+=his->GetAxis(dIndex)->GetBinWidth(index1[dIndex]);
  derivative/=dx;
  value=derivative;
  for (Int_t i=0;i<nDim; i++) {
    deltaXYZ[i] = (his->GetAxis(i)->GetBinCenter(index0[i]) - xyz[i])/his->GetAxis(i)->GetBinWidth(index0[i]);
    Int_t dBin = (deltaXYZ[i] > 0) ? 1 : -1;
    if (verbose>2) dBin*=-1;
    index0[i] += dBin;
    if (index0[i]<1) {index0[i]=1; dBin=0;};
    if (index0[i]>=his->GetAxis(i)->GetNbins()) {index0[i]=his->GetAxis(i)->GetNbins(); dBin=0;};
    //
    memcpy(index1,index0,nDim*sizeof(Int_t));
    if (index1[dIndex]<=1)  index1[dIndex]=2;
    if (index1[dIndex]>=his->GetAxis(dIndex)->GetNbins()-1)  index1[dIndex]=his->GetAxis(dIndex)->GetNbins();
    //
    index1[dIndex]+=1;
    Double_t derivativeL=his->GetBinContent(index1);
    index1[dIndex]-=2;
    derivativeL-=his->GetBinContent(index1);
    derivativeL/=dx;
    ddXYZ[i] = (derivativeL-derivative)*dBin;
    index0[i] -= dBin;
    value -= ddXYZ[i]*deltaXYZ[i];
    if (verbose>1){ printf("%d\t%f\t%f\t%f\n", i, derivative, deltaXYZ[i], ddXYZ[i]); }
  }
  return value;
}

/// FitMVAClassification - TODO add documentation and test
/// \param output
/// \param inputTrees
/// \param cuts
/// \param variables
/// \param methods
/// \param factoryString
/// \return
Int_t  AliNDFunctionInterface::FitMVAClassification(const char * output, const char *inputTrees, const char *cuts, const char * variables, const char *methods, const char * factoryString){

  TObjArray *treeList = TString(inputTrees).Tokenize(":");
  TObjArray *cutList = TString(cuts).Tokenize(":");

  TString inputTreesString=TString(inputTrees);
  if(treeList->GetEntries()==1){
    for(Int_t i=1;i<cutList->GetEntries() ; i++){
      inputTreesString+=":"+TString(inputTrees);
    }
  }
  delete treeList;
  treeList = inputTreesString.Tokenize(":");

  if (treeList->GetEntries()<=1){
    ::Error("AliNDFunctionInterface::FitMVAClassification","More than one input tree required.");
    return -1;
  }

  if (treeList->GetEntries()<cutList->GetEntries()){
    ::Error("AliNDFunctionInterface::FitMVAClassification","Entries in cut list exceed number of trees.");
    return -1;
  }
  
  TTree * tree_in[treeList->GetEntries()];
  TFile * input[treeList->GetEntries()];
  TCut cut[treeList->GetEntries()];
  TString tree_name[treeList->GetEntries()];
  
  for (Int_t i=0; i<treeList->GetEntries(); i++) cut[i]="";
  for (Int_t i=0; i<cutList->GetEntries(); i++) cut[i]=cutList->At(i)->GetName();
  for (Int_t i=0; i<treeList->GetEntries(); i++){
  	 TObjArray * inTree = TString(treeList->At(i)->GetName()).Tokenize("#");
 	 input[i] = TFile::Open( TString(inTree->At(0)->GetName() ));
    //::Info("AliNDFunctionInterface::FitMVAClassification","Load file %s",inTree->At(0)->GetName());
 	 tree_in[i] = (TTree*) input[i]->Get(inTree->At(1)->GetName());
    ::Info("AliNDFunctionInterface::FitMVAClassification","Load tree %s",inTree->At(1)->GetName());
    tree_name[i]=inTree->At(1)->GetName();
    delete inTree;
  }

  TObjArray *dirFile=TString(output).Tokenize("#");
  if(dirFile->GetEntries()<=1){
    ::Error("AliNDFunctionInterface::FitMVAClassification","Empty directory list");
    return -1;
  }

  auto outputFile = TFile::Open(dirFile->At(0)->GetName(), "update");
  if(outputFile==NULL || outputFile->GetListOfKeys()->Contains(dirFile->At(1)->GetName())){
    ::Error("AliNDFunctionInterface::FitMVAClassification","Output directory already exists");
    return -1;
  }
  TString directory=TString(dirFile->At(1)->GetName());
  TString fstring;
  if(treeList->GetEntries()==2){
    fstring="!V:!Silent:Color:DrawProgressBar:AnalysisType=classification";
  } 
  else {
    fstring="!V:!Silent:Color:DrawProgressBar:AnalysisType=multiclass";
  }
  TMVA::Factory factory("TMVAMultiClass", outputFile, fstring );

  // Declare DataLoader
  TMVA::DataLoader loader(directory);
  TObjArray *variableList = TString(variables).Tokenize(":");
  if (variableList->GetEntries()<=0){
    ::Error("AliNDFunctionInterface::FitMVAClassification","Empty variable list");
    return -1;
  }
  for (Int_t i=0; i<variableList->GetEntries(); i++){
    loader.AddVariable(variableList->At(i)->GetName());
  }
  delete variableList;

  for (Int_t i=0; i<treeList->GetEntries(); i++){
    loader.AddTree(tree_in[i],tree_name[i]+Form("%d",i),1.0,cut[i]);
  }

  loader.PrepareTrainingAndTestTree("", FactorySetting[factoryString]);

  TObjArray *methodList = TString(methods).Tokenize(":");
  if (methodList->GetEntries()<=0){
    ::Error("AliNDFunctionInterface::FitMVAClassification","Empty method list");
    return -1;
  }
  for (Int_t i=0; i<methodList->GetEntries(); i++){
    std::string methodName=methodList->At(i)->GetName();
    printf("Booking method %s\t%d\t%s\n",methodName.data(), regressionMethodID[methodName], regressionMethodSetting[methodName].data());
    factory.BookMethod(&loader,regressionMethodID[methodName], methodName.data(), regressionMethodSetting[methodName].data());
  }
  factory.TrainAllMethods();

  /// Write ascii weight files to root file
  outputFile->cd(directory);
  for (Int_t i=0; i<methodList->GetEntries(); i++) {
    std::string methodName = methodList->At(i)->GetName();
    TString weightFile = TString::Format(directory+"/weights/TMVAMultiClass_%s.weights.xml",methodName.data());
    printf("Weight file\t%s\n", weightFile.Data());
    TObjString weight(gSystem->GetFromPipe(TString::Format("cat %s",weightFile.Data())));
    weight.Write(methodName.data());
  }
  TObjString varList(variables);
  varList.Write("variables");
  delete methodList;
  factory.TestAllMethods();
  factory.EvaluateAllMethods();
  outputFile->Close();

  return 0;
}

///  FitMVARegression wrapper = Do TMVA regression and save weights into output root file
/// \param output              - output path <file>#directory/
///                            - output file is opened in update mode - in case key (directory) exist - fit should failed - data are not written (currently it failed)
/// \param tree                - input tree
/// \param varFit              - variable to fit (for the moment 1)
/// \param cut                 - cut formula
/// \param variables           - : separated list of explanatory variables
/// \param methods             - : separated list of the regression methods (have to be registered before)
/// \param factoryString      - factory string (e.g specifying number of training and test samples). Factory string to be registered before. If empty default used
/// \return                    - Fit parameters will be saved in the output destination - return value 0 in case of success

/// Example usage in test macro ( QAtrendingFitExample.C function makeMVAFits() and  makeMVABootstrapMI() ), e.g:
/// -----------
/// \code
/// TString output="TMVA_RegressionOutput.root#";
///  for (Int_t iBoot=0; iBoot<nRegression; iBoot++) {
///    AliNDFunctionInterface::FitMVARegression(output+"resolutionMIP"+iBoot,treeCache, "resolutionMIP", "interactionRate>0", "interactionRate:bz0:qmaxQASum:qmaxQASumR", "BDTRF25_8:BDTRF12_16:KNN", "");
///    AliNDFunctionInterface::FitMVARegression(output+"meanMIPeleR"+iBoot,treeCache, "meanMIPeleR", "interactionRate>0", "interactionRate:bz0:qmaxQASum:qmaxQASumR", "BDTRF25_8:BDTRF12_16:KNN","");
///    AliNDFunctionInterface::FitMVARegression(output+"tpcItsMatchA"+iBoot,treeCache, "tpcItsMatchA", "interactionRate>0", "interactionRate:bz0:qmaxQASum:qmaxQASumR", "BDTRF25_8:BDTRF12_16:KNN","");
///  }
/// \endcode
/// Algorithm:
/// ---------
Int_t  AliNDFunctionInterface::FitMVARegression(const char * output, TTree *tree, const char *varFit, TCut cut, const char * variables, const char *methods,const char * factoryString){
  /// * 0. Declare Factory
  TObjArray *dirFile=TString(output).Tokenize("#");
  if(dirFile->GetEntries()<=1){
	 ::Error("AliNDFunctionInterface::FitMVARegression","Empty directory list");
	 return -1;
  }
  auto outputFile = TFile::Open(dirFile->At(0)->GetName(), "update");
  if(outputFile->GetListOfKeys()->Contains(dirFile->At(1)->GetName())){
    ::Error("AliNDFunctionInterface::FitMVA","Output directory already exists");
    return -1;
  } 
  TString directory=TString(dirFile->At(1)->GetName());
  TMVA::Factory factory("TMVARegression", outputFile, "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression" );
  TString varNameFit=varFit;
  /// * 1.) Declare DataLoader
  TMVA::DataLoader loader(directory);
  //  * 2.) Add the feature variables from the variables list
  TObjArray *weight_split = TString(variables).Tokenize(";");
  TObjArray *variableList = TString(weight_split->At(0)->GetName()).Tokenize(":");
  TString weights = "";
  if(weight_split->GetEntries()==2){
    weights = TString(weight_split->At(1)->GetName());
  }
  if (variableList->GetEntries()<=0){
    ::Error("AliNDFunctionInterface::FitMVARegression","Empty variable list");
    return -1;
  }
  for (Int_t i=0; i<variableList->GetEntries(); i++){
    loader.AddVariable(variableList->At(i)->GetName());
    /// TODO check existence - validity of variable - in case of error - exit with error message
  }
  delete variableList;
  loader.AddTarget(varFit);
  /// * 3.) Setup DataSet
  loader.AddRegressionTree(tree, 1.0);   // link the TTree to the loader, weight for each event  = 1
  
  Int_t entries = tree->Draw("1",cut,"goff");
  loader.PrepareTrainingAndTestTree(cut, FactorySetting[factoryString]);
  if (weights!="") loader.SetWeightExpression(weights,"Regression");
  /// * 4.) Book regression methods from the methods list
  TObjArray *methodList = TString(methods).Tokenize(":");
  if (methodList->GetEntries()<=0){
    ::Error("AliNDFunctionInterface::FitMVARegression","Empty method list");
    return -1;
  }
  for (Int_t i=0; i<methodList->GetEntries(); i++){
    std::string methodName=methodList->At(i)->GetName();
    printf("Booking method %s\t%d\t%s\n",methodName.data(), regressionMethodID[methodName], regressionMethodSetting[methodName].data());
    factory.BookMethod(&loader,regressionMethodID[methodName], methodName.data(), regressionMethodSetting[methodName].data());
  }
  /// * 5.) Train all methods
  factory.TrainAllMethods();
  /// * 6.) Write ascii weight files to root file
  outputFile->cd(directory);
  for (Int_t i=0; i<methodList->GetEntries(); i++) {
    std::string methodName = methodList->At(i)->GetName();
    TString weightFile = TString::Format(directory+"/weights/TMVARegression_%s.weights.xml",methodName.data());
    printf("Weight file\t%s\n", weightFile.Data());
    TObjString weight(gSystem->GetFromPipe(TString::Format("cat %s",weightFile.Data())));
    weight.Write(methodName.data());
  }
  
  TObjString varList(weight_split->At(0)->GetName());
  varList.Write("variables");
  delete methodList;
  /// * 7.) Evaluate all MVAs using the set of test events             /// TODO - make an option
  factory.TestAllMethods();
  /// * 8.) Evaluate and compare performance of all configured MVAs    /// TODO - make an option
  factory.EvaluateAllMethods();
  /// * Save the output
  outputFile->Close();
  return 0;
}

/// Load MVA regression object and register it in the map of available readers fr evaluation in TFormula
/// Current ION of TMVA does not allow standard persistence in root file
///       Way around:
///          Writing -  ASCII files are stored as TSting in the root file
///          Reading -  String written as ASCII file and read back by TMVA::Reader
/// \param id            - regression method ID  (using e.g. hash value)
/// \param inputFile     - input file  (e.g TMVA_RegressionOutput.root)
/// \param method        - method name (e.g MLP)
/// \param dir           - directory usually coding regression variable with dir_ prefix (e.g. dir_meanMIPele)
/// \return              - error code
/// TODOs:
/// -----------------------------------------------------------------------------------------------------------
/// * TODO - ASK TMVA team to implement standard IO
/// * TODO - ASK TMVA team to disable requirement to set variable addresses or find switch to make it
///
/// Example usage (see also QAtrendingFitExample.C loadMVAReaders()):
/// -----------------------------------------------------------------------------------------------------------
/// \code
/// void loadMVAReaders(){
///  AliNDFunctionInterface::LoadMVAReader(0,"TMVA_RegressionOutput.root","BDTRF25_8","resolutionMIP0");
///  AliNDFunctionInterface::LoadMVAReader(1,"TMVA_RegressionOutput.root","BDTRF12_16","resolutionMIP0");
///  AliNDFunctionInterface::LoadMVAReader(2,"TMVA_RegressionOutput.root","KNN","resolutionMIP0");
///}
/// \endcode
TMVA::MethodBase * AliNDFunctionInterface::LoadMVAReader(Int_t id/*=0*/, const char * inputFile/*="TMVA_RegressionOutput.root"*/, const char *method/*="MLP"*/, const char *dir/*="dir_resolutionMIP"*/){
  TMVA::Reader *reader;
  reader = new TMVA::Reader( "!Color:!Silent" );
  auto fin = TFile::Open(inputFile);
  if (fin==NULL){
    ::Error("LoadMVAReader","Invalid input file\t%s", inputFile);
    return NULL;
  }
  TDirectoryFile *dir0 = 0,*dirVar=0;
  fin->GetObject(dir,dir0);
  if (dir0==NULL) {
    ::Error("LoadMVAReader","Invalid input directory\t%s\t%s", dir, inputFile);
    return NULL;
  }
  TObjString *varList;
  dir0->GetObject("variables",varList);
  TObjArray *variableList = varList->GetString().Tokenize(":");
  for (Int_t i=0; i<variableList->GetEntries(); i++){
       Float_t value;
       reader->AddVariable(variableList->At(i)->GetName(),&value);
       /// TODO check existence - validity of variable - in case of error - exit with error message
  }
  delete varList;
  delete variableList;

  /// write weight from the root file to txt files as it is expected by reader
  FILE * pFile;
  pFile = fopen ("weights.xml","w");
  TObjString *methodWeights= NULL;

  if (dir0==NULL) {
    ::Error("LoadMVAReader","Object %s not present in  input directory\t%s/%s",  method, inputFile,dir);
    delete pFile;
    return NULL;
  }
  dir0->GetObject(method,methodWeights);
  if (methodWeights==NULL){
    ::Error("LoadMVAReader","Object %s not present in  input directory\t%s/%s",  method, inputFile,dir);
    delete pFile;
    return NULL;
  }
  fprintf (pFile, "%s",methodWeights->GetName());
  fclose (pFile);
  reader->BookMVA( method, "weights.xml");
  if (id>=0) AliNDFunctionInterface::readerMethodBase[id]=dynamic_cast<TMVA::MethodBase *>(reader->FindMVA(method));
  return dynamic_cast<TMVA::MethodBase *>(reader->FindMVA(method));
}

/// Load array of the MVA reader and register it in the AliNDFunctionInterface method array maps ( readerMethodBaseArray)
/// Method arrays could be later use  for the TMVA array evaluation ()
/// \param id
/// \param inputFile
/// \param methodMask
/// \param dirMask
/// \return  error code - 0 mean OK
///
/// Example usage in  (see also QAtrendingFitExample.C loadMVAReadersBootstrap())
/// -----------------------------------------------------------------------------------------------------------
/// \code
/// void loadMVAReadersBootstrap() {
///  AliNDFunctionInterface::LoadMVAReaderArray(0,"TMVA_RegressionOutput.root","BDTRF12_16",".*resolutionMIP");
///  AliNDFunctionInterface::LoadMVAReaderArray(1,"TMVA_RegressionOutput.root","BDTRF25_8",".*resolutionMIP");
///  AliNDFunctionInterface::LoadMVAReaderArray(2,"TMVA_RegressionOutput.root","KNN",".*resolutionMIP");
/// }
/// \endcode
Int_t  AliNDFunctionInterface::LoadMVAReaderArray(Int_t id, const char * inputFile, const char *methodMask, const char *dirMask){
  auto fin = TFile::Open(inputFile);
  if (fin==NULL){
    ::Error("LoadMVAReader","Invalid input file\t%s", inputFile);
    return 1;
  }
  TList * keyArrayDir = fin->GetListOfKeys();
  TPRegexp regDir(dirMask);
  TPRegexp regMask(methodMask);
  for (Int_t iDir=0; iDir<keyArrayDir->GetEntries(); iDir++){
    if (regDir.Match(keyArrayDir->At(iDir)->GetName())==0) continue;
    printf("%s\n",keyArrayDir->At(iDir)->GetName());
    if (fin->cd(keyArrayDir->At(iDir)->GetName()) != 0) {
      TList *keyArrayMethod = fin->GetDirectory(keyArrayDir->At(iDir)->GetName())->GetListOfKeys();
      for (Int_t iMethod = 0; iMethod < keyArrayDir->GetEntries(); iMethod++) {
        TKey *key = (TKey *) keyArrayMethod->At(iMethod);
        if (key == 0) continue;
        if (regMask.Match(keyArrayMethod->At(iMethod)->GetName()) == 0) continue;
        TString className = key->GetClassName();
        if (!className.Contains("TObjString")) continue;
        printf("%s/%s\n", keyArrayDir->At(iDir)->GetName(), keyArrayMethod->At(iMethod)->GetName());
        TMVA::MethodBase *method = AliNDFunctionInterface::LoadMVAReader(-1, inputFile, keyArrayMethod->At(iMethod)->GetName(), keyArrayDir->At(iDir)->GetName());
        AliNDFunctionInterface::AppendMethodToArray(id, method);
      }
    }
  }
  return 0;
}


/// \brief Register TMVA method to the array at index
/// Not assumed to be used by users
/// \param index     - registered array index
/// \param method    - pointer to the method
/// \return
Int_t  AliNDFunctionInterface::AppendMethodToArray(Int_t index, TMVA::MethodBase * method){
  TObjArray * array = readerMethodBaseArray[index];
  if (array==NULL){
    array=new TObjArray(100);
    readerMethodBaseArray[index]=array;
  }
  array->AddLast(method);
  return 0;
};

/// Return statistic variable using array of MVA methods
/// \param id           - id of the TMVA method array
/// \param statType     - type of the statistic (mean, median, rms, TODO cumulant)
/// \param point        - point to evaluate
/// \return             - requested statistic
/// TODO:
/// ------------------------------
/// TODO - for the moment in evaluation we assume only one variable
///       - using the DNN or MLP we should start to support vectors
/// TODO  - extend statistic list
Double_t   AliNDFunctionInterface::EvalMVAStatArray(int id, int statType, vector<float> point){
  TObjArray * array = readerMethodBaseArray[id];
  if (array==NULL) return 0;
  Int_t size = array->GetEntries();
  if (size<=0) return 0;
  TMVA::Event  event = TMVA::Event(point,point.size());
  vector<float> value(size);
  for (Int_t i=0;i<size; i++){
    TMVA::MethodBase * method= (TMVA::MethodBase *)array->UncheckedAt(i);
    if (method==NULL) continue;                     /// some optional verbosity needed
    value[i]=method->GetRegressionValues(&event)[0];
  }
  if (statType==0) return TMath::Mean(value.size(),value.data());
  if (statType==1) return TMath::Median(value.size(),value.data());
  if (statType==2) return TMath::RMS(value.size(),value.data());
  return 0;
}

/// Register example MVA methods
void AliNDFunctionInterface::registerDefaultMVAMethods(){
  AliNDFunctionInterface::registerMethod("BDTRF25_8","!H:!V:NTrees=25:Shrinkage=0.1:UseRandomisedTrees:nCuts=20:MaxDepth=8",TMVA::Types::kBDT);
  AliNDFunctionInterface::registerMethod( "BDTRF12_16", "!H:!V:NTrees=12:Shrinkage=0.1:UseRandomisedTrees:nCuts=20:MaxDepth=16", TMVA::Types::kBDT);
  AliNDFunctionInterface::registerMethod("KNN","nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim", TMVA::Types::kKNN);
  AliNDFunctionInterface::registerMethod("MLP", "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator",TMVA::Types::kMLP);
}
