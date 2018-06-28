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


#include "TH1.h"
#include "TGraph.h"
#include "TTreeStream.h"
#include "TVirtualFFT.h"
#include "TMath.h"
#include "TVector.h"
#include "TStatToolkit.h"
#include <stdio.h>

#include "THn.h"
#include "TMVA/Factory.h"
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/DataLoader.h"
#include "TMVA/MethodBase.h"
#include "AliNDFunctionInterface.h"


///
/// Linear interpolation of the bin content
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

///
/// Linear interpolation of the numerical derivative
/// dV/dx \approx (V(i+1)-V(i-1))/(2(delta)
///
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

/// Fit MVA regression
/// \param tree            - input tree (or chain)
/// \param varFit          - colon separated variable string
/// \param cut             - selection (TTree cut)
/// \param variables       - colon separated explanatory variable list
/// \param methods         - colon separated method list - only registered methods can be used (using registerMethod(std::string method,  std::string content, TMVA::Types::EMVA id))
/// \return                - status
/// weights are saved in the root file
/// #### Example usage in test macro ( QAtrendingFitExample.C function makeMVAFits() )
///
/// #### TODOs:
/// * TODO - add weights  -equivalent to 1/error 2 - similar to AliTMinuit toolkit?
/// * TODO - user defined splitting training/test
/// * TODO  - enable support for entry list ???
/// * TODO  - Error handling: file output  currently in UPDATE mode - in case of existing entries - regression fail
/// * TODO  - implement BOOTSTRAP regression ()
/// #### Algorithm:
///  * 1.) Declare DataLoader
///  * 2.) Add the feature variables from the variables list and target variable(s)
///  * 4.) Book regression methods from the methods list
///  * 5.) Train all methods
///  * 6.) Write ascii weight files to root file
Int_t  AliNDFunctionInterface::FitMVA(TTree *tree, const char *varFit, TCut cut, const char * variables, const char *methods){
  ///  0. Declare Factory
  auto outputFile = TFile::Open("TMVA_RegressionOutput.root", "UPDATE");
  TMVA::Factory factory("TMVARegression", outputFile, "!V:!Silent:Color:DrawProgressBar:AnalysisType=Regression" );
  /// 1.) Declare DataLoader
  TMVA::DataLoader loader(TString::Format("dir_%s",varFit).Data());
  //  2.) Add the feature variables from the variables list
  TObjArray *variableList = TString(variables).Tokenize(":");
  if (variableList->GetEntries()<=0){
    ::Error("AliNDFunctionInterface::FitMVA","Empty variable list");
    return -1;
  }
  for (Int_t i=0; i<variableList->GetEntries(); i++){
    loader.AddVariable(variableList->At(i)->GetName());
    /// TODO check existence - validity of variable - in case of errror - exit with error message
  }
  delete variableList;
  loader.AddTarget(varFit);
  /// 3.) Setup DataSet
  loader.AddRegressionTree(tree, 1.0);   // link the TTree to the loader, weight for each event  = 1
  Int_t entries = tree->Draw("1",cut,"goff");
  loader.PrepareTrainingAndTestTree(cut, TString::Format("nTrain_Regression=%d:nTest_Regression=%d:SplitMode=Random:NormMode=NumEvents:!V",entries/2,entries/2));
  /// 4.) Book regression methods from the methods list
  TObjArray *methodList = TString(methods).Tokenize(":");
  if (methodList->GetEntries()<=0){
    ::Error("AliNDFunctionInterface::FitMVA","Empty method list");
    return -1;
  }
  for (Int_t i=0; i<methodList->GetEntries(); i++){
    std::string methodName=methodList->At(i)->GetName();
    printf("Booking method %s\t%d\t%s\n",methodName.data(), regressionMethodID[methodName], regressionMethodSetting[methodName].data());
    factory.BookMethod(&loader,regressionMethodID[methodName], methodName.data(), regressionMethodSetting[methodName].data());
  }
  /// 5.) Train all methods
  factory.TrainAllMethods();
  /// 6.) Write ascii weight files to root file
  for (Int_t i=0; i<methodList->GetEntries(); i++) {
    std::string methodName = methodList->At(i)->GetName();
    TString weightFile = TString::Format("dir_%s/weights/TMVARegression_%s.weights.xml", varFit,methodName.data());
    printf("Weight file\t%s\n", weightFile.Data());
    TObjString weight(gSystem->GetFromPipe(TString::Format("cat %s",weightFile.Data())));
    outputFile->cd(TString::Format("dir_%s", varFit).Data());
    weight.Write(methodName.data());
  }
  delete methodList;
  // 7.) Evaluate all MVAs using the set of test events             /// TODO - make an option
  factory.TestAllMethods();
  // 8.) Evaluate and compare performance of all configured MVAs    /// TODO - make an option
  factory.EvaluateAllMethods();
  // Save the output
  outputFile->Close();
  return 0;
}

// Load MVA regression object and register it in the map of available readers fr evaluation in TFormula
/// Current ION of TMVA does not allow standard persistency in root file
///       Way around:
///          Writing -  ASCII files are stored as TSting in the root file
///          Reading -  String written as ASCII file and read back by TMVA::Reader
/// \param id            - regression method ID  (using e.g. hash value)
/// \param inputFile     - input file  (e.g TMVAOutut.root)
/// \param method        - method name (e.g MLP)
/// \param dir           - directory usually coding regression variable with dir_ prefix (e.g. dir_meanMIPele)
/// \return              - error code
///
/// #### TODOs:
/// * TODO - ASK TMVA team to implement stadard IO
/// * TODO - ASK TMVA team to disable requirenment to set variable adresses or find switch to make it
/// #### Example parameters (see also QAtrendingFitExample.C loadMVAreaders())
/// \code
/// Int_t id=0;  const char * inputFile="TMVA_RegressionOutput.root"; const char *method="MLP"; const char *dir="dir_resolutionMIP"
/// \endcode
Int_t AliNDFunctionInterface::LoadMVAReader(Int_t id/*=0*/, const char * inputFile/*="TMVA_RegressionOutput.root"*/, const char *method/*="MLP"*/, const char *dir/*="dir_resolutionMIP"*/){
  TMVA::Reader *reader;
  reader = new TMVA::Reader( "!Color:!Silent" );
  auto fin = TFile::Open(inputFile);
  if (fin==NULL){
    ::Error("LoadMVAReader","Invalid input file\t%s", inputFile);
    return 2;
  }
  TDirectoryFile *dir0 = 0,*dirVar=0;
  fin->GetObject(dir,dir0);
  if (dir0==NULL) {
    ::Error("LoadMVAReader","Invalid input directory\t%s\t%s", dir, inputFile);
    return 4;
  }
  dir0->GetObject("InputVariables_Id",dirVar);
  for (Int_t i=0; i<dirVar->GetNkeys()-2; i++){
    TString varHis = dirVar->GetListOfKeys()->At(i)->GetName();
    TString var(varHis(0,varHis.First("__")));
    printf("%s\n",var.Data());
    Float_t value;
    reader->AddVariable(var.Data(),&value);   /// dummy booking  it is not used - TODO - ask TMVA to remove requirement
  }
  /// write weight from the root file to txt files as it is expected by reader
  FILE * pFile;
  pFile = fopen ("weights.xml","w");
  TObjString *methodWeights= NULL;
  dir0->GetObject(method,methodWeights);
  if (dir0==NULL) {
    ::Error("LoadMVAReader","Object %s not present in  input directory\t%s/%s",  method, inputFile,dir);
    delete pFile;
    return 8;
  }
  fprintf (pFile, "%s",methodWeights->GetName());
  fclose (pFile);
  reader->BookMVA( method, "weights.xml");
  AliNDFunctionInterface::readerMethodBase[id]=dynamic_cast<TMVA::MethodBase *>(reader->FindMVA(method));
  return 0;
}


/// Example MVA methods
void AliNDFunctionInterface::registerDefaultMVAMethods(){
  AliNDFunctionInterface::registerMethod("BDTRF25_8","!H:!V:NTrees=25:Shrinkage=0.1:UseRandomisedTrees:nCuts=20:MaxDepth=8",TMVA::Types::kBDT);
  AliNDFunctionInterface::registerMethod( "BDTRF12_16", "!H:!V:NTrees=12:Shrinkage=0.1:UseRandomisedTrees:nCuts=20:MaxDepth=16", TMVA::Types::kBDT);
  AliNDFunctionInterface::registerMethod("KNN","nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim", TMVA::Types::kKNN);
  AliNDFunctionInterface::registerMethod("MLP", "!H:!V:VarTransform=Norm:NeuronType=tanh:NCycles=20000:HiddenLayers=N+20:TestRate=6:TrainingMethod=BFGS:Sampling=0.3:SamplingEpoch=0.8:ConvergenceImprove=1e-6:ConvergenceTests=15:!UseRegulator",TMVA::Types::kMLP);
}

