/*************************************************************************
* Copyright(c) 1998-2009, ALICE Experiment at CERN, All rights reserved. *
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
//                Dielectron TMVACuts                                  //
//                                                                       //
//                                                                       //
/*
Detailed description


*/
//                                                                       //
///////////////////////////////////////////////////////////////////////////

#include "TSystem.h"

#include "TMVA/Tools.h"
#include "TMVA/Reader.h"

#include "AliLog.h"
#include "AliDielectronVarManager.h"
#include "AliDielectronTMVACuts.h"

ClassImp(AliDielectronTMVACuts)

AliDielectronTMVACuts::AliDielectronTMVACuts() :
AliAnalysisCuts(),
  TMVAReader(0),
  TMVAReaderName(""),
  TMVAWeightPathName(""),
  TMVAWeightFileName(""),
  fUsedVars(new TBits(AliDielectronVarManager::kNMaxValues)),
  fIsSpectator(new TBits(nInputFeatureMax)),
  nInputFeatureActive(0),
  mvaCutValue(0.),
  isInitialized(kFALSE)
{
  //
  // Default Constructor
  //

  // --- Create the Reader object
  TMVAReader = new TMVA::Reader( "!Color:!Silent" );

  // --- initialize input feature array
  for(Int_t i = 0; i < nInputFeatureMax; i++){
    inputFeatureNumber[i] = AliDielectronVarManager::kPx;
    inputFeature[i]       = 0.;
  }
}

//______________________________________________
AliDielectronTMVACuts::AliDielectronTMVACuts(const char* name, const char* title) :
	     AliAnalysisCuts(name, title),
	     TMVAReader(0),
	     TMVAReaderName(""),
	     TMVAWeightPathName(""),
	     TMVAWeightFileName(""),
	     fUsedVars(new TBits(AliDielectronVarManager::kNMaxValues)),
	     fIsSpectator(new TBits(nInputFeatureMax)),
	     nInputFeatureActive(0),
	     mvaCutValue(0.),
	     isInitialized(kFALSE)
{
  //
  // Named Constructor
  //

  // --- Create the Reader object
  TMVAReader = new TMVA::Reader( "!Color:!Silent" );

  // --- initialize input feature array
  for(Int_t i = 0; i < nInputFeatureMax; i++){
    inputFeatureNumber[i] = AliDielectronVarManager::kPx;
    inputFeature[i]       = 0.;
  }
}

//______________________________________________
AliDielectronTMVACuts::~AliDielectronTMVACuts()
{
  //
  // Default Destructor
  //
  if (TMVAReader) delete TMVAReader;
  if (fUsedVars)  delete fUsedVars;

}

//______________________________________________
void AliDielectronTMVACuts::InitTMVAReader()
{
  //
  // TMVA reader initialization (needed for GRID running)
  //

  AliInfo("Initialize TMVA Reader");

  if(!TMVAReader){
    AliFatal("TMVA reader not available");
    return;
  }

  for(Int_t i = 0; i < nInputFeatureActive; i++){

    // adding spectator value to TMVA reader
    if(fIsSpectator->TestBitNumber(i)){
      AliInfo(Form("Input spectator %d = %s (%s)",i,inputFeatureName[i].Data(),AliDielectronVarManager::GetValueName(inputFeatureNumber[i])));
      TMVAReader->AddSpectator(inputFeatureName[i].Data(), &inputFeature[i]);
    }
    // adding feature value to TMVA reader
    else{
      AliInfo(Form("Input feature %d = %s (%s)",i,inputFeatureName[i].Data(),AliDielectronVarManager::GetValueName(inputFeatureNumber[i])));
      TMVAReader->AddVariable(inputFeatureName[i].Data(), &inputFeature[i]);
    }
  }

  // initialize reader (copy weight file and add weight file name)
  AliInfo(Form("Initialize TMVA reader %s with weight file %s from path %s",TMVAReaderName.Data(),TMVAWeightFileName.Data(),TMVAWeightPathName.Data()));
  gSystem->Exec(Form("alien_cp %s/%s .",TMVAWeightPathName.Data(),TMVAWeightFileName.Data()));
  TMVAReader->BookMVA(TMVAReaderName.Data(),TMVAWeightFileName.Data());
  
  // set to initialized
  isInitialized = kTRUE;
}


//______________________________________________
void AliDielectronTMVACuts::AddTMVAInput(TString featureName, AliDielectronVarManager::ValueTypes dielectronVar)
{
  //
  // Adding an input feature for the TMVA reader
  // - featureName: name of the feature during training/in the weight file
  // - dielectronVar: variable number in the AliDielectronVarManager, e.g. AliDielectronVarManager::kTPCnSigmaEle
  //

  if(!TMVAReader){
    AliFatal("TMVA reader not available");
    return;
  }
    
  if(featureName == "" || dielectronVar < 0 || dielectronVar > AliDielectronVarManager::kNMaxValues){
    AliWarning(Form("Missing feature name (%s) or invalid dielectron variable number (%d) --> Do not add input",featureName.Data(),dielectronVar));
    return;
  }

  if(nInputFeatureActive >= nInputFeatureMax){
    AliWarning(Form("Too many input features defined --> Do not add input"));
    return;
  }

  // specifying name of input feature in weight file
  inputFeatureName[nInputFeatureActive]   = featureName;
  
  // specifying location of input feature in AliDielectronVarManager::ValueTypes
  inputFeatureNumber[nInputFeatureActive] = dielectronVar;

  // setting bit map of used features
  fUsedVars->SetBitNumber((UInt_t)dielectronVar,kTRUE);

  // count input features
  nInputFeatureActive++;
  
}

//______________________________________________
void AliDielectronTMVACuts::AddTMVASpectator(TString featureName, AliDielectronVarManager::ValueTypes dielectronVar)
{
  //
  // Adding an input feature spectator (written in weight file, but not used for decision) for the TMVA reader
  // - featureName: name of the feature during training/in the weight file
  // - dielectronVar: variable number in the AliDielectronVarManager, e.g. AliDielectronVarManager::kTPCnSigmaEle
  //

  if(!TMVAReader){
    AliFatal("TMVA reader not available");
    return;
  }
    
  if(featureName == "" || dielectronVar < 0 || dielectronVar > AliDielectronVarManager::kNMaxValues){
    AliWarning(Form("Missing feature name (%s) or invalid dielectron variable number (%d) --> Do not add input",featureName.Data(),dielectronVar));
    return;
  }

  if(nInputFeatureActive >= nInputFeatureMax){
    AliWarning(Form("Too many input features defined --> Do not add input"));
    return;
  }

  // specifying name of input feature in weight file
  inputFeatureName[nInputFeatureActive]   = featureName;

  // specifying location of input feature in AliDielectronVarManager::ValueTypes
  inputFeatureNumber[nInputFeatureActive] = dielectronVar;
  
  // setting bit map of used features
  fUsedVars->SetBitNumber((UInt_t)dielectronVar,kTRUE);

  // setting bit map of spectators in input feature array
  fIsSpectator->SetBitNumber(nInputFeatureActive,kTRUE);
  
  // count input features
  nInputFeatureActive++;
  
}

//______________________________________________
void AliDielectronTMVACuts::SetTMVAWeights(TString TMVAName, TString weightName)
{
  //
  // Setting the weight file for the TMVA reader
  // - TMVAName: name of the TMVA object
  // - weightName: name of the weight xml file
  //

  if(!TMVAReader){
    AliFatal("TMVA reader not available");
    return;
  }

  if(TMVAName == "" || weightName == ""){
    AliWarning(Form("Missing TMVA name (%s) or weight input file name (%s) --> Do not add input",TMVAName.Data(),weightName.Data()));
    return;
  }

  if(weightName.Contains("alien://")){
    AliInfo(Form("Use TMVA weight input file from Alien: %s",weightName.Data()));
  }
  else{
    AliWarning(Form("TMVA weight input file not from Alien? Not supported at the moment: %s",weightName.Data()));
    return;
  }

  // set name of the reader
  TMVAReaderName = TMVAName;

  // set location of input weight file
  TObjArray* strings = weightName.Tokenize("/");  
  if(strings->GetEntriesFast()) {
    TIter iString(strings);
    TObjString* oString = NULL;
    Int_t i = 0;
    TMVAWeightPathName = "alien:///";
    
    while ((oString = (TObjString*)iString()) && i < strings->GetEntriesFast()-1) {

      // need this due to special treatment of first token
      if(!oString->GetString().Contains("alien"))
	TMVAWeightPathName.Append(Form("%s/",oString->GetString().Data()));
      i++;
    }
    
    TMVAWeightFileName = oString->GetString().Data();

  }
}

//______________________________________________
Bool_t AliDielectronTMVACuts::IsSelected(TObject* track)
{
  //
  // Apply configured cuts
  //

  // first check if TMVA reader is initialized (needed for running on GRID)
  if(!isInitialized)
    InitTMVAReader();

  AliVTrack *vtrack=dynamic_cast<AliVTrack*>(track);
  if (!vtrack) return kFALSE;
  
  Bool_t accept=kTRUE;

  if(!TMVAReader){
    AliFatal("TMVA reader not available");
    return kFALSE;
  }

  if(TMVAReaderName == ""){
    AliWarning("TMVA Reader not set up properly");
    return kFALSE;
  }

  // // Just testing....
  // TMVAReader->BookMVA(TMVAReaderName.Data(),"TMVA_e_0_05_05_05_electron0555first.weights.xml");


  // set input features
  Double_t values[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::SetFillMap(fUsedVars);
  AliDielectronVarManager::Fill(track, values);
  
  for(Int_t i = 0; i < nInputFeatureActive; i++){
    inputFeature[i] = (Float_t) values[inputFeatureNumber[i]];
  }

  // evaluate MVA output value
  Float_t mvaOutput = (Float_t)TMVAReader->EvaluateMVA(TMVAReaderName.Data());

  // check if above cut value
  if(mvaOutput < mvaCutValue){
    accept = kFALSE;
  }

  return accept;
}
