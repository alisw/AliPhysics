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
  nInputFeatureActive(0),
  fUsedVars(new TBits(AliDielectronVarManager::kNMaxValues)),
  mvaCutValue(0.)
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
	     nInputFeatureActive(0),
	     fUsedVars(new TBits(AliDielectronVarManager::kNMaxValues)),
	     mvaCutValue(0.)
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

  // specifying location of input feature in AliDielectronVarManager::ValueTypes
  inputFeatureNumber[nInputFeatureActive] = dielectronVar;

  // setting bit map of used features
  fUsedVars->SetBitNumber((UInt_t)dielectronVar,kTRUE);

  // adding feature value to TMVA reader
  TMVAReader->AddVariable(featureName.Data(), &inputFeature[nInputFeatureActive]);

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

  // specifying location of input feature in AliDielectronVarManager::ValueTypes
  inputFeatureNumber[nInputFeatureActive] = dielectronVar;
  
  // setting bit map of used features
  fUsedVars->SetBitNumber((UInt_t)dielectronVar,kTRUE);

  // adding spectator value to TMVA reader
  TMVAReader->AddSpectator(featureName.Data(), &inputFeature[nInputFeatureActive]);

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
  
  TMVAReader->BookMVA(TMVAName.Data(),weightName.Data());
  TMVAReaderName = TMVAName;
}

//______________________________________________
Bool_t AliDielectronTMVACuts::IsSelected(TObject* track)
{
  //
  // Apply configured cuts
  //

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

  // set input features
  Double_t values[AliDielectronVarManager::kNMaxValues];
  AliDielectronVarManager::SetFillMap(fUsedVars);
  AliDielectronVarManager::Fill(track, values);
  
  for(Int_t i = 0; i < nInputFeatureActive; i++){
    inputFeature[i] = values[inputFeatureNumber[i]];
  }

  // evaluate MVA output value
  Float_t mvaOutput = (Float_t)TMVAReader->EvaluateMVA(TMVAReaderName.Data());

  // check if above cut value
  if(mvaOutput < mvaCutValue){
    accept = kFALSE;
  }

  return accept;
}
