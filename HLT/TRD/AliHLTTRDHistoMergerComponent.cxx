// $Id$
//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors:                                                       *
//*                  for The ALICE HLT Project.                            *
//*                                                                        *
//* Permission to use, copy, modify and distribute this software and its   *
//* documentation strictly for non-commercial purposes is hereby granted   *
//* without fee, provided that the above copyright notice appears in all   *
//* copies and that both the copyright notice and this permission notice   *
//* appear in the supporting documentation. The authors make no claims     *
//* about the suitability of this software for any purpose. It is          *
//* provided "as is" without express or implied warranty.                  *
//**************************************************************************

/** @file   AliHLTTRDHistoMergerComponent.cxx
    @author Theodor Rascanu
    @brief  Component for adding histos of partition wise working histo components
*/

#if __GNUC__>= 3
using namespace std;
#endif

#include "TFile.h"
#include "TString.h"
#include "TObjString.h"
#include "TObjArray.h"
#include "TH1F.h"

#include "AliHLTTRDHistoMergerComponent.h"
#include "AliHLTTRDDefinitions.h"
#include "AliHLTTRDUtils.h"

#define stringCompare(a,b) !(strcmp(a,b))

//#include "AliHLTTRD.h"
//#include <stdlib.h>
//#include <cerrno>

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTRDHistoMergerComponent)

AliHLTTRDHistoMergerComponent::AliHLTTRDHistoMergerComponent()
: AliHLTProcessor(),
  fOutputSize(500000)
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

}

AliHLTTRDHistoMergerComponent::~AliHLTTRDHistoMergerComponent()
{
  // see header file for class documentation
}

// Public functions to implement AliHLTComponent's interface.
// These functions are required for the registration process

const char* AliHLTTRDHistoMergerComponent::GetComponentID()
{
  // see header file for class documentation
  
  return "TRDHistoMerger";
}

void AliHLTTRDHistoMergerComponent::GetInputDataTypes(AliHLTComponentDataTypeList& list)
{
  // see header file for class documentation
  list.clear();
  list.push_back(kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD);
}

AliHLTComponentDataType AliHLTTRDHistoMergerComponent::GetOutputDataType()
{
  // see header file for class documentation
  return kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD;

}

void AliHLTTRDHistoMergerComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
  // see header file for class documentation
  constBase = fOutputSize;
  inputMultiplier = 0;
}

AliHLTComponent* AliHLTTRDHistoMergerComponent::Spawn()
{
  // see header file for class documentation
  return new AliHLTTRDHistoMergerComponent;
}

int AliHLTTRDHistoMergerComponent::DoInit(int argc, const char** argv)
{

  int iResult=0;
  
  TString configuration="";
  TString argument="";
  for (int i=0; i<argc && iResult>=0; i++) {
    argument=argv[i];
    if (!configuration.IsNull()) configuration+=" ";
    configuration+=argument;
  }

  if (!configuration.IsNull()) {
    iResult=Configure(configuration.Data());
  }

  for(int i=0; i<9; i++)
    fHistoArr[i]=NULL;

  for(int i=0; i<18; i++)
    fIncSM[i]=kFALSE;
 
  return 0;
}
  
int AliHLTTRDHistoMergerComponent::DoDeinit()
{
  // see header file for class documentation

  return 0;
}

int AliHLTTRDHistoMergerComponent::DoEvent(const AliHLTComponentEventData& /*evtData*/,
					    AliHLTComponentTriggerData& /*trigData*/)
{
  if(!IsDataEvent())return 0;

  int histNr = 0;
  int lastSM = -1;

  for(const TObject* iter = GetFirstInputObject(kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD);
	iter != NULL; iter = GetNextInputObject() ) {

    if(!dynamic_cast<const TH1*>(iter))
      continue;

    AliHLTUInt32_t spec = GetSpecification(iter);
    int SM = AliHLTTRDUtils::GetSM(spec);

    if(SM!=lastSM){
      if(fIncSM[SM]){
	for(int i = 0; fHistoArr[i]; i++){
	  PushBack((TObject*)fHistoArr[i], kAliHLTDataTypeHistogram | kAliHLTDataOriginTRD, 0);
	  delete fHistoArr[i];
	  fHistoArr[i] = NULL;
	}
	for(int i=0; i<18; i++)
	  fIncSM[i]=kFALSE;
      }
      lastSM = SM;
      histNr = 0;
      fIncSM[SM]=kTRUE;
    }

    if(histNr>9){
      HLTError("Got more histogramms than expected.");
      return 0;
    }

    if(!fHistoArr[histNr]) fHistoArr[histNr] = (TH1*)iter->Clone();
    else if(stringCompare(fHistoArr[histNr]->GetName(),iter->GetName()))
      fHistoArr[histNr]->Add((TH1*)iter);

    histNr++;
  }

  return 0;
}

int AliHLTTRDHistoMergerComponent::Configure(const char* arguments){
  int iResult=0;
  if (!arguments) return iResult;
  
  TString allArgs=arguments;
  TString argument;
  int bMissingParam=0;

  TObjArray* pTokens=allArgs.Tokenize(" ");
  if (pTokens) {
    for (int i=0; i<pTokens->GetEntries() && iResult>=0; i++) {
      argument=((TObjString*)pTokens->At(i))->GetString();
      if (argument.IsNull()) continue;
      
      if (argument.CompareTo("output_size")==0) {
	if ((bMissingParam=(++i>=pTokens->GetEntries()))) break;
	HLTInfo("Setting output size to: %s", ((TObjString*)pTokens->At(i))->GetString().Data());
	fOutputSize=((TObjString*)pTokens->At(i))->GetString().Atoi();
	continue;
      } 
      else {
	HLTError("unknown argument: %s", argument.Data());
	iResult=-EINVAL;
	break;
      }
    }
    delete pTokens;
  }
  if (bMissingParam) {
    HLTError("missing parameter for argument %s", argument.Data());
    iResult=-EINVAL;
  }
  return iResult;
}
