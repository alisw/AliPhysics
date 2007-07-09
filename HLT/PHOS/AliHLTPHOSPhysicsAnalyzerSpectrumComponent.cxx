/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Øystein Djuvsland <oysteind@ift.uib.no>                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSPhysicsAnalyzerSpectrumComponent.h"
#include "AliHLTPHOSPhysicsAnalyzerPeakFitter.h"
#include "AliHLTPHOSPhysicsDefinitions.h"
#include "AliHLTPHOSPhysicsAnalyzerSpectrum.h"
#include "AliHLTPHOSPhysicsAnalyzerSpectrumComponent.h"
#include "Rtypes.h"

class AliHLTPHOSDefinitions;

const AliHLTComponentDataType AliHLTPHOSPhysicsAnalyzerSpectrumComponent::fgkInputDataTypes[]={kAliHLTVoidDataType,{0,"",""}};
UInt_t AliHLTPHOSPhysicsAnalyzerSpectrumComponent::fgCount = 0; 

AliHLTPHOSPhysicsAnalyzerSpectrumComponent gAliHLTPHOSPhysicsAnalyzerSpectrumComponent;

AliHLTPHOSPhysicsAnalyzerSpectrumComponent::AliHLTPHOSPhysicsAnalyzerSpectrumComponent():AliHLTProcessor(), fAnalyzerPtr(0), 
											       fRootHistPtr(0)
{
  //Constructor
}

AliHLTPHOSPhysicsAnalyzerSpectrumComponent::~AliHLTPHOSPhysicsAnalyzerSpectrumComponent()
{
  //Destructor
  if(fPeakFitter)
    {
      delete fPeakFitter;
      fPeakFitter = 0;
    }

  if(fAnalyzerPtr)
    {
      delete fAnalyzerPtr;
      fAnalyzerPtr = 0;
    }
  
  if(fRootHistPtr)
    {
      delete fRootHistPtr;
      fRootHistPtr = 0;
    }
      
}

AliHLTPHOSPhysicsAnalyzerSpectrumComponent::AliHLTPHOSPhysicsAnalyzerSpectrumComponent(const AliHLTPHOSPhysicsAnalyzerSpectrumComponent &):AliHLTProcessor(),
																		    fAnalyzerPtr(0),
																		    fRootHistPtr(0)
{
  //Copy constructor not implemented 
}

Int_t
AliHLTPHOSPhysicsAnalyzerSpectrumComponent::Deinit()
{
  //Deinitialize the component
  if(fPeakFitter)
    {
      fPeakFitter->SetHistogram(fRootHistPtr);
      fPeakFitter->FitLorentzian();
      delete fPeakFitter;
      fPeakFitter = 0;
    }

  if(fAnalyzerPtr)
    {
      delete fAnalyzerPtr;
      fAnalyzerPtr = 0;
    }
  
  if(fRootHistPtr)
    {
      delete fRootHistPtr;
      fRootHistPtr = 0;
    }
      
  return 0;
}

Int_t
AliHLTPHOSPhysicsAnalyzerSpectrumComponent::DoDeinit()
{
  //Deinitialize the component
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSPhysicsAnalyzerSpectrumComponent DoDeinit");

  return 0;
}


const Char_t* 
AliHLTPHOSPhysicsAnalyzerSpectrumComponent::GetComponentID()
{
  return "AliHltPhosPhysicsAnalyzerSpectrum";
}

void
AliHLTPHOSPhysicsAnalyzerSpectrumComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  //Get the input data types
  const AliHLTComponentDataType* pType=fgkInputDataTypes;
  while (pType->fID!=0) {
    list.push_back(*pType);
    pType++;
  }
}

AliHLTComponentDataType 
AliHLTPHOSPhysicsAnalyzerSpectrumComponent::GetOutputDataType()
{
  return AliHLTPHOSPhysicsDefinitions::fgkAliHLTSpectrumDataType;
}

void
AliHLTPHOSPhysicsAnalyzerSpectrumComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )
  
{
  //Get the data size of the output
  constBase = 30;
  inputMultiplier = 1;
}


Int_t 
AliHLTPHOSPhysicsAnalyzerSpectrumComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
						    AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* /*outputPtr*/, AliHLTUInt32_t& /*size*/,
						    std::vector<AliHLTComponentBlockData>& /*outputBlocks*/)
{
  //Do event
  const AliHLTComponentBlockData* iter = NULL; 
  unsigned long ndx; 
  
  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      
      if(iter->fDataType != AliHLTPHOSPhysicsDefinitions::fgkAliHLTClusterDataType)
	{
	  cout << "Warning: data type is not fgkAliHLTClusterDataType " << endl;
	  continue;
	}
      
      fClusterArrayPtr[ndx] = reinterpret_cast<AliHLTPHOSClusterDataStruct*>(iter->fPtr);
      
    } 
  
  fAnalyzerPtr->Analyze(fClusterArrayPtr, ndx);

  if(fgCount%fWriteInterval == 0 && fgCount != 0)
    {
      PushBack(fRootHistPtr, kAliHLTAnyDataType, (AliHLTUInt32_t)0);
    }

  fgCount++; 
  
  return 0;
  
}

Int_t
AliHLTPHOSPhysicsAnalyzerSpectrumComponent::DoInit(int argc, const char** argv )
{
  //Initialize the component
  Float_t firstThreshold = atof(argv[0]);
  Float_t secondThreshold = atof(argv[1]);
  fWriteInterval = atoi(argv[2]);
  Int_t nBins = atoi(argv[3]);
  Float_t lowLimit  = atof(argv[4]);
  Float_t highLimit = atof(argv[5]);

  fPeakFitter = new AliHLTPHOSPhysicsAnalyzerPeakFitter();
  fRootHistPtr = new TH1F("hist", "hist", nBins, lowLimit, highLimit);
  fAnalyzerPtr = new AliHLTPHOSPhysicsAnalyzerSpectrum();
  fAnalyzerPtr->SetThreshold(firstThreshold,secondThreshold);
  fAnalyzerPtr->SetHistogram(fRootHistPtr);

  if (argc==0 && argv==NULL) {
    // this is currently just to get rid of the warning "unused parameter"
  }
  
  return 0;
}

AliHLTComponent*
AliHLTPHOSPhysicsAnalyzerSpectrumComponent::Spawn()
{
  //Spawn a new AliHLTPHOSPhysicsAnalyzerSpectrumComponent, for the HLT framework
  return new AliHLTPHOSPhysicsAnalyzerSpectrumComponent();
}
