
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
#include "AliHLTPHOSPhysicsAnalyzerSpectrum.h"
#include "AliHLTPHOSPhysicsAnalyzerPeakFitter.h"
#include "AliHLTPHOSPhysicsDefinitions.h"
#include "AliHLTPHOSDefinitions.h"
#include "TCanvas.h"
#include "TFile.h"
#include <iostream>
#include "stdio.h"			    
#include <cstdlib>

const AliHLTComponentDataType AliHLTPHOSPhysicsAnalyzerSpectrumComponent::inputDataTypes[]={kAliHLTVoidDataType,{0,"",""}};
int AliHLTPHOSPhysicsAnalyzerSpectrumComponent::fEventCount = 0; 

AliHLTPHOSPhysicsAnalyzerSpectrumComponent gAliHLTPHOSPhysicsAnalyzerSpectrumComponent;

AliHLTPHOSPhysicsAnalyzerSpectrumComponent::AliHLTPHOSPhysicsAnalyzerSpectrumComponent():AliHLTProcessor(), fAnalyzerPtr(0), 
											       fRootHistPtr(0)
{
}

AliHLTPHOSPhysicsAnalyzerSpectrumComponent::~AliHLTPHOSPhysicsAnalyzerSpectrumComponent()
{
}

AliHLTPHOSPhysicsAnalyzerSpectrumComponent::AliHLTPHOSPhysicsAnalyzerSpectrumComponent(const AliHLTPHOSPhysicsAnalyzerSpectrumComponent &):AliHLTProcessor(),
																		    fAnalyzerPtr(0),
																		    fRootHistPtr(0)
{
  
  cout << "AliHLTPHOSPhysicsAnalyzerSpectrumComponent: Copy constructor not implemented yet!" << endl;
  
}

Int_t
AliHLTPHOSPhysicsAnalyzerSpectrumComponent::Deinit()
{

 
  fPeakFitter->SetHistogram(fRootHistPtr);
  cout << "Fitting..." << endl;
  //pf->FitGaussian();
  fPeakFitter->FitLorentzian();

  return 0;
}

Int_t
AliHLTPHOSPhysicsAnalyzerSpectrumComponent::DoDeinit()
{
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
  const AliHLTComponentDataType* pType=inputDataTypes;
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
  constBase = 30;
  inputMultiplier = 1;
}


Int_t 
AliHLTPHOSPhysicsAnalyzerSpectrumComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
					AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
					std::vector<AliHLTComponentBlockData>& outputBlocks)
{

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

  if(fEventCount%fWriteInterval == 0 && fEventCount != 0)
    {
      PushBack(fRootHistPtr, kAliHLTAnyDataType, (AliHLTUInt32_t)0);
    }

  fEventCount++; 
  
  return 0;
  
}

Int_t
AliHLTPHOSPhysicsAnalyzerSpectrumComponent::DoInit(Int_t argc, const Char_t** argv )
{

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
  return new AliHLTPHOSPhysicsAnalyzerSpectrumComponent();
}
