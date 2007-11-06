/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Ãystein Djuvsland <oysteind@ift.uib.no>                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/




#include "AliHLTPHOSClusterizerComponent.h"
#include "AliHLTPHOSClusterizer.h"
//#include "AliHLTPHOSPhysicsDefinitions.h"
//#include "AliHLTPHOSDefinitions.h"
#include "AliHLTPHOSRecPointDataStruct.h"
#include "AliHLTPHOSClusterDataStruct.h"
#include "AliHLTPHOSRecPointListDataStruct.h"
#include "AliHLTPHOSDigitContainerDataStruct.h"

using namespace std;


const AliHLTComponentDataType AliHLTPHOSClusterizerComponent::fgkInputDataTypes[]=
  {
    kAliHLTVoidDataType,{0,"",""}
  };

AliHLTPHOSClusterizerComponent gAliHLTPHOSClusterizerComponent;

//AliHLTPHOSClusterizerComponent::AliHLTPHOSClusterizerComponent(): AliHLTPHOSBase(), AliHLTProcessor(), fClusterizerPtr(0), fOutPtr(0),
//								 fRecPointStructArrayPtr(0), fRecPointListPtr(0)
AliHLTPHOSClusterizerComponent::AliHLTPHOSClusterizerComponent(): AliHLTPHOSProcessor(), fClusterizerPtr(0), fOutPtr(0),
    fRecPointStructArrayPtr(0), fRecPointListPtr(0)
{
  //Constructor
}

AliHLTPHOSClusterizerComponent::~AliHLTPHOSClusterizerComponent()
{
  //Destructor

  if (fClusterizerPtr)
    {
      delete fClusterizerPtr;
      fClusterizerPtr = 0;
    }

  if (fRecPointListPtr)
    {
      delete fRecPointListPtr;
      fRecPointListPtr = 0;
    }

  if (fRecPointStructArrayPtr)
    {
      for (int i = 0; i < 1000; i++)
        {
          //	  fRecPointStructArrayPtr[i].Del();
        }
      delete fRecPointStructArrayPtr;
      fRecPointStructArrayPtr = 0;
    }

}

/*
int
AliHLTPHOSClusterizerComponent::AliHLTPHOSClusterizerComponent::Deinit()
{
  ////////// PTH WARNING you should Define a class AliHLTPHOSModuleProcessor
}
*/

// PTH AliHLTPHOSClusterizerComponent::AliHLTPHOSClusterizerComponent(const AliHLTPHOSClusterizerComponent &):AliHLTProcessor(),
//												       fClusterizerPtr(0),
//												       fOutPtr(0),
//												       fRecPointStructArrayPtr(0),
//												       fRecPointListPtr(0)
//{
//Copy constructor, not implemented
//}

int
AliHLTPHOSClusterizerComponent::Deinit()
{
  //Deinitialization

  if (fClusterizerPtr)
    {
      delete fClusterizerPtr;
      fClusterizerPtr = 0;
    }

  if (fRecPointListPtr)
    {
      delete fRecPointListPtr;
      fRecPointListPtr = 0;
    }

  for (int i = 0; i < 1000; i++)
    {
      //    fRecPointStructArrayPtr[i].Del();
    }

  if (fRecPointStructArrayPtr)
    {
      for (int i = 0; i < 1000; i++)
        {
          //	  fRecPointStructArrayPtr[i].Del();
        }
      delete fRecPointStructArrayPtr;
      fRecPointStructArrayPtr = 0;
    }

  return 0;
}

int
AliHLTPHOSClusterizerComponent::DoDeinit()
{
  //Do deinitialization
  Logging(kHLTLogInfo, "HLT", "PHOS", ",AliHLTPHOSClusterizerComponent DoDeinit");

  return 0;
}


const Char_t*
AliHLTPHOSClusterizerComponent::GetComponentID()
{
  //comment
  return "AliHltPhosClusterizer";
}

void
AliHLTPHOSClusterizerComponent::GetInputDataTypes( vector<AliHLTComponentDataType>& list)
{
  //Get datatypes for input
  const AliHLTComponentDataType* pType=fgkInputDataTypes;
  while (pType->fID!=0)
    {
      list.push_back(*pType);
      pType++;
    }
}

AliHLTComponentDataType
AliHLTPHOSClusterizerComponent::GetOutputDataType()
{
  //  return AliHLTPHOSPhysicsDefinitions::fgkAliHLTClusterDataType;
  return AliHLTPHOSDefinitions::fgkAliHLTClusterDataType;
}

void
AliHLTPHOSClusterizerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier )

{
  constBase = 30;
  inputMultiplier = 0.2;
}

int
AliHLTPHOSClusterizerComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
                                        AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
                                        std::vector<AliHLTComponentBlockData>& outputBlocks)
{
  //Do event

  UInt_t tSize            = 0;
  UInt_t offset           = 0;
  UInt_t mysize           = 0;
  Int_t nRecPoints        = 0;
  Int_t nDigits           = 0;
  Int_t index             = 0;
  Int_t j =0;

  AliHLTUInt8_t* outBPtr;
  outBPtr = outputPtr;
  const AliHLTComponentBlockData* iter = 0;
  unsigned long ndx;

  AliHLTPHOSDigitContainerDataStruct *digitContainerPtr = 0;
  
  //AliHLTPHOSRecPointContainerStruct *recPointContainerPtr = (AliHLTPHOSRecPointContainerStruct*)outBPtr;
  fClusterizerPtr->SetRecPointContainer((AliHLTPHOSRecPointContainerStruct*)outBPtr);
  for ( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      digitContainerPtr = reinterpret_cast<AliHLTPHOSDigitContainerDataStruct*>(iter->fPtr);
      if (iter->fDataType != AliHLTPHOSDefinitions::fgkAliHLTDigitDataType)
        {
          //	  cout << "Warning: data type is not fgkAliHLTDigitDataType " << endl;
          continue;
        }
      for (Int_t i = 0; i < digitContainerPtr->fNDigits; i++)
        {
	  if(fNoCrazyness && digitContainerPtr->fDigitDataStruct[i].fCrazyness)
	    continue;
	    
          fAllDigitsPtr->fDigitDataStruct[j+nDigits].fX = digitContainerPtr->fDigitDataStruct[i].fX;
          fAllDigitsPtr->fDigitDataStruct[j+nDigits].fZ = digitContainerPtr->fDigitDataStruct[i].fZ;
          fAllDigitsPtr->fDigitDataStruct[j+nDigits].fAmplitude = digitContainerPtr->fDigitDataStruct[i].fAmplitude;
          fAllDigitsPtr->fDigitDataStruct[j+nDigits].fTime = digitContainerPtr->fDigitDataStruct[i].fTime;
         // fAllDigitsPtr->fDigitDataStruct[i+nDigits].fCrazyness = digitContainerPtr->fDigitDataStruct[i].fCrazyness;
	  j++;
        }
      nDigits++;
    }

// nRecPoints = fClusterizerPtr->CreateRecPointStructArray(fRecPointStructArrayPtr, fRecPointListPtr, index);

  fOutPtr =  (AliHLTPHOSRecPointContainerStruct*)outBPtr;
  nRecPoints = fClusterizerPtr->ClusterizeEvent();
  //nRecPoints = fClusterizerPtr->CalculateCenterOfGravity(&fRecPointStructArrayPtr[i]);
  cout << "Number of clusters found: " << nRecPoints << " extracted from " << nDigits << " digits" << endl;

  mysize = 0;
  offset = tSize;

  //      fClusterizerPtr->CalculateMoments(&fRecPointStructArrayPtr[i], 0);
  //     fClusterizerPtr->ClusterizeStruct(&fRecPointStructArrayPtr[i], fOutPtr);

  mysize += sizeof(AliHLTPHOSClusterDataStruct);

  AliHLTComponentBlockData bd;
  FillBlockData( bd );
  bd.fOffset = offset;
  bd.fSize = mysize;
  // PTH      bd.fDataType = AliHLTPHOSPhysicsDefinitions::fgkAliHLTClusterDataType;
  bd.fDataType = AliHLTPHOSDefinitions::fgkAliHLTClusterDataType;
  bd.fSpecification = 0xFFFFFFFF;
  outputBlocks.push_back( bd );

  tSize += mysize;
  outBPtr += mysize;

  if ( tSize > size )
    {
      Logging( kHLTLogFatal, "HLT::AliHLTPHOSClusterizerComponent::DoEvent", "Too much data",
               "Data written over allowed buffer. Amount written: %lu, allowed amount: %lu."
               , tSize, size );
      return EMSGSIZE;
    }


  size = tSize;
// fClusterizerPtr->ResetCellEnergyArray();

  return 0;

}

int
AliHLTPHOSClusterizerComponent::DoInit(int argc, const char** argv )
{
  //Do initialization
  fAllDigitsPtr = new AliHLTPHOSDigitContainerDataStruct();
  fClusterizerPtr = new AliHLTPHOSClusterizer();
  //fClusterizerPtr->SetNoCrazyness(true);
  //
  for (int i = 0; i < argc; i++)
    {
      /*
      if(!strcmp("-threshold", argv[i]))
      fClusterizerPtr->SetThreshold(atof(argv[i+1]));
      if(!strcmp("-clusterthreshold", argv[i]))
      fClusterizerPtr->SetClusterThreshold(atof(argv[i+1]));
      if(!strcmp("-highgain", argv[i]))
      fClusterizerPtr->SetHighGainFactor(atof(argv[i+1]));
      if(!strcmp("-lowgain", argv[i]))
      fClusterizerPtr->SetLowGainFactor(atof(argv[i+1]));
      if(!strcmp("-arraysize", argv[i]))
      fClusterizerPtr->SetArraySize(atoi(argv[i+1]));*/
    }
  //  fClusterizerPtr->ResetCellEnergyArray();
  fRecPointListPtr = new AliHLTPHOSRecPointListDataStruct[N_ZROWS_MOD*N_XCOLUMNS_MOD];
  fRecPointStructArrayPtr = new AliHLTPHOSRecPointDataStruct[1000];
  for (int i = 0; i < 1000; i++)
    {
      fRecPointStructArrayPtr[i].fMultiplicity = atoi(argv[4])* atoi(argv[4]);
      // fRecPointStructArrayPtr[i].New();
    }
  /*
  printf("Clusterizer component started with:\n");
  printf(" Cell threshold:     %f\n", fClusterizerPtr->GetThreshold());
  printf(" Cluster threshold:  %f\n", fClusterizerPtr->GetClusterThreshold());
  printf(" High gain factor:   %f\n", fClusterizerPtr->GetHighGainFactor());
  printf(" Low gain factor:    %f\n", fClusterizerPtr->GetLowGainFactor());
  printf(" Cluster array size: %d\n\n", fClusterizerPtr->GetArraySize());
  */
  return 0;
}

AliHLTComponent*
AliHLTPHOSClusterizerComponent::Spawn()
{
  //Spawn a new AliHLTPHOSClusterizerComponent, for HLT framework
  return new AliHLTPHOSClusterizerComponent();
}
