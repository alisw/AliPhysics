 // $Id$

 /**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTEMCALDigitMakerComponent.h"
#include "AliHLTCaloDigitMaker.h"
#include "AliHLTCaloDigitDataStruct.h"
#include "AliHLTCaloChannelDataHeaderStruct.h"
#include "AliHLTCaloChannelDataStruct.h"
#include "AliHLTEMCALMapper.h"
#include "AliHLTEMCALDefinitions.h"
#include "AliCaloCalibPedestal.h"
#include "AliEMCALCalibData.h"
#include "AliCDBEntry.h"
#include "AliCDBPath.h"
#include "AliCDBManager.h"
#include "TFile.h"
#include <sys/stat.h>
#include <sys/types.h>


/** 
 * @file   AliHLTEMCALDigitMakerComponent.cxx
 * @author Oystein Djuvsland
 * @date   
 * @brief  A digit maker component for EMCAL HLT
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


ClassImp(AliHLTEMCALDigitMakerComponent)

AliHLTEMCALDigitMakerComponent gAliHLTEMCALDigitMakerComponent;

AliHLTEMCALDigitMakerComponent::AliHLTEMCALDigitMakerComponent() :
  AliHLTCaloProcessor(),
  AliHLTCaloConstantsHandler("EMCAL"),
  fDigitMakerPtr(0),
  fDigitContainerPtr(0),
  fPedestalData(0),
  fCalibData(0),
  fBCMInitialised(true),
  fGainsInitialised(true)
{
  //see header file for documentation
}


AliHLTEMCALDigitMakerComponent::~AliHLTEMCALDigitMakerComponent()
{
  //see header file for documentation
}

int 
AliHLTEMCALDigitMakerComponent::Deinit()
{ 
  //see header file for documentation
  if(fDigitMakerPtr)
    {
      delete fDigitMakerPtr;
      fDigitMakerPtr = 0;
    }
  return 0;
}

const char*
AliHLTEMCALDigitMakerComponent::GetComponentID()
{
  //see header file for documentation
  return "EmcalDigitMaker";
}


void
AliHLTEMCALDigitMakerComponent::GetInputDataTypes(vector<AliHLTComponentDataType>& list)
{ 
  //see header file for documentation
  list.clear();
  list.push_back(AliHLTEMCALDefinitions::fgkChannelDataType);
}

AliHLTComponentDataType 
AliHLTEMCALDigitMakerComponent::GetOutputDataType()
{
  //see header file for documentation
//  return AliHLTCaloDefinitions::fgkDigitDataType|kAliHLTDataOriginEMCAL;
  return AliHLTEMCALDefinitions::fgkDigitDataType;
}


void 
AliHLTEMCALDigitMakerComponent::GetOutputDataSize(unsigned long& constBase, double& inputMultiplier)
{
  //see header file for documentation
  constBase = 0;
  inputMultiplier = (float)sizeof(AliHLTCaloDigitDataStruct)/sizeof(AliHLTCaloChannelDataStruct) + 1;
}

int 
AliHLTEMCALDigitMakerComponent::DoEvent(const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
					AliHLTComponentTriggerData& /*trigData*/, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size,
					std::vector<AliHLTComponentBlockData>& outputBlocks)
{
  //see header file for documentation
  UInt_t offset           = 0; 
  UInt_t mysize           = 0;
  Int_t digitCount        = 0;
  Int_t ret               = 0;

  AliHLTUInt8_t* outBPtr;
  outBPtr = outputPtr;
  const AliHLTComponentBlockData* iter = 0; 
  unsigned long ndx; 

  UInt_t specification = 0;
  AliHLTCaloChannelDataHeaderStruct* tmpChannelData = 0;
  
  //  fDigitMakerPtr->SetDigitHeaderPtr(reinterpret_cast<AliHLTCaloDigitHeaderStruct*>(outputPtr));

  fDigitMakerPtr->SetDigitDataPtr(reinterpret_cast<AliHLTCaloDigitDataStruct*>(outputPtr));

  for( ndx = 0; ndx < evtData.fBlockCnt; ndx++ )
    {
      iter = blocks+ndx;
      
      if(iter->fDataType != AliHLTEMCALDefinitions::fgkChannelDataType)
	{
	  continue;
	}
      if(!fBCMInitialised)
      {
	 AliHLTEMCALMapper mapper(iter->fSpecification);
	 Int_t module = mapper.GetModuleFromSpec(iter->fSpecification);
	 for(Int_t x = 0; x < fCaloConstants->GetNXCOLUMNSMOD(); x++)
	 {
	    for(Int_t z = 0; z < fCaloConstants->GetNZROWSMOD(); z++)
	    {
	       fDigitMakerPtr->SetBadChannel(x, z, fPedestalData->IsBadChannel(module, z+1, x+1));
	    }
	 }
	 //delete fBadChannelMap;
	 fBCMInitialised = true;
      }
      if(!fGainsInitialised)
      {
	 AliHLTEMCALMapper mapper(iter->fSpecification);;
	 Int_t module = mapper.GetModuleFromSpec(iter->fSpecification);
	 for(Int_t x = 0; x < fCaloConstants->GetNXCOLUMNSMOD(); x++)
	 {
	    for(Int_t z = 0; z < fCaloConstants->GetNZROWSMOD(); z++)
	    {
	       module = 0; //removing warning
		//fDigitMakerPtr->SetGain(x, z, fCalibData->GE(module, z+1, x+1), fCalibData->GetADCchannelEmc(module, z+1, x+1));
	    }
	 }
	 fGainsInitialised = true;
      }
      specification |= iter->fSpecification;
      tmpChannelData = reinterpret_cast<AliHLTCaloChannelDataHeaderStruct*>(iter->fPtr);
    
      ret = fDigitMakerPtr->MakeDigits(tmpChannelData, size-(digitCount*sizeof(AliHLTCaloDigitDataStruct)));
      if(ret == -1) 
	{
	  HLTError("Trying to write over buffer size");
	  return -ENOBUFS;
	}
      digitCount += ret; 
    }
  
  mysize += digitCount*sizeof(AliHLTCaloDigitDataStruct);

  HLTDebug("# of digits: %d, used memory size: %d, available size: %d", digitCount, mysize, size);

  if(mysize > 0) 
    {
      AliHLTComponentBlockData bd;
      FillBlockData( bd );
      bd.fOffset = offset;
      bd.fSize = mysize;
      bd.fDataType = AliHLTEMCALDefinitions::fgkDigitDataType;
      bd.fSpecification = specification;
      outputBlocks.push_back(bd);
    }

  fDigitMakerPtr->Reset();

  size = mysize; 

  return 0;
}


int
AliHLTEMCALDigitMakerComponent::DoInit(int argc, const char** argv )
{
  //see header file for documentation

  fDigitMakerPtr = new AliHLTCaloDigitMaker("EMCAL");

  AliHLTCaloMapper *mapper = new AliHLTEMCALMapper(2);
  fDigitMakerPtr->SetMapper(mapper);
  
  for(int i = 0; i < argc; i++)
    {
      if(!strcmp("-lowgainfactor", argv[i]))
	{
	  fDigitMakerPtr->SetGlobalLowGainFactor(atof(argv[i+1]));
	}
      if(!strcmp("-highgainfactor", argv[i]))
	{
	  fDigitMakerPtr->SetGlobalHighGainFactor(atof(argv[i+1]));
	}
    }
 GetBCMFromCDB();
  //fDigitMakerPtr->SetDigitThreshold(2);

  return 0;
}

int AliHLTEMCALDigitMakerComponent::GetBCMFromCDB()
{
   // See header file for class documentation
   
     fBCMInitialised = false;
   
//   HLTInfo("Getting bad channel map...");

  AliCDBPath path("EMCAL","Calib","Pedestals");
  if(path.GetPath())
    {
      //      HLTInfo("configure from entry %s", path.GetPath());
      AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
	if (pEntry) 
	{
	    fPedestalData = (AliCaloCalibPedestal*)pEntry->GetObject();
	}
      else
	{
//	    HLTError("can not fetch object \"%s\" from CDB", path);
	    return -1;
	}
    }
   if(!fPedestalData) 
   {
	return -1;
   }

   return 0;
}

int AliHLTEMCALDigitMakerComponent::GetGainsFromCDB()
{
   // See header file for class documentation
   
     fGainsInitialised = false;
   
//   HLTInfo("Getting bad channel map...");

  AliCDBPath path("EMCAL","Calib","Data");
  if(path.GetPath())
    {
      //      HLTInfo("configure from entry %s", path.GetPath());
      AliCDBEntry *pEntry = AliCDBManager::Instance()->Get(path/*,GetRunNo()*/);
	if (pEntry) 
	{
	    fCalibData = (AliEMCALCalibData*)pEntry->GetObject();
	}
      else
	{
//	    HLTError("can not fetch object \"%s\" from CDB", path);
	    return -1;
	}
    }
   if(!fCalibData) return -1;



   return 0;
}



AliHLTComponent*
AliHLTEMCALDigitMakerComponent::Spawn()
{
  //see header file for documentation
  return new AliHLTEMCALDigitMakerComponent();
}
