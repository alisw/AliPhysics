// @(#) $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland <oysteind@ift.uib.no>               *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** @file   AliHLTCaloDigitPublisherComponent.cxx
    @author Oystein Djuvsland
    @date
    @brief  Calo digit publisher component (input from offline).
*/

// see header file for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTCaloDigitPublisherComponent.h"
#include "AliRunLoader.h"
#include "AliLog.h"
#include "TTree.h"
#include "AliHLTTPCDigitData.h"
#include "AliHLTTPCDefinitions.h"
#include "AliHLTTPCFileHandler.h"
#include "AliHLTCaloDigitDataStruct.h"
#include "AliHLTCaloDefinitions.h"
#include "offline/AliHLTEMCALDigitHandler.h"
#include "offline/AliHLTPHOSDigitHandler.h"
#include "AliHLTEMCALDefinitions.h"
#include "AliHLTPHOSDefinitions.h"


/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTCaloDigitPublisherComponent)

AliHLTCaloDigitPublisherComponent gAliHLTCaloDigitPublisherComponent;

AliHLTCaloDigitPublisherComponent::AliHLTCaloDigitPublisherComponent() : AliHLTOfflineDataSource()
        ,fDigitHandler(0)
        ,fSpecification(0)
	,fDataType(kAliHLTAnyDataType)
        ,fModule(-1)
{
    // see header file for class documentation
    // or
    // refer to README to build package
    // or
    // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTCaloDigitPublisherComponent::~AliHLTCaloDigitPublisherComponent()
{
    // see header file for class documentation
    DoDeinit();
}

const char* AliHLTCaloDigitPublisherComponent::GetComponentID()
{
    // see header file for class documentation
    return "CaloDigitPublisher";
}

AliHLTComponentDataType AliHLTCaloDigitPublisherComponent::GetOutputDataType()
{
    return fDataType;
}

void AliHLTCaloDigitPublisherComponent::GetOutputDataSize( unsigned long& constBase, double& inputMultiplier )
{
    constBase=1000;
    inputMultiplier=10;
}

AliHLTComponent* AliHLTCaloDigitPublisherComponent::Spawn()
{
    // see header file for class documentation
    return new AliHLTCaloDigitPublisherComponent;
}

int AliHLTCaloDigitPublisherComponent::DoInit( int argc, const char** argv )
{
    // see header file for class documentation
    // scan arguments
    for (Int_t i = 0; i < argc; i++)
    {
        if (!strcmp(argv[i], "-detector"))
        {
            if (!strcmp(argv[i + 1], "EMCAL"))
            {

                fDigitHandler = AliHLTEMCALDigitHandler::Instance();
                if (fDigitHandler) fDigitHandler->Init(AliRunLoader::Instance());
                else
                {
                    HLTFatal("Could not create EMCAL digit handler");
                    return -1;
                }
                fDataType = fDigitHandler->GetDataType();
            }
            else if (!strcmp(argv[i + 1], "PHOS"))
            {
                fDigitHandler = AliHLTPHOSDigitHandler::Instance();
                if (fDigitHandler) fDigitHandler->Init(AliRunLoader::Instance());
                else
                {
                    HLTFatal("Could not create PHOS digit handler");
                    return -1;
                }
                fDataType = fDigitHandler->GetDataType();
            }
        }
        if(!strcmp(argv[i], "-module"))
	{
	  fModule = atoi(argv[i+1]);
	}
	if(!strcmp(argv[i], "-spec"))
	{
	  fSpecification = atoi(argv[i+1]);
	  HLTDebug("Data specification: 0x%x", fSpecification);
	}
    }
    return 0;
}

int AliHLTCaloDigitPublisherComponent::DoDeinit()
{
    // see header file for class documentation
    int iResult=0;
    return iResult;
}

int AliHLTCaloDigitPublisherComponent::GetEvent(const AliHLTComponentEventData& /*evtData*/,
        AliHLTComponentTriggerData& /*trigData*/,
        AliHLTUInt8_t* outputPtr,
        AliHLTUInt32_t& osize,
        vector<AliHLTComponentBlockData>& outputBlocks)
{
  
    // see header file for class documentation
    int iResult=0;
    AliHLTUInt32_t capacity=osize;
    osize=0;

    AliHLTCaloDigitDataStruct *digOut = reinterpret_cast<AliHLTCaloDigitDataStruct*>(outputPtr);
    
    // process data events only
    if (!IsDataEvent()) return 0;

    if (outputPtr==NULL || capacity==0)
    {
        HLTError("no target buffer provided");
        return -EFAULT;
    }
    
    int event=GetEventCount();
    
    fDigitHandler->ProcessEvent(event);

    Int_t nDigs = fDigitHandler->GetDigits(fModule, digOut);
    if (nDigs > 0)
    {
        AliHLTComponentBlockData bd;
        FillBlockData( bd );
        bd.fOffset = 0;
        bd.fSize = nDigs*sizeof(AliHLTCaloDigitDataStruct);
        bd.fDataType = fDataType;
        bd.fSpecification = fSpecification;
        outputBlocks.push_back(bd);
    }

    return iResult;
}
