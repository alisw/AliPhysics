
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Oystein Djuvsland <oysteind@ift.uib.no>                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include "AliHLTPHOSClusterizerComponentNbyN.h"
#include "AliHLTCaloClusterizerNbyN.h"

ClassImp(AliHLTPHOSClusterizerComponentNbyN);

AliHLTPHOSClusterizerComponentNbyN gAliHLTPHOSClusterizerComponentNbyN;

AliHLTPHOSClusterizerComponentNbyN::AliHLTPHOSClusterizerComponentNbyN() : AliHLTPHOSClusterizerComponent()
{
// Constructor
}

AliHLTPHOSClusterizerComponentNbyN::~AliHLTPHOSClusterizerComponentNbyN()
{
// Destructor
}


const char* AliHLTPHOSClusterizerComponentNbyN::GetComponentID()
{
    // See header file for class documentation
    return "PhosClusterizerNbyN";
}


int AliHLTPHOSClusterizerComponentNbyN::DoInit(int argc, const char** argv)
{
    // See header file for class documentation
    fClusterizerPtr = new AliHLTCaloClusterizerNbyN("PHOS");

    return AliHLTPHOSClusterizerComponent::DoInit(argc, argv);
}

int AliHLTPHOSClusterizerComponentNbyN::ScanConfigurationArgument(int argc, const char** argv)
{
// See header file for class documentation
   
   if (argc <= 0) return 0;

    int i=0;

    TString argument=argv[i];

    if (argument.CompareTo("-gridsize") == 0)
    {
        if (++i >= argc) return -EPROTO;
        argument = argv[i];
	AliHLTCaloClusterizerNbyN* tmp = dynamic_cast<AliHLTCaloClusterizerNbyN*>(fClusterizerPtr);
        tmp->SetGridDimension(argument.Atoi());
        return 1;
    }
    
    return AliHLTCaloClusterizerComponent::ScanConfigurationArgument(argc, argv);
    
}

AliHLTComponent*
AliHLTPHOSClusterizerComponentNbyN::Spawn()
{
  //See headerfile for documentation

  return new AliHLTPHOSClusterizerComponentNbyN();
}
