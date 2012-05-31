
/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Anders Knospe <anders.knospe@cern.ch>                         *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include "AliHLTEMCALClusterizerComponentNbyN.h"
#include "AliHLTCaloClusterizerNbyN.h"
#include "AliHLTEMCALRecoParamHandler.h"

ClassImp(AliHLTEMCALClusterizerComponentNbyN);

AliHLTEMCALClusterizerComponentNbyN gAliHLTEMCALClusterizerComponentNbyN;

AliHLTEMCALClusterizerComponentNbyN::AliHLTEMCALClusterizerComponentNbyN() : AliHLTEMCALClusterizerComponent()
{
// Constructor
}

AliHLTEMCALClusterizerComponentNbyN::~AliHLTEMCALClusterizerComponentNbyN()
{
// Destructor
}


const char* AliHLTEMCALClusterizerComponentNbyN::GetComponentID()
{
    // See header file for class documentation
    return "EmcalClusterizerNbyN";
}


int AliHLTEMCALClusterizerComponentNbyN::DoInit(int argc, const char** argv)
{
    // See header file for class documentation
    fClusterizerPtr = new AliHLTCaloClusterizerNbyN("EMCAL");

    //return AliHLTEMCALClusterizerComponent::DoInit(argc, argv);
   
    fRecoParamsPtr = new AliHLTEMCALRecoParamHandler();

    return AliHLTCaloClusterizerComponent::DoInit(argc, argv);
}

int AliHLTEMCALClusterizerComponentNbyN::ScanConfigurationArgument(int argc, const char** argv)
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
	if(tmp) 
	  {
	    tmp->SetGridDimension(argument.Atoi());
	  }
	else 
	  {
	    return -1;
	  }
        return 1;
    }
    
    return AliHLTCaloClusterizerComponent::ScanConfigurationArgument(argc, argv);
    
}

AliHLTComponent*
AliHLTEMCALClusterizerComponentNbyN::Spawn()
{
  //See headerfile for documentation

  return new AliHLTEMCALClusterizerComponentNbyN();
}
