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

#include "AliHLTCaloDigitHandler.h"
#include "AliHLTCaloDigitDataStruct.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "AliHLTCaloGeometry.h"
#include "AliDigitNew.h"

ClassImp(AliHLTCaloDigitHandler);

AliHLTCaloDigitHandler::AliHLTCaloDigitHandler(TString detName) : AliHLTLogging(), AliHLTCaloConstantsHandler(detName.Data())
        ,fRunLoader(0)
	,fDetLoader(0)
        ,fNumberOfEvents(0)
        ,fDigitTree(0)
        ,fDigits(0)
        ,fDigitsInModule(0)
        ,fGeometry(0)
        ,fCurrentEvent(999999)
{
    // See header file for class documentation
}

AliHLTCaloDigitHandler::~AliHLTCaloDigitHandler()
{
    // See header file for class documentation
    if (fDigits)
    {
        for (Int_t m = 0; m < fCaloConstants->GetNMODULES(); m++)
        {
            delete [] fDigits[m];
        }
    }
    if (fDigitsInModule)
    {
        delete [] fDigitsInModule;
    }
}

Int_t AliHLTCaloDigitHandler::Init(AliRunLoader *runLoader)
{
    // See header file for class documentation

    fGeometry->InitialiseGeometry();

    fRunLoader = runLoader;

    if (!fRunLoader)
    {
        HLTFatal("Run loader is NULL pointer");
        return -1;
    }

    fNumberOfEvents = fRunLoader->GetNumberOfEvents();
    HLTInfo("Found %d events", fNumberOfEvents);
    
    InitDigitArray();

    return fNumberOfEvents;

}


Int_t AliHLTCaloDigitHandler::GetDigits(Int_t module, AliHLTCaloDigitDataStruct* buffer)
{
    for (Int_t iDig = 0; iDig < fDigitsInModule[module]; iDig++)
    {
        buffer[iDig] = fDigits[module][iDig];
    }
    return fDigitsInModule[module];
}


void AliHLTCaloDigitHandler::InitDigitArray()
{

    // See header file for class documentation

    fDigitsInModule = new Int_t[fCaloConstants->GetNMODULES()];

    fDigits = new AliHLTCaloDigitDataStruct*[fCaloConstants->GetNMODULES()];
    for (Int_t m = 0; m < fCaloConstants->GetNMODULES(); m++)
    {
        fDigits[m] = new AliHLTCaloDigitDataStruct[fCaloConstants->GetNXCOLUMNSMOD()*fCaloConstants->GetNZROWSMOD()];
    }
}

void AliHLTCaloDigitHandler::ResetDigitArray()
{
    // See header file for class documentation

    for (Int_t m = 0; m < fCaloConstants->GetNMODULES(); m++)
    {
        fDigitsInModule[m] = 0;
    }
}


Int_t AliHLTCaloDigitHandler::ProcessEvent(UInt_t ev)
{
    if (fCurrentEvent != ev)
    {
        ResetDigitArray();
        fCurrentEvent = ev;
        if (fRunLoader)
        {
            fRunLoader->GetEvent(ev);
            fDetLoader->LoadDigits("read");
            TTree *tree = fDetLoader->TreeD();
            if (!tree)
            {
                HLTFatal("Could not get digit tree");
                return -1;
            }
            TClonesArray *digits = 0;
            tree->SetBranchAddress(fCaloConstants->GetDETNAME(), &digits);
            tree->GetEvent(0);
            Int_t nDigits = digits->GetEntries();
	    HLTDebug("Found %d digits for branch: %s", nDigits, fCaloConstants->GetDETNAME().Data());
            for (Int_t iDig = 0; iDig < nDigits; iDig++)
            {
	      ConvertDigit(dynamic_cast<AliDigitNew*>(digits->At(iDig)));
	    }
        }
    }
    return 0;
}
