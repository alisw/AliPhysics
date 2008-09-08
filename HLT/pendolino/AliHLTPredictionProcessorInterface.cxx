// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Sebastian Bablok <Sebastian.Bablok@ift.uib.no>        *
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

/** @file   AliHLTPredictionProcessorInterface.cxx
    @author Sebastian Bablok
    @date   
    @brief  
*/

#include "AliHLTPredictionProcessorInterface.h"
#include "AliHLTPendolino.h"
//#include "AliShuttleInterface.h"


ClassImp(AliHLTPredictionProcessorInterface)

AliHLTPredictionProcessorInterface::AliHLTPredictionProcessorInterface(
			const char* detector, AliHLTPendolino* pendolino) : 
			AliPreprocessor(detector, reinterpret_cast<AliShuttleInterface*>
					(pendolino)), fpPend(pendolino) {
}


AliHLTPredictionProcessorInterface::~AliHLTPredictionProcessorInterface() {

}

Int_t AliHLTPredictionProcessorInterface::GetRunNumber() {
	return fpPend->GetRunNumber();
}


Bool_t AliHLTPredictionProcessorInterface::includeAliCDBEntryInList(
            const TString& entryPath) {

    return fpPend->includeAliCDBEntryInList(entryPath);
}



