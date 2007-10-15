
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

ClassDef(AliHLTPHOSEnergyEqualizer);

using namespace PhosHLTConst;

AliHLTPHOSEnergyEqualizer::AliHLTPHOSEnergyEqualizer() : 
  AliHLTPHOSBase(),
  fDigitContainerPtr(0),
  fGlobalHighGainFactor(0),
  fGlobalLowGainFactor(0),
  fGlobalConversion(true)
{

}

Int_t
AliHLTPHOSEnergyEqualizer::MakeDigits(AliHLTPHOSRcuCellEnergyDataStruct *cellEnergies)
{
  
  Int_t nDigits = fDigitDataStruct-fNDigits;
  if(fGlobalConversion)
    {
      for(Int_t i = 0; i < fCnt; i++)
	{
	  if(
	  fDigitContainerPtr->fDigitDataStruct[nDigits].fModule = cellEnergies->fValidData.fModuleID;
	  fDigitContainerPtr->fDigitDataStruct[nDigits].fX = cellEnergies->fValidData.fX + (cellEnergies->fRcuX)*N_XCOLUMNS_RCU;
	  fDigitContainerPtr->fDigitDataStruct[nDigits].fZ = cellEnergies->fValidData.fZ + (cellEnergies->fRcuZ)*N_ZROWS_RCU;
	}
    }
  else
    //not implemented yet
