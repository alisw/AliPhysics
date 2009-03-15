// $Id$

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Dep. of Physics                                         *
 * University of Oslo, Norway, 2007                                       *
 *                                                                        *
 * Author: Per Thomas Hille <perthi@fys.uio.no> for the ALICE HLT Project.*
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthi@fys.uio.no                                *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliHLTPHOSRcuProperties.h"
#include "AliHLTPHOSUtilities.h"



AliHLTPHOSRcuProperties::AliHLTPHOSRcuProperties() :fkEquippmentID(0),
						    fModID(2), 
						    fRcuID(0),
						    fRcuX(0), 
						    fRcuZ(0),
						    fRcuZOffset(0),
						    fRcuXOffset(0),
						    fPrintInfo(false),
						    fPrintInfoFrequncy(10000),
						    fIsSetEquippmentID(kFALSE),
						    fLog(),
						    fIsInitialized(false),
						    fUtilitiesPtr(0)
{
 
  fUtilitiesPtr = new AliHLTPHOSUtilities();

}


AliHLTPHOSRcuProperties::~AliHLTPHOSRcuProperties()
{

}


int
AliHLTPHOSRcuProperties::GetEquippmentID() const
{
  return fkEquippmentID;
}


int 
AliHLTPHOSRcuProperties::GetRCUID() const
{
  return fRcuID;
}


int
AliHLTPHOSRcuProperties::ScanArguments(int argc, const char** argv)
{
  int iResult=0; 
  //       -equipmentID
  fIsSetEquippmentID = fUtilitiesPtr->ScanSingleIntArgument(argc, argv, "-equipmentID", &fkEquippmentID);

  if(fIsSetEquippmentID == true)
    {
      //     cout << "AliHLTPHOSRcuProperties::ScanArguments fIsSetEquippmentID  == true" << endl;
    }
  else
    {
      //     cout << "AliHLTPHOSRcuProperties::ScanArguments  fIsSetEquippmentID  == false" << endl;
    }
  

  InitializeCoordinates(fkEquippmentID);
  //  fPrintInfo = fIsSetEquippmentID = fUtilitiesPtr->ScanSingleIntArgument(argc, argv, "-printinfo", &fPrintInfoFrequncy);
  fPrintInfo  = fUtilitiesPtr->ScanSingleIntArgument(argc, argv, "-printinfo", &fPrintInfoFrequncy);
  
  if(fIsSetEquippmentID == kFALSE)
    {
      fLog.Logging( kHLTLogFatal, "HLT::AliHLTPHOSRcuProperties::DoInt( int argc, const char** argv )", "Missing argument",
	       "The argument equippmentID is not set: set it with a component argumet like this: -equippmentID  <number>");
      iResult = -3; 
    }
  
  return iResult;
}


void 
AliHLTPHOSRcuProperties::InitializeCoordinates(AliHLTUInt16_t /*equippmentID*/)
{
  fRcuID =  (fkEquippmentID - 1792)%NRCUSPERMODULE;
  fModID = (fkEquippmentID  -1792 - fRcuID)/NRCUSPERMODULE;
 
  if( fRcuID  == 0)
    {
      fRcuX = 0; 
      fRcuZ = 0;
    }

  if( fRcuID  == 1)
    {
      fRcuX = 0; 
      fRcuZ = 1;
    }
 
  if( fRcuID == 2)
    {
      fRcuX = 1; 
      fRcuZ = 0;
    }

  if( fRcuID == 3)
    {
      fRcuX = 1; 
      fRcuZ = 1;
    }

  fRcuZOffset =  NZROWSRCU*fRcuZ;
  fRcuXOffset =  NXCOLUMNSRCU*fRcuX;
  fIsInitialized = true;
}

