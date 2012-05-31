// $Id$

//**************************************************************************
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//*                                                                        *
//* Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
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

/** @file   AliHLTOUTHandlerDetectorDDL.cxx
    @author Matthias Richter
    @date   2008-09-09
    @brief  HLTOUT handler returning equipment id from data type and spec.
*/

#include "AliHLTOUTHandlerDetectorDDL.h"
#include "AliHLTOUT.h"
#include "AliHLTDAQ.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTOUTHandlerDetectorDDL)

AliHLTOUTHandlerDetectorDDL::AliHLTOUTHandlerDetectorDDL(const char* detector, AliHLTComponentDataType dt)
  :
  fDDLOffset(-1),
  fNumberOfDDLs(-1),
  fDt(dt)
{ 
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
  
  fDDLOffset=AliHLTDAQ::DdlIDOffset(detector);
  fNumberOfDDLs=AliHLTDAQ::NumberOfDdls(detector);
}

AliHLTOUTHandlerDetectorDDL::~AliHLTOUTHandlerDetectorDDL()
{
  // see header file for class documentation
}

int AliHLTOUTHandlerDetectorDDL::ProcessData(AliHLTOUT* pData)
{
  // see header file for class documentation
  if (!pData) return -EINVAL;
  if (fDDLOffset<0 || fNumberOfDDLs<0) return -ENODEV;

  static int errorCount=0;
  const int maxErrorCount=10;
  AliHLTComponentDataType dt=kAliHLTVoidDataType;
  AliHLTUInt32_t spec=kAliHLTVoidDataSpec;
  int iResult=pData->GetDataBlockDescription(dt, spec);
  if (iResult>=0 && dt!=kAliHLTVoidDataType && spec!=kAliHLTVoidDataSpec) {
    if (dt==fDt) {
      int ddlNo=0;
      for (;ddlNo<32 && ddlNo<fNumberOfDDLs; ddlNo++) {
	if (spec&(0x1<<ddlNo)) break;
      }
      if (ddlNo>=32 || ddlNo>=fNumberOfDDLs) {
	HLTError("invalid specification 0x%08x: can not extract DDL id for data block %s", spec, AliHLTComponent::DataType2Text(dt).c_str());
	iResult=-ENODEV;
      } else if (spec^(0x1<<ddlNo)) {
	iResult=-EEXIST;
	HLTError("multiple links set in specification 0x%08x: can not extract DDL id for data block %s", spec, AliHLTComponent::DataType2Text(dt).c_str());
      } else {
	iResult=fDDLOffset+ddlNo;
      }
    } else {
      if (errorCount++<10) {
	HLTError("wrong data type: expecting %s, got %s%s",
		 AliHLTComponent::DataType2Text(fDt).c_str(),
		 AliHLTComponent::DataType2Text(dt).c_str(),
		 errorCount==maxErrorCount?"; suppressing further error messages":"");
      }
      iResult=-EBADF;
    }
  } else {
    if (errorCount++<10) {
      HLTError("can not get a valid data type and specification from HLTOUT: type %s, specification 0x%08x%s",
	       AliHLTComponent::DataType2Text(dt).c_str(), spec,
	       errorCount==maxErrorCount?"; suppressing further error messages":"");
    }
    iResult=-ENODATA;
  }
  return iResult;
}
