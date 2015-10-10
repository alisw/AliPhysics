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

/** @file   AliHLTDAQInterfaceImplementation.cxx
    @author Matthias Richter
    @date   
    @brief  Interface to the AliDAQ class
*/

#include "AliHLTDAQInterfaceImplementation.h"
#include "AliDAQ.h"

/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTDAQInterfaceImplementation)

AliHLTDAQInterfaceImplementation::AliHLTDAQInterfaceImplementation()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTDAQInterfaceImplementation::~AliHLTDAQInterfaceImplementation()
{
  // see header file for class documentation
}

Int_t       AliHLTDAQInterfaceImplementation::VirtNumberOfDetectors()
{
  // see header file for class documentation
  return AliDAQ::kNDetectors;
}

Int_t       AliHLTDAQInterfaceImplementation::VirtHLTId()
{
  // see header file for class documentation
  return AliDAQ::kHLTId;
}

Int_t       AliHLTDAQInterfaceImplementation::VirtDetectorID(const char *detectorName)
{
  // see header file for class documentation
  return AliDAQ::DetectorID(detectorName);
}

const char *AliHLTDAQInterfaceImplementation::VirtDetectorName(Int_t detectorID)
{
  // see header file for class documentation
  return AliDAQ::DetectorName(detectorID);
}


Int_t       AliHLTDAQInterfaceImplementation::VirtDdlIDOffset(const char *detectorName)
{
  // see header file for class documentation
  return AliDAQ::DdlIDOffset(detectorName);
}

Int_t       AliHLTDAQInterfaceImplementation::VirtDdlIDOffset(Int_t detectorID)
{
  // see header file for class documentation
  return AliDAQ::DdlIDOffset(detectorID);
}

const char *AliHLTDAQInterfaceImplementation::VirtDetectorNameFromDdlID(Int_t ddlID, Int_t &ddlIndex)
{
  // see header file for class documentation
  return AliDAQ::DetectorNameFromDdlID(ddlID, ddlIndex);
}

Int_t       AliHLTDAQInterfaceImplementation::VirtDetectorIDFromDdlID(Int_t ddlID, Int_t &ddlIndex)
{
  // see header file for class documentation
  return AliDAQ::DetectorIDFromDdlID(ddlID, ddlIndex);
}

Int_t       AliHLTDAQInterfaceImplementation::VirtDdlID(const char *detectorName, Int_t ddlIndex)
{
  // see header file for class documentation
  return AliDAQ::DdlID(detectorName, ddlIndex);
}

Int_t       AliHLTDAQInterfaceImplementation::VirtDdlID(Int_t detectorID, Int_t ddlIndex)
{
  // see header file for class documentation
  return AliDAQ::DdlID(detectorID, ddlIndex);
}

const char *AliHLTDAQInterfaceImplementation::VirtDdlFileName(const char *detectorName, Int_t ddlIndex)
{
  // see header file for class documentation
  return AliDAQ::DdlFileName(detectorName, ddlIndex);
}

const char *AliHLTDAQInterfaceImplementation::VirtDdlFileName(Int_t detectorID, Int_t ddlIndex)
{
  // see header file for class documentation
  return AliDAQ::DdlFileName(detectorID, ddlIndex);
}

Int_t       AliHLTDAQInterfaceImplementation::VirtNumberOfDdls(const char *detectorName)
{
  // see header file for class documentation
  return AliDAQ::NumberOfDdls(detectorName);
}

Int_t       AliHLTDAQInterfaceImplementation::VirtNumberOfDdls(Int_t detectorID)
{
  // see header file for class documentation
  return AliDAQ::NumberOfDdls(detectorID);
}

const char *AliHLTDAQInterfaceImplementation::VirtListOfTriggeredDetectors(UInt_t detectorPattern)
{
  // see header file for class documentation
  return AliDAQ::ListOfTriggeredDetectors(detectorPattern);
}

UInt_t      AliHLTDAQInterfaceImplementation::VirtDetectorPattern(const char *detectorList)
{
  // see header file for class documentation
  return AliDAQ::DetectorPattern(detectorList);
}

const char *AliHLTDAQInterfaceImplementation::VirtOfflineModuleName(const char *detectorName)
{
  // see header file for class documentation
  return AliDAQ::OfflineModuleName(detectorName);
}

const char *AliHLTDAQInterfaceImplementation::VirtOfflineModuleName(Int_t detectorID)
{
  // see header file for class documentation
  return AliDAQ::OfflineModuleName(detectorID);
}

const char *AliHLTDAQInterfaceImplementation::VirtOnlineName(const char *detectorName)
{
  // see header file for class documentation
  return AliDAQ::OnlineName(detectorName);
}

const char *AliHLTDAQInterfaceImplementation::VirtOnlineName(Int_t detectorID)
{
  // see header file for class documentation
  return AliDAQ::OnlineName(detectorID);
}
