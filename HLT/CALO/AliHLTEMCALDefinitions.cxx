/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Matthias Richter <Matthias.Richter@ift.uib.no>                *
 *          Timm Steinbeck <timm@kip.uni-heidelberg.de>                   *
 *          for The ALICE Off-line Project.                               *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Definitions for the HLT EMCAL components                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliHLTEMCALDefinitions.h"


ClassImp(AliHLTEMCALDefinitions)

const AliHLTComponentDataType AliHLTEMCALDefinitions::fgkDDLRawDataType = 
AliHLTComponentDataTypeInitializer("DDL_RAW ", kAliHLTDataOriginEMCAL);

const AliHLTComponentDataType AliHLTEMCALDefinitions::fgkDigitDataType =
AliHLTComponentDataTypeInitializer("DIGITTYP", kAliHLTDataOriginEMCAL);

const AliHLTComponentDataType AliHLTEMCALDefinitions::fgkClusterDataType =
AliHLTComponentDataTypeInitializer("CLUSTERS", kAliHLTDataOriginEMCAL);

const AliHLTComponentDataType AliHLTEMCALDefinitions::fgkESDDataType = 
AliHLTComponentDataTypeInitializer("GLOBALESD", kAliHLTDataOriginEMCAL);

const AliHLTComponentDataType AliHLTEMCALDefinitions::fgkEMCALESDDataType =
AliHLTComponentDataTypeInitializer("EMCALESD", kAliHLTDataOriginEMCAL);

const AliHLTComponentDataType AliHLTEMCALDefinitions::fgkCalibrationDataType = 
AliHLTComponentDataTypeInitializer("CALIBRAH", kAliHLTDataOriginEMCAL);

const AliHLTComponentDataType AliHLTEMCALDefinitions::fgkChannelDataType = 
AliHLTComponentDataTypeInitializer("CHANNELT", kAliHLTDataOriginEMCAL);

const AliHLTComponentDataType AliHLTEMCALDefinitions::fgkTriggerRawDigitDataType =
AliHLTComponentDataTypeInitializer("TDIGT   ", kAliHLTDataOriginEMCAL);

const AliHLTComponentDataType AliHLTEMCALDefinitions::fgkTriggerSTUDataType =
AliHLTComponentDataTypeInitializer("STUT    ", kAliHLTDataOriginEMCAL);

const AliHLTComponentDataType AliHLTEMCALDefinitions::fgkFastorDataType =
AliHLTComponentDataTypeInitializer("FASTORT ", kAliHLTDataOriginEMCAL);

const AliHLTComponentDataType AliHLTEMCALDefinitions::fgkTriggerPatchDataType =
AliHLTComponentDataTypeInitializer("TRIGGERT", kAliHLTDataOriginEMCAL);

AliHLTEMCALDefinitions::AliHLTEMCALDefinitions()
{
  // see header file for class documentation
  // or
  // refer to README to build package
}

AliHLTEMCALDefinitions::~AliHLTEMCALDefinitions()
{
  // see header file for class documentation
}
