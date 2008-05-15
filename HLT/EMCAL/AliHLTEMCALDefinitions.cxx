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

const AliHLTComponentDataType AliHLTEMCALDefinitions::fgkDDLRawDataType = { sizeof(AliHLTComponentDataType), {'D','D','L','_','R','A','W',' '},{'E','M','C','L'}};;

const AliHLTComponentDataType AliHLTEMCALDefinitions::fgkClusterDataType = { sizeof(AliHLTComponentDataType), {'C','L','U','S','T','E','R','S'},{'E','M','C','L'}};;

const AliHLTComponentDataType AliHLTEMCALDefinitions::fgkESDDataType = { sizeof(AliHLTComponentDataType), {'G','L','O','B','L','E','S','D'},{'E','M','C','L'}};;

const AliHLTComponentDataType AliHLTEMCALDefinitions::fgkEMCALESDDataType = { sizeof(AliHLTComponentDataType), {'E','M','C','A','L','E','S','D'},{'E','M','C','L'}};;

const AliHLTComponentDataType AliHLTEMCALDefinitions::fgkCalibrationDataType = { sizeof(AliHLTComponentDataType), {'C','A','L','I','B','R','A','H'},{'E','M','C','L'}};;

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
