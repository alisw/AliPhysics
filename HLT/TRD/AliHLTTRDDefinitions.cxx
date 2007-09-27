// $Id$

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
// Definitions for the HLT TRD components                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliHLTTRDDefinitions.h"


ClassImp(AliHLTTRDDefinitions)

const AliHLTComponentDataType AliHLTTRDDefinitions::fgkDDLRawDataType = { sizeof(AliHLTComponentDataType), {'D','D','L','_','R','A','W',' '},{'T','R','D',' '}};;

const AliHLTComponentDataType AliHLTTRDDefinitions::fgkClusterDataType = { sizeof(AliHLTComponentDataType), {'C','L','U','S','T','E','R','S'},{'T','R','D',' '}};;

const AliHLTComponentDataType AliHLTTRDDefinitions::fgkTRDSATracksDataType = { sizeof(AliHLTComponentDataType), {'T','R','A','C','K','S','S','A'},{'T','R','D',' '}};;

const AliHLTComponentDataType AliHLTTRDDefinitions::fgkTRDSAEsdDataType = { sizeof(AliHLTComponentDataType), {'S','A','E','S','D','T','R','D'},{'T','R','D',' '}};;

const AliHLTComponentDataType AliHLTTRDDefinitions::fgkMCMtrackletDataType = { sizeof(AliHLTComponentDataType), {'M','C','M','T','R','L','E','T'},{'T','R','D',' '}};;

const AliHLTComponentDataType AliHLTTRDDefinitions::fgkMCMcalibrationDataType = { sizeof(AliHLTComponentDataType), {'M','C','M','C','A','L','I','H'},{'T','R','D',' '}};;

AliHLTTRDDefinitions::AliHLTTRDDefinitions()
{
  // see header file for class documentation
  // or
  // refer to README to build package
}

AliHLTTRDDefinitions::~AliHLTTRDDefinitions()
{
  // see header file for class documentation
}
