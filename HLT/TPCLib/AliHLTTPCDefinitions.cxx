// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
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

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
// Definitions for the HLT TPC components                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliHLTTPCDefinitions.h"


/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTTPCDefinitions)

const AliHLTComponentDataType AliHLTTPCDefinitions::fgkDDLPackedRawDataType = { sizeof(AliHLTComponentDataType),         {'D','D','L','_','R','W','P','K'},kAliHLTDataOriginTPC};;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkDDLEncodedEntropyRawDataType = { sizeof(AliHLTComponentDataType), {'D','D','L','E','N','C','E','N'},kAliHLTDataOriginTPC};;

const AliHLTComponentDataType AliHLTTPCDefinitions::fgkPackedRawDataType = { sizeof(AliHLTComponentDataType), {'R','A','W','P','A','K','E','D'},kAliHLTDataOriginTPC};;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkUnpackedRawDataType = { sizeof(AliHLTComponentDataType), {'R','A','W','U','N','P','A','K'},kAliHLTDataOriginTPC};;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkClustersDataType = { sizeof(AliHLTComponentDataType), {'C','L','U','S','T','E','R','S'},kAliHLTDataOriginTPC};;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkVertexDataType = { sizeof(AliHLTComponentDataType), {'V','E','R','T','E','X',' ',' '},kAliHLTDataOriginTPC};;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkTrackSegmentsDataType = { sizeof(AliHLTComponentDataType), {'T','R','A','K','S','E','G','S'},kAliHLTDataOriginTPC};;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkTracksDataType = { sizeof(AliHLTComponentDataType), {'T','R','A','C','K','S',' ',' '},kAliHLTDataOriginTPC};;

const AliHLTComponentDataType AliHLTTPCDefinitions::fgkCalibPedestalDataType = { sizeof(AliHLTComponentDataType), {'C','A','L','_','P','E','D',' '},kAliHLTDataOriginTPC};;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkCalibSignalDataType = { sizeof(AliHLTComponentDataType), {'C','A','L','_','S','I','G',' '},kAliHLTDataOriginTPC};;


AliHLTTPCDefinitions::AliHLTTPCDefinitions()
{
  // see header file for class documentation
  // or
  // refer to README to build package
  // or
  // visit http://web.ift.uib.no/~kjeks/doc/alice-hlt
}

AliHLTTPCDefinitions::~AliHLTTPCDefinitions()
{
  // see header file for class documentation
}

    
