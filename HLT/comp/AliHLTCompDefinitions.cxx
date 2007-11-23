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
// Definitions for the HLT COMP components                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliHLTCompDefinitions.h"


/** ROOT macro for the implementation of ROOT specific class methods */
ClassImp(AliHLTCompDefinitions)

/*
general concept for data types is needed: how do we define data types
across different libraries. Currently, all definitions have data origin
TPC but this has to change.
 */
//const AliHLTComponentDataType AliHLTCompDefinitions::fgkDDLRawDataType = { sizeof(AliHLTComponentDataType),         {'D','D','L','_','R','A','W',' '},kAliHLTDataOriginTPC};;
const AliHLTComponentDataType AliHLTCompDefinitions::fgkDDLEncodedHuffmanAltroDataType = { sizeof(AliHLTComponentDataType), {'E','N','C','_','H','U','F','F'},kAliHLTDataOriginAny};;

//const AliHLTComponentDataType AliHLTCompDefinitions::fgkPackedRawDataType = { sizeof(AliHLTComponentDataType), {'R','A','W','P','A','K','E','D'},kAliHLTDataOriginTPC};;
//const AliHLTComponentDataType AliHLTCompDefinitions::fgkUnpackedRawDataType = { sizeof(AliHLTComponentDataType), {'R','A','W','U','N','P','A','K'},kAliHLTDataOriginTPC};;
//const AliHLTComponentDataType AliHLTCompDefinitions::fgkClustersDataType = { sizeof(AliHLTComponentDataType), {'C','L','U','S','T','E','R','S'},kAliHLTDataOriginTPC};;
//const AliHLTComponentDataType AliHLTCompDefinitions::fgkTrackSegmentsDataType = { sizeof(AliHLTComponentDataType), {'T','R','A','K','S','E','G','S'},kAliHLTDataOriginTPC};;
//const AliHLTComponentDataType AliHLTCompDefinitions::fgkTracksDataType = { sizeof(AliHLTComponentDataType), {'T','R','A','C','K','S',' ',' '},kAliHLTDataOriginTPC};;

const AliHLTComponentDataType AliHLTCompDefinitions::fgkHuffmanAltroCalDataType = { sizeof(AliHLTComponentDataType), {'C','A','L','_','H','U','F','F'},kAliHLTDataOriginAny};;
