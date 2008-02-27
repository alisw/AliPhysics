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

const AliHLTComponentDataType AliHLTTPCDefinitions::fgkDDLPackedRawDataType = 
  (AliHLTComponentDataType){sizeof(AliHLTComponentDataType),                          {'D','D','L','_','R','W','P','K'},  kAliHLTDataOriginAny} | kAliHLTDataOriginTPC;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkDDLEncodedEntropyRawDataType = 								  	      
  (AliHLTComponentDataType){sizeof(AliHLTComponentDataType),                          {'D','D','L','E','N','C','E','N'},  kAliHLTDataOriginAny} | kAliHLTDataOriginTPC;
										      								  	      
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkPackedRawDataType = 	      								  	      
  (AliHLTComponentDataType){sizeof(AliHLTComponentDataType),                          {'R','A','W','P','A','K','E','D'},  kAliHLTDataOriginAny} | kAliHLTDataOriginTPC;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkUnpackedRawDataType =	      								  	      
  (AliHLTComponentDataType){sizeof(AliHLTComponentDataType),                          {'R','A','W','U','N','P','A','K'},  kAliHLTDataOriginAny} | kAliHLTDataOriginTPC;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkClustersDataType = 	      								  	      
  (AliHLTComponentDataType){sizeof(AliHLTComponentDataType),                          {'C','L','U','S','T','E','R','S'},  kAliHLTDataOriginAny} | kAliHLTDataOriginTPC;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkVertexDataType =		      								  	      
  (AliHLTComponentDataType){sizeof(AliHLTComponentDataType),                          {'V','E','R','T','E','X',' ',' '},  kAliHLTDataOriginAny} | kAliHLTDataOriginTPC;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkTrackSegmentsDataType =	      								  	      
  (AliHLTComponentDataType){sizeof(AliHLTComponentDataType),                          {'T','R','A','K','S','E','G','S'},  kAliHLTDataOriginAny} | kAliHLTDataOriginTPC;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkTracksDataType =		      								  	      
  (AliHLTComponentDataType){sizeof(AliHLTComponentDataType),                          {'T','R','A','C','K','S',' ',' '},  kAliHLTDataOriginAny} | kAliHLTDataOriginTPC;

const AliHLTComponentDataType AliHLTTPCDefinitions::fgkClusterTracksModelDataType =
  (AliHLTComponentDataType){ sizeof(AliHLTComponentDataType),                         {'C','L','S','T','R','K','M','D'},  kAliHLTDataOriginAny} | kAliHLTDataOriginTPC;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkRemainingClustersModelDataType =
  (AliHLTComponentDataType) { sizeof(AliHLTComponentDataType),                        {'R','E','M','C','L','S','M','D'},  kAliHLTDataOriginAny} | kAliHLTDataOriginTPC;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkClusterTracksCompressedDataType =
  (AliHLTComponentDataType) { sizeof(AliHLTComponentDataType),                        {'C','L','S','T','R','K','C','M'},  kAliHLTDataOriginAny} | kAliHLTDataOriginTPC;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkRemainingClustersCompressedDataType =
  (AliHLTComponentDataType) { sizeof(AliHLTComponentDataType),                        {'R','E','M','C','L','S','C','M'},  kAliHLTDataOriginAny} | kAliHLTDataOriginTPC;
										      								  	      
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkCalibPedestalDataType =	      								  	      
  (AliHLTComponentDataType){sizeof(AliHLTComponentDataType),                          {'C','A','L','_','P','E','D',' '},  kAliHLTDataOriginAny} | kAliHLTDataOriginTPC;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkCalibPulserDataType =	      								  	      
  (AliHLTComponentDataType){sizeof(AliHLTComponentDataType),                          {'C','A','L','_','P','U','L','S'},  kAliHLTDataOriginAny} | kAliHLTDataOriginTPC;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkCalibCEDataType =	      								  	      
  (AliHLTComponentDataType){sizeof(AliHLTComponentDataType),                          {'C','A','L','_','C','E',' ',' '},  kAliHLTDataOriginAny} | kAliHLTDataOriginTPC;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkActivePadsDataType =	      								  	      
  (AliHLTComponentDataType){sizeof(AliHLTComponentDataType),                          {'A','C','T','I','V','P','A','D'},  kAliHLTDataOriginAny} | kAliHLTDataOriginTPC;
const AliHLTComponentDataType AliHLTTPCDefinitions::fgkNoiseHistoDataType =
  (AliHLTComponentDataType){sizeof(AliHLTComponentDataType),                          {'N','O','I','S','E','M','A','P'},  kAliHLTDataOriginAny} | kAliHLTDataOriginTPC;

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

    
