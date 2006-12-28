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
// Definitions for the HLT TPC components                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliHLTTPCDefinitions.h"


ClassImp(AliHLTTPCDefinitions)

const AliHLTComponentDataType AliHLTTPCDefinitions::gkDDLPackedRawDataType = { sizeof(AliHLTComponentDataType), {'D','D','L','_','R','W','P','K'},{'T','P','C',' '}};;
const AliHLTComponentDataType AliHLTTPCDefinitions::gkPackedRawDataType = { sizeof(AliHLTComponentDataType), {'R','A','W','P','A','K','E','D'},{'T','P','C',' '}};;
const AliHLTComponentDataType AliHLTTPCDefinitions::gkUnpackedRawDataType = { sizeof(AliHLTComponentDataType), {'R','A','W','U','N','P','A','K'},{'T','P','C',' '}};;
const AliHLTComponentDataType AliHLTTPCDefinitions::gkClustersDataType = { sizeof(AliHLTComponentDataType), {'C','L','U','S','T','E','R','S'},{'T','P','C',' '}};;
const AliHLTComponentDataType AliHLTTPCDefinitions::gkVertexDataType = { sizeof(AliHLTComponentDataType), {'V','E','R','T','E','X',' ',' '},{'T','P','C',' '}};;
const AliHLTComponentDataType AliHLTTPCDefinitions::gkTrackSegmentsDataType = { sizeof(AliHLTComponentDataType), {'T','R','A','K','S','E','G','S'},{'T','P','C',' '}};;

    
