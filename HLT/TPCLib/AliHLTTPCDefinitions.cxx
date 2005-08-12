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

const AliHLTComponent_DataType AliHLTTPCDefinitions::gkDDLPackedRawDataType = { sizeof(AliHLTComponent_DataType), {'D','D','L','_','R','W','P','K'},{'T','P','C',' '}};;
const AliHLTComponent_DataType AliHLTTPCDefinitions::gkPackedRawDataType = { sizeof(AliHLTComponent_DataType), {'R','A','W','P','A','K','E','D'},{'T','P','C',' '}};;
const AliHLTComponent_DataType AliHLTTPCDefinitions::gkUnpackedRawDataType = { sizeof(AliHLTComponent_DataType), {'R','A','W','U','N','P','A','K'},{'T','P','C',' '}};;
const AliHLTComponent_DataType AliHLTTPCDefinitions::gkClustersDataType = { sizeof(AliHLTComponent_DataType), {'C','L','U','S','T','E','R','S'},{'T','P','C',' '}};;
const AliHLTComponent_DataType AliHLTTPCDefinitions::gkVertexDataType = { sizeof(AliHLTComponent_DataType), {'V','E','R','T','E','X',' ',' '},{'T','P','C',' '}};;
const AliHLTComponent_DataType AliHLTTPCDefinitions::gkTrackSegmentsDataType = { sizeof(AliHLTComponent_DataType), {'T','R','A','K','S','E','G','S'},{'T','P','C',' '}};;

    
