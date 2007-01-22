// $Id$

/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Author: Per Thomas Hille for the ALICE HLT Project.                    *
 * Contributors are mentioned in the code where appropriate.              *
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
// Definitions for the HLT PHOS components                                    //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include "AliHLTPHOSDefinitions.h"


//ClassImp(AliHLTPHOSDefinitions)

const AliHLTComponentDataType AliHLTPHOSDefinitions::gkDDLPackedRawDataType = { sizeof(AliHLTComponentDataType), {'D','D','L','_','R','W','P','K'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::gkPackedRawDataType = { sizeof(AliHLTComponentDataType), {'R','A','W','P','A','K','E','D'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::gkUnpackedRawDataType = { sizeof(AliHLTComponentDataType), {'R','A','W','U','N','P','A','K'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::gkClustersDataType = { sizeof(AliHLTComponentDataType), {'C','L','U','S','T','E','R','S'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::gkVertexDataType = { sizeof(AliHLTComponentDataType), {'V','E','R','T','E','X',' ',' '},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::gkTrackSegmentsDataType = { sizeof(AliHLTComponentDataType), {'T','R','A','K','S','E','G','S'},{'P','H','O','S'}};;

    
