
// $Id$

/**************************************************************************
 * Copyright(c) 2006, ALICE Experiment at CERN, All rights reserved.      *
 *                                                                        *
 * Authors: Per Thomas Hille <perthi@fys.uio.no>,                         *
 *          �ystein Djuvsland <oysteind@ift.uib.no>after                  * 
 *          Matthias Richter <Matthias.Richter@ift.uib.no>                *
 *          Timm Steinbeck <timm@kip.uni-heidelberg.de>                   *
 *          for the ALICE Offline Project.                                *
 *	                                                                  *
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


const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkDDLPackedRawDataType          = { sizeof(AliHLTComponentDataType),       {'D','D','L','_','R','W','P','K'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkCellEnergyDataType            = { sizeof(AliHLTComponentDataType),       {'C','E','L','L','E','N','E','R'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkCellEnergyHistogramDataType   = { sizeof(AliHLTComponentDataType),       {'E','N','E','R','H','I','S','T'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkCellAverageEnergyDataType     = { sizeof(AliHLTComponentDataType),       {'E','N','E','R','A','V','E','R'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkCellAccumulatedEnergyDataType = { sizeof(AliHLTComponentDataType),       {'E','N','E','R','A','C','C','U'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkCellTimingHistogramDataType   = { sizeof(AliHLTComponentDataType),       {'T','I','M','E','H','I','S','T'},{'P','H','O','S'}};;    
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkCellTimingAverageDataType     = { sizeof(AliHLTComponentDataType),       {'T','I','M','E','A','V','E','R'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkCellChannelDataDataType       = { sizeof(AliHLTComponentDataType),       {'C','H','A','N','D','A','T','A'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkAliHLTClusterDataType         = { sizeof(AliHLTComponentDataType),       {'C','L','U','S','T','R','T','Y'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkAliHLTHistDataType            = { sizeof(AliHLTComponentDataType),       {'H','I','S','T','T','Y','P','E'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkAliHLTSpectrumDataType        = { sizeof(AliHLTComponentDataType),       {'S','P','E','C','T','Y','P','E'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkAliHLTRootTreeDataType        = { sizeof(AliHLTComponentDataType),       {'R','T','R','E','T','Y','P','E'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkAliHLTBaselineDataType        = { sizeof(AliHLTComponentDataType),       {'B','A','S','L','T','Y','P','E'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkAliHLTDigitDataType           = { sizeof(AliHLTComponentDataType),       {'D','I','G','T','Y','P','E'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkAliHLTNoiseMapDataType        = { sizeof(AliHLTComponentDataType),       {'N','O','M','T','Y','P','E'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkAliHLTMIPDataType             = { sizeof(AliHLTComponentDataType),       {'M','I','P','T','Y','P','E'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkAliHLTSandboxDataType         = { sizeof(AliHLTComponentDataType),       {'S','B','X','T','Y','P','E'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkEmcCalibDataType         = { sizeof(AliHLTComponentDataType),       {'C','A','L','T','Y','P','E'},{'P','H','O','S'}};;
