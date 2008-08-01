
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
//CRAP PTH
//const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkDDLPackedRawDataType          = { sizeof(AliHLTComponentDataType),       {'D','D','L','_','R','W','P','K'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkDDLPackedRawDataType          = { sizeof(AliHLTComponentDataType),       {'D','D','L','_','R','A','W',' '},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkCellEnergyDataType            = { sizeof(AliHLTComponentDataType),       {'C','E','L','L','E','N','E','R'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkCellEnergyHistogramDataType   = { sizeof(AliHLTComponentDataType),       {'E','N','E','R','H','I','S','T'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkCellAverageEnergyDataType     = { sizeof(AliHLTComponentDataType),       {'E','N','E','R','A','V','E','R'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkCellAccumulatedEnergyDataType = { sizeof(AliHLTComponentDataType),       {'E','N','E','R','A','C','C','U'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkCellTimingHistogramDataType   = { sizeof(AliHLTComponentDataType),       {'T','I','M','E','H','I','S','T'},{'P','H','O','S'}};;    
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkCellTimingAverageDataType     = { sizeof(AliHLTComponentDataType),       {'T','I','M','E','A','V','E','R'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkCellChannelDataDataType       = { sizeof(AliHLTComponentDataType),       {'C','H','A','N','D','A','T','A'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkClusterDataType               = { sizeof(AliHLTComponentDataType),       {'C','L','U','S','T','R','T','Y'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkRecPointDataType              = { sizeof(AliHLTComponentDataType),       {'R','E','C','P','O','I','N','T'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkHistDataType                  = { sizeof(AliHLTComponentDataType),       {'H','I','S','T','O','G','R','A'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkSpectrumDataType              = { sizeof(AliHLTComponentDataType),       {'S','P','E','C','T','R','U','M'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkRootTreeDataType              = { sizeof(AliHLTComponentDataType),       {'R','O','O','T','T','R','E','E'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkBaselineDataType              = { sizeof(AliHLTComponentDataType),       {'B','A','S','E','L','I','N','E'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkDigitDataType                 = { sizeof(AliHLTComponentDataType),       {'D','I','G','I','T','T','Y','P'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkNoiseMapDataType              = { sizeof(AliHLTComponentDataType),       {'N','O','I','S','E','M','A','P'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkMIPDataType                   = { sizeof(AliHLTComponentDataType),       {'M','I','P','D','T','Y','P','E'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkSandboxDataType               = { sizeof(AliHLTComponentDataType),       {'S','A','N','D','B','O','X','T'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkEmcCalibDataType              = { sizeof(AliHLTComponentDataType),       {'C','A','L','I','T','Y','P','E'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTPHOSDefinitions::fgkCaloClusterDataType           = { sizeof(AliHLTComponentDataType),       {'C','A','L','C','L','U','S','T'},{'P','H','O','S'}};;


