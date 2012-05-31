

// $Id: AliHLTPHOSDefinitions.cxx 34264 2009-08-14 18:29:23Z odjuvsla $

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

#include "AliHLTCaloDefinitions.h"
 const AliHLTComponentDataType AliHLTCaloDefinitions::fgkPhosHistDataType              = { sizeof(AliHLTComponentDataType),       {'H','I','S','T','O','G','R','A'},{'P','H','O','S'}};;
//CRAP PTH
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkFourierTransform              = { sizeof(AliHLTComponentDataType),       {'F','O','U','R','I','E','R','T'},{'P','H','O','S'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkChannelDataType               = { sizeof(AliHLTComponentDataType),       {'C','H','A','N','N','E','L','T'},{'C','A','L','O'}};;
//const AliHLTComponentDataType AliHLTCaloDefinitions::fgkDDLPackedRawDataType          = { sizeof(AliHLTComponentDataType),       {'D','D','L','_','R','W','P','K'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkDDLPackedRawDataType          = { sizeof(AliHLTComponentDataType),       {'D','D','L','_','R','A','W',' '},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkCellEnergyDataType            = { sizeof(AliHLTComponentDataType),       {'C','E','L','L','E','N','E','R'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkCellEnergyHistogramDataType   = { sizeof(AliHLTComponentDataType),       {'E','N','E','R','H','I','S','T'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkCellAverageEnergyDataType     = { sizeof(AliHLTComponentDataType),       {'E','N','E','R','A','V','E','R'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkCellAccumulatedEnergyDataType = { sizeof(AliHLTComponentDataType),       {'E','N','E','R','A','C','C','U'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkCellTimingHistogramDataType   = { sizeof(AliHLTComponentDataType),       {'T','I','M','E','H','I','S','T'},{'C','A','L','O'}};;    
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkCellTimingAverageDataType     = { sizeof(AliHLTComponentDataType),       {'T','I','M','E','A','V','E','R'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkCellChannelDataDataType       = { sizeof(AliHLTComponentDataType),       {'C','H','A','N','D','A','T','A'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkClusterDataType               = { sizeof(AliHLTComponentDataType),       {'C','L','U','S','T','R','T','Y'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkRecPointDataType              = { sizeof(AliHLTComponentDataType),       {'R','E','C','P','O','I','N','T'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkHistDataType                  = { sizeof(AliHLTComponentDataType),       {'H','I','S','T','O','G','R','A'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkSpectrumDataType              = { sizeof(AliHLTComponentDataType),       {'S','P','E','C','T','R','U','M'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkRootTreeDataType              = { sizeof(AliHLTComponentDataType),       {'R','O','O','T','T','R','E','E'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkBaselineDataType              = { sizeof(AliHLTComponentDataType),       {'B','A','S','E','L','I','N','E'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkNoiseMapDataType              = { sizeof(AliHLTComponentDataType),       {'N','O','I','S','E','M','A','P'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkMIPDataType                   = { sizeof(AliHLTComponentDataType),       {'M','I','P','D','T','Y','P','E'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkSandboxDataType               = { sizeof(AliHLTComponentDataType),       {'S','A','N','D','B','O','X','T'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkEmcCalibDataType              = { sizeof(AliHLTComponentDataType),       {'C','A','L','I','T','Y','P','E'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkCaloClusterDataType           = { sizeof(AliHLTComponentDataType),       {'C','A','L','C','L','U','S','T'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkESDCaloClusterDataType        = { sizeof(AliHLTComponentDataType),       {'E','S','D','C','L','U','S','T'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkESDCaloCellsDataType          = { sizeof(AliHLTComponentDataType),       {'E','S','D','C','E','L','L','S'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkPhysicsHistogramsDataType     = { sizeof(AliHLTComponentDataType),       {'P','H','Y','S','H','I','S','T'},{'C','A','L','O'}};;
const AliHLTComponentDataType AliHLTCaloDefinitions::fgkDigitDataType     = { sizeof(AliHLTComponentDataType),       {'D','I','G','I','T','T','Y','P'},{'*','*','*','*'}};;


