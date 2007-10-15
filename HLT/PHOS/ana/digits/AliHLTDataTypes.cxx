// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Matthias Richter <Matthias.Richter@ift.uib.no>        *
 *                  Timm Steinbeck <timm@kip.uni-heidelberg.de>           *
 *                  Jochen Thaeder <thaeder@kip.uni-heidelberg.de>        *
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

/** @file   AliHLTDataTypes.cxx
    @author Matthias Richter, Timm Steinbeck, Jochen Thaeder
    @date   
    @brief  Implementation of data types. */

#include "AliHLTDataTypes.h"

// those types can not be implemented in the header files as rootcint
// can not cope with the type id and origin defines.

/** multiple output data types */
const AliHLTComponentDataType kAliHLTMultipleDataType = {
  sizeof(AliHLTComponentDataType),
  {'M','U','L','T','I','P','L','E'},
  kAliHLTDataOriginPrivate
};

/** data to file exchange subscriber */
const AliHLTComponentDataType kAliHLTDataTypeFXSCalib = {
  sizeof(AliHLTComponentDataType),
  kAliHLTFXSCalibDataTypeID,
  kAliHLTDataOriginOut
};

/** DDL list data type */
const AliHLTComponentDataType kAliHLTDataTypeDDL  = {
  sizeof(AliHLTComponentDataType),
  kAliHLTDDLDataTypeID,
  kAliHLTDataOriginOut
};

/** SOR data type */
const AliHLTComponentDataType kAliHLTDataTypeSOR  = {
  sizeof(AliHLTComponentDataType),
  kAliHLTSORDataTypeID,
  kAliHLTDataOriginPrivate
};

/** EOR data type */
const AliHLTComponentDataType kAliHLTDataTypeEOR  = {
  sizeof(AliHLTComponentDataType),
  kAliHLTEORDataTypeID,
  kAliHLTDataOriginPrivate
};

/** Event type specification */
const AliHLTComponentDataType kAliHLTDataTypeEvent  = {
  sizeof(AliHLTComponentDataType),
  kAliHLTEventDataTypeID,
  kAliHLTDataOriginPrivate
};
