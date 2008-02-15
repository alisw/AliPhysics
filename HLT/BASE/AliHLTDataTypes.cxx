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
const AliHLTComponentDataType kAliHLTMultipleDataType =  (AliHLTComponentDataType) {
  sizeof(AliHLTComponentDataType),
  {'M','U','L','T','I','P','L','E'},
  kAliHLTDataOriginAny
}|kAliHLTDataOriginPrivate;

/** data to file exchange subscriber */
const AliHLTComponentDataType kAliHLTDataTypeFXSCalib =  (AliHLTComponentDataType) {
  sizeof(AliHLTComponentDataType),
  kAliHLTFXSCalibDataTypeID,
  kAliHLTDataOriginAny
}|kAliHLTDataOriginOut;

/** DDL list data type */
const AliHLTComponentDataType kAliHLTDataTypeDDL  =  (AliHLTComponentDataType) {
  sizeof(AliHLTComponentDataType),
  kAliHLTDDLDataTypeID,
  kAliHLTDataOriginAny
}|kAliHLTDataOriginOut;

/** SOR data type */
const AliHLTComponentDataType kAliHLTDataTypeSOR  =  (AliHLTComponentDataType) {
  sizeof(AliHLTComponentDataType),
  kAliHLTSORDataTypeID,
  kAliHLTDataOriginAny
}|kAliHLTDataOriginPrivate;

/** EOR data type */
const AliHLTComponentDataType kAliHLTDataTypeEOR  =  (AliHLTComponentDataType) {
  sizeof(AliHLTComponentDataType),
  kAliHLTEORDataTypeID,
  kAliHLTDataOriginAny
}|kAliHLTDataOriginPrivate;

/** Event type specification */
const AliHLTComponentDataType kAliHLTDataTypeEvent  =  (AliHLTComponentDataType) {
  sizeof(AliHLTComponentDataType),
  kAliHLTEventDataTypeID,
  kAliHLTDataOriginAny
}|kAliHLTDataOriginPrivate;

/** Configuration event data type */
const AliHLTComponentDataType kAliHLTDataTypeComConf  =  (AliHLTComponentDataType) {
  sizeof(AliHLTComponentDataType),
  kAliHLTComConfDataTypeID,
  kAliHLTDataOriginAny
}|kAliHLTDataOriginPrivate;

/** DCS value update event */
const AliHLTComponentDataType kAliHLTDataTypeUpdtDCS  =  (AliHLTComponentDataType) {
  sizeof(AliHLTComponentDataType),
  kAliHLTUpdtDCSDataTypeID,
  kAliHLTDataOriginAny
}|kAliHLTDataOriginPrivate;

/** RAW DDL data specification, data publisher will set type id and origin correctly */
const AliHLTComponentDataType kAliHLTDataTypeDDLRaw =  (AliHLTComponentDataType) {
  sizeof(AliHLTComponentDataType),
  kAliHLTDDLRawDataTypeID,
  kAliHLTDataOriginAny
};

/** ESD data specification */
const AliHLTComponentDataType kAliHLTDataTypeESDObject =  (AliHLTComponentDataType) {
  sizeof(AliHLTComponentDataType),
  kAliHLTESDObjectDataTypeID,
  kAliHLTDataOriginAny
};
/** ESD tree data specification */
const AliHLTComponentDataType kAliHLTDataTypeESDTree =  (AliHLTComponentDataType) {
  sizeof(AliHLTComponentDataType),
  kAliHLTESDTreeDataTypeID,
  kAliHLTDataOriginAny
};

/** 16 bit Hardware address selection data specification, origin is 'any' */
const AliHLTComponentDataType kAliHLTDataTypeHwAddr16 = (AliHLTComponentDataType) {
  sizeof(AliHLTComponentDataType),
  kAliHLTHwAddr16DataTypeID,
  kAliHLTDataOriginAny
};



//////////////////////////////////////////////////////////////////////////
//
// Data urigin variables, to be used with the operator|
//
// AliHLTComponentDataType dt;
// dt = kAliHLTDataTypeDDLRaw | gkAliHLTDataOriginTPC;
//
//////////////////////////////////////////////////////////////////////////

/** HLT out */
const char kAliHLTDataOriginOut[kAliHLTComponentDataTypefOriginSize]    = {'H','L','T',' '};

/** HLT/PubSub private internal */
const char kAliHLTDataOriginPrivate[kAliHLTComponentDataTypefOriginSize]= {'P','R','I','V'};

/** TPC */
const char kAliHLTDataOriginTPC[kAliHLTComponentDataTypefOriginSize]    = {'T','P','C',' '};

/** PHOS */
const char kAliHLTDataOriginPHOS[kAliHLTComponentDataTypefOriginSize]   = {'P','H','O','S'};

/** MUON */
const char kAliHLTDataOriginMUON[kAliHLTComponentDataTypefOriginSize]   = {'M','U','O','N'};

/** TRD */
const char kAliHLTDataOriginTRD[kAliHLTComponentDataTypefOriginSize]    = {'T','R','D',' '};

/** ITS */
const char kAliHLTDataOriginITS[kAliHLTComponentDataTypefOriginSize]    = {'I','T','S',' '};
