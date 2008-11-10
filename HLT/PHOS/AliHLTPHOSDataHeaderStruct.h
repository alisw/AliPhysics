//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTPHOSDATAHEADERSTRUCT_H
#define ALIHLTPHOSDATAHEADERSTRUCT_H
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Author:  Per Thomas Hille  <perthi@fys.uio.no>                 *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


#include "AliHLTDataTypes.h"

struct AliHLTPHOSDataHeaderStruct
{
  AliHLTUInt32_t fSize;              /**<Total size of datablock in bytes, incuding the header*/
  AliHLTComponentDataType fDataType; /**<Data type stored in this file */
  AliHLTUInt32_t fEventID;           /**<The HLT internal event ID for this event */
  AliHLTUInt32_t fAlgorithm;         /**<Wich algorithm was uses estimate cellenergies*/
  AliHLTUInt32_t fFormatVersion;     /**<Header format version, currently 1*/
  AliHLTUInt32_t fFutureUse0;
  AliHLTUInt32_t fFutureUse1;
  AliHLTUInt32_t fFutureUse2;
};


#endif
