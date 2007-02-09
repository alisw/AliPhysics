#ifndef ALIHLTPHOSDATAHEADERSTRUCT_H
#define ALIHLTPHOSDATAHEADERSTRUCT_H

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
