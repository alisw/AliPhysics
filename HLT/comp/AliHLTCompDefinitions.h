// XEmacs -*-C++-*-
// @(#) $Id$

#ifndef ALIHLTCOMPDEFINITIONS_H
#define ALIHLTCOMPDEFINITIONS_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/* AliHLTCompDefinitions
 */

#include "AliHLTDataTypes.h"
#include "Rtypes.h"

class AliHLTCompDefinitions
{
public:
  AliHLTCompDefinitions() {}
  virtual ~AliHLTCompDefinitions() {}

//   /** DDL packed RAW data */
//   static const AliHLTComponentDataType fgkDDLRawDataType;    // see above
  /** DDL entropy encoded data */
  static const AliHLTComponentDataType fgkDDLEncodedHuffmanAltroDataType; // see above
//   /** packed RAW data */
//   static const AliHLTComponentDataType fgkPackedRawDataType;       // see above
//   /** unpacked RAW data */
//   static const AliHLTComponentDataType fgkUnpackedRawDataType;     // see above
//   /** cluster data */
//   static const AliHLTComponentDataType fgkClustersDataType;        // see above
//   /** track segments in local coordinates */
//   static const AliHLTComponentDataType fgkTrackSegmentsDataType;   // see above
//   /** tracks in global koordinates */
//   static const AliHLTComponentDataType fgkTracksDataType;          // see above
  /** altro huffman compression table */
  static const AliHLTComponentDataType fgkHuffmanAltroCalDataType; // see above

  ClassDef(AliHLTCompDefinitions, 0);

};

#endif
