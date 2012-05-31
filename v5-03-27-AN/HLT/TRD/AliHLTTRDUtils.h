// $Id$
#ifndef ALIHLTTRDUTILS_H
#define ALIHLTTRDUTILS_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

////////////////////////////////////////////////////////////////////////////
//                                                                        //
//  HLT TRD Utillities Class                                              //
//                                                                        //
////////////////////////////////////////////////////////////////////////////


#include "AliHLTDataTypes.h"
#include "TObject.h"
//#include "AliHLTProcessor.h"

class TClonesArray;
class AliESDEvent;
class AliTRDtransform;
class AliHLTTRDUtils
{
public:
  virtual ~AliHLTTRDUtils(){}
  static AliHLTUInt32_t AddClustersToOutput(const TClonesArray *const inClusterArray, AliHLTUInt8_t *const outBlockPtr, Int_t nTimeBins=24);
  static AliHLTUInt32_t AddTracksToOutput(const TClonesArray *const inTrackArray, AliHLTUInt8_t *const output, Int_t nTimeBins=24);
  static AliHLTUInt32_t ReadClusters(TClonesArray *const outArray, const void *const inputPtr, AliHLTUInt32_t size, Int_t* nTimeBins=0x0);
  static AliHLTUInt32_t ReadTracks(TClonesArray *const outArray, const void *const inputPtr, AliHLTUInt32_t size, Int_t* nTimeBins=0x0);
  static AliHLTUInt32_t AddESDToOutput(const AliESDEvent* const esd, AliHLTUInt8_t* const outBlockPtr);
  static void EmulateHLTClusters(TClonesArray *clusterArray);
  static void EmulateHLTTracks(TClonesArray *trackArray);
  static AliHLTUInt32_t GetSM(AliHLTUInt32_t spec);
  static AliHLTUInt32_t AddTracksToOutputAlt(const TClonesArray *const inTrackArray, AliHLTUInt8_t *const output, Int_t nTimeBins=24);
  static AliHLTUInt32_t ReadTracksAlt(TClonesArray *const outArray, const void *const inputPtr, AliHLTUInt32_t size, Int_t* nTimeBins=0x0);

  ClassDef(AliHLTTRDUtils, 0)

};

#endif
