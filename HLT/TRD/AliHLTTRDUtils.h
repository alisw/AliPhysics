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
class AliHLTTRDUtils
{
public:
  virtual ~AliHLTTRDUtils(){}
  static AliHLTUInt32_t AddClustersToOutput(TClonesArray* inClusterArray, AliHLTUInt8_t* outBlockPtr, Int_t nTimeBins=24);
  static AliHLTUInt32_t AddTracksToOutput(TClonesArray* inTrackArray, AliHLTUInt8_t* output, Int_t nTimeBins=24);
  static AliHLTUInt32_t ReadClusters(TClonesArray *outArray, void* inputPtr, AliHLTUInt32_t size, Int_t* nTimeBins=0x0);
  static AliHLTUInt32_t ReadTracks(TClonesArray *outArray, void* inputPtr, AliHLTUInt32_t size, Int_t* nTimeBins=0x0);
  static AliHLTUInt32_t AddESDToOutput(const AliESDEvent* const esd, AliHLTUInt8_t* const outBlockPtr);
  static void EmulateHLTClusters(TClonesArray *clusterArray);
  static void EmulateHLTTracks(TClonesArray *trackArray);

  ClassDef(AliHLTTRDUtils, 0)


};

#endif
