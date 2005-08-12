// @(#) $Id$

#ifndef ALIHLTPROCESSOR_H
#define ALIHLTPROCESSOR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTProcessor
   base class for HLT processing components
 */

#include "AliHLTComponent.h"

class AliHLTProcessor : public AliHLTComponent {
 public:
  AliHLTProcessor();
  virtual ~AliHLTProcessor();

  int Init( AliHLTComponentEnvironment* environ, void* environ_param, int argc, const char** argv );
  int Deinit();
  int ProcessEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
		    AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		    AliHLTUInt32_t& size, AliHLTUInt32_t& outputBlockCnt, 
		    AliHLTComponent_BlockData*& outputBlocks,
		    AliHLTComponent_EventDoneData*& edd );

  // Information member functions for registration.
  TComponentType GetComponentType() { return AliHLTComponent::kProcessor;}

 private:
  virtual int DoInit( int argc, const char** argv ){
    return 0;
  }

  virtual int DoDeinit(){
    return 0;
  }
  virtual int DoEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
		       AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		       AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks ) = 0;


  ClassDef(AliHLTProcessor, 0)
};
#endif
