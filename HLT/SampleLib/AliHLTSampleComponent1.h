// @(#) $Id$

#ifndef ALIHLTSAMPLECOMPONENT1_H
#define ALIHLTSAMPLECOMPONENT1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTSampleComponent1.h
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  A sample processing component for the HLT. */

#include "AliHLTProcessor.h"

/**
 * @class AliHLTSampleComponent1
 * @brief An HLT sample component
 * @ingroup alihlt_tutorial
 */
class AliHLTSampleComponent1 : public AliHLTProcessor {
public:
  AliHLTSampleComponent1();
  virtual ~AliHLTSampleComponent1();

  const char* GetComponentID() { return "Sample-component1";}
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
    const AliHLTComponentDataType* pType=fgInputDataTypes;
    while (pType->fID!=0) {
      list.push_back(*pType);
      pType++;
    }
  }
  AliHLTComponentDataType GetOutputDataType() {return kAliHLTVoidDataType;}
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {constBase = 0;inputMultiplier = 0;};

  // Spawn function, return new class instance
  AliHLTComponent* Spawn() {return new AliHLTSampleComponent1;};

 protected:
  
  int DoInit( int argc, const char** argv );
  int DoDeinit();
  int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		       AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		       AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );

private:
  /** array of input data types */
  static const AliHLTComponentDataType fgInputDataTypes[]; // see above

  ClassDef(AliHLTSampleComponent1, 0)
};
#endif
