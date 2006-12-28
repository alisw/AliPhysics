// @(#) $Id$

#ifndef ALIHLTSAMPLECOMPONENT2_H
#define ALIHLTSAMPLECOMPONENT2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/** @file   AliHLTSampleComponent2.h
    @author Matthias Richter, Timm Steinbeck
    @date   
    @brief  Another sample processing component for the HLT. */

#include "AliHLTProcessor.h"

/**
 * @class AliHLTSampleComponent2
 * @brief An HLT sample component
 * @ingroup alihlt_tutorial
 */
class AliHLTSampleComponent2 : public AliHLTProcessor {
public:
  AliHLTSampleComponent2();
  virtual ~AliHLTSampleComponent2();

  const char* GetComponentID() { return "Sample-component2";}
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
    const AliHLTComponentDataType* pType=inputDataTypes;
    while (pType->fID!=0) {
      list.push_back(*pType);
      pType++;
    }
  }
  AliHLTComponentDataType GetOutputDataType() {return outputDataType;}
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {constBase = 0;inputMultiplier = 0;};

  // Spawn function, return new class instance
  AliHLTComponent* Spawn() {return new AliHLTSampleComponent2;};

 protected:
  
  int DoInit( int argc, const char** argv );
  int DoDeinit();
  int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		       AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		       AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );

private:
  static const AliHLTComponentDataType inputDataTypes[];
  static const AliHLTComponentDataType outputDataType;

  ClassDef(AliHLTSampleComponent2, 0)
};
#endif
