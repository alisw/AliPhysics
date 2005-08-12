// @(#) $Id$

#ifndef ALIHLTSAMPLECOMPONENT2_H
#define ALIHLTSAMPLECOMPONENT2_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* AliHLTSampleComponent2
 */

#include "AliHLTProcessor.h"

class AliHLTSampleComponent2 : public AliHLTProcessor {
public:
  AliHLTSampleComponent2();
  virtual ~AliHLTSampleComponent2();

  const char* GetComponentID() { return "Sample-component2";}
  void GetInputDataTypes( vector<AliHLTComponent_DataType>& list) {
    const AliHLTComponent_DataType* pType=inputDataTypes;
    while (pType->fID!=0) {
      list.push_back(*pType);
      pType++;
    }
  }
  AliHLTComponent_DataType GetOutputDataType() {return outputDataType;}
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {constBase = 0;inputMultiplier = 0;};

  // Spawn function, return new class instance
  AliHLTComponent* Spawn() {return new AliHLTSampleComponent2;};

 protected:
  
  int DoInit( int argc, const char** argv );
  int DoDeinit();
  int DoEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
		       AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		       AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks );

private:
  static const AliHLTComponent_DataType inputDataTypes[];
  static const AliHLTComponent_DataType outputDataType;

  ClassDef(AliHLTSampleComponent2, 0)
};
#endif
