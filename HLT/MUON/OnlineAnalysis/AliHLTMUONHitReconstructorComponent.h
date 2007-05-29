#ifndef ALIHLTMUONHITRECONSTRUCTORCOMPONENT_H
#define ALIHLTMUONHITRECONSTRUCTORCOMPONENT_H
/* Copyright(c) 1998-2007, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/*  @file   AliHLTMUONHitReconstructorComponent.h
 *  @author Indranil Das <indra.das@saha.ac.in> | <indra.ehep@gmail.com>
 *  @date   
 *  @brief  Hit Reconstruction processing component for the dimuon HLT. 
 */

#include "AliHLTProcessor.h"
#include "AliHLTMUONConstants.h"
#include "AliHLTMUONHitReconstructor.h"

//class AliHLTMUONHitReconstructor;

class AliHLTMUONHitReconstructorComponent : public AliHLTProcessor {

public:
  AliHLTMUONHitReconstructorComponent();
  virtual ~AliHLTMUONHitReconstructorComponent();

  const char* GetComponentID() { return "MUONHitRec";}

  void GetInputDataTypes( vector<AliHLTComponentDataType>& list) {
    list.clear();
    list.push_back( AliHLTMUONConstants::TrackingDDLRawDataType() );
  }
  
  AliHLTComponentDataType GetOutputDataType() {return AliHLTMUONConstants::RecHitsBlockDataType();}
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) {constBase = 0;inputMultiplier = 0;};

  // Spawn function, return new class instance
  AliHLTComponent* Spawn() {return new AliHLTMUONHitReconstructorComponent;};

 protected:
  
  int DoInit( int argc, const char** argv );
  int DoDeinit();
  int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		       AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		       AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );

private:
  /** array of input data types */

  AliHLTMUONHitReconstructor* fHitRec;
  bool ReadLookUpTable(DHLTLut* lookupTable, const char* lutpath, int iDDL);
  bool ReadBusPatchToDetElemFile(BusToDetElem& busToDetElem, const char* buspatchmappath);


  ClassDef(AliHLTMUONHitReconstructorComponent, 0)
};

#endif // ALIHLTMUONHITRECONSTRUCTORCOMPONENT_H
