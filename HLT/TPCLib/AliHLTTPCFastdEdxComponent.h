//-*- Mode: C++ -*-
// ************************************************************************
// This file is property of and copyright by the ALICE HLT Project        *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************

///  @file   AliHLTTPCFastdEdxComponent.h
///  @author Sergey Gorbunov <sergey.gorbunov@kip.uni-heidelberg.de>
///  @date   June 2009
///  @brief  An ITS tracker processing component for the HLT

#ifndef ALIHLTTPCFASTDEDXCOMPONENT_H
#define ALIHLTTPCFASTDEDXCOMPONENT_H

#include "AliHLTProcessor.h"
#include "AliHLTDataTypes.h"
class AliHLTTPCRawClusterData;

class AliHLTTPCFastdEdxComponent : public AliHLTProcessor
{
public:
  AliHLTTPCFastdEdxComponent();
  AliHLTTPCFastdEdxComponent( const AliHLTTPCFastdEdxComponent& );
  AliHLTTPCFastdEdxComponent& operator=( const AliHLTTPCFastdEdxComponent& );
  
  virtual ~AliHLTTPCFastdEdxComponent();

  const char* GetComponentID() ;
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list )  ;
  AliHLTComponentDataType GetOutputDataType() ;
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) ;
  AliHLTComponent* Spawn() ;

protected:
  int DoInit( int argc, const char** argv );
  int DoDeinit();
  int Reconfigure( const char* cdbEntry, const char* chainId );
  int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );

private:
  float GetSortTruncMean(float* array, int count, int trunclow, int trunchigh);

  Int_t fNPatchClusters[36][6]; //!
  AliHLTTPCRawClusterData* fPatchClusters[36][6];//!
  
  float fBz;
  float* fBufMax;
  float* fBufTot;
  int fMaxClusterCount;

  void SetDefaultConfiguration();
  int ReadConfigurationString(  const char* arguments );
  int Configure( const char* cdbEntry, const char* chainId, const char *commandLine );

  ClassDef( AliHLTTPCFastdEdxComponent, 0 );
};
#endif
