//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTRDCLUSTERIZERCOMPONENT_H
#define ALIHLTTRDCLUSTERIZERCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTRDClusterizerComponent.h
    @author Theodor Rascanu
    @date   
    @brief  Declaration of a TRDClusterizer component.
*/


#include "AliHLTProcessor.h"
class AliCDBManager;
class AliHLTTRDClusterizer;
class AliRawReaderMemory;
class TFile;
class AliTRDrecoParam;
class AliTRDReconstructor;

/**
 * @class AliHLTTRDClusterizerComponent
 * @brief A TRDClusterizer HLT processing component. 
 *
 * An implementiation of a TRDClusterizer component that just copies its input data
 * as a test, demonstration, and example of the HLT component scheme.
 * @ingroup alihlt_tutorial
 */
class AliHLTTRDClusterizerComponent : public AliHLTProcessor
{
public:
  AliHLTTRDClusterizerComponent();
  virtual ~AliHLTTRDClusterizerComponent();

  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process

  const char* GetComponentID();
  void GetInputDataTypes( vector<AliHLTComponentDataType>& list);
  AliHLTComponentDataType GetOutputDataType();
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  AliHLTComponent* Spawn();
  void GetOCDBObjectDescription(TMap* const targetMap);
	
protected:
	
  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing
  // capabilities of the component. 

  int DoInit( int argc, const char** argv );
  int DoDeinit();
  int DoEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
	       AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
	       AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks );
  int Reconfigure(const char* cdbEntry, const char* chainId);
  void PrintObject( TClonesArray* inClustersArray);

  using AliHLTProcessor::DoEvent;

  int Configure(const char* arguments);
  int SetParams();
	
protected:
  /** copy constructor prohibited */
  AliHLTTRDClusterizerComponent(const AliHLTTRDClusterizerComponent&);
  /** assignment operator prohibited */
  AliHLTTRDClusterizerComponent& operator=(const AliHLTTRDClusterizerComponent&);

  // The size of the output data produced, as a percentage of the input data's size.
  // Can be greater than 100 (%)

  unsigned int fOutputPercentage; // Output volume in percentage of the input
  unsigned int fOutputConst;

  AliHLTTRDClusterizer *fClusterizer; //! Offline derived HLT clusterizer
  AliTRDrecoParam *fRecoParam; //! Offline reco params
  AliRawReaderMemory *fMemReader; //! Input raw data reader
  AliTRDReconstructor *fReconstructor;

  Int_t fRecoParamType;        // default will be the low flux
  Int_t fRecoDataType;         // default will be simulation
  Int_t fRawDataVersion;       // depreceated ?
  Int_t fyPosMethod;           // 0=COG 1=LUT 2=Gauss 
  TString fgeometryFileName;
  Bool_t fProcessTracklets;    // write the L1 tracklets to output
  Bool_t fHLTstreamer;         // use FastStreamer
  Bool_t fTC;                  // using tail cancellation
  Bool_t fHLTflag;             // use HLT flag in reconstructor
  Bool_t fHighLevelOutput;     // do we what to have high level output (only for debuging)
  Bool_t fEmulateHLTClusters;  // for debugging data containers

  ClassDef(AliHLTTRDClusterizerComponent, 5)

};
#endif
