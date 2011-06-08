// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCCOMPMODELCONVERTERCOMPONENT_H
#define ALIHLTTPCCOMPMODELCONVERTERCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCCompModelConverterComponent.h
    @author Timm Steinbeck
    @author changed by J. Wagner
    @date   17-11-2007
    @brief  Declaration of a copy component. */


#include "AliHLTProcessor.h"
#include "AliHLTTPCCompModelConverter.h"
#include "AliHLTStdIncludes.h"

class AliHLTSpacePointContainer;
class AliHLTComponentBenchmark;

/**
 * @class AliHLTTPCCompModelConverterComponent
 * @brief A dummy HLT processing component. 
 *
 * An implementiation of a converter component that
 * takes in clusters and tracks in the standard HLT format
 * and converts them into the Vestbo-format
 * such that the Vestbo compression can then be 
 * applied to these tracks and clusters
 * @ingroup alihlt_tpc
 */
class AliHLTTPCCompModelConverterComponent : public AliHLTProcessor
{
public:
  /** standard constructor */
  AliHLTTPCCompModelConverterComponent();
  /** standard destructor */
  virtual ~AliHLTTPCCompModelConverterComponent();
  
  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process
  
  /** function to get component id 
   * @return const char* pointer to componentid
   */
  const char* GetComponentID();

  /** function to get input data types
   * @param list vecotr of AliHLTComponent_DataType
   */
  void GetInputDataTypes( vector<AliHLTComponent_DataType>& list);

  /** function to get output data type
   * @return AliHLTComponent_DataType
   */
  AliHLTComponent_DataType GetOutputDataType();

  /** function to get output data size
   * @param constBase address of an unsigned long
   * @param inputMultiplier address of a double
   */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );

  /** spawn function
   * @return AliHLTComponent* pointer to instance
   */
  AliHLTComponent* Spawn();
  
protected:
  
  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing
  // capabilities of the component. 
  
  /** initialisation function
   * @param argc integer counting number of input arguments
   * @param argv const char** for parameter values
   * @return zero upon success
   */
  int DoInit( int argc, const char** argv );

  /** deinitialisation function
   * @return zero upon success
   */
  int DoDeinit();

  /** do event function
   * @param evtData      const AliHLTComponent_EventData& to event data
   * @param blocks       const AliHLTComponent_BlockData* to blocks of event data
   * @param trigData     AliHLTComponent_TriggerData& of trigger data
   * @param outputPtr    AliHLTUInt8_t* pointer to output data
   * @param size         AliHLTUInt32_t& of output size
   * @param outputBlocks vector<AliHLTComponent_BlockData>& of output block data
   * @return zero upon success
   */
  int DoEvent( const AliHLTComponent_EventData& evtData, const AliHLTComponent_BlockData* blocks, 
	       AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr, 
	       AliHLTUInt32_t& size, vector<AliHLTComponent_BlockData>& outputBlocks );
  
  /** pointer to model converter instance */
  AliHLTTPCCompModelConverter* fConverter; // pointer to model converter instance 
  /** pointer to model analysis instance */
  AliHLTTPCCompModelAnalysis* fModelAnalysisInstance; // pointer to model analysis instance
  /** filename (including path) where to write model analysis results 
      if set by command line argument -dumptofile result-folder/results */
  TString fDumpFileName; // filename (including path) 
  /** root file name for graphs from track comparison */
  TString fGraphFileName; // root file name for graphs
  
private:
  /** copy constructor prohibited */
  AliHLTTPCCompModelConverterComponent(const AliHLTTPCCompModelConverterComponent&);
  /** assignment operator prohibited */
  AliHLTTPCCompModelConverterComponent& operator=(const AliHLTTPCCompModelConverterComponent&);

  AliHLTComponentBenchmark* GetBenchmarkInstance() const {return fpBenchmark;}

  /** flags to decide wheter to do track or model loss analysis */
  /** switch on model analysis */
  Bool_t fModelAnalysis;            // switch on model analysis
  /** switch on track analysis */
  Bool_t fTrackAnalysis;           // switch on track analysis
  /** flag to check if first or second array is filled */
  Bool_t fFillingFirstTrackArray;  // 1 if first array is processed, 0 for second filling process

  /// input cluster handler
  AliHLTSpacePointContainer* fInputClusters; //! input cluster handler

  /// benchmark
  AliHLTComponentBenchmark* fpBenchmark; //! benchmark instance

  ClassDef(AliHLTTPCCompModelConverterComponent, 0)
    
    };
#endif
