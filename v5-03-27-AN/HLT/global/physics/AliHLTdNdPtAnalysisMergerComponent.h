#ifndef ALIHLTDNDPTANALYSISMERGERCOMPONENT_H
#define ALIHLTDNDPTANALYSISMERGERCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTdNdPtAnalysisMergerComponent.h
    @author Sergey Gorbunov
    @brief  Component for merge output from several dNdPt analysis components 
*/

#include "AliHLTProcessor.h"
#include "AliHLTComponentBenchmark.h"

class AlidNdPtEventCuts;
class AlidNdPtAcceptanceCuts;
class AliESDtrackCuts;
class AlidNdPtAnalysisPbPb;


/**
 * @class AliHLTdNdPtAnalysisMergerComponent
 * Component for dNdPt analysis
 */
class AliHLTdNdPtAnalysisMergerComponent : public AliHLTProcessor
{
public:
  /** default constructor */
  AliHLTdNdPtAnalysisMergerComponent();
  /** destructor */
  virtual ~AliHLTdNdPtAnalysisMergerComponent();

  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process

  /** interface function, see AliHLTComponent for description */
  const char* GetComponentID();
  /** interface function, see AliHLTComponent for description */
  void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  /** interface function, see AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();
  /** interface function, see AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  /** interface function, see AliHLTComponent for description */
  AliHLTComponent* Spawn();

protected:

  // Protected functions to implement AliHLTComponent's interface.
  // These functions provide initialization as well as the actual processing
  // capabilities of the component. 

  /** interface function, see AliHLTComponent for description */
  int DoInit( int argc, const char** argv );
  /** interface function, see AliHLTComponent for description */
  int DoDeinit();
  /** interface function, see AliHLTComponent for description */
  int DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& trigData );

  int Reconfigure(const char* cdbEntry, const char* chainId);

  using AliHLTProcessor::DoEvent;
  
private:
  /** copy constructor prohibited */
  AliHLTdNdPtAnalysisMergerComponent(const AliHLTdNdPtAnalysisMergerComponent&);
  /** assignment operator prohibited */
  AliHLTdNdPtAnalysisMergerComponent& operator=(const AliHLTdNdPtAnalysisMergerComponent&);
  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */ 
  void SetDefaultConfiguration();
  int ReadConfigurationString(  const char* arguments );
  int ReadCDBEntry( const char* cdbEntry, const char* chainId );
  int Configure( const char* cdbEntry, const char* chainId, const char *commandLine );

  
  AliHLTUInt32_t fUID;// uID of the component
  AliHLTComponentBenchmark fBenchmark; // benchmarks
  AlidNdPtAnalysisPbPb *fAnalysis; // analysis object

};
#endif
