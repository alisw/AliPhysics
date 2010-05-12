// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCCFCOMPARISONCOMPONENT_H
#define ALIHLTTPCCFCOMPARISONCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCCFComparisonComponent.h
    @author Kalliopi Kanaki
    @date   
    @brief  Comparison component for FCF and SCF
*/

#include "AliHLTProcessor.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCTrackSegmentData.h"

class TNtuple;
class TH1F;

/**
 * @class AliHLTTPCCFComparisonComponent
 * Component for plotting proparties of Tracks. 
 * The component gives out 2 NTuples. One for cluster and one for track properties
 * 
 * <h2>General properties:</h2> 
 *
 * Component ID: \b TPCCFComparison <br>
 * Library: \b libAliHLTTPC.so <br>
 * Input Data Types: AliHLTTPCDefinitions::fgkClustersDataType,
 *                   AliHLTTPCDefinitions::fgkTracksDataType <br>
 * Output Data Types: ::kAliHLTDataTypeTNtuple <br> 
 *
 * <h2> Mandatory arguments: </h2>
 * \li No mandaroty arguments. 
 * 
 * <h2> Optional arguments: </h2>
 * 
 * <h2>Configuration:</h2>
 * 
 *
 * <h2>Default CDB entries:</h2>
 * The component has for now no CDB entry
 *
 * <h2>Performance:</h2>
 * Not Tested 
 *
 * <h2>Memory consumption:</h2>
 * Not Tested
 *
 * <h2>Output size:</h2>
 * Size variables in Ntuple
 *
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCCFComparisonComponent : public AliHLTProcessor
{
public:
  /** default constructor */
  AliHLTTPCCFComparisonComponent();
  /** destructor */
  virtual ~AliHLTTPCCFComparisonComponent();

  // Public functions to implement AliHLTComponent's interface.
  // These functions are required for the registration process

  /** interface function, see AliHLTComponent for description */
  const char* GetComponentID();
  /** interface function, see AliHLTComponent for description */
  void GetInputDataTypes(AliHLTComponentDataTypeList& list);
  /** interface function, see AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();
  /** interface function, see AliHLTComponent for description */
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);
  /** interface function, see AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );
  /** interface function, see AliHLTComponent for description */
  AliHLTComponent* Spawn();
  /** interface function, see @ref AliHLTComponent for description */
  void GetOCDBObjectDescription( TMap* const targetMap);

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
  /** inherited from AliHLTComponent: handle re-configuration event */
  int Reconfigure(const char* cdbEntry, const char* chainId);
  /** inherited from AliHLTComponent, scan one argument and its parameters */
  int ScanConfigurationArgument(int argc, const char** argv);

  using AliHLTProcessor::DoEvent;
  
private:
  /** copy constructor prohibited */
  AliHLTTPCCFComparisonComponent(const AliHLTTPCCFComparisonComponent&);
  /** assignment operator prohibited */
  AliHLTTPCCFComparisonComponent& operator=(const AliHLTTPCCFComparisonComponent&);
  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */ 
  
  void ReadTracks(const AliHLTComponentBlockData* iter,Int_t &tt);

  void PushHisto();
 
  Int_t fEvtMod;     //! number of events reached to reset the counter
  Int_t fBufferSize; //! size of circular buffer (number of entries) for the ntuples
 
  TH1F *fMultiplicity;     //! transient (track multiplicity by Z.Y.)

  TNtuple *fClusters;                             //! transient  
  TNtuple *fTracks;                               //! transient

  AliHLTTPCSpacePointData *fFCFClustersArray[36][6]; //! transient
  UInt_t                   fFCFNSpacePoints[36][6];  //! transient

  AliHLTTPCSpacePointData *fSCFClustersArray[36][6]; //! transient
  UInt_t                   fSCFNSpacePoints[36][6];  //! transient
  
  /** the default configuration entry for this component */
  static const char* fgkOCDBEntry; //!transient

  ClassDef(AliHLTTPCCFComparisonComponent, 0);

};
#endif
