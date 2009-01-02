// XEmacs -*-C++-*-
// $Id$

#ifndef ALIHLTTPCTRACKHISTOCOMPONENT_H
#define ALIHLTTPCTRACKHISTOCOMPONENT_H
//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTTPCTrackHistoComponent.h
    @author Gaute Ovrebekk
    @date   
    @brief  Component for track histo
*/

#include "AliHLTProcessor.h"
#include "AliHLTTPCSpacePointData.h"
#include "AliHLTTPCTrackSegmentData.h"

class AliHLTTPCConfMapper;
class TNtuple;
class AliHLTTPCTrackArray;

/**
 * @class AliHLTTPCTrackHistoComponent
 * Component for ploting proparties of Tracks. 
 * The component gives out two NTuples. One for cluster- and one for track proprties
 * 
 * <h2>General properties:</h2> 
 *
 * Component ID: \b TPCTrackHisto <br>
 * Library: \b libAliHLTTPC.so <br>
 * Input Data Types: AliHLTTPCDefinitions::fgkClustersDataType,
 *                   AliHLTTPCDefinitions::fgkTrackSegmentsDataType or
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
 * Size varibles in Ntuple
 *
 * 
 *
 * @ingroup alihlt_tpc_components
 */
class AliHLTTPCTrackHistoComponent : public AliHLTProcessor
{
public:
  /** default constructor */
  AliHLTTPCTrackHistoComponent();
  /** destructor */
  virtual ~AliHLTTPCTrackHistoComponent();

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

  using AliHLTProcessor::DoEvent;
  
private:
  /** copy constructor prohibited */
  AliHLTTPCTrackHistoComponent(const AliHLTTPCTrackHistoComponent&);
  /** assignment operator prohibited */
  AliHLTTPCTrackHistoComponent& operator=(const AliHLTTPCTrackHistoComponent&);
  /**
   * Configure the component.
   * Parse a string for the configuration arguments and set the component
   * properties.
   */ 
  int Configure(const char* arguments);
 
  void ReadTracks(const AliHLTComponentBlockData* iter,Int_t &tt);

  void PushHisto();
  void FillResidual( UInt_t pos,AliHLTUInt8_t slice,AliHLTUInt8_t patch,Float_t& resy,Float_t& resz);
 
  TNtuple *fClusters;                                              //! transient  
  TNtuple *fTracks;                                                //! transient

  vector<UInt_t> fTrackClusterID[36][6];                           //! transient

  AliHLTTPCTrackArray *fTracksArray;                               //! transient
  AliHLTTPCSpacePointData *fClustersArray[36][6];                  //! transient
  UInt_t fNcl[36][6];                                              //! transient
  
  ClassDef(AliHLTTPCTrackHistoComponent, 1);

};
#endif
