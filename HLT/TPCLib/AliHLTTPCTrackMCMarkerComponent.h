// ************************************************************************
// This file is property of and copyright by the ALICE HLT Project        *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************


#ifndef ALIHLTTPTRACKMCMARKERCOMPONENT_H
#define ALIHLTTPTRACKMCMARKERCOMPONENT_H

/** @file   AliHLTTPCTrackMCMarkerComponent.h
    @author Matthias Kretz
    @date
    @brief  HLT TPC CA global merger component.
*/

#include "AliHLTProcessor.h"
class AliHLTTPCClusterMCData;

/**
 * @class AliHLTTPCTrackMCMarkerComponent
 * Assigns MC labels to the TPC tracks 
 */
class AliHLTTPCTrackMCMarkerComponent : public AliHLTProcessor
{
  public:
    /**
     * Constructs a AliHLTTPCTrackMCMarkerComponent.
     */
    AliHLTTPCTrackMCMarkerComponent();

    /**
     * Destructs the AliHLTTPCTrackMCMarkerComponent
     */
    virtual ~AliHLTTPCTrackMCMarkerComponent() {};

    // Public functions to implement AliHLTComponent's interface.
    // These functions are required for the registration process

    /**
     * @copydoc AliHLTComponent::GetComponentID
     */
    const char *GetComponentID();

    /**
     * @copydoc AliHLTComponent::GetInputDataTypes
     */
    void GetInputDataTypes( AliHLTComponentDataTypeList &list );

    /**
     * @copydoc AliHLTComponent::GetOutputDataType
     */
    AliHLTComponentDataType GetOutputDataType();

    /**
     * @copydoc AliHLTComponent::GetOutputDataSize
     */
    virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );

    /**
     * @copydoc AliHLTComponent::Spawn
     */
    AliHLTComponent *Spawn();

  protected:

    // Protected functions to implement AliHLTComponent's interface.
    // These functions provide initialization as well as the actual processing
    // capabilities of the component.

    /**
     * @copydoc AliHLTComponent::DoInit
     */
    int DoInit( int argc, const char **argv );

    /**
     * @copydoc AliHLTComponent::DoDeinit
     */
    int DoDeinit();

    /** reconfigure **/
    int Reconfigure( const char* cdbEntry, const char* chainId );

    /**
     * @copydoc @ref AliHLTProcessor::DoEvent
     */
    int DoEvent( const AliHLTComponentEventData &evtData, const AliHLTComponentBlockData *blocks,
                 AliHLTComponentTriggerData &trigData, AliHLTUInt8_t *outputPtr,
                 AliHLTUInt32_t &size, AliHLTComponentBlockDataList &outputBlocks );

    using AliHLTProcessor::DoEvent;

  private:

    static AliHLTTPCTrackMCMarkerComponent fgAliHLTTPCTrackMCMarkerComponent;

    // disable copy
    AliHLTTPCTrackMCMarkerComponent( const AliHLTTPCTrackMCMarkerComponent & );
    AliHLTTPCTrackMCMarkerComponent &operator=( const AliHLTTPCTrackMCMarkerComponent & );

    /** set configuration parameters **/
    void SetDefaultConfiguration();
    int ReadConfigurationString(  const char* arguments );
    int ReadCDBEntry( const char* cdbEntry, const char* chainId );
    int Configure( const char* cdbEntry, const char* chainId, const char *commandLine );

    /**
     * Get MC label for a track
     * @param hits array of track hit Id's
     * @param nHits n of hits
     * @return neg. -1 if failed
     */
    int GetTrackMCLabel( unsigned int *hits, int nHits );

    /** array of pointers to cluster MC labels **/

    AliHLTTPCClusterMCData *fClusterLabels[36*6]; //! cluster MC labels for each TPC patch

    ClassDef( AliHLTTPCTrackMCMarkerComponent, 0 )
};

#endif
