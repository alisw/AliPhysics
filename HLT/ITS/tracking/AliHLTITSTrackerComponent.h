//-*- Mode: C++ -*-
// $Id$
// ************************************************************************
// This file is property of and copyright by the ALICE HLT Project        *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************

///  @file   AliHLTITSTrackerComponent.h
///  @author Sergey Gorbunov <sergey.gorbunov@kip.uni-heidelberg.de>
///  @date   June 2009
///  @brief  An ITS tracker processing component for the HLT

#ifndef ALIHLTITSTRACKERCOMPONENT_H
#define ALIHLTITSTRACKERCOMPONENT_H

#include "AliHLTProcessor.h"
#include "AliHLTDataTypes.h"
#include "AliHLTComponentBenchmark.h"
class AliITStrackerHLT;


/**
 * @class AliHLTITSTrackerComponent
 * The HL ITS tracker component.
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b ITSTracker                              <br>
 * Library: \b libAliHLTITS.so                              <br>
 * Input Data Types:                                        <br> 
 *    kAliHLTDataTypeTrack|kAliHLTDataOriginTPC             <br>
 *    kAliHLTDataTypeClusters|kAliHLTDataOriginITSSSD       <br>
 *    kAliHLTDataTypeClusters|kAliHLTDataOriginITSSPD       <br>
 *    kAliHLTDataTypeClusters|kAliHLTDataOriginITSSDD       <br>
 *      
 * Output Data Types:                                       <br>
 *    kAliHLTDataTypeTrack|kAliHLTDataOriginITS             <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li -config1      <i> teststring   </i> <br>
 *      a configuration argument with one parameter
 * \li -config2                            <br>
 *      a configuration argument without parameters
 *
 * <h2>Default CDB entries:</h2>
 *
 * ITS/Align/Data
 * ITS/Calib/SPDNoisy
 * ITS/Calib/SPDDead
 * ITS/Calib/PITConditions
 * ITS/Calib/CalibSDD
 * ITS/Calib/RespSDD
 * ITS/Calib/DriftSpeedSDD
 * ITS/Calib/DDLMapSDD
 * ITS/Calib/MapsTimeSDD
 * ITS/Calib/NoiseSSD
 * ITS/Calib/GainSSD
 * ITS/Calib/BadChannelsSSD
 *
 * <h2>Performance:</h2>
 * TODO
 *
 * <h2>Memory consumption:</h2>
 * TODO
 *
 * <h2>Output size:</h2>
 * TODO
 * 
 * @ingroup alihlt_its_components
 */
class AliHLTITSTrackerComponent : public AliHLTProcessor
{
  public:
    /** standard constructor */
    AliHLTITSTrackerComponent();

    /** dummy copy constructor, defined according to effective C++ style */
    AliHLTITSTrackerComponent( const AliHLTITSTrackerComponent& );

    /** dummy assignment op, but defined according to effective C++ style */
    AliHLTITSTrackerComponent& operator=( const AliHLTITSTrackerComponent& );

    /** standard destructor */
    virtual ~AliHLTITSTrackerComponent();

    // Public functions to implement AliHLTComponent's interface.
    // These functions are required for the registration process

    /** @see component interface @ref AliHLTComponent::GetComponentID */
    const char* GetComponentID() ;

    /** @see component interface @ref AliHLTComponent::GetInputDataTypes */
    void GetInputDataTypes( vector<AliHLTComponentDataType>& list )  ;

    /** @see component interface @ref AliHLTComponent::GetOutputDataType */
    AliHLTComponentDataType GetOutputDataType() ;

    /** @see component interface @ref AliHLTComponent::GetOutputDataType */
    int  GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);

    /** @see component interface @ref AliHLTComponent::GetOutputDataSize */
    virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier ) ;

    /** @see component interface @ref AliHLTComponent::Spawn */
    AliHLTComponent* Spawn() ;

  protected:

    // Protected functions to implement AliHLTComponent's interface.
    // These functions provide initialization as well as the actual processing
    // capabilities of the component.

    /** @see component interface @ref AliHLTComponent::DoInit */
    int DoInit( int argc, const char** argv );

    /** @see component interface @ref AliHLTComponent::DoDeinit */
    int DoDeinit();

    /** reconfigure **/
    int Reconfigure( const char* cdbEntry, const char* chainId );

    /** @see component interface @ref AliHLTProcessor::DoEvent */
    int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
                 AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr,
                 AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );

  private:

    /** magnetic field */
    double fSolenoidBz;                                            // see above
    AliHLTComponentBenchmark fBenchmark;// benchmark
    AliITStrackerHLT *fTracker; // the tracker itself

    /** set configuration parameters **/
    void SetDefaultConfiguration();
    int ReadConfigurationString(  const char* arguments );
    int ReadCDBEntry( const char* cdbEntry, const char* chainId );
    int Configure( const char* cdbEntry, const char* chainId, const char *commandLine  );

    ClassDef( AliHLTITSTrackerComponent, 0 );

};
#endif
