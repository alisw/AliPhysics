// ************************************************************************
// This file is property of and copyright by the ALICE HLT Project        *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************

#ifndef ALIHLTITSTRACKERCOMPONENT_H
#define ALIHLTITSTRACKERCOMPONENT_H

#include "AliHLTProcessor.h"
#include "AliHLTDataTypes.h"
class AliITStrackerHLT;


/**
 * @class AliHLTITSTrackerComponent
 * The ITS tracker component.
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
    double fFullTime; //* total time for DoEvent() [s]
    double fRecoTime; //* total reconstruction time [s]
    Long_t    fNEvents;  //* number of reconstructed events
    AliITStrackerHLT *fTracker; // the tracker itself

    /** set configuration parameters **/
    void SetDefaultConfiguration();
    int ReadConfigurationString(  const char* arguments );
    int ReadCDBEntry( const char* cdbEntry, const char* chainId );
    int Configure( const char* cdbEntry, const char* chainId, const char *commandLine  );

    ClassDef( AliHLTITSTrackerComponent, 0 );

};
#endif
