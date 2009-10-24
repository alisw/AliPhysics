//-*- Mode: C++ -*-
// ************************************************************************
// This file is property of and copyright by the ALICE HLT Project        *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************

///  @file   AliHLTTPCdEdxComponent.h
///  @author Sergey Gorbunov <sergey.gorbunov@kip.uni-heidelberg.de>
///  @date   June 2009
///  @brief  An ITS tracker processing component for the HLT

#ifndef ALIHLTTPCDEDXCOMPONENT_H
#define ALIHLTTPCDEDXCOMPONENT_H

#include "AliHLTProcessor.h"
#include "AliHLTDataTypes.h"
class AliTPCclusterMI;

/**
 * @class AliHLTTPCdEdxComponent
 * The dEdx calculator component for the HLT TPC
 *
 */
class AliHLTTPCdEdxComponent : public AliHLTProcessor
{
public:
  /** standard constructor */
  AliHLTTPCdEdxComponent();
  
  /** dummy copy constructor, defined according to effective C++ style */
  AliHLTTPCdEdxComponent( const AliHLTTPCdEdxComponent& );
  
  /** dummy assignment op, but defined according to effective C++ style */
  AliHLTTPCdEdxComponent& operator=( const AliHLTTPCdEdxComponent& );
  
  /** standard destructor */
  virtual ~AliHLTTPCdEdxComponent();

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
    double fStatTime; //* total time for DoEvent() [s]   
    Long_t    fStatNEvents;  //* number of reconstructed events
  static const Int_t fkNPatches = 36*6; // number of patches in TPC
  AliTPCclusterMI *fPatchClusters[fkNPatches]; //! arrays of cluster data for each TPC patch
  Int_t fNPatchClusters[fkNPatches]; //! N of clusters for each TPC patch

    /** set configuration parameters **/
    void SetDefaultConfiguration();
    int ReadConfigurationString(  const char* arguments );
    int ReadCDBEntry( const char* cdbEntry, const char* chainId );
    int Configure( const char* cdbEntry, const char* chainId, const char *commandLine  );

    ClassDef( AliHLTTPCdEdxComponent, 0 );

};
#endif
