//-*- Mode: C++ -*-
// $Id$
// ************************************************************************
// This file is property of and copyright by the ALICE HLT Project        *
// ALICE Experiment at CERN, All rights reserved.                         *
// See cxx source for full Copyright notice                               *
//                                                                        *
//*************************************************************************

///  @file   AliHLTITSVertexerSPDComponent.h
///  @author Sergey Gorbunov <sergey.gorbunov@kip.uni-heidelberg.de>
///  @date   Nov 2009
///  @brief  An ITS pixel vertexer component for the HLT

#ifndef ALIHLTITSVERTEXERSPDCOMPONENT_H
#define ALIHLTITSVERTEXERSPDCOMPONENT_H

#include "AliHLTProcessor.h"
#include "AliHLTDataTypes.h"


/**
 * @class AliHLTITSVertexerSPDComponent
 * The HLT ITS SPD z-vertexer component.
 * The vertexer uses approximate XY position of the beam
 * It can treat initial beam offset up to 1.5 cm in XY
 * The component does a calibration of the beam position 
 * every n==fAutoCalibration events
 * 
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b ITSVertexerSPD                            <br>
 * Library: \b libAliHLTITS.so                              <br>
 * Input Data Types:                                        <br> 
 *    kAliHLTDataTypeClusters|kAliHLTDataOriginITSSPD       <br>
 *      
 * Output Data Types:                                       <br>
 *    kAliHLTDataTypeESDVertex|kAliHLTDataOriginITS            <br>
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
 * TODO
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
class AliHLTITSVertexerSPDComponent : public AliHLTProcessor
{
  public:
    /** standard constructor */
    AliHLTITSVertexerSPDComponent();

    /** dummy copy constructor, defined according to effective C++ style */
    AliHLTITSVertexerSPDComponent( const AliHLTITSVertexerSPDComponent& );

    /** dummy assignment op, but defined according to effective C++ style */
    AliHLTITSVertexerSPDComponent& operator=( const AliHLTITSVertexerSPDComponent& );

    /** standard destructor */
    virtual ~AliHLTITSVertexerSPDComponent();

    // Public functions to implement AliHLTComponent's interface.
    // These functions are required for the registration process

    /** @see component interface @ref AliHLTComponent::GetComponentID */
    const char* GetComponentID() ;

    /** @see component interface @ref AliHLTComponent::GetInputDataTypes */
    void GetInputDataTypes( vector<AliHLTComponentDataType>& list )  ;

    /** @see component interface @ref AliHLTComponent::GetOutputDataType */
    AliHLTComponentDataType GetOutputDataType() ;

    /** @see component interface @ref AliHLTComponent::GetOutputDataTypes */
    int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);

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

  struct AliHLTITSVZCluster{
    float fX, fY, fZ;
  };
  
  Double_t fZRange;// Z range for the vertex seearch
  Double_t fZBinSize; // size of the Z bin [cm] 
  double fRunVtx[3];// default vertex position
  double fFullTime; //* total time for DoEvent() [s]
  double fRecoTime; //* total reconstruction time [s]
  Long_t fNEvents;  //* number of reconstructed events
  double *fSum[9]; // coefficients for the LSM method
  double *fSumW; // sum of weights per Z bin 
  int *fSumN; // N entries per Z bin
  double fZMin; // Z of the first bin ( == -fZRange)
  int fNZBins; // N of bins

  /** set configuration parameters **/
  void SetDefaultConfiguration();
  int ReadConfigurationString(  const char* arguments );
  int ReadCDBEntry( const char* cdbEntry, const char* chainId );
  int Configure( const char* cdbEntry, const char* chainId, const char *commandLine  );
  
  ClassDef( AliHLTITSVertexerSPDComponent, 0 );
  
};
#endif
