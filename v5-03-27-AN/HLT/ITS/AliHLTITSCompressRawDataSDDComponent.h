//-*- Mode: 
// $Id$

#ifndef ALIHLTITSCOMPRESSRAWDATASDDCOMPONENT_H
#define ALIHLTITSCOMPRESSRAWDATASDDCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTITSCompressRawDataSDDComponent.cxx
    @author Jochen Thaeder <thaeder@kip.uni-heidelberg.de>
    @date   
    @brief  Component to run data compression for SDD
*/

#include "AliHLTProcessor.h"

#include "AliRawReaderMemory.h"

#include "AliITSCompressRawDataSDD.h"

/**
 * @class AliHLTITSCompressRawDataSDDComponent
 * Component to run data compression for SDD
 *
 * @ingroup alihlt_its_components
 */

class AliHLTITSCompressRawDataSDDComponent : public AliHLTProcessor
{
 public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */
  
  /** constructor */
  AliHLTITSCompressRawDataSDDComponent();

  /** destructor */
  virtual ~AliHLTITSCompressRawDataSDDComponent();

  /*
   * ---------------------------------------------------------------------------------
   * Public functions to implement AliHLTComponent's interface.
   * These functions are required for the registration process
   * ---------------------------------------------------------------------------------
   */

  /** interface function, see @ref AliHLTComponent for description */
  const char* GetComponentID();

  /** interface function, see @ref AliHLTComponent for description */
   void GetInputDataTypes( vector<AliHLTComponentDataType>& list);

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponentDataType GetOutputDataType();

  /** interface function, see @ref AliHLTComponent for description */
  virtual void GetOutputDataSize( unsigned long& constBase, double& inputMultiplier );

  /** interface function, see @ref AliHLTComponent for description */
  AliHLTComponent* Spawn();

 protected:
	
  /*
   * ---------------------------------------------------------------------------------
   * Protected functions to implement AliHLTComponent's interface.
   * These functions provide initialization as well as the actual processing
   * capabilities of the component. 
   * ---------------------------------------------------------------------------------
   */
	
  /** Initialization */
  Int_t DoInit( int argc, const char** argv );

  /** DeInitialization */
  Int_t DoDeinit();
  
  /** EventLoop */
  Int_t DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
		 AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr, 
		 AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );

  using AliHLTProcessor::DoEvent;
  
  ///////////////////////////////////////////////////////////////////////////////////
    
    private:
  
  /** copy constructor prohibited */
  AliHLTITSCompressRawDataSDDComponent(const AliHLTITSCompressRawDataSDDComponent&);
  /** assignment operator prohibited */
  AliHLTITSCompressRawDataSDDComponent& operator=(const AliHLTITSCompressRawDataSDDComponent&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  /** the cluster finder object */
  AliITSCompressRawDataSDD* fDataCompressor;                 //!transient

  /** the reader object for data decoding */
  AliRawReaderMemory* fRawReader;                            //!transient
  
  ClassDef(AliHLTITSCompressRawDataSDDComponent, 0)
    
    };
#endif
