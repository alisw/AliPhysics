//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTITSCLUSTERFINDERSDDCOMPONENT_H
#define ALIHLTITSCLUSTERFINDERSDDCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTITSClusterFinderSDDComponent.cxx
    @author Gaute Øvrebekk <st05886@alf.uib.no>
    @date   
    @brief  Component to run the offline clusterfinder.
*/

#include "AliHLTProcessor.h"
#include "AliRawReaderMemory.h"
#include "AliITSClusterFinderV2SDD.h"
#include "TClonesArray.h"
#include "AliITSDetTypeRec.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSgeom.h"
#include "AliITSInitGeometry.h"
#include "AliITSsegmentationSDD.h"

/**
 * @class AliHLTITSClusterFinderSDDComponent
 * Component to run the offline clusterfinder.
 *
 * @ingroup alihlt_its_components
 */

class AliHLTITSClusterFinderSDDComponent : public AliHLTProcessor
{
 public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */
  
  /** constructor */
  AliHLTITSClusterFinderSDDComponent();

  /** destructor */
  virtual ~AliHLTITSClusterFinderSDDComponent();

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
  
  Int_t fNModules;             // total number of modules

  /** EventLoop */
  //Int_t DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks, 
  //	 AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr
  //	 ,AliHLTUInt32_t& size, AliHLTComponentBlockList& outputBlocks);

  Int_t DoEvent( const AliHLTComponentEventData& /*evtData*/, AliHLTComponentTriggerData& /*trigData*/);
		
  using AliHLTProcessor::DoEvent;
 
  ///////////////////////////////////////////////////////////////////////////////////
    
    private:
  
  /** copy constructor prohibited */
  AliHLTITSClusterFinderSDDComponent(const AliHLTITSClusterFinderSDDComponent&);
  /** assignment operator prohibited */
  AliHLTITSClusterFinderSDDComponent& operator=(const AliHLTITSClusterFinderSDDComponent&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  /** the cluster finder object */
  AliITSClusterFinderV2SDD* fClusterFinder;                   //!transient

  /** the reader object for data decoding */
  AliRawReaderMemory* fRawReader;                             //!transient

  AliITSDetTypeRec* fDettype;                                 //!transient

  TClonesArray** fClusters;                                   //!transient

  AliITSgeom* fgeom;                                          //!transient

  AliITSInitGeometry* fgeomInit;                              //!transient

  AliITSsegmentationSDD* fSeg;                                //!transient
  
  ClassDef(AliHLTITSClusterFinderSDDComponent, 1)
    
    };
#endif
