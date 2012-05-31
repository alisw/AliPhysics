//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTITSDIGITPUBLISHERCOMPONENT_H
#define ALIHLTITSDIGITPUBLISHERCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTITSDigitPublisherComponent.cxx
    @author Kenneth Aamodt, Sergey Gorbunov
    @date   
    @brief  Component to run the offline clusterfinder.
*/

#include "AliHLTOfflineDataSource.h"
class AliRunLoader;
class AliITSLoader;
class TClonesArray;
class TTree;
/**
 * @class AliHLTITSDigitPublisherComponent
 *
 */
class AliHLTITSDigitPublisherComponent : public AliHLTOfflineDataSource
{
 public:
  
  
  /** 
   * constructor
   */
  AliHLTITSDigitPublisherComponent();

  /** destructor */
  virtual ~AliHLTITSDigitPublisherComponent();

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
  Int_t GetEvent(const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);

  ///////////////////////////////////////////////////////////////////////////////////
    
 private:
  /** copy constructor prohibited */
  AliHLTITSDigitPublisherComponent(const AliHLTITSDigitPublisherComponent&);
  /** assignment operator prohibited */
  AliHLTITSDigitPublisherComponent& operator=(const AliHLTITSDigitPublisherComponent&);

  AliRunLoader * fRunLoader;
  AliITSLoader * fITSLoader;
  Int_t fNumberOfEvents;
  Int_t fEventNumber;
  TTree *tD;

  ClassDef(AliHLTITSDigitPublisherComponent, 0)
    
};
#endif
