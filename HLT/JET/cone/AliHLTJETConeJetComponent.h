//-*- Mode: C++ -*-

// $Id: AliHLTJETConeJetComponent.h  $

#ifndef ALIHLTJETCONEJETCOMPONENT_H
#define ALIHLTJETCONEJETCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTJETConeJetComponent.h
    @author Jochen Thaeder <thaeder@kip.uni-heidelberg.de>
    @date   
    @brief  Component to run the ConeJet jetfinder
*/

#include "AliHLTProcessor.h"

#include "AliHLTJETReader.h"
#include "AliHLTJETReaderHeader.h"

#include "AliHLTJETTrackCuts.h"

/**
 * @class AliHLTJETConeJetComponent
 * Component to run the ConeJet jetfinder
 *
 * @ingroup alihlt_jet
 */

class AliHLTJETConeJetComponent : public AliHLTProcessor {
public:
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */
  
  /** constructor */
  AliHLTJETConeJetComponent();

  /** destructor */
  virtual ~AliHLTJETConeJetComponent();

  /*
   * ---------------------------------------------------------------------------------
   * Public functions to implement AliHLTComponent's interface.
   * These functions are required for the registration process
   * ---------------------------------------------------------------------------------
   */

  /** interface function, see @ref AliHLTComponent for description */
  const Char_t* GetComponentID();

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
	
  /** Initialization
   * Overwrites the AliHLTProcessor::DoInit() method.  
   * @param argc           size of the argument array
   * @param argv           agument array for component initialization
   * @return number of processed members of the argv <br>
   *         -EINVAL unknown argument <br>
   *         -EPROTO parameter for argument missing
   */
  Int_t DoInit( Int_t argc, const Char_t** argv );

  /** DeInitialization 
   * Calls also the one of AliHLTProcessor.
   */
  Int_t DoDeinit();
  
  /** EventLoop 
   * Data processing method for the component.
   * The component uses the @ref alihltcomponent-high-level-interface
   * to retrieve and put serialized Root object into the output stream.
   * @param evtData       event data structure
   * @param trigData	  trigger data structure
   * @return
   */
  Int_t DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData );

  using AliHLTProcessor::DoEvent;
  
  ///////////////////////////////////////////////////////////////////////////////////
  
private:
  /*
   * ---------------------------------------------------------------------------------
   * Private functions to implement AliHLTComponent's interface.
   * These functions provide initialization as well as the actual processing
   * capabilities of the component. 
   * ---------------------------------------------------------------------------------
   */
  
  /** copy constructor prohibited */
  AliHLTJETConeJetComponent(const AliHLTJETConeJetComponent&);

  /** assignment operator prohibited */
  AliHLTJETConeJetComponent& operator=(const AliHLTJETConeJetComponent&);
  
  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */
#if 0  
  /** Ptr to the jet finder */
  AliHLTJETConeJetFinder    *fJetFinder;                      //!transient

  /** Ptr to the jet finder header */
  AliHLTJETConeJetHeader    *fJetHeader;                      //!transient
#endif
  /** Ptr to jet reader */ 
  AliHLTJETReader           *fJetReader;                      //!transient
  
  /** Ptr to jet reader header */ 
  AliHLTJETReaderHeader     *fJetReaderHeader;                //!transient

  /** Ptr to track cuts */ 
  AliHLTJETTrackCuts        *fJetTrackCuts;                   //!transient

  ClassDef(AliHLTJETConeJetComponent, 0)
    
};
#endif
