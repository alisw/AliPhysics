//-*- Mode: C++ -*-

// $Id: AliHLTJETFastJetComponent.h  $

#ifndef ALIHLTJETFASTJETCOMPONENT_H
#define ALIHLTJETFASTJETCOMPONENT_H

//* This file is property of and copyright by the ALICE HLT Project        * 
//* ALICE Experiment at CERN, All rights reserved.                         *
//* See cxx source for full Copyright notice                               *

/** @file   AliHLTJETFastJetComponent.h
    @author Jochen Thaeder <thaeder@kip.uni-heidelberg.de>
    @date   
    @brief  Component to run the FastJet jetfinder
*/

#include "AliHLTProcessor.h"

#include "AliHLTJETReader.h"
#include "AliHLTJETReaderHeader.h"

#include "AliHLTJETTrackCuts.h"
#include "AliHLTJETJetCuts.h"

#include "AliHLTJETFastJetFinder.h"
#include "AliHLTJETFastJetHeader.h"

/**
 * @class AliHLTJETFastJetComponent
 * Component to run the FastJet jetfinder
 *
 * <h2>General properties:</h2>
 *
 * Component ID: \b JETFaseJetFinder <br>
 * Library: \b libAliHLTJET.so     <br>
 * Input Data Types: <br>
 *  -  kAliHLTDataTypeMCObject|kAliHLTDataOriginHLT --> class AliHLTMCEvent<br>
 *  -  kAliHLTDataTypeESDObject|kAliHLTDataOriginOffline --> class AliHLTESDEvent<br>
 *  -  kAliHLTDataTypeESDObject|kAliHLTDataOriginHLT --> class AliHLTESDEvent<br>
 * Output Data Types: <br>
 *  - kAliHLTDataTypeJet|kAliHLTDataOriginHLT --> class AliHLTJets<br>
 *
 * <h2>Mandatory arguments:</h2>
 * There are no mandatrory arguments <br>
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 * \li  -algorithm   <i> JetAlgorithm to be run </u><br>
 *       - Possible values : Kt, AntiKt <br>
 *       - Default : Kt <br>
 *
 * \li  -coneRadius    <i> Cone radius for cone finder </i> <br>
 *       - Default : 0.4 <br>
 *
 * \li  -trackCutMinPt <i> min pt for cut on tracks, in GeV/c </i> <br>
 *       - Default : 1.0 <br>
 *
 * \li  -jetCutMinPt   <i> min Et for cut on found jets, in GeV/c </i> <br>
 *       - Default : 15.0 <br>
 *
 * @ingroup alihlt_jet
 * @ingroup alihlt_jet_fastjet
 */

class AliHLTJETFastJetComponent : public AliHLTProcessor {
public:
  
  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */
  
  /** constructor */
  AliHLTJETFastJetComponent();

  /** destructor */
  virtual ~AliHLTJETFastJetComponent();

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
  AliHLTJETFastJetComponent(const AliHLTJETFastJetComponent&);

  /** assignment operator prohibited */
  AliHLTJETFastJetComponent& operator=(const AliHLTJETFastJetComponent&);
  
  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */
  
   /** Ptr to the jet finder */
  AliHLTJETFastJetFinder    *fJetFinder;                      //!transient

  /** Ptr to the jet finder header */
  AliHLTJETFastJetHeader    *fJetHeader;                      //!transient

  /** Ptr to jet reader */ 
  AliHLTJETReader           *fJetReader;                      //!transient
  
  /** Ptr to jet reader header */ 
  AliHLTJETReaderHeader     *fJetReaderHeader;                //!transient

  /** Ptr to track cuts */ 
  AliHLTJETTrackCuts        *fTrackCuts;                      //!transient

  /** Ptr to jet cuts */ 
  AliHLTJETJetCuts          *fJetCuts;                        //!transient

  /** Ptr to jet container holding AliAODJets */
  AliHLTJets                *fJets;       

  ClassDef(AliHLTJETFastJetComponent, 1)
    
};
#endif
