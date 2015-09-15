//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTVZERORECOCOMPONENT_H
#define ALIHLTVZERORECOCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file    AliHLTVZERORecoComponent.h
    @author  Jochen Thaeder <jochen@thaeder.de>
    @brief   VZERO reconstruction component
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


#include "AliHLTProcessor.h"

class TTree;

class AliRunInfo;
class AliESDVZERO;
class AliRawReaderMemory;
class AliVZERORecoParam;
class AliVZEROReconstructor;

/**
 * @class AliHLTVZERORecoComponent
 * Reconstruction of VZERO data
 * 
 * <h2>General properties:</h2>
 *
 * Component ID: \b VZEROReconstruction <br>
 * Library: \b libAliHLTVZERO.so     <br>
 * Input Data Types:  @ref kAliHLTDataTypeDDLRaw <br>
 * Output Data Types: @ref kAliHLTDataTypeESDContent|kAliHLTDataOriginVZERO <br>
 *
 * <h2>Mandatory arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting --> 
 *
 * <h2>Optional arguments:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Configuration:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <h2>Default CDB entries:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <tt>HLT/ConfigVZERO/VZEROReconstruction</tt>
 * \li -TObjString object holding a string with the configuration parameters
 *      currently empty 
 *
 * <tt>GRP/GRP/Data</tt>
 * \li -GRP object - run information
 *
 *  <tt>GRP/CTP/CTPtiming</tt>
 * \li -GRP object - CTP information
 *
 * <tt>GRP/CTP/TimeAlign</tt>
 * \li -GRP object - CTP information
 * 
 * <tt>GRP/Calib/LHCClockPhase</tt>
 * \li -GRP object - time calibration
 *
 * <tt>VZERO/Calib/Data</tt>
 * \li -VZERO calibration object
 *
 * <tt>VZERO/Calib/TimeDelays</tt>
 * \li -VZERO calibration object
 *
 * <tt>VZERO/Calib/TimeSlewing</tt>
 * \li -VZERO calibration object
 *
 * <h2>Performance:</h2>
 *
 * <h2>Memory consumption:</h2>
 *
 * <h2>Input size:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * \li pp: 5968 Byte
 *
 * <h2>Output size:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * \li pp: Average : 1.8 kByte
 *
 * <h2>Macros Tests</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <tt>macros/makeConfigurationObjectVZEROReconstruction.C</tt>
 * \li - Create configuration TObjString
 *
 * <tt>macros/HLTVZEROTest.C</tt>
 * \li - Test macro for VZERO test in off-line environment
 *
 * <tt>macros/runVZEROTest.sh</tt>
 * \li - Run Test macro HLTVZEROTest.C
 *
 * @ingroup alihlt_vzero
 */
class AliHLTVZERORecoComponent : public AliHLTProcessor {
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** constructor */
  AliHLTVZERORecoComponent();
  
  /** destructor */
  virtual ~AliHLTVZERORecoComponent();

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
  int GetOutputDataTypes(AliHLTComponentDataTypeList& tgtList);

  /** interface function, see @ref AliHLTComponent for description */
  void GetOutputDataSize( ULong_t& constBase, Double_t& inputMultiplier );

  /** interface function, see @ref AliHLTComponent for description */
  void GetOCDBObjectDescription( TMap* const targetMap);

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

  // AliHLTComponent interface functions

  /** interface function, see @ref AliHLTComponent for description */
  Int_t DoInit( Int_t argc, const Char_t** argv );

  /** interface function, see @ref AliHLTComponent for description */
  Int_t DoDeinit();

  /** interface function, see @ref AliHLTComponent for description */
  Int_t DoEvent( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);

  using AliHLTProcessor::DoEvent;

  /** interface function, see @ref AliHLTComponent for description */
  Int_t ScanConfigurationArgument(Int_t argc, const Char_t** argv);

  /** interface function, see @ref AliHLTComponent for description */
  Int_t Reconfigure(const Char_t* cdbEntry, const Char_t* chainId);

  /** interface function, see @ref AliHLTComponent for description */
  Int_t ReadPreprocessorValues(const Char_t* modules);
 
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
  AliHLTVZERORecoComponent(const AliHLTVZERORecoComponent&);

  /** assignment operator prohibited */
  AliHLTVZERORecoComponent& operator=(const AliHLTVZERORecoComponent&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */
  
  /** runInfo Object */
  AliRunInfo            *fRunInfo;            // see above

  /** VZERO reco param instance */
  AliVZERORecoParam     *fVZERORecoParam;     //! transient

  /** VZERO reconstructor instance */
  AliVZEROReconstructor *fVZEROReconstructor; //! transient

  /** Rawreader instance */
  AliRawReaderMemory    *fRawReader;          //! transient
  
  ClassDef(AliHLTVZERORecoComponent, 0)
};
#endif
