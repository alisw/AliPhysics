//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTTOFRECOCOMPONENT_H
#define ALIHLTTOFRECOCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file    AliHLTTOFRecoComponent.h
    @author  Jochen Thaeder <jochen@thaeder.de>
    @brief   TOF reconstruction component
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


#include "AliHLTProcessor.h"

class TTree;

class AliRunInfo;
class AliESDTOF;
class AliRawReaderMemory;
class AliTOFRecoParam;
class AliTOFReconstructor;
class AliTOFRawStream;

/**
 * @class AliHLTTOFRecoComponent
 * Reconstruction of TOF data
 * 
 * <h2>General properties:</h2>
 *
 * Component ID: \b TOFReconstruction <br>
 * Library: \b libAliHLTTOF.so     <br>
 * Input Data Types:  @ref kAliHLTDataTypeDDLRaw <br>
 * Output Data Types: @ref kAliHLTDataTypeESDContent|kAliHLTDataOriginTOF <br>
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
 * <tt>HLT/ConfigTOF/TOFReconstruction</tt>
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
 * <tt>TOF/Calib/Data</tt>
 * \li -TOF calibration object
 *
 * <tt>TOF/Calib/TimeDelays</tt>
 * \li -TOF calibration object
 *
 * <tt>TOF/Calib/TimeSlewing</tt>
 * \li -TOF calibration object
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
 * <tt>macros/makeConfigurationObjectTOFReconstruction.C</tt>
 * \li - Create configuration TObjString
 *
 * <tt>macros/HLTTOFTest.C</tt>
 * \li - Test macro for TOF test in off-line environment
 *
 * <tt>macros/runTOFTest.sh</tt>
 * \li - Run Test macro HLTTOFTest.C
 *
 * @ingroup alihlt_vzero
 */
class AliHLTTOFRecoComponent : public AliHLTProcessor {
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** constructor */
  AliHLTTOFRecoComponent();
  
  /** destructor */
  virtual ~AliHLTTOFRecoComponent();

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
  Int_t DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
		 AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr,
		 AliHLTUInt32_t& size, vector<AliHLTComponentBlockData>& outputBlocks );

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
  AliHLTTOFRecoComponent(const AliHLTTOFRecoComponent&);

  /** assignment operator prohibited */
  AliHLTTOFRecoComponent& operator=(const AliHLTTOFRecoComponent&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */

  UShort_t fDebugLevel;                              //! set debug checks/output level, 0: debug off
  
  /** runInfo Object */
  AliRunInfo            *fRunInfo;            // see above

  /** TOF reco param instance */
  AliTOFRecoParam     *fTOFRecoParam;     //! transient

  /** TOF reconstructor instance */
  AliTOFReconstructor *fTOFReconstructor; //! transient

  /** Rawreader instance */
  AliRawReaderMemory    *fRawReaderMem;          //! transient

  /** TOF Rawreader instance */
  AliTOFRawStream    *fTOFRawStream;          //! transient
  
  ClassDef(AliHLTTOFRecoComponent, 0)
};
#endif
