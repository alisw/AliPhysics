//-*- Mode: C++ -*-
// $Id: AliHLTAnaManagerComponent $

#ifndef ALIHLTANAMANAGERCOMPONENT_H
#define ALIHLTANAMANAGERCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file    AliHLTAnaManagerComponent.h
    @author  Jochen Thaeder <jochen@thaeder.de>
    @brief   Component for Multiplicty Correlations
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "AliHLTProcessor.h"

class TH1F;
class TList;

class AliESDEvent;
class AliESDVZERO;
class AliESDtrackCuts;
class AliHLTCTPData;
class AliHLTMultiplicityCorrelations;
class AliHLTGlobalTriggerDecision;
class AliAnalysisManager;
class AliHLTTestInputHandler;

/**
 * @class AliHLTAnaManagerComponent
 * Create Correlations for Multiplicities
 * 
 * <h2>General properties:</h2>
 *
 * Component ID: \b MultiplicityCorrelations <br>
 * Library: \b libAliHLTGlobal.so     <br>
 * Input Data Types:  @ref kAliHLTDataTypeESDObject <br>
 * Output Data Types: @ref kAliHLTDataTypeTObject|kAliHLTDataOriginHLT <br>
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
 * \li -minpt    <i> pt  </i> <br>
 *      minimum pt - pt range
 * \li -maxpt    <i> pt  </i> <br>
 *      maximum pt - pt range
 * \li -min-ldca    <i> dca  </i> <br>
 *      minimum longitudinal dca to reference point
 * \li -max-ldca    <i> dca  </i> <br>
 *      maximum longitudinal dca to reference point
 * \li -min-tdca    <i> dca  </i> <br>
 *      minimum transverse dca to reference point
 * \li -max-tdca    <i> dca  </i> <br>
 *      maximum transverse dca to reference point
 * \li -etarange    <i> eta  </i> <br>
 *      +/- eta range
 *
 * \li -binningVzero    <i> bins min max  </i> <br>
 *       bins (Int_t), minBin (Float_t), maxBin (Float_t)
 * \li -binningTpc      <i> bins min max  </i> <br>
 *       bins (Int_t), minBin (Float_t), maxBin (Float_t)
 * \li -binningZdc      <i> bins min max  </i> <br>
 *       bins (Int_t), minBin (Float_t), maxBin (Float_t)
 * \li -binningZnp     <i> bins min max  </i> <br>
 *       bins (Int_t), minBin (Float_t), maxBin (Float_t)
 * \li -binningZem     <i> bins min max  </i> <br>
 *       bins (Int_t), minBin (Float_t), maxBin (Float_t)
 * \li -binningCalo    <i> bins min max  </i> <br>
 *       bins (Int_t), minBin (Float_t), maxBin (Float_t)
 *
 * \li -addTrigger     <i> TriggerClass beginning (eg CPBI1)  </i> <br>

 * <h2>Default CDB entries:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <tt>HLT/ConfigGlobal/MultiplicityCorrelations</tt>
 * \li -TObjString object holding a string with the configuration parameters
 *      currently empty 
 *
 * <h2>Performance:</h2>
 *
 * <h2>Memory consumption:</h2>
 *
 * <h2>Input size:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * \li pp: 
 *
 * <h2>Output size:</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * \li pp: Average : 
 *
 * <h2>Macros Tests</h2>
 * <!-- NOTE: ignore the \li. <i> and </i>: it's just doxygen formatting -->
 *
 * <tt>macros/makeConfigurationObjectMultiplicityCorrelations.C</tt>
 * \li - Create configuration TObjString
 *
 * <tt>macros/HLTMultiplicityCorrelationsTest.C</tt>
 * \li - Test macro for test in off-line environment
 *
 * <tt>macros/runMultiplicityCorrelationsTest.sh</tt>
 * \li - Run Test macro HLTMultiplicityCorrelationsTest.C
 *
 * @ingroup alihlt_physics
 */
class AliHLTAnaManagerComponent : public AliHLTProcessor {
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** constructor */
  AliHLTAnaManagerComponent();
  
  /** destructor */
  virtual ~AliHLTAnaManagerComponent();

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
  AliHLTAnaManagerComponent(const AliHLTAnaManagerComponent&);

  /** assignment operator prohibited */
  AliHLTAnaManagerComponent& operator=(const AliHLTAnaManagerComponent&);


  /*
   * ---------------------------------------------------------------------------------
   *                              Helper
   * ---------------------------------------------------------------------------------
   */

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */
  
  /** UID for merging */
  AliHLTUInt32_t fUID;                        // see above

  AliAnalysisManager *fAnalysisManager;        // Manger

  AliHLTTestInputHandler *fInputHandler;    // input handler

  ClassDef(AliHLTAnaManagerComponent, 0)
};
#endif
