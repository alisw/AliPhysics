//-*- Mode: C++ -*-
// $Id: AliHLTMultiplicityCorrelationsComponent $

#ifndef ALIHLTMULTIPLICITYCORRELATIONSCOMPONENT_H
#define ALIHLTMULTIPLICITYCORRELATIONSCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file    AliHLTMultiplicityCorrelationsComponent.h
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

/**
 * @class AliHLTMultiplicityCorrelationsComponent
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
class AliHLTMultiplicityCorrelationsComponent : public AliHLTProcessor {
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** constructor */
  AliHLTMultiplicityCorrelationsComponent();
  
  /** destructor */
  virtual ~AliHLTMultiplicityCorrelationsComponent();

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
  AliHLTMultiplicityCorrelationsComponent(const AliHLTMultiplicityCorrelationsComponent&);

  /** assignment operator prohibited */
  AliHLTMultiplicityCorrelationsComponent& operator=(const AliHLTMultiplicityCorrelationsComponent&);


  /*
   * ---------------------------------------------------------------------------------
   *                              Helper
   * ---------------------------------------------------------------------------------
   */

  /** Set Default Configuartion for track cuts */
  void SetDefaultConfiguration();

  /** Select event as PhysicsSelection */
  Bool_t SelectEvent(AliESDEvent *esdEvent, AliESDVZERO* esdV0);

  /** Checks if event was triggered by the selected trigger classes */
  Bool_t IsEventTriggered();
  
  /** Create selection list of triggers to be checked */
  void CreateTriggerList();

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */
  
  /** ESD track cuts */
  AliESDtrackCuts *fESDTrackCuts;             //! transient

  /** correlations object */
  AliHLTMultiplicityCorrelations *fCorrObj;   //! transient

  /** UID for merging */
  AliHLTUInt32_t fUID;                        // see above

  /** Centrality estimation histogram */
  TH1F *fCentHistV0Mpercentile;               // see above

  /** List of possible trigger descriptor */
  TList *fListTriggerDescriptor;              //! transient
  
  /** List of trigger classnames */
  TList *fListTrigger;                        //! transient

  /** Ptr to CTP data */
  const AliHLTCTPData *fCTPData;              //! transient

  ClassDef(AliHLTMultiplicityCorrelationsComponent, 0)
};
#endif
