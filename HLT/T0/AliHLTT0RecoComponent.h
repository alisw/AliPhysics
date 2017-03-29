//-*- Mode: C++ -*-
// $Id$

#ifndef ALIHLTT0RECOCOMPONENT_H
#define ALIHLTT0RECOCOMPONENT_H

/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/** @file    AliHLTT0RecoComponent.h
    @author  Alla Maevskaya <Alla.Maevskaya@cern.ch>
    @brief   T0 reconstruction component
*/

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt


#include "AliHLTProcessor.h"

class TTree;
class TH1F;
class TObjArray;

class AliRunInfo;
class AliESDTZERO;
class AliRawReaderMemory;
class AliT0RecoParam;
class AliT0Reconstructor;
class AliRawReader;
class AliT0RawReader;

/**
 * @class AliHLTT0RecoComponent
 * Reconstruction of T0 data
 * 
 * <h2>General properties:</h2>
 *
 * Component ID: \b T0Reconstruction <br>
 * Library: \b libAliHLTT0.so     <br>
 * Input Data Types:  @ref kAliHLTDataTypeDDLRaw <br>
 * Output Data Types: @ref kAliHLTDataTypeESDContent|kAliHLTDataOriginT0 <br>
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
 * <tt>HLT/ConfigT0/T0Reconstruction</tt>
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
 * <tt>T0/Calib/Data</tt>
 * \li -T0 calibration object
 *
 * <tt>T0/Calib/TimeDelays</tt>
 * \li -T0 calibration object
 *
 * <tt>T0/Calib/TimeSlewing</tt>
 * \li -T0 calibration object
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
 * <tt>macros/makeConfigurationObjectT0Reconstruction.C</tt>
 * \li - Create configuration TObjString
 *
 * <tt>macros/HLTT0Test.C</tt>
 * \li - Test macro for T0 test in off-line environment
 *
 * <tt>macros/runT0Test.sh</tt>
 * \li - Run Test macro HLTT0Test.C
 *
 * @ingroup alihlt_vzero
 */
class AliHLTT0RecoComponent : public AliHLTProcessor {
public:

  /*
   * ---------------------------------------------------------------------------------
   *                            Constructor / Destructor
   * ---------------------------------------------------------------------------------
   */

  /** constructor */
  AliHLTT0RecoComponent();
  
  /** destructor */
  virtual ~AliHLTT0RecoComponent();

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
  AliHLTT0RecoComponent(const AliHLTT0RecoComponent&);

  /** assignment operator prohibited */
  AliHLTT0RecoComponent& operator=(const AliHLTT0RecoComponent&);

  /*
   * ---------------------------------------------------------------------------------
   *                             Members - private
   * ---------------------------------------------------------------------------------
   */
  
  /** runInfo Object */
  AliRunInfo            *fRunInfo;            // see above

  /** T0 reco param instance */
  AliT0RecoParam     *fT0RecoParam;     //! transient

  /** T0 reconstructor instance */
  AliT0Reconstructor *fT0Reconstructor; //! transient

  /** Rawreader instance */
  AliRawReaderMemory    *fRawReader;          //! transient

  AliESDTZERO *fESDTZERO;
  void RecT0Raw(AliRawReader*rawReader);
  void GetMeanAndSigma(TH1F* hist,  Float_t &mean, Float_t &sigma);
 
  
  // my members
  Int_t fNevent;
  TH1F* fhTimeDiff[24];
  TH1F* fhCFD[24];
  TH1F* fhT0[4];
  Double_t fVertexSPDz;
  TObjArray  *fWalk;   //walk correction function
  TObjArray  *fT0CalibHisto;
  Float_t fMeanCFD[24];
  Float_t fDiffCFD[24];
  Float_t fT0shift[4];
  
  ClassDef(AliHLTT0RecoComponent, 0)
};
#endif
