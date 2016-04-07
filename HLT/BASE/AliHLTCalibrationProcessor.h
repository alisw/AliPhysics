//-*- Mode: C++ -*-
// $Id$
#ifndef ALIHLTCALIBRATIONPROCESSOR_H
#define ALIHLTCALIBRATIONPROCESSOR_H
/* This file is property of and copyright by the ALICE HLT Project        * 
 * ALICE Experiment at CERN, All rights reserved.                         *
 * See cxx source for full Copyright notice                               */

/**
 * Class takes care of handling and shipping of freshly produced calibration
 * data. this data will be shipped to the FXS of the HLT
 *
 * @file   AliHLTCalibrationProcessor.h
 * @author Jochen Thaeder, Sebastian Bablok
 * @date 
 * @brief  Base class of HLT calibration components.
 */

#include "AliHLTProcessor.h"
#include "AliHLTMessage.h"
#include "AliHLTDataTypes.h"

class AliHLTReadoutList;

/**
 * @class AliHLTCalibrationProcessor
 * Base class of HLT calibration components.
 * The class provides a common interface for the implementation of HLT
 * calibration components. It takes care of handling and shipping of 
 * produced calibration data to the FXS of the HLT. The child class must 
 * implement the functions:
 * - @ref InitCalibration (optional)
 * - @ref ScanArgument (optional)
 * - @ref DeinitCalibration (optional)
 * - @ref ProcessCalibration
 * - @ref ShipDataToFXS
 * - @ref GetComponentID
 * - @ref GetInputDataTypes
 * - @ref GetOutputDataType
 * - @ref GetOutputDataSize
 * - @ref Spawn

 * @ingroup alihlt_component
 */
class AliHLTCalibrationProcessor : public AliHLTProcessor {
  
 public: 
  /** standard constructor */
  AliHLTCalibrationProcessor();

  /** standard destructor */
  virtual ~AliHLTCalibrationProcessor();

 protected:

  /*
   * ######################## PushToFXS #####################
   */

  /**
   * Insert an object into the output. FXS header will be inserted before the root object.
   * @param pObject     pointer to root object
   * @param pDetector   4 byte Detector identifier
   * @param pFileID     name of the file to which the data shall be stored
   * @param pDDLList    pointer to ReadoutList ( class @ref AliHLTReadoutList) of
   *                    participating DDLs. Have to be set or unset with @ref
   *                    AliHLTComponent::EnableDDLBit() and @ref AliHLTComponent::DisableDDLBit(). 
   *			Will be filled automatically if not supplied by the component.
   * @return neg. error code if failed 
   */
  Int_t PushToFXS(const TObject* pObject, const char* pDetector, const char* pFileID, const AliHLTReadoutList* pDDLList);

  /**
   * Insert an object into the output. FXS header will be inserted before the root object.
   * @param pBuffer     pointer to buffer
   * @param iSize       size of the buffer
   * @param pDetector   4 byte Detector identifier
   * @param pFileID     name of the file to which the data shall be stored
   * @param pDDLList    pointer to ReadoutList ( class @ref AliHLTReadoutList) of
   *                    participating DDLs. Have to be set or unset with @ref
   *                    AliHLTComponent::EnableDDLBit() and @ref AliHLTComponent::DisableDDLBit().
   *			Will be filled automatically if not supplied by the component.
   * @return neg. error code if failed 
   */
   Int_t PushToFXS(void* pBuffer, int iSize, const char* pDetector, const char* pFileID, const AliHLTReadoutList* pDDLList);

  /** Constants  */ 
  static const AliHLTUInt32_t fgkFXSProtocolHeaderSize;
  static const AliHLTUInt32_t fgkFXSProtocolHeaderVersion;

  /*
   * ######################## ProcessCalibaration #####################
   */

  /**
   * The low-level data processing method for the component,
   * called for every data event. This is the custom processing 
   * method and can be overloaded by the component.
   * @param [in] evtData       event data structure
   * @param [in] blocks        input data block descriptors
   * @param [in] trigData	  trigger data structure
   * @param [in] outputPtr	  pointer to target buffer
   * @param [in,out] size	  <i>input</i>: size of target buffer
   *            	  <i>output</i>:size of produced data
   * @param [in] outputBlocks  list to receive output block descriptors
   * @return neg. error code if failed
   */
  virtual Int_t ProcessCalibration(const AliHLTComponent_EventData& evtData,
			  const AliHLTComponent_BlockData* blocks,
			  AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr,
			  AliHLTUInt32_t& size,
			  vector<AliHLTComponent_BlockData>& outputBlocks);
 
  /**  
   * The high-level data processing method for the component,
   * called for every data event. This is the custom processing 
   * method and can be overloaded by the component.
   * This is the default processing method; the method is called
   * if no low level @ref ProcessCalibration method is overloaded by the component.
   * @param evtData       event data structure
   * @param trigData	  trigger data structure
   * @return neg. error code if failed
   */
  virtual Int_t ProcessCalibration( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);

  /*
   * ######################## ShipDataToFXS #####################
   */

  /**
   * The low-level data processing method for the component,
   * called for the END_OF_RUN event. This is the custom processing 
   * method and can be overloaded by the component.
   * @param [in] evtData       event data structure
   * @param [in] blocks        input data block descriptors
   * @param [in] trigData	  trigger data structure
   * @param [in] outputPtr	  pointer to target buffer
   * @param [in,out] size	  <i>input</i>: size of target buffer
   *            	  <i>output</i>:size of produced data
   * @param [in] outputBlocks  list to receive output block descriptors
   * @return neg. error code if failed
   */
  virtual Int_t ShipDataToFXS(const AliHLTComponent_EventData& evtData,
			  const AliHLTComponent_BlockData* blocks,
			  AliHLTComponent_TriggerData& trigData, AliHLTUInt8_t* outputPtr,
			  AliHLTUInt32_t& size,
			  vector<AliHLTComponent_BlockData>& outputBlocks);
 
  /**  
   * The high-level data processing method for the component,
   * called for the END_OF_RUN event. This is the custom processing 
   * method and can be overloaded by the component.
   * This is the default processing method; the method is called
   * if no low level @ref ShipDataToFXS method is overloaded by the component.
   * @param evtData       event data structure
   * @param trigData	  trigger data structure
   * @return neg. error code if failed
   */
  virtual Int_t ShipDataToFXS( const AliHLTComponentEventData& evtData, AliHLTComponentTriggerData& trigData);

 private:

  /*
   * ######################## InitCalibration #####################
   */

  /**
   * Internal initialization method, which is not available for child classes. InitCalibration is the
   * corresponding function for classes derived from AliHLTCalibrationProcessor.
   */
  Int_t DoInit( int argc, const char** argv );

  /**
   * Default method for the internal initialization.
   * The method is called by @ref DoInit.
   * This class can be overridden by the child class.
   */
  virtual Int_t InitCalibration();

  /**
   * Scan one argument and adjacent parameters.
   * Can be overloaded by child classes in order to add additional arguments
   * beyond the standard arguments of the calibration processor. The method is called
   * whenever a non-standard argument is recognized. Make sure to return 
   * <tt> -EPROTO </tt> if the argument is not recognized be the child.
   * @param argc           size of the argument array
   * @param argv           agument array for component initialization
   * @return number of processed members of the argv <br>
   *         -EINVAL unknown argument <br>
   *         -EPROTO parameter for argument missing <br>
   */
  virtual Int_t ScanArgument(int argc, const char** argv);

  /*
   * ######################## DeinitCalibration #####################
   */

  /**
   * Internal deinitialization method, which is not available for child classes. DeinitCalibration is the
   * corresponding function for classes derived from AliHLTCalibrationProcessor.
   */
  Int_t DoDeinit();

  /**
   * Default method for the internal clean-up.
   * The method is called by @ref DoDeinit.
   * This class can be overridden by the child class.
   */
  virtual Int_t DeinitCalibration();

  /*
   * ######################## DoEvent #####################
   */
  
  /**
   * The low-level data processing method for the component.
   * It decides wether to call @ref ProcessCalibration or @ref ShipDataToFXS
   * according to the event type - END_OF_RUN / DATA / CALIBRATION
   * If commandline paramater  "-eventmodulo x"is given and x > 0, data will 
   * be also shipped to the FXS with this modulo.
   * @param [in] evtData       event data structure
   * @param [in] blocks        input data block descriptors
   * @param [in] trigData	  trigger data structure
   * @param [in] outputPtr	  pointer to target buffer
   * @param [in,out] size	  <i>input</i>: size of target buffer
   *            	  <i>output</i>:size of produced data
   * @param [in] outputBlocks  list to receive output block descriptors
   * @return neg. error code if failed
   */
  Int_t DoEvent( const AliHLTComponentEventData& evtData,
		 const AliHLTComponentBlockData* blocks, 
		       AliHLTComponentTriggerData& trigData,
		       AliHLTUInt8_t* outputPtr, 
		       AliHLTUInt32_t& size,
		       vector<AliHLTComponentBlockData>& outputBlocks );

  using AliHLTProcessor::DoEvent;
  
  /*
   * ######################## CreateFXSHeader #####################
   */

  /**
   * Insert an object into the output. FXS header will be inserted before the root object.
   * @param pHeader     pointer to AliHLTFXSHeader
   * @param pDetector   4 byte Detector identifier
   * @param pFileID     name of the file to which the data shall be stored
   * @param pDDLList    pointer to ReadoutList ( class @ref AliHLTReadoutList) of
   *                    participating DDLs. Will be filled automatically if not 
   *                    supplied by the component.
   * @return neg. error code if failed 
   */
  Int_t CreateFXSHeader( AliHLTFXSHeader &pHeader, const char* pDetector, const char* pFileID, const AliHLTReadoutList* pDDLList );

  /*
   * ######################## Members #####################
   */

  /** 
   * Event modulo, when Data should shipped to FXS additionally.
   * Default is 0, eventmodulo is deactivated.
   */
  Int_t fEventModulo;                                  // see above

  /** 
   * if kTrue corrupt events will be passed through,
   * if kFalse (default) they will be filtered
   */
  Bool_t fUseCorruptEvents;                            // see above

  /** Event counter */
  Int_t fEventCounter;                                 // see above

  /** 
   * Bit Array of participating DDL Numbers.
   * Scheme: every 4 Bits of each Byte convert digit 
   * (0 - F), before writing to FXS Header
   */
  //Char_t fDDLNumber[gkAliHLTFXSHeaderfDDLNumberSize];// see above
  Char_t fDDLNumber[64];                               // see above

  /** Dummy in order to cope with RC 11 */
  Int_t fDummy;                                        // see above

  ClassDef(AliHLTCalibrationProcessor, 0)

};

#endif // ALIHLTCALIBRATIONPROCESSOR_H
