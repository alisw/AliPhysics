/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Markus Fasel                                          *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#ifndef AliHLTEMCALRawAnalyzerComponentSTU_H
#define AliHLTEMCALRawAnalyzerComponentSTU_H

/**
 * @file AliHLTEMCALTriggerRawDigitMaker.h
 * @brief  Component converting raw data to EMCAL digits
 * @author Markus Fasel <markus.fasel@cern.ch>
 * @date Nov 3, 2015
 */

#include "AliHLTCaloConstants.h"
#include "AliHLTCaloProcessor.h"
#include <vector>

class AliEMCALTriggerData;
class AliHLTCaloMapper;
class AliHLTEMCALGeometry;
class AliHLTEMCALSTURawDigitMaker;
class AliRawReaderMemory;

/**
 * @class AliHTLEMCALTriggerRawDigitMaker
 * @brief Creating EMCAL Raw Digits from the raw data
 */
class AliHLTEMCALRawAnalyzerComponentSTU : public AliHLTCaloProcessor, protected AliHLTCaloConstantsHandler {
public:

  /**
   * Constructor
   */
  AliHLTEMCALRawAnalyzerComponentSTU();

  /**
   * Destructor
   */
  virtual ~AliHLTEMCALRawAnalyzerComponentSTU();

  /**
   * Initialize component
   * @param argc
   * @param argv
   * @return status code
   */
  virtual int DoInit(int argc =0, const char** argv  = 0) ;

  /**
   * Deinit component
   * @return status code
   */
  virtual int DoDeinit();

  /**
   * Get the ID (name) of the component
   * @return name of the component
   */
  virtual const char* GetComponentID();

  /**
   * Get the list of input data types used by this component
   * @param list Input data type list
   */
  virtual void GetInputDataTypes( std::vector <AliHLTComponentDataType>& list);

  /**
   * Get the output data type
   * @return Return output data type
   */
  virtual AliHLTComponentDataType GetOutputDataType();

  /**
   * Get the size of the output data
   * @param constBase Offset
   * @param inputMultiplier Size mulitplier
   */
  virtual void GetOutputDataSize(unsigned long& constBase, double& inputMultiplier);

  /**
   * Create new trigger raw digit maker component
   * @return The new trigger raw digit maker component
   */
  virtual AliHLTComponent* Spawn();

protected:

  /**
   * Process the event
   * @param evtData
   * @param blocks
   * @param trigData
   * @param outputPtr
   * @param size
   * @param outputBlocks
   * @return
   */
  virtual int DoEvent( const AliHLTComponentEventData& evtData, const AliHLTComponentBlockData* blocks,
          AliHLTComponentTriggerData& trigData, AliHLTUInt8_t* outputPtr,
          AliHLTUInt32_t& size, std::vector<AliHLTComponentBlockData>& outputBlocks );

  /**
   * Init EMCAL mapper
   * @param specification Specification
   */
  void InitMapping( const int specification );

  bool CheckInputDataType(const AliHLTComponentDataType &datatype);

private:
  AliHLTEMCALRawAnalyzerComponentSTU(const AliHLTEMCALRawAnalyzerComponentSTU &ref);
  AliHLTEMCALRawAnalyzerComponentSTU &operator=(const AliHLTEMCALRawAnalyzerComponentSTU &ref);

  /** Pointer to the raw data reader which reads from memory */
  AliRawReaderMemory                  *fRawReaderMemoryPtr;            //!transient

  /** STU Raw Digit Maker */
  AliHLTEMCALSTURawDigitMaker         *fSTURawDigitMaker;

  /** Pointer to the geometry class */
  AliHLTEMCALGeometry                 *fGeometry;

  ClassDef(AliHLTEMCALRawAnalyzerComponentSTU, 1);
};

#endif
