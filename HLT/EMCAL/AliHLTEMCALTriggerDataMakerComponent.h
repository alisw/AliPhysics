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
#ifndef ALIHLTEMCALTRIGGERDATAMAKERCOMPONENT_H
#define ALIHLTEMCALTRIGGERDATAMAKERCOMPONENT_H

/**
 * @file AliHLTEMCALTriggerRawDigitMaker.h
 * @brief  Component converting raw data to EMCAL digits
 * @author Markus Fasel <markus.fasel@cern.ch>
 * @date Nov 3, 2015
 */

#include <vector>
#include "AliHLTCaloConstants.h"
#include "AliHLTCaloProcessor.h"
#include "AliHLTCaloTriggerRawDigitDataStruct.h"
#include "AliHLTEMCALSTUHeaderStruct.h"

struct AliHLTCaloTriggerDataStruct;
class AliHLTEMCALGeometry;

/**
 * @class AliHTLEMCALTriggerRawDigitMaker
 * @brief Combining data from STU and TRU in EMCAL trigger data
 */
class AliHLTEMCALTriggerDataMakerComponent : public AliHLTCaloProcessor, protected AliHLTCaloConstantsHandler {
public:

  /**
   * Constructor
   */
  AliHLTEMCALTriggerDataMakerComponent();

  /**
   * Destructor
   */
  virtual ~AliHLTEMCALTriggerDataMakerComponent();

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

  bool CheckInputDataType(const AliHLTComponentDataType &datatype);

  void ReadSTUData(AliHLTEMCALSTUHeaderStruct *headerptr, AliHLTCaloTriggerRawDigitDataStruct *dataptr);
  void ReadTRUData(UShort_t ndigits, AliHLTCaloTriggerRawDigitDataStruct * triggerdata);
  Int_t MakeTriggerData(AliHLTCaloTriggerDataStruct *outputptr, AliHLTUInt32_t &availableSize);
  void Reset();

  void CombineTRUSTUDigit(AliHLTCaloTriggerRawDigitDataStruct &target, const AliHLTCaloTriggerRawDigitDataStruct &trudigit, const AliHLTCaloTriggerRawDigitDataStruct &studigit);
  void ConvertRawDigit(AliHLTCaloTriggerDataStruct *target, const AliHLTCaloTriggerRawDigitDataStruct * source, Int_t col, Int_t row);

private:
  enum {
    kMaxChannels = 5962
  };
  AliHLTEMCALTriggerDataMakerComponent(const AliHLTEMCALTriggerDataMakerComponent &ref);
  AliHLTEMCALTriggerDataMakerComponent &operator=(const AliHLTEMCALTriggerDataMakerComponent &ref);

  /** Pointer to the geometry class */
  AliHLTEMCALGeometry                       *fGeometry;

  AliHLTEMCALSTUHeaderStruct                fSTUHeader;
  Short_t                                   fNRawDigitsTRU;
  Short_t                                   fRawIndexesTRU[kMaxChannels];
  Short_t                                   fNRawDigitsSTU;
  Short_t                                   fRawIndexesSTU[kMaxChannels];
  Short_t                                   fMaxChannel;
  AliHLTCaloTriggerRawDigitDataStruct       fTRURawDigitBuffer[kMaxChannels];
  AliHLTCaloTriggerRawDigitDataStruct       fSTURawDigitBuffer[kMaxChannels];

  ClassDef(AliHLTEMCALTriggerDataMakerComponent, 1);
};

#endif
