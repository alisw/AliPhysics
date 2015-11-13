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
#ifndef AliHLTEMCALSTURawDigitMaker_H
#define AliHLTEMCALSTURawDigitMaker_H

#include <AliHLTCaloTriggerRawDigitDataStruct.h>
#include <TObject.h>
#include <vector>
#include "AliCaloBunchInfo.h"
#include "AliHLTLogging.h"

class AliCaloRawStreamV3;
class AliEMCALTriggerData;
class AliEMCALTriggerSTUDCSConfig;
class AliEMCALTriggerSTURawStream;
class AliRawReader;

class AliHLTEMCALGeometry;

class AliHLTEMCALSTURawDigitMaker : public TObject, public AliHLTLogging {
public:
  /**
   * Constructor
   */
  AliHLTEMCALSTURawDigitMaker();

  /**
   * Destructor
   */
  virtual ~AliHLTEMCALSTURawDigitMaker();

  /**
   * Post process digits
   */
  void ProcessSTUStream(AliEMCALTriggerSTURawStream *stustream);

  /**
   * Set the geometry ptr for this class
   * @param geo Geometry ptr
   */
  void SetGeometry(const AliHLTEMCALGeometry *geo) { fkGeometryPtr = geo; }

  /**
   * Reset indices
   */
  void Reset();

  Int_t WriteRawDigitsBuffer(AliHLTCaloTriggerRawDigitDataStruct *bufferptr) const;

  Int_t GetNumberOfRawDigits() const { return fNRawDigits; }

  const AliEMCALTriggerData *GetTriggerData() const { return fTriggerData; }

protected:
  enum{
    fgkNRawDigits = 5952,
    fgkMaxFastorModule = 144
  };


  /**
   * Get the raw digit for a given index. If not existing, create it.
   * @param index
   * @return
   */
  AliHLTCaloTriggerRawDigitDataStruct &GetRawDigit(Int_t index);

  /** Ptr to the EMCAL geometry */
  const AliHLTEMCALGeometry                               *fkGeometryPtr;
	/** DCS Config */
	AliEMCALTriggerSTUDCSConfig                             *fDCSConfigSTU;
	/** Raw digit container */
	/** Trigger data */
	AliEMCALTriggerData                                     *fTriggerData;

  AliHLTCaloTriggerRawDigitDataStruct                       fRawDigitBuffer[fgkNRawDigits];
  /** Raw digit indexes */
  Short_t                                                   fNRawDigits;
  Short_t                                                   fRawDigitIndex[fgkNRawDigits];

private:
  AliHLTEMCALSTURawDigitMaker(const AliHLTEMCALSTURawDigitMaker &);
  AliHLTEMCALSTURawDigitMaker &operator=(const AliHLTEMCALSTURawDigitMaker &);

	ClassDef(AliHLTEMCALSTURawDigitMaker, 1);
};

#endif
