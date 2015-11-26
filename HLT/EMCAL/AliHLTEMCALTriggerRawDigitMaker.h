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
#ifndef ALIHLTEMCALTRIGGERRAWDIGITMAKER_H
#define ALIHLTEMCALTRIGGERRAWDIGITMAKER_H

#include <TObject.h>
#include <vector>
#include "AliCaloBunchInfo.h"
#include "AliHLTEMCALTriggerRawDigitDataStruct.h"
#include "AliHLTLogging.h"

class AliCaloRawAnalyzerFakeALTRO;
class AliCaloRawStreamV3;
class AliEMCALTriggerData;
class AliEMCALTriggerSTUDCSConfig;
class AliEMCALTriggerSTURawStream;
class AliRawReader;

class AliHLTEMCALGeometry;

class AliHLTEMCALTriggerRawDigitMaker : public TObject, public AliHLTLogging {
public:
  /**
   * Constructor
   */
  AliHLTEMCALTriggerRawDigitMaker();

  /**
   * Destructor
   */
  virtual ~AliHLTEMCALTriggerRawDigitMaker();

  /**
   * Add bunch list of type AliCaloBunchInfo
   * @param bunchlist
   */
  void Add(const std::vector<AliCaloBunchInfo> &bunchlist);

  /**
   * Connect I/O
   * @param reader ALICE Raw Reader
   * @param in Calorimeter Raw Stream
   * @param inSTU EMCAL STU Raw Stream
   * @param digitcont container for raw digits
   * @param data
   */
  void SetIO(AliRawReader* reader, AliCaloRawStreamV3& in, AliEMCALTriggerSTURawStream& inSTU, AliEMCALTriggerData* data);

  const std::vector<AliHLTEMCALTriggerRawDigitDataStruct> &GetRawDigits()  const { return fRawDigitVector; }

  /**
   * Post process digits
   */
  void PostProcess();

  /**
   * Set the geometry ptr for this class
   * @param geo Geometry ptr
   */
  void SetGeometry(const AliHLTEMCALGeometry *geo) { fkGeometryPtr = geo; }

  /**
   * Reset indices
   */
  void Reset();

protected:
  static const Int_t fgkSTUEqId;
  static const Int_t fgkNRawDigits;

  /**
   * Get the raw digit for a given index. If not existing, create it.
   * @param index
   * @return
   */
  AliHLTEMCALTriggerRawDigitDataStruct &GetRawDigit(Int_t index);

  /** Ptr to the EMCAL geometry */
  const AliHLTEMCALGeometry                               *fkGeometryPtr;
  /** Raw reader */
	AliRawReader                                            *fRawReader;
  /** Calo raw stream */
  AliCaloRawStreamV3                                      *fCaloRawStream;
  /** STU raw stream */
  AliEMCALTriggerSTURawStream                             *fSTURawStream;
  /** Raw analyzer */
	AliCaloRawAnalyzerFakeALTRO                             *fRawAnalyzer;
	/** DCS Config */
	AliEMCALTriggerSTUDCSConfig                             *fDCSConfigSTU;
	/** Raw digit container */
	std::vector<AliHLTEMCALTriggerRawDigitDataStruct>       fRawDigitVector;
	/** Trigger data */
	AliEMCALTriggerData                                     *fTriggerData;
	/** Raw digit indexes */
  Int_t                                                   fRawDigitIndex[5952];

private:
  AliHLTEMCALTriggerRawDigitMaker(const AliHLTEMCALTriggerRawDigitMaker &);
  AliHLTEMCALTriggerRawDigitMaker &operator=(const AliHLTEMCALTriggerRawDigitMaker &);

	ClassDef(AliHLTEMCALTriggerRawDigitMaker, 1);
};

#endif
