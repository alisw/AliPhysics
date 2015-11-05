#ifndef ALIHLTEMCALTRURAWDIGITMAKER_H
#define ALIHLTEMCALTRURAWDIGITMAKER_H

#include "AliHLTCaloTriggerRawDigitDataStruct.h"
#include "AliHLTLogging.h"
#include <vector>

class AliCaloRawStreamV3;
class AliCaloBunchInfo;
class AliHLTEMCALGeometry;

class AliHLTEMCALTRURawDigitMaker : protected AliHLTLogging {
public:
  AliHLTEMCALTRURawDigitMaker();
  virtual ~AliHLTEMCALTRURawDigitMaker();

  void Initialize(Int_t runno);
  void Add(const std::vector<AliCaloBunchInfo> &bunchlist);
  Int_t WriteRawDigitsBuffer(AliHLTCaloTriggerRawDigitDataStruct *bufferptr) const;

  Int_t GetNumberOfRawDigits() const { return fNRawDigits; }

  /**
   * Connect I/O
   * @param in Calorimeter Raw Stream
   */
  void SetRawReader(AliCaloRawStreamV3* in);

  const std::vector<AliHLTCaloTriggerRawDigitDataStruct> GetRawDigits()  const;

  void Reset();

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

  AliCaloRawStreamV3                                        *fCaloRawStream;
  const AliHLTEMCALGeometry                                 *fGeometryPtr;
  AliHLTCaloTriggerRawDigitDataStruct                       fRawDigitBuffer[fgkNRawDigits];
  /** Raw digit indexes */
  Short_t                                                   fNRawDigits;
  Short_t                                                   fRawDigitIndex[fgkNRawDigits];

  ClassDef(AliHLTEMCALTRURawDigitMaker, 1);
};

#endif
