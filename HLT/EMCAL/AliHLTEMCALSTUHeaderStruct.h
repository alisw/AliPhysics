#ifndef ALIHLTEMCALSTUHEADERSTRUCT_H
#define ALIHLTEMCALSTUHEADERSTRUCT_H

struct AliHLTEMCALSTUHeaderStruct {
  /** L1 thresholds from raw data */
  Int_t    fL1Threshold[4];
  /** L1 threshold components */
  Int_t    fL1V0[2];
  /** Validation flag for L1 data */
  Int_t    fL1FrameMask;
  /** Number of consecutive STU Raw Digits */
  Int_t     fNRawDigits;
};

#endif
