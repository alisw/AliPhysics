/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Salvatore Aiola                                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/
#ifndef ALIEMCALTRIGGERFASTOR_H
#define ALIEMCALTRIGGERFASTOR_H
/**
 * Trigger FastOR data struct
 *
 * @file   AliEMCALTriggerFastOR.h
 * @author Salvatore Aiola
 * @date
 * @brief  Trigger FastOR data struct
 */

#include <Rtypes.h>

class AliEMCALGeometry;

/**
 * @class AliEMCALTriggerFastOR
 * Trigger FastOR data struct
 *
 */
class AliEMCALTriggerFastOR {

 public:

  AliEMCALTriggerFastOR();
  AliEMCALTriggerFastOR(UInt_t L0amp, UInt_t L1amp, Int_t absId, const AliEMCALGeometry* geom);
  AliEMCALTriggerFastOR(UInt_t L0amp, UInt_t L1amp, Int_t globalRow, Int_t glocalCol, const AliEMCALGeometry* geom);
  
  void Initialize(UInt_t L0amp, UInt_t L1amp, Int_t absId, const AliEMCALGeometry* geom);
  void Initialize(UInt_t L0amp, UInt_t L1amp, Int_t globalRow, Int_t glocalCol, Int_t L0time, const AliEMCALGeometry* geom);

  void Initialize(Int_t absId, const AliEMCALGeometry* geom);
  void Initialize(Int_t globalRow, Int_t glocalCol, const AliEMCALGeometry* geom);

  UInt_t                   GetAbsId()        const { return fAbsId     ; }
  UChar_t                  GetGlobalCol()    const { return fGlobalCol ; }
  UChar_t                  GetGlobalRow()    const { return fGlobalRow ; }
  UChar_t                  GetSM()           const { return fSM        ; }
  UChar_t                  GetCol()          const { return fCol       ; }
  UChar_t                  GetRow()          const { return fRow       ; }
  UInt_t                   GetL0Amp()        const { return fL0Amp     ; }
  UInt_t                   GetL1Amp()        const { return fL1Amp     ; }
  Int_t                    GetL0Time()       const { return fL0Time    ; }

  void                     SetL0Amp(UInt_t amp)    { fL0Amp     = amp  ; }
  void                     SetL1Amp(UInt_t amp)    { fL1Amp     = amp  ; }
  void                     SetL0Time(Int_t t)      { fL0Time    = t    ; }

 private:
  /**Abs ID of the trigger FastOR */
  UInt_t                    fAbsId;
  /** Global column of the trigger FastOR */
  UChar_t                   fGlobalCol;
  /** Global row of the trigger FastOR */
  UChar_t                   fGlobalRow;
  /** SM of the trigger FastOR */
  UChar_t                   fSM;
  /** Column of the trigger FastOR  within the SM */
  UChar_t                   fCol;
  /** Row of the trigger FastOR within the SM*/
  UChar_t                   fRow;
  /** ADC counts in the trigger FastOR */
  UInt_t                    fL0Amp;
  /** L1 time sum in the trigger FastOR */
  UInt_t                    fL1Amp;
  /** Time of the trigger FastOR */
  Int_t                     fL0Time;
};

#endif
