/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        *
 * ALICE Experiment at CERN, All rights reserved.                         *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland <oysteind@ift.uib.no>               *
 *                  for The ALICE HLT Project.                            *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTEMCALDIGITHANDLER_H
#define ALIHLTEMCALDIGITHANDLER_H

#include "offline/AliHLTCaloDigitHandler.h"
#include "AliHLTEMCALDefinitions.h"

class AliEMCALCalibData;
class AliDigitNew;
class AliHLTEMCALDigitHandler : public AliHLTCaloDigitHandler
{

public:
      
    virtual ~AliHLTEMCALDigitHandler();
    
    static AliHLTEMCALDigitHandler* Instance();
      
    virtual Int_t Init(AliRunLoader* runLoader);
    
    virtual AliHLTComponentDataType GetDataType() { return AliHLTEMCALDefinitions::fgkDigitDataType; }
    

protected:
  
    virtual Int_t ConvertDigit(AliDigitNew *digit);    
    
    int GetGainsFromCDB();

    
private:
  
  /** Constructor, private */
  AliHLTEMCALDigitHandler();

  /** The one and only instance of class */
  static AliHLTEMCALDigitHandler *fgkInstance;
  
  /** Calibration data */
  AliEMCALCalibData *fCalibData;
  
  /** Prohibited */
  AliHLTEMCALDigitHandler(const AliHLTEMCALDigitHandler& );
    
  /** Prohibited */
  AliHLTEMCALDigitHandler& operator=(const AliHLTEMCALDigitHandler& );
};

#endif // ALIHLTEMCALDIGITHANDLER_H
