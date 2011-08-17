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

#ifndef ALIHLTPHOSDIGITHANDLER_H
#define ALIHLTPHOSDIGITHANDLER_H

#include "offline/AliHLTCaloDigitHandler.h"
#include "AliHLTPHOSDefinitions.h"

class AliHLTPHOSDigitHandler : public AliHLTCaloDigitHandler
{

public:
      
    virtual ~AliHLTPHOSDigitHandler();
    
    static AliHLTPHOSDigitHandler* Instance();
      
    virtual Int_t Init(AliRunLoader* runLoader);

    virtual AliHLTComponentDataType GetDataType() { return AliHLTPHOSDefinitions::fgkDigitDataType; }

protected:
  
    virtual Int_t ConvertDigit(AliDigitNew *digit);    
    
private:
  /** Constructor, private */
  AliHLTPHOSDigitHandler();

  /** The one and only instance of class */
  static AliHLTPHOSDigitHandler *fgkInstance;

  /** Prohibited */
  AliHLTPHOSDigitHandler(const AliHLTPHOSDigitHandler& );
    
  /** Prohibited */
  AliHLTPHOSDigitHandler& operator=(const AliHLTPHOSDigitHandler& );
};

#endif // ALIHLTPHOSDIGITHANDLER_H
