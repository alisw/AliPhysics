
/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Per Thomas Hille, Oystein Djuvsland                   *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTPHOSRAWANALYZERCRUDECOMPONENTV3_H
#define ALIHLTPHOSRAWANALYZERCRUDECOMPONENTV3_H

#include "AliHLTPHOSRawAnalyzerComponentv3.h"


class AliHLTPHOSRawAnalyzerCrudeComponentv3: public AliHLTPHOSRawAnalyzerComponentv3
{
 public:
  AliHLTPHOSRawAnalyzerCrudeComponentv3();
  virtual ~AliHLTPHOSRawAnalyzerCrudeComponentv3();
  AliHLTPHOSRawAnalyzerCrudeComponentv3(const AliHLTPHOSRawAnalyzerCrudeComponentv3 & );
  AliHLTPHOSRawAnalyzerCrudeComponentv3 & operator = (const AliHLTPHOSRawAnalyzerCrudeComponentv3)
  {
    return *this;
  };

  virtual int Deinit();
  virtual const char* GetComponentID();
  virtual AliHLTComponent* Spawn();
};

#endif
