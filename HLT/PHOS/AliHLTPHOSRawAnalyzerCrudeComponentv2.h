
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

#ifndef ALIHLTPHOSRAWANALYZERCRUDECOMPONENTV2_H
#define ALIHLTPHOSRAWANALYZERCRUDECOMPONENTV2_H

#include "AliHLTPHOSRawAnalyzerComponentv2.h"


class AliHLTPHOSRawAnalyzerCrudeComponentv2: public AliHLTPHOSRawAnalyzerComponentv2
{
 public:
  AliHLTPHOSRawAnalyzerCrudeComponentv2();
  virtual ~AliHLTPHOSRawAnalyzerCrudeComponentv2();
  AliHLTPHOSRawAnalyzerCrudeComponentv2(const AliHLTPHOSRawAnalyzerCrudeComponentv2 & );
  AliHLTPHOSRawAnalyzerCrudeComponentv2 & operator = (const AliHLTPHOSRawAnalyzerCrudeComponentv2)
  {
    return *this;
  };

  virtual int Deinit();
  virtual const char* GetComponentID();
  virtual AliHLTComponent* Spawn();
};

#endif
