#ifndef ALICALORAWANALYZERLMS_H
#define ALICALORAWANALYZERLMS_H
/**************************************************************************
 * This file is property of and copyright by                              *
 * the Relativistic Heavy Ion Group (RHIG), Yale University, US, 2009     *
 *                                                                        *
 * Primary Author: Per Thomas Hille <p.t.hille@fys.uio.no>                *
 *                                                                        *
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to p.t.hille@fys.uio.no                             *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/


// Extraction of amplitude and peak position
// FRom CALO raw data using
// Chi square fit

#include "AliCaloRawAnalyzerFitter.h"
#include "AliCaloConstants.h"

using namespace ALTRO;

class  AliCaloRawAnalyzerLMS : public AliCaloRawAnalyzerFitter
{
  friend class AliCaloRawAnalyzerFactory;

 public:
  virtual ~AliCaloRawAnalyzerLMS();
  virtual AliCaloFitResults  Evaluate( const std::vector<AliCaloBunchInfo> &bunchvector, const UInt_t altrocfg1,  const UInt_t altrocfg2 );
  
 private:
  AliCaloRawAnalyzerLMS();
  AliCaloRawAnalyzerLMS(const AliCaloRawAnalyzerLMS & );
  AliCaloRawAnalyzerLMS  & operator = (const AliCaloRawAnalyzerLMS  &);
  ClassDef(AliCaloRawAnalyzerLMS, 2)

};

#endif
