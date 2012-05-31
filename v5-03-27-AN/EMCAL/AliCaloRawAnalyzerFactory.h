// -*- mode: c++ -*-

#ifndef ALICALORAWANALYZERFACTORY_H
#define ALICALORAWANALYZERFACTORY_H

/**************************************************************************
 * This file is property of and copyright by the Experimental Nuclear     *
 * Physics Group, Yale University, US 2011                                *
 *                                                                        *
 * Author: Per Thomas Hille <perthomas.hille@yale.edu> for the ALICE      *
 * experiment. Contributors are mentioned in the code where appropriate.  *
 * Please report bugs to  perthomas.hille@yale.edu                        *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#include "AliCaloConstants.h"
using namespace Algo;

class AliCaloRawAnalyzer;

class  AliCaloRawAnalyzerFactory
{
public:
  virtual ~AliCaloRawAnalyzerFactory();
  static AliCaloRawAnalyzer* CreateAnalyzer( const int algo ); 

private:
  AliCaloRawAnalyzerFactory();
};

#endif
