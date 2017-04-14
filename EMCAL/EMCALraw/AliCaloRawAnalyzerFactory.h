// -*- mode: c++ -*-

#ifndef ALICALORAWANALYZERFACTORY_H
#define ALICALORAWANALYZERFACTORY_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliCaloRawAnalyzerFactory
/// \ingroup EMCALraw
/// \brief  Raw data fitters handling
///
/// Raw data fitters handling. 
/// The selected fitter is executed here.
///
/// \author Per Thomas Hille <p.t.hille@fys.uio.no>, Yale. 
//_________________________________________________________________________

#include "AliCaloConstants.h"
using namespace Algo;

class AliCaloRawAnalyzer;

class  AliCaloRawAnalyzerFactory
{
  
public:
  
  virtual ~AliCaloRawAnalyzerFactory() { ; }
  
  static AliCaloRawAnalyzer* CreateAnalyzer( const int algo ); 

private:
  
  AliCaloRawAnalyzerFactory();
  
};

#endif // ALICALORAWANALYZERFACTORY_H
