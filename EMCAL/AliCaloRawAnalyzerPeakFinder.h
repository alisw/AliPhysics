// -*- mode: c++ -*-
#ifndef ALICALORAWANALYZERPEAKFINDER_H
#define ALICALORAWANALYZERPEAKFINDER_H
/**************************************************************************
 * This file is property of and copyright by                              *
 * the Relativistic Heavy Ion Group (RHIG), Yale University, US, 2009     *
 *                                                                        *
 * Primary Author: Per Thomas Hille  <perthomas.hille@yale.edu>           *
 *                                                                        *
 * Contributors are mentioned in the code where appropriate.              *
 * Please report bugs to perthomas.hille@yale.edu                         *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

// The Peak-Finder algorithm
// The amplitude is extracted  as a
// weighted sum of the samples using the 
// best possible weights.

#include "AliCaloRawAnalyzer.h"
#include "AliCaloConstants.h"

class AliCaloBunchInfo;
class AliCaloPeakFinderVectors;

class  AliCaloRawAnalyzerPeakFinder : public AliCaloRawAnalyzer
{
  friend class AliCaloRawAnalyzerFactory; // shutting up the rule checker
 public:
  virtual ~AliCaloRawAnalyzerPeakFinder();
  virtual AliCaloFitResults Evaluate( const std::vector<AliCaloBunchInfo> &bunchvector, 
				      const UInt_t altrocfg1,  const UInt_t altrocfg2 );

 private:
  AliCaloRawAnalyzerPeakFinder();
  AliCaloRawAnalyzerPeakFinder( const AliCaloRawAnalyzerPeakFinder   & );
  AliCaloRawAnalyzerPeakFinder   & operator = ( const  AliCaloRawAnalyzerPeakFinder  & );
  void     LoadVectorsOCDB();
  void     CopyVectors(const AliCaloPeakFinderVectors *const pfvectors );
  void     ResetVectors();
  void     InitOCDB(bool alien) const;
  void     WriteRootFile() const;
  void PrintVectors();
  Double_t ScanCoarse(const Double_t *const array, const int length ) const ; // Find a rough estimate of peak position and t0
  Double_t fPFAmpVectorsCoarse[PF::MAXSTART][PF::SAMPLERANGE][100];  // Vectors for Amplitude extraction, first iteration
  Double_t fPFTofVectorsCoarse[PF::MAXSTART][PF::SAMPLERANGE][100];  // Vectors for TOF extraction, first iteration
  Double_t fPFAmpVectors[PF::MAXSTART][PF::SAMPLERANGE][100];        // Vectors for Amplitude extraction, second iteration
  Double_t fPFTofVectors[PF::MAXSTART][PF::SAMPLERANGE][100];        // Vectors for TOF extraction, second iteration
  AliCaloPeakFinderVectors  *fPeakFinderVectors; // Collection of Peak-Fincer vectors
  bool fRunOnAlien; // Wether or not we are running on the GRID
  bool fIsInitialized;
  ClassDef( AliCaloRawAnalyzerPeakFinder, 1 )
};


#endif
