#ifndef ALICALORAWANALYZERPEAKFINDER_H
#define ALICALORAWANALYZERPEAKFINDER_H


/**************************************************************************
 * This file is property of and copyright by                              *
 * the Relativistic Heavy Ion Group (RHIG), Yale University, US, 2009     *
 *                                                                        *
 * Primary Author: Per Thomas Hille  <perthomas.hille@yale.edu>           *
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


// The Peak-Finder algorithm
// The amplitude is extracted  as a
// weighted sum of the samples using the 
// best possible weights.


#include "AliCaloRawAnalyzer.h"

#define MAXSTART 3
#define SAMPLERANGE 15
#define SHIF 0.5

class AliCaloBunchInfo;


class  AliCaloRawAnalyzerPeakFinder : public AliCaloRawAnalyzer
{
 public:
  AliCaloRawAnalyzerPeakFinder();
  virtual ~AliCaloRawAnalyzerPeakFinder();
  virtual AliCaloFitResults Evaluate( const vector<AliCaloBunchInfo> &bunchvector, const UInt_t altrocfg1,  const UInt_t altrocfg2 );

 private:
  AliCaloRawAnalyzerPeakFinder( const AliCaloRawAnalyzerPeakFinder   & );
  AliCaloRawAnalyzerPeakFinder   & operator = ( const  AliCaloRawAnalyzerPeakFinder  & );

  void LoadVectors();
  double *fPFAmpVectors[MAXSTART][SAMPLERANGE]; // Vectors for Amplitude extraction 
  double *fPFTofVectors[MAXSTART][SAMPLERANGE]; // Vectors for TOF extraction

  //  double fTof; 
  //  double fAmp;

  double fAmpA[3]; // The  amplitude of the signal (eveluate 3 times using 3 differtnt phase shifts of the input samples )
  // double fAmp2;
  // double fAmp3;

  ClassDef( AliCaloRawAnalyzerPeakFinder, 1 )

};

#endif
