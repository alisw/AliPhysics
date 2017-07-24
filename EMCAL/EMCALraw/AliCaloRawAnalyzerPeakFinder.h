// -*- mode: c++ -*-
#ifndef ALICALORAWANALYZERPEAKFINDER_H
#define ALICALORAWANALYZERPEAKFINDER_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliCaloRawAnalyzerPeakFinder
/// \ingroup EMCALraw
/// \brief  Raw data fitting: Peak Finder
///
/// The Peak-Finder algorithm
/// The amplitude is extracted  as a
/// weighted sum of the samples using the 
/// best possible weights.
/// The weights are calculated only once and the
/// actual extraction of amplitude and peak position
/// is done with a simple vector multiplication, allowing
/// extremely fast computations.
///
/// \author Per Thomas Hille <p.t.hille@fys.uio.no>, Yale. 
//_________________________________________________________________________


#include "AliCaloRawAnalyzer.h"
#include "AliCaloConstants.h"

class AliCaloBunchInfo;
class AliCaloPeakFinderVectors;

class  AliCaloRawAnalyzerPeakFinder : public AliCaloRawAnalyzer
{
  friend class AliCaloRawAnalyzerFactory; // rule checker request
  
 public:
  
  virtual ~AliCaloRawAnalyzerPeakFinder() { ; }
  
  virtual AliCaloFitResults Evaluate( const std::vector<AliCaloBunchInfo> &bunchvector, 
                                      UInt_t altrocfg1, UInt_t altrocfg2 );

 private:
  
  AliCaloRawAnalyzerPeakFinder();
  AliCaloRawAnalyzerPeakFinder(                const AliCaloRawAnalyzerPeakFinder & );
  AliCaloRawAnalyzerPeakFinder  & operator = ( const AliCaloRawAnalyzerPeakFinder & );
  
  void     LoadVectorsOCDB();
  void     CopyVectors(const AliCaloPeakFinderVectors * pfvectors );
  void     ResetVectors();
  void     WriteRootFile() const;
  void     PrintVectors();
  Double_t ScanCoarse( const Double_t *array, Int_t  length ) const ; // Find a rough estimate of peak position and t0
  
  Double_t fPFAmpVectorsCoarse[PF::MAXSTART][PF::SAMPLERANGE][100];   ///< Vectors for Amplitude extraction, first iteration
  Double_t fPFTofVectorsCoarse[PF::MAXSTART][PF::SAMPLERANGE][100];   ///< Vectors for TOF extraction, first iteration
  Double_t fPFAmpVectors      [PF::MAXSTART][PF::SAMPLERANGE][100];   ///< Vectors for Amplitude extraction, second iteration
  Double_t fPFTofVectors      [PF::MAXSTART][PF::SAMPLERANGE][100];   ///< Vectors for TOF extraction, second iteration
  
  AliCaloPeakFinderVectors * fPeakFinderVectors; ///< Collection of Peak-Fincer vectors
  
  bool fRunOnAlien;    ///< Wether or not we are running on the GRID
  bool fIsInitialized; ///< init flag
 
  /// \cond CLASSIMP
  ClassDef( AliCaloRawAnalyzerPeakFinder, 1 ) ;
  /// \endcond
  
};

#endif //ALICALORAWANALYZERPEAKFINDER_H
