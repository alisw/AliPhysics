// -*- mode: c++ -*-
#ifndef ALICALOPEAKFINDERVECTORS_H
#define ALICALOPEAKFINDERVECTORS_H

/* Copyright(c) 1998-2010, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice     */

//_________________________________________________________________________
/// \class AliCaloPeakFinderVectors
/// \ingroup EMCALraw
/// \brief  Container class for Peak Finder vectors
///
/// Container class for Peak Finder vectors
///
/// \author Per Thomas Hille <p.t.hille@fys.uio.no>, Yale. 
//_________________________________________________________________________

#include "TObject.h"

#include "AliCaloConstants.h"

class AliCaloPeakFinderVector;

class  AliCaloPeakFinderVectors : public TObject
{
 
public: 
  
  AliCaloPeakFinderVectors();
  
  virtual ~AliCaloPeakFinderVectors();
  
  void SetVector(const int i, const int j, 
                 const Double_t *const  a, const Double_t *const t,  
                 const Double_t *const ac, const Double_t *const tc );
  
  void GetVector(const int i, const int j, 
                 Double_t *const  a, Double_t *const  t,  
                 Double_t *const ac, Double_t *const tc ) const;
  
  void PrintVectors() const;
  
  void ResetVectors();

 private:
  
  Double_t fPFAmpVC[PF::MAXSTART][PF::SAMPLERANGE][100]; ///< Vectors for Amplitude extraction, first iteration 
  Double_t fPFTofVC[PF::MAXSTART][PF::SAMPLERANGE][100]; ///< Vectors for TOF extraction, first iteration
  Double_t fPFAmpV [PF::MAXSTART][PF::SAMPLERANGE][100]; ///< Vectors for Amplitude extraction, second iteration 
  Double_t fPFTofV [PF::MAXSTART][PF::SAMPLERANGE][100]; ///< Vectors for TOF extraction, second iteration  
  
  /// \cond CLASSIMP
  ClassDef( AliCaloPeakFinderVectors, 2 ) ;
  /// \endcond

};

#endif //ALICALOPEAKFINDERVECTORS_H

