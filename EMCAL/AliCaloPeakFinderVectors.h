#ifndef ALICALOPEAKFINDERVECTORS_H
#define ALICALOPEAKFINDERVECTORS_H

/**************************************************************************
 * This file is property of and copyright by the Relativistic Heavy Ion   *
 * Group (RHIG), Department of Physics Yale University, US, 2010          *
 *                                                                        *
 * Author: Per Thomas Hille <perthomas.hille@yale.edu> for the ALICE EMCAL*
 * project. Contributors are mentioned in the code where appropriate.     *
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

//Container class for Peak Finder vectors

#include "TObject.h"
#include "AliCaloPeakFinderConstants.h"
using namespace PeakFinderConstants;


class AliCaloPeakFinderVector;


class  AliCaloPeakFinderVectors : public TObject
{
 public: 
  AliCaloPeakFinderVectors();
  virtual ~AliCaloPeakFinderVectors();
  void SetVector(const int i, const int j, const Double_t  *const a, const Double_t *const t,  
		 const Double_t *const ac, const Double_t *const tc );
  void GetVector(const int i, const int j, Double_t *const a, Double_t *const  t,  
		 Double_t *const ac, Double_t *const tc ) const;
  void PrintVectors() const;
  void ResetVectors();

 private:
  Double_t fPFAmpVC[MAXSTART][SAMPLERANGE][100]; // Vectors for Amplitude extraction, first iteration 
  Double_t fPFTofVC[MAXSTART][SAMPLERANGE][100]; // Vectors for TOF extraction, first iteration
  Double_t fPFAmpV[MAXSTART][SAMPLERANGE][100];  // Vectors for Amplitude extraction, second iteration 
  Double_t fPFTofV[MAXSTART][SAMPLERANGE][100];  // Vectors for TOF extraction, second iteration  
 
  ClassDef( AliCaloPeakFinderVectors, 1 )
    
};

#endif

