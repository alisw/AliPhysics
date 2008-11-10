//-*- Mode: C++ -*-
// $Id$

/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Authors: Oystein Djuvsland <oysteind@ift.uib.no>                       *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALIHLTPHOSPHYSICSANALYZERPEAKFITTER_H
#define ALIHLTPHOSPHYSICSANALYZERPEAKFITTER_H

/**
 * Class does fitting on an histogram
 *
 * @file   AliHLTPHOSPhysicsAnalyzerPeakFitter.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Fitter for PHOS HLT
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "Rtypes.h"
//#include "AliHLTPHOSClusterDataStruct.h"
#include "AliHLTDataTypes.h"
#include "TH1F.h"
#include "TMath.h"

/** 
 * @class AliHLTPHOSPhysicsAnalyzerPeakFitter
 * Fitter for PHOS HLT
 * Makes a fit on a histogram, either a Lorentzian fit or
 * Gaussian fit
 *
 * @ingroup alihlt_phos
 */
class AliHLTPHOSPhysicsAnalyzerPeakFitter
{

 public:

  /** Constructor */
  AliHLTPHOSPhysicsAnalyzerPeakFitter();

  /** Destructor */
  virtual ~AliHLTPHOSPhysicsAnalyzerPeakFitter();

  /** Copy constructor */
  AliHLTPHOSPhysicsAnalyzerPeakFitter(const AliHLTPHOSPhysicsAnalyzerPeakFitter &);

  /** Assignment operator */
  AliHLTPHOSPhysicsAnalyzerPeakFitter & operator = (const AliHLTPHOSPhysicsAnalyzerPeakFitter &) {return *this;}

  /** Set histogram */
  void    SetHistogram(TH1F* histPtr)            { fRootHistPtr = histPtr; }

  /** Fit with a gaussian */
  Int_t   FitGaussian();

  /** Fit with a lorentzian */
  Int_t   FitLorentzian();
  

 private:
  
  /** Factor for low gain */
  Float_t fGainLow;             //COMMENT

  /** Factor for high gain */
  Float_t fGainHigh;            //COMMENT

  /** Pointer to the ROOT histogram */
  TH1F* fRootHistPtr;           //! transient
 
  ClassDef(AliHLTPHOSPhysicsAnalyzerPeakFitter, 1);

};






#endif
