//-*- Mode: C++ -*-
// $Id$

/**************************************************************************
 * This file is property of and copyright by the ALICE HLT Project        * 
 * All rights reserved.                                                   *
 *                                                                        *
 * Primary Authors: Oystein Djuvsland                                     *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          * 
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/** 
 * @file   AliHLTPHOSPhysicsAnalyzerSpectrum.h
 * @author Oystein Djuvsland
 * @date 
 * @brief  Invariant mass spectrum from 2 gammas  */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#ifndef ALIHLTPHOSPHYSICSANALYZERSPECTRUM_H
#define ALIHLTPHOSPHYSICSANALYZERSPECTRUM_H

#include "AliHLTPHOSPhysicsAnalyzer.h"
#include "Rtypes.h"

class AliHLTPHOSRecPointContainerStruct;


/**
 * @class AliHLTPHOSPhysicsAnalyzerSpectrum
 * Invariant mass spectrum from 2 gammas
 * @ingroup alihlt_phos
 */

class AliHLTPHOSPhysicsAnalyzerSpectrum : public AliHLTPHOSPhysicsAnalyzer
{
 public:

  /** Constructor */
  AliHLTPHOSPhysicsAnalyzerSpectrum();

  /** Copy constructor */
  AliHLTPHOSPhysicsAnalyzerSpectrum(const AliHLTPHOSPhysicsAnalyzerSpectrum &);

  /** Assignment */
  AliHLTPHOSPhysicsAnalyzerSpectrum & operator = (const AliHLTPHOSPhysicsAnalyzerSpectrum)
    {
      //Assignment
      return *this; 
    }


  /** Destructor */
  virtual ~AliHLTPHOSPhysicsAnalyzerSpectrum();

  /** Set the threshold for the 2 photon energies */
  Int_t SetThreshold(Float_t photonEnergy0, Float_t photonEnergy1);

  /** 
   * Evaluate the distance between 2 rec points
   * @return the distance
   */
  Float_t EvalDistance();
 

  /** 
   * Analyze the event 
   * @param recPointsArrayPtr is an array of rec points in the event
   * @param nRecPoints is the number of rec points in the event
   */
    virtual void Analyze(AliHLTPHOSRecPointContainerStruct* recPointsArrayPtr, Int_t nRecPoints);

 private:
  /** Position of the first cluster */
  Float_t* fPos0Ptr;                        //! transient

  /** Position of the second cluster */
  Float_t* fPos1Ptr;                        //! transient

  /** Cut thresholds */
  Float_t* fThresholdPtr;                   //! transient

  /** Energy of the clusters */
  Float_t* fEnergyPtr;                      //! transient

  ClassDef(AliHLTPHOSPhysicsAnalyzerSpectrum, 1);
  
};

#endif
 
