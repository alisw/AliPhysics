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

#ifndef ALIHLTPHOSPHYSICSANALYZER_H
#define ALIHLTPHOSPHYSICSANALYZER_H

/**
 * Class is intended to be a base class for physics analysis for 
 * PHOS in HLT
 *
 * @file   AliHLTPHOSPhysicsAnalyzer.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Physics analysis base class
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "Rtypes.h"
#include "AliHLTPHOSConstants.h"

using namespace PhosHLTConst;

class TObjArray;
class TH1F;
class AliHLTPHOSRecPointDataStruct;
class AliHLTPHOSRecPointContainerStruct;

const Float_t kCRYSTALSIZE = 2.25;

/**
 * @class AliHLTPHOSPhysicsAnalyzer
 * Base class for physics analysis for PHOS in HLT
 * @ingroup alihlt_phos
 */

class AliHLTPHOSPhysicsAnalyzer
{ 
 public:

  /** Constructor */
  AliHLTPHOSPhysicsAnalyzer();

  /** Destructor */
  virtual ~AliHLTPHOSPhysicsAnalyzer();

  /** Copy constructor */  
  AliHLTPHOSPhysicsAnalyzer(const AliHLTPHOSPhysicsAnalyzer &):
  fRecPointsPtr(0), 
  fRootHistPtr(0), 
  fPHOSRadius(0)
  {
    //Copy constructor not implemented
  }
  
  /** Assignment */
  AliHLTPHOSPhysicsAnalyzer & operator = (const AliHLTPHOSPhysicsAnalyzer &) 
  {
    //Assignement
    return *this;
  }
  
  void SetHistogram(TH1F* histPtr) {fRootHistPtr = histPtr;}

  /** 
   * Get the local position of a reconstruction point in the PHOS module 
   * @param recPointPtr is a pointer to the reconstruction point
   * @param locPositionPtr is a pointer to the array of local coordinates
   */
  void LocalPosition(AliHLTPHOSRecPointDataStruct* recPointPtr, Float_t* locPositionPtr);


  /** 
   * Get the global position of a reconstruction point in ALICE
   * @param recPointPtr is a pointer to the reconstruction point
   * @param positionPtr is a pointer to the array where to fill the coord
   */
  void GlobalPosition(AliHLTPHOSRecPointDataStruct* recPointPtr, Float_t* positionPtr);
  
  /** 
   * Get the global position of local coordinates in ALICE
   * @param locPositionPtr is a pointer to the array of local coordinates
   * @param positionPtr is a pointer to the array of global coordinates
   * @param module is the module number (0 - 4)
   */
  void GlobalPosition(Float_t* locPositionPtr , Float_t* positionPtr, Int_t module);

  virtual void WriteHistogram(const Char_t* fileName = "histogram.root");

  /** 
   * Abstract function, for doing analysis
   * @param recPointsArrayPtr is an array of pointers to recPoints
   * @param nRecPoints is the numbers of recPoints
   */
  virtual void Analyze(AliHLTPHOSRecPointContainerStruct* recPointsArrayPtr, Int_t nRecPoints);

 protected:

  /** Pointer to the clusters to be analyzed */
  TObjArray* fRecPointsPtr;                         //! transient

  /** Pointer to the histograms which is to be filled */
  TH1F* fRootHistPtr;                              //! transient 

 private:

  /** Parameters for calculating global position */
  Float_t fRotParametersCos[5];                        

  /** Parameters for calculating global position */
  Float_t fRotParametersSin[5];                   
    
  /** Distance from the IP to the crystals */
  Float_t fPHOSRadius;                                

  ClassDef(AliHLTPHOSPhysicsAnalyzer,1);
};

#endif
