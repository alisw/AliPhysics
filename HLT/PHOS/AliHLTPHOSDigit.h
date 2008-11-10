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

#ifndef ALIHLTPHOSDIGIT_H
#define ALIHLTPHOSDIGIT_H

/**
 * Digit class for PHOS HLT
 *
 * @file   AliHLTPHOSDigit.h
 * @author Oystein Djuvsland
 * @date
 * @brief  Digit class for PHOS HLT
 */

// see below for class documentation
// or
// refer to README to build package
// or
// visit http://web.ift.uib.no/~kjeks/doc/alice-hlt

#include "TObject.h"
//#include "AliHLTPHOSAltroConfig.h"
#include "AliHLTPHOSBase.h"

//class AliHLTPHOSDigit : public TObject, public AliHLTPHOSAltroConfig

/**
 * @class AliHLTPHOSDigit
 * Digit class for PHOS HLT, used for writing ROOT files
 *
 * @ingroup alihlt_phos
 */
class AliHLTPHOSDigit : public TObject, public AliHLTPHOSBase
{
   
public: 

  /** Constructor */
  AliHLTPHOSDigit();

  /** Destructor */
  virtual ~AliHLTPHOSDigit();

   /** Copy constructor */  
  AliHLTPHOSDigit(const AliHLTPHOSDigit &) : 
    TObject(),
    AliHLTPHOSBase(),
    fX(-1),
    fZ(-1),
    fAmplitude(-1),
    fTime(-1),
    fEnergy(-1),
    fGain(-1),
    fSamples(55),
    fPreSamples(15),
    fTotalSamples(70),
    fDebugVar(-1),
    fData(0),
    fCrazyness(0),
    fBaseline(0)
  {
    //Copy constructor not implemented
  }
  
  /** Assignment */
  AliHLTPHOSDigit & operator = (const AliHLTPHOSDigit)
  {
    //Assignment
    return *this; 
  }

  /** Set x */
  void SetX(Int_t x) { fX = x; }

  /** Set z */
  void SetZ(Int_t z) { fZ = z; }

  /** Set the amplitude in ADC counts */
  void SetAmplitude(Float_t amp) { fAmplitude = amp; }

  /** Set the time in sample counts ? */
  void SetTime(Float_t time) { fTime = time; }

  /** Set the energy in GeV */
  void SetEnergy(Float_t energy) { fEnergy = energy; }

  /** Set the gain */
  void SetGain(Int_t gain) { fGain = gain; }

  /** 
   * Set the raw data 
   * @param rawData is a pointer to an array of raw data 
   */
  void SetRawData(Int_t* rawData);

  /** Set the crazyness */
  void SetCrazyness(Int_t crazyness) { fCrazyness = crazyness; }

  /** Set the baseline value */
  void SetBaseline(Float_t baseline) { fBaseline = baseline; }
  
  /** Set the number of samples */
  void SetSamples(Int_t samples) { fSamples = samples; }

  /** Set the number of pre samples */
  void SetPreSamples(Int_t presamples) { fPreSamples = presamples; }

  /** Reset the digit */
  void ResetDigit();
   
  /** Set the debug variable */
  void SetDebugVar(Int_t val) { fDebugVar = val; }
  
  /** Get x */
  Int_t GetX() { return fX; }

  /** Get z */ 
  Int_t GetZ() { return fZ; }

  /** Get the amplitude */
  Float_t GetAmplitude() { return fAmplitude; }

  /** Get the time */ 
  Float_t GetTime() { return fTime; }

  /** Get the energy */
  Float_t GetEnergy() { return fEnergy; }

  /** Get the gain */
  Int_t GetGain() { return fGain; }

  /** 
   * Get the raw data
   * @return a pointer to an array of raw data 
   */
  Int_t* GetRawData() { return fData; }
 
  /** Get the crazyness */
  Int_t GetCrazyness() {return fCrazyness; }

  /** Get the baseline */
  Float_t GetBaseline() { return fBaseline; }
  
  /** Get number of samples */ 
  Int_t GetSamples() { return fSamples; }

  /** Get number of pre samples */ 
  Int_t GetPreSamples() { return  fPreSamples; }

  /** Get the total number of samples */
  Int_t GetTotalSamples(){ return fNTotalSamples;}
  
  /** Get the debug variable */
  Int_t GetDebugVar() { return fDebugVar; }
  

private:
  
  /** The x coordinate */
  Int_t fX;                     //COMMENT

  /** The z coordinate */
  Int_t fZ;                     //COMMENT

  /** The amplitude in ADC counts*/
  Float_t fAmplitude;           //COMMENT

  /** The time */
  Float_t fTime;                //COMMENT

  /** The energy in GeV */
  Float_t fEnergy;              //COMMENT

  /** The gain */ 
  Int_t fGain;                  //COMMENT

  /** The number of samples */
  Int_t fSamples;               //COMMENT

  /** The number of pre samples */ 
  Int_t fPreSamples;            //COMMENT

  /** The total number of samples */ 
  Int_t fTotalSamples;          //COMMENT
  
  /** A debug variable */
  Int_t fDebugVar;              //COMMENT
  
  /** Pointer to the digit raw data */
  Int_t *fData;                 //[fTotalSamples]

  /** The crazyness */
  Int_t fCrazyness;             //COMMENT

  /** The baseline */
  Float_t fBaseline;            //COMMENT

  ClassDef(AliHLTPHOSDigit, 1);
  
};

#endif
