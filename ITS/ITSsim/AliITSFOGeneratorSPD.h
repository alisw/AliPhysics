#ifndef ALIITS_FOGENERATORSPD_H
#define ALIITS_FOGENERATORSPD_H

/////////////////////////////////////////////////////////////////////
// Author: Henrik Tydesjo                                          //
//                                                                 //
// This class is used to generate Fast-OR signals from SPD chips.  //
//                                                                 //
// This procedure is meant to be used during the digitization,     //
// and will be based on the number of pixels firing in each chip.  //
// The method 'ProcessPixelHit' should be used for each fired      //
// pixel. An efficiency value on Fast-Or signal creation upon a    //
// single fired pixel will then be used. Optionally, there may be  //
// one value per chip or even one value per column. These values   //
// are taken from the class AliITSFOEfficiencySPD, normally placed //
// in OCDB.                                                        //
//                                                                 //
// Through a similar class, AliITSFONoiseSPD, there is a           //
// possibility to apply random noise to the generation of fast-or  //
// signals. This will then be performed by method 'ProcessNoise',  //
// normally called after the processing of the fired pixels.       //
//                                                                 //
// The output signals are represented by the AliITSFOsignalsSPD    //
// class. Basically, it contains a bit map with all the 1200 pixel //
// chips.                                                          //
/////////////////////////////////////////////////////////////////////

#include "AliITSFOEfficiencySPD.h"
#include "AliITSFONoiseSPD.h"
#include "AliITSFOSignalsSPD.h"

class AliITSFOGeneratorSPD {

 public:
  AliITSFOGeneratorSPD();
  AliITSFOGeneratorSPD(AliITSFOEfficiencySPD* ocdbEff, AliITSFONoiseSPD* ocdbNoise);
  AliITSFOGeneratorSPD(const AliITSFOGeneratorSPD& handle);
  virtual ~AliITSFOGeneratorSPD();
  AliITSFOGeneratorSPD& operator=(const AliITSFOGeneratorSPD& handle);

  virtual void   SetEfficiencyAndNoise(AliITSFOEfficiencySPD* ocdbEff, AliITSFONoiseSPD* ocdbNoise);
  virtual void   SetEfficiency(AliITSFOEfficiencySPD* ocdbEff);
  virtual void   SetNoise(AliITSFONoiseSPD* ocdbNoise);
  virtual Bool_t EfficiencyAndNoiseAreSet() {return fOCDBEff!=NULL && fOCDBNoise!=NULL;}

  virtual void   ResetSignals() {fSignals.ResetSignals();}

  virtual void   ProcessPixelHitM(UInt_t module, UInt_t colM, UInt_t rowM);
  virtual void   ProcessPixelHit(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row);
  virtual void   ProcessNoise();

  virtual AliITSFOSignalsSPD* GetFOSignals() {return &fSignals;}


 protected:
  AliITSFOSignalsSPD     fSignals;   // Fast-OR signals object

  AliITSFOEfficiencySPD *fOCDBEff;   // link to FO efficiency obj
  AliITSFONoiseSPD      *fOCDBNoise; // link to FO noise obj

};

#endif
