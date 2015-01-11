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

#include "AliITSFOGeneratorSPD.h"
#include "AliITSRawStreamSPD.h"
#include <TRandom.h>

AliITSFOGeneratorSPD::AliITSFOGeneratorSPD() :
  fSignals(), fOCDBEff(NULL), fOCDBNoise(NULL)
{
  // default constructor
}
//______________________________________________________________________
AliITSFOGeneratorSPD::AliITSFOGeneratorSPD(AliITSFOEfficiencySPD* ocdbEff, AliITSFONoiseSPD* ocdbNoise) :
    fSignals(), fOCDBEff(ocdbEff), fOCDBNoise(ocdbNoise)
{
  // constructor
}
//______________________________________________________________________
AliITSFOGeneratorSPD::AliITSFOGeneratorSPD(const AliITSFOGeneratorSPD& handle): 
  fSignals(handle.fSignals), fOCDBEff(handle.fOCDBEff), fOCDBNoise(handle.fOCDBNoise)
{
  // copy constructor
}
//______________________________________________________________________
AliITSFOGeneratorSPD::~AliITSFOGeneratorSPD() {
  // destructor
}
//______________________________________________________________________
AliITSFOGeneratorSPD& AliITSFOGeneratorSPD::operator=(const AliITSFOGeneratorSPD& handle) {
  // assignment operator
  if (this!=&handle) {
    fSignals = handle.fSignals;
    fOCDBEff = handle.fOCDBEff;
    fOCDBNoise = handle.fOCDBNoise;
  }
  return *this;
}
//______________________________________________________________________
void AliITSFOGeneratorSPD::SetEfficiencyAndNoise(AliITSFOEfficiencySPD* ocdbEff, AliITSFONoiseSPD* ocdbNoise) {
  // Method to give pointers to the OCDB entries, needed by methods ProcessPixelHit and ProcessNoise
  SetEfficiency(ocdbEff);
  SetNoise(ocdbNoise);
}
//______________________________________________________________________
void AliITSFOGeneratorSPD::SetEfficiency(AliITSFOEfficiencySPD* ocdbEff) {
  // Method to give pointer to the OCDB efficiency entry
  fOCDBEff = ocdbEff;
}
//______________________________________________________________________
void AliITSFOGeneratorSPD::SetNoise(AliITSFONoiseSPD* ocdbNoise) {
  // Method to give pointer to the OCDB noise entry
  fOCDBNoise = ocdbNoise;
}
//______________________________________________________________________
void AliITSFOGeneratorSPD::ProcessPixelHit(UInt_t eq, UInt_t hs, UInt_t chip, UInt_t col, UInt_t row) {
  // Here it will be decided wether a fired pixel will give rise to a fast-or signal or not
  if (eq>=20 || hs>=6 || chip>=10 || col>=32 || row>=256) {
    Error("AliITSFOGeneratorSPD::ProcessPixelHit", "eq,hs,chip,col,row (%d,%d,%d,%d,%d) out of bounds.",
	  eq,hs,chip,col,row);
    return;
  }
  if (fOCDBEff==NULL) {
    Error("AliITSFOGeneratorSPD::ProcessPixelHit", "No AliITSFOEfficiencySPD entry has been provided.");
    return;
  }
  // simulate if this fired pixel gives rise to a fast-or signal:
  if (gRandom->Rndm() < fOCDBEff->GetColumnEfficiency(eq,hs,chip,col)) {
    fSignals.SetSignal(eq,hs,chip);
  }
}
//______________________________________________________________________
void AliITSFOGeneratorSPD::ProcessPixelHitM(UInt_t module, UInt_t colM, UInt_t rowM) {
  // Converts offline coordinates to online, and calls ProcessPixelHit
  if (module>=240 || colM>=160 || rowM>=256) {
    Error("AliITSFOGeneratorSPD::ProcessPixelHitM", "module,colM,rowM (%d,%d,%d) out of bounds.",
	  module,colM,rowM);
    return;
  }
  UInt_t eq,hs,chip,col,row;
  if (AliITSRawStreamSPD::OfflineToOnline(module,colM,rowM,eq,hs,chip,col,row)) {
    ProcessPixelHit(eq,hs,chip,col,row);
  }
}
//______________________________________________________________________
void AliITSFOGeneratorSPD::ProcessNoise() {
  // 
  if (fOCDBNoise==NULL) {
    Error("AliITSFOGeneratorSPD::ProcessNoise", "No AliITSFONoiseSPD entry has been provided.");
    return;
  }
  // simulate if each pixel chip will give rise to a random noise induced fast-or signal:
  for (UInt_t eq=0; eq<20; eq++) {
    for (UInt_t hs=0; hs<6; hs++) {
      for (UInt_t chip=0; chip<10; chip++) {
	if (gRandom->Rndm() < fOCDBNoise->GetChipNoise(eq,hs,chip)) {
	  fSignals.SetSignal(eq,hs,chip);
	}
      }
    }
  }
}
