#ifndef ALIITSUSIMULATIONPIX_H
#define ALIITSUSIMULATIONPIX_H

/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
* See cxx source for full Copyright notice                               */

////////////////////////////////////////////////////////////
// Simulation class for upgrade pixels                    //
////////////////////////////////////////////////////////////

#include <TObjArray.h>
#include "AliITSUSimulation.h"
#include "AliITSUSegmentationPix.h"

class TH1F;
class AliITSUModule;
class AliITSUSimuParam;

//-------------------------------------------------------------------

class AliITSUSimulationPix : public AliITSUSimulation {
public:
  AliITSUSimulationPix();
  AliITSUSimulationPix(AliITSUSimuParam* sim,AliITSUSensMap* map);
  virtual ~AliITSUSimulationPix();
  AliITSUSimulationPix(const AliITSUSimulationPix &source); 
  AliITSUSimulationPix& operator=(const AliITSUSimulationPix &s);
  void Init();
  //
  void FinishSDigitiseModule();
  void DigitiseModule(AliITSUModule *mod,Int_t mask, Int_t event, AliITSsegmentation* seg);
  //
  void SDigitiseModule(AliITSUModule *mod, Int_t mask, Int_t event, AliITSsegmentation* seg);
  void WriteSDigits();
  void Hits2SDigits(AliITSUModule *mod);
  void Hits2SDigitsFast(AliITSUModule *mod);
  void AddNoisyPixels();   
  void RemoveDeadPixels();
  void FrompListToDigits();
  Int_t CreateNoisyDigits(Int_t minID,Int_t maxID,double probNoisy, double noise, double base);
  Bool_t SetTanLorAngle(Double_t WeightHole=1.0);
  Double_t GetTanLorAngle() const {return fTanLorAng;};
  //
  // For backwards compatibility
  void SDigitsToDigits(){ FinishSDigitiseModule();}
  
  // This sets fStrobe flag and allows generating the strobe and applying it to select hits 
  void SetStrobeGeneration(Bool_t b=kFALSE) {fStrobe=b;};
  void GenerateStrobePhase();
  //
  
 private:
  void SpreadCharge(Double_t x0,Double_t z0,Int_t ix0,Int_t iz0,Double_t el,Double_t sig,Double_t ld,Int_t t,Int_t hi);
  void SpreadChargeAsym(Double_t x0,Double_t z0,Int_t ix0,Int_t iz0,Double_t el,Double_t sigx,Double_t sigz,Double_t ld,Int_t t,Int_t hi);
  //
  void SetCoupling(AliITSUSDigit* old,Int_t ntrack,Int_t idhit);     // "New" coupling routine  Tiziano Virgili
  void SetCouplingOld(AliITSUSDigit* old,Int_t ntrack,Int_t idhit);  // "Old" coupling routine  Rocco Caliandro
  //   
 protected:
   Double_t      fTanLorAng;    //! Tangent of the Lorentz Angle (weighted average for hole and electrons)
   Bool_t        fStrobe;       // kTRUE if readout strobe with proper phase applied to select hits
   Int_t         fStrobeLenght; // Strobe signal lenght in units of 25 ns
   Double_t      fStrobePhase;  // The phase of the strobe signal with respect to the trigger
   ClassDef(AliITSUSimulationPix,1)  // Simulation of pixel clusters
 };
#endif 
