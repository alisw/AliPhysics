#ifndef ALIZDCRECOPARAM_H
#define ALIZDCRECOPARAM_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////
//                                                      //
//   Class with ZDC reconstruction parameters           //
//   Origin: Chiara.Oppedisano@to.infn.it               //
//                                                      //
//////////////////////////////////////////////////////////

#include <TH1D.h>
#include <TF1.h>
#include "AliDetectorRecoParam.h"

//class TF1;

class AliZDCRecoParam : public AliDetectorRecoParam {
 public:
  AliZDCRecoParam();
  virtual ~AliZDCRecoParam();
  
  virtual Float_t GetBeamEnergy() {return fBeamEnergy;}
  virtual TH1D*   GethNpartDist() const {return fhNpartDist;}	      
  virtual TH1D*   GethbDist()     const {return fhbDist;}
  virtual Float_t GetClkCenter()  const {return fClkCenter;}
      
  virtual void PrintParameters() const {;} 
  
  virtual void SetGlauberMCDist(Float_t beamEnergy);
  virtual void SetBeamEnergy(Float_t beamEnergy) {fBeamEnergy = beamEnergy;}
  
  virtual void SetNpartDist(TH1D *hDist) {fhNpartDist = hDist;}    
  virtual void SetbDist(TH1D *hbDist) {fhbDist = hbDist;}    
  virtual void SetClkCenter(Float_t xValue) {fClkCenter = xValue;}    
  
 protected:
  
  AliZDCRecoParam(const AliZDCRecoParam&);
  AliZDCRecoParam& operator =(const AliZDCRecoParam&);
  
  Float_t fBeamEnergy;    // beam energy
  
  // *** PARAMETERS FOR Pb-Pb 
  TH1D *  fhNpartDist;    // Npart distribution from Glauber MC
  TH1D *  fhbDist;	  // b distribution from Glauber MC
  Float_t fClkCenter;     // clock center: value of x-axis 
   
 ClassDef(AliZDCRecoParam, 4)

};

#endif
