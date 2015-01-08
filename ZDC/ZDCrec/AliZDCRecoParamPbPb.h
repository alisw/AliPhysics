#ifndef ALIZDCRECOPARAMPBPB_H
#define ALIZDCRECOPARAMPBPB_H
/* Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//////////////////////////////////////////////////////////
//                                                      //
//   Class with ZDC reconstruction parameters           //
//   		Pb - Pb collisions	                //
//   Origin: Chiara.Oppedisano@to.infn.it               //
//                                                      //
//////////////////////////////////////////////////////////

#include "AliLog.h"
#include "AliCDBEntry.h"
#include "AliCDBManager.h"
#include "AliZDCRecoParam.h"

class AliZDCRecoParamPbPb : public AliZDCRecoParam {
 public:
  AliZDCRecoParamPbPb();
  AliZDCRecoParamPbPb(TH1D * hNpart, TH1D * hb, Float_t clkCent);
  AliZDCRecoParamPbPb(const AliZDCRecoParamPbPb &oldrecopar);
  AliZDCRecoParamPbPb& operator= (const AliZDCRecoParamPbPb &recpar);
  virtual ~AliZDCRecoParamPbPb();


  // make reco parameters for A-A collisions
  static AliZDCRecoParamPbPb *GetHighFluxParam(Float_t beamEnergy);
  
  TH1D* GethNpartDist()  const {return fhNpartDist;} 
  TH1D* GethbDist() 	 const {return fhbDist;}
  Float_t GetClkCenter() const {return fClkCenter;}
  
  void SetNpartDist(TH1D *hDist) {fhNpartDist = hDist;}    
  void SetbDist(TH1D *hbDist) {fhbDist = hbDist;}    
  void SetClkCenter(Float_t xValue) {fClkCenter = xValue;}    
  void SetGlauberMCDist(Float_t beamEnergy); 
    
  //void Print(Option_t *) const; 
  
 protected:
  
  // *** PARAMETERS FOR Pb-Pb COLLISIONS
  // --- Correlation E_ZDC vs. E_ZEM
  TH1D *  fhNpartDist;    // Npart distribution from Glauber MC
  TH1D *  fhbDist;	  // b distribution from Glauber MC
  Float_t fClkCenter;     // clock center: value of x-axis 
 
 ClassDef(AliZDCRecoParamPbPb, 3)

};

#endif
