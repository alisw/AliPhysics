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
  AliZDCRecoParamPbPb(TH2F * hZDCvsZEM, TH2F * hZDCCvsZEM, TH2F * hZDCAvsZEM);
  AliZDCRecoParamPbPb(TH2F * hZDCvsZEM, TH2F * hZDCCvsZEM, TH2F * hZDCAvsZEM,
  		      TH1D * hNpart, TH1D * hb, Float_t clkCent);
  AliZDCRecoParamPbPb(const AliZDCRecoParamPbPb &oldrecopar);
  AliZDCRecoParamPbPb& operator= (const AliZDCRecoParamPbPb &recpar);
  virtual ~AliZDCRecoParamPbPb();


  // make reco parameters for A-A collisions
  static AliZDCRecoParamPbPb *GetHighFluxParam();
  
  TH2F* GethZDCvsZEM()   const {return fhZDCvsZEM;}  
  TH2F* GethZDCCvsZEM()  const {return fhZDCCvsZEM;}  
  TH2F* GethZDCAvsZEM()  const {return fhZDCAvsZEM;}  
  TH1D* GethNpartDist()  const {return fhNpartDist;} 
  TH1D* GethbDist() 	 const {return fhbDist;}
  Float_t GetClkCenter() const {return fClkCenter;}
 
  void SetZDCvsZEM(TH2F *hCorr)  {fhZDCvsZEM  = hCorr;}    
  void SetZDCCvsZEM(TH2F *hCorr) {fhZDCCvsZEM = hCorr;}    
  void SetZDCAvsZEM(TH2F *hCorr) {fhZDCAvsZEM = hCorr;}   
  void SetNpartDist(TH1D *hDist) {fhNpartDist = hDist;}    
  void SetbDist(TH1D *hbDist) {fhbDist = hbDist;}    
  void SetClkCenter(Float_t xValue) {fClkCenter = xValue;}    
  void SetGlauberMCDist(); 
    
  //void Print(Option_t *) const; 
  
 protected:
  
  // *** PARAMETERS FOR Pb-Pb COLLISIONS
  // --- Correlation E_ZDC vs. E_ZEM
  TH2F *  fhZDCvsZEM;	// E_ZDC (total) vs. E_ZEM 
  TH2F *  fhZDCCvsZEM;  // E_ZDC vs. E_ZEM sideC
  TH2F *  fhZDCAvsZEM;  // E_ZDC vs. E_ZEM sideA
  TH1D *  fhNpartDist;  // Npart distribution from Glauber MC
  TH1D *  fhbDist;	// b distribution from Glauber MC
  Float_t fClkCenter;   // clock center: value of x-axis 
 
 ClassDef(AliZDCRecoParamPbPb, 2)

};

#endif
