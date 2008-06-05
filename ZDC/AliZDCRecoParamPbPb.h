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

#include <TF1.h>
#include "AliZDCRecoParam.h"
#include "AliZDCRecoParamPbPb.h"

//class TF1;

class AliZDCRecoParamPbPb : public AliZDCRecoParam {
 public:
  AliZDCRecoParamPbPb();
  virtual ~AliZDCRecoParamPbPb();

  // make reco parameters for Pb-Pb collisions
  static AliZDCRecoParamPbPb *GetPbPbRecoParam();
  
  TF1* GetfZNCen()  const {return fZNCen;}
  TF1* GetfZNPer()  const {return fZNPer;} 
  TF1* GetfZPCen()  const {return fZPCen;} 
  TF1* GetfZPPer()  const {return fZPPer;} 
  TF1* GetfZDCCen() const {return fZDCCen;}
  TF1* GetfZDCPer() const {return fZDCPer;}
  TF1* GetfbCen()   const {return fbCen;}   
  TF1* GetfbPer()   const {return fbPer;}   
  TF1* GetfZEMn()   const {return fZEMn;}   
  TF1* GetfZEMp()   const {return fZEMp;}   
  TF1* GetfZEMsp()  const {return fZEMsp;}
  TF1* GetfZEMb()   const {return fZEMb;}   
  
  void SetfZNCen(const char* formula, Double_t xmin, Double_t xmax)
       {fZNCen = new TF1("fZNCen", formula, xmin, xmax);}
  void SetfZNPer(const char* formula, Double_t xmin, Double_t xmax)
       {fZNPer = new TF1("fZNPer", formula, xmin, xmax);}
  void SetfZPCen(const char* formula, Double_t xmin, Double_t xmax)
       {fZPCen = new TF1("fZPCen", formula, xmin, xmax);}
  void SetfZPPer(const char* formula, Double_t xmin, Double_t xmax)
       {fZPPer = new TF1("fZPPer", formula, xmin, xmax);}
  void SetfZDCCen(const char* formula, Double_t xmin, Double_t xmax)
       {fZDCCen = new TF1("fZDCCen", formula, xmin, xmax);}
  void SetfZDCPer(const char* formula, Double_t xmin, Double_t xmax)
       {fZDCPer = new TF1("fZDCPer", formula, xmin, xmax);}
  void SetfbCen(const char* formula, Double_t xmin, Double_t xmax)
       {fbCen = new TF1("fbCen", formula, xmin, xmax);}
  void SetfbPer(const char* formula, Double_t xmin, Double_t xmax)
       {fbPer = new TF1("fbPer", formula, xmin, xmax);}
  void SetfZEMn(const char* formula, Double_t xmin, Double_t xmax)
       {fZEMn = new TF1("fZEMn", formula, xmin, xmax);}
  void SetfZEMp(const char* formula, Double_t xmin, Double_t xmax)
       {fZEMp = new TF1("fZEMp", formula, xmin, xmax);}
  void SetfZEMsp(const char* formula, Double_t xmin, Double_t xmax)
       {fZEMsp = new TF1("fZEMsp", formula, xmin, xmax);}
  void SetfZEMb(const char* formula, Double_t xmin, Double_t xmax)
       {fZEMb = new TF1("fZEMb", formula, xmin, xmax);}
  
  Float_t GetZEMEndValue()     const {return fZEMEndValue;}
  Float_t GetZEMCutFraction()  const {return fZEMCutFraction;}
  Float_t GetDZEMSup()	       const {return fDZEMSup;}
  Float_t GetDZEMInf()	       const {return fDZEMInf;}
  //
  Float_t GetEZN1MaxValue()  const {return fEZN1MaxValue;}
  Float_t GetEZP1MaxValue()  const {return fEZP1MaxValue;}
  Float_t GetEZDC1MaxValue() const {return fEZDC1MaxValue;}
  Float_t GetEZN2MaxValue()  const {return fEZN2MaxValue;}
  Float_t GetEZP2MaxValue()  const {return fEZP2MaxValue;}
  Float_t GetEZDC2MaxValue() const {return fEZDC2MaxValue;}
    
  void  SetZEMEndValue(Float_t ZEMEndValue) {fZEMEndValue = ZEMEndValue;}
  void  SetZEMCutFraction(Float_t ZEMCutFraction) {fZEMCutFraction = ZEMCutFraction;}
  void  SetDZEMSup(Float_t DZEMSup) {fDZEMSup = DZEMSup;}
  void  SetDZEMInf(Float_t DZEMInf) {fDZEMInf = DZEMInf;}
  //
  void	SetEZN1MaxValue(Float_t value)  {fEZN1MaxValue = value;}
  void	SetEZP1MaxValue(Float_t value)  {fEZP1MaxValue = value;}
  void	SetEZDC1MaxValue(Float_t value) {fEZDC1MaxValue = value;}
  void	SetEZN2MaxValue(Float_t value)  {fEZN2MaxValue = value;}
  void	SetEZP2MaxValue(Float_t value)  {fEZP2MaxValue = value;}
  void	SetEZDC2MaxValue(Float_t value) {fEZDC2MaxValue = value;}
  
  void PrintParameters() const; 
  
 protected:
  
  AliZDCRecoParamPbPb(const AliZDCRecoParamPbPb&);
  AliZDCRecoParamPbPb& operator =(const AliZDCRecoParamPbPb&);
 
  // *** PARAMETERS FOR Pb-Pb COLLISIONS
  // --- Functions to evaluate centrality variables from defined functions
  TF1* fZNCen;   //! Nspectator n true vs. EZN
  TF1* fZNPer;   //! Nspectator n true vs. EZN
  TF1* fZPCen;   //! Nspectator p true vs. EZP
  TF1* fZPPer;   //! Nspectator p true vs. EZP
  TF1* fZDCCen;  //! Nspectators true vs. EZDC
  TF1* fZDCPer;  //! Nspectators true vs. EZDC
  TF1* fbCen;	 //! b vs. EZDC
  TF1* fbPer;	 //! b vs. EZDC
  TF1* fZEMn;	 //! Nspectators n from ZEM energy
  TF1* fZEMp;	 //! Nspectators p from ZEM energy
  TF1* fZEMsp;   //! Nspectators from ZEM energy
  TF1* fZEMb;	 //! b from ZEM energy
  // --- Coefficients for centrality selection from ZEM signal
  Float_t  fZEMEndValue;    	 // End point value of ZEM energy spectrum
  Float_t  fZEMCutFraction; 	 // Fraction of ZEM energy spectrum used to cut
  Float_t  fDZEMSup;// Upper value of EZDCvs.ZEM correlation where ZEM signal is used
  Float_t  fDZEMInf;// Lower value of EZDCvs.ZEM correlation where ZEM signal is used
  // --- Parameters from EZDC vs. Nspec correlation
  Float_t  fEZN1MaxValue;	 // Max value of ZN1 vs. Nspec n correlation
  Float_t  fEZP1MaxValue;	 // Max value of ZP1 vs. Nspec p correlation
  Float_t  fEZDC1MaxValue;	 // Max value of ZDC1 vs. Nspec n+p correlation
  Float_t  fEZN2MaxValue;	 // Max value of ZN2 vs. Nspec n correlation
  Float_t  fEZP2MaxValue;	 // Max value of ZP2 vs. Nspec p correlation
  Float_t  fEZDC2MaxValue;	 // Max value of ZDC2 vs. Nspec n+p correlation
 
 ClassDef(AliZDCRecoParamPbPb, 1)

};

#endif
