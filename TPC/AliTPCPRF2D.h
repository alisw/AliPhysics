#ifndef ALITPCPRF2D_H
#define ALITPCPRF2D_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
////////////////////////////////////////////////
//  Manager class for AliTPCPRF2D             //
////////////////////////////////////////////////
  

// include files and class forward declarations
//DSTEP in cm
//NPRF in number of interpolation points


#include "TObject.h"
#include "TMath.h"
class TF2;
class TArrayF;

class AliTPCPRF2D : public TObject {
public : 
  AliTPCPRF2D();
   AliTPCPRF2D(const AliTPCPRF2D &prf);
  AliTPCPRF2D & operator = (const AliTPCPRF2D &prf);
  ~AliTPCPRF2D();
  void Update();  //recalculate tables for charge calculation
  Float_t GetGRF(Float_t xin, Float_t yin); 
  //return generic response function  in xin
  Float_t GetPRF(Float_t xin, Float_t yin, Bool_t inter=kFALSE); 
  //return PRF in point xin,yin
  void DrawX(Float_t x1, Float_t x2,Float_t y, Bool_t inter=kFALSE);  
  //draw one dimensional response for
  //fixed y
  // void DrawY(Float_t y1, Float_t y2,Float_t x);
  //draw one dimensional response for fixed x
  void DrawPRF(Float_t x1, Float_t x2, Float_t y1, Float_t y2,
	    Bool_t inter=kFALSE, Int_t Nx=20, Int_t Ny=20);
  //draw two dimensional PRF

  void DrawDist(Float_t x1, Float_t x2, Float_t y1, Float_t y2,
		Bool_t inter=kFALSE, Int_t Nx=20, Int_t Ny=20, 
		Float_t  thr=0);
  //draw distortion of COG method
  //we suppose threshold equal to thr

  void SetPad(Float_t width, Float_t height);
  //set base chevron parameters
  void SetChevron(Float_t hstep, Float_t shifty, Float_t fac);
  //set chevron parameters   
  void SetChParam(Float_t width, Float_t height,
		  Float_t hstep, Float_t shifty, Float_t fac);
  //set all geometrical parameters     
  void SetY(Float_t y1, Float_t y2, Int_t nYdiv) ;
  void SetChargeAngle(Float_t angle){fChargeAngle = angle;} //set angle of pad and charge distribution
                                                            //axes
  void SetGauss(Float_t sigmaX,Float_t sigmaY , Float_t kNorm=1);
  //adjust PRF with GAUSIAN as generic GRF 
  //if  direct = kTRUE then it does't convolute distribution
  void SetCosh(Float_t sigmaX,Float_t sigmaY , Float_t kNorm=1);
  //adjust PRF with 1/Cosh  as generic GRF
  void  SetGati(Float_t K3X, Float_t K3Y,
		     Float_t padDistance,
		     Float_t kNorm=1);
  void SetParam(TF2 * GRF,Float_t kNorm, 
		Float_t sigmaX=0, Float_t sigmaY=0);
  void SetNdiv(Int_t Ndiv){fNdiv=Ndiv;}
  Float_t GetSigmaX() const {return fSigmaX;}
  Float_t GetSigmaY() const {return fSigmaY;}
  
  
protected:
  void Update1(); 
  void UpdateSigma();  //recalculate sigma of PRF
  Float_t GetPRFActiv(Float_t xin); //return PRF in point xin and actual y
  Float_t  * fcharge; //field with PRF 
  Float_t fY1;        //position of first "virtual" vire 
  Float_t fY2;        //position of last virtual vire
  Int_t fNYdiv;       //number of wires
  Float_t * ffcharge;  //pointer to array of arrays

private: 

  //chevron parameters
  Float_t fHeightFull;  //height of the full pad
  Float_t fHeightS;     //height of the one step
  Float_t fShiftY;      //shift of the step
  Float_t fWidth;       //width of the pad
  Float_t fK;           //k factor of the chewron

  Double_t funParam[5];//parameters of used charge function
  Int_t  fNPRF;      //number of interpolations point
  Int_t  fNdiv;      //number of division to calculate integral
  Float_t fDStep;    //element step for point 
  Float_t fkNorm;     //normalisation factor of the charge integral
  Float_t fInteg;     //integral of GRF on +- infinity
  TF2 *  fGRF;        //charge distribution function

  Float_t fK3X;       //KX parameter (only for Gati parametrization)
  Float_t fK3Y;       //KY parameter (only for Gati parametrisation)
  Float_t fPadDistance; //pad anode distnce (only for Gati parametrisation)

  Float_t  fOrigSigmaX; //sigma of original distribution;  
  Float_t  fOrigSigmaY; //sigma of original distribution;  
  Float_t  fChargeAngle;//'angle' of charge distribution refernce system to pad reference system
  Float_t  fCosAngle;   //'angle' of the pad assymetry

  Float_t  fSigmaX;    //sigma X of PAD response function
  Float_t  fSigmaY;    //sigma Y of PAD response function
  Float_t  fMeanX;     //mean X value
  Float_t  fMeanY;     //mean Y value
  //calculated during update

   

  char  fType[5];       //charge type
  Float_t fCurrentY;    //in reality we calculate PRF only for one fixed y 
  Float_t fDYtoWire;    //! used to make PRF calculation faster in GetPRF
  Float_t fDStepM1;     //! used in GetPRFActiv to make calculation faster
  
  static const Float_t fgSQRT12; //numeric constant
  static const Int_t   fgNPRF;   //default number of division

  ClassDef(AliTPCPRF2D,1) 
}; 

#endif /* ALITPCPRF2D_H */
  
