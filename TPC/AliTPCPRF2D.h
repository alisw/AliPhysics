#ifndef ALITPCPRF2D_H
#define ALITPCPRF2D_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////
//  Manager class for AliTPCPRF2D                  //
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
  ~AliTPCPRF2D();

  Float_t GetGRF(Float_t xin, Float_t yin); 
  //return generic response function  in xin
  Float_t GetPRF(Float_t xin, Float_t yin, Bool_t inter=kFALSE); 
  //return PRF in point xin,yin
  void SetY(Float_t y1, Float_t y2, Int_t nYdiv) ;
  void DrawX(Float_t x1, Float_t x2,Float_t y, Bool_t inter=kFALSE);  
  //draw one dimensional response for
  //fixed y
  // void DrawY(Float_t y1, Float_t y2,Float_t x);
  //draw one dimensional response for fixed x
  void Draw(Option_t *) {}
  void Draw(Float_t x1, Float_t x2, Float_t y1, Float_t y2,
	    Bool_t inter=kFALSE, Int_t Nx=20, Int_t Ny=20);
  //draw two dimensional PRF

  void DrawDist(Float_t x1, Float_t x2, Float_t y1, Float_t y2,
		Bool_t inter=kFALSE, Int_t Nx=20, Int_t Ny=20, 
		Float_t  thr=0);
  //draw distortion of COG method
  //we suppose threshold equal to thr
   
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
  void SetPad(Float_t width, Float_t height);
  //set base chevron parameters
  void SetChevron(Float_t hstep, Float_t shifty, Float_t fac);
  //set chevron parameters   
  void SetChParam(Float_t width, Float_t height,
		  Float_t hstep, Float_t shifty, Float_t fac);
  //set all geometrical parameters  
  void SetNdiv(Int_t Ndiv){fNdiv=Ndiv;}
  void Update();  
protected:
  void Update1();  
  Float_t GetPRFActiv(Float_t xin); //return PRF in point xin and actual y
  Float_t  * fcharge; // field with PRF 
  Float_t fY1;
  Float_t fY2;
  Int_t fNYdiv;  
  Float_t * ffcharge;  //pointer to array of arrays

private: 
  Double_t funParam[5];//parameters of used charge function
  Int_t  fNPRF;      //number of interpolations point
  Int_t  fNdiv;      //number of division to calculate integral
  Float_t fDStep;    //element step for point 
  Float_t fkNorm;     //normalisation factor of the charge integral
  Float_t fInteg;     //integral of GRF on +- infinity
  TF2 *  fGRF;        //charge distribution function

  Float_t fK3X;
  Float_t fK3Y;
  Float_t fPadDistance;

  Float_t  forigsigmaX;//sigma of original distribution;  
  Float_t  forigsigmaY;//sigma of original distribution;  
  Float_t fSigmaX;     //sigma of PAD response function
  //calculated during update

  //chewron parameters
  Float_t fHeightFull;  //height of the full pad
  Float_t fHeightS;     //height of the one step
  Float_t fShiftY;  //shift of the step
  Float_t fWidth;       //width of the pad
  Float_t fK;           //k factor of the chewron
  Float_t fActualY;          //in reality we calculate PRF only for 
  //one fixed y
  char  fType[5];  //charge type
  //to make calculation faster we reduce  division
  Float_t fDYtoWire;    //used to make PRF calculation faster in GetPRF
  Float_t fDStepM1;     //used in GetPRFActiv
  ClassDef(AliTPCPRF2D,1) 
}; 

#endif /* ALITPCPRF2D_H */
  
