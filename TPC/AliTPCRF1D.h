#ifndef AliTPCRF1D_H
#define AliTPCRF1D_H
////////////////////////////////////////////////
//  Manager class for AliTPCRF1D                  //
////////////////////////////////////////////////
  

// include files and class forward declarations
//DSTEP in cm
//NPRF in number of interpolation points
const Int_t   NRF=100;
const Float_t RFDSTEP=0.01;

#include "TObject.h"
#include "TMath.h"
class TF1;


class AliTPCRF1D : public TObject {
public : 
  AliTPCRF1D(Bool_t direct=kFALSE,Int_t np=NRF,Float_t step=RFDSTEP );
  ~AliTPCRF1D();  
  Float_t GetRF(Float_t xin); //return RF in point xin
  Float_t GetGRF(Float_t xin); //return generic response function  in xin
  void SetGauss(Float_t sigma,Float_t padWidth, Float_t kNorm);
  //adjust RF with GAUSIAN as generic GRF 
  //if  direct = kTRUE then it does't convolute distribution
  void SetCosh(Float_t sigma,Float_t padWidth, Float_t kNorm);
  void SetGati(Float_t K3, Float_t padDistance, Float_t padWidth,
	       Float_t kNorm);
  //adjust RF with 1/Cosh  as generic GRF
  void SetParam(TF1 * GRF,Float_t padwidth,Float_t kNorm, 
		Float_t sigma=0);
  //adjust RF with general function 
  void SetOffset(Float_t xoff) {fOffset=xoff;}
  //set offset value 
  Float_t GetPadWidth(){ return fpadWidth;};       
  //return  pad width 
  Float_t  GetSigma(){return fSigma;}
  //return estimated sigma of RF
  void Draw(Option_t*) {}
  void Draw(Float_t x1=-3 ,Float_t x2 =3.,Int_t N = 200);
  //draw RF it don't delete histograms after drawing
  /// it's on user !!!!
  void Update();  
private: 
  Double_t funParam[5];//parameters of used charge function
  Int_t  fNRF;      //number of interpolations point
  Float_t fDSTEPM1;    //element step for point
  Float_t* fcharge; // field with RF
  Float_t  forigsigma;//sigma of original distribution;
  Float_t fpadWidth;  //width of pad
  Float_t fkNorm;     //normalisation factor of the charge integral
  Float_t fInteg;     //integral of GRF on +- infinity
  TF1 *  fGRF;        //charge distribution function
  Float_t fSigma;     //sigma of PAD response function

  Float_t fOffset;    //offset of response function (for time reponse we 
  //have for expample shifted gauss)
  //calculated during update
 
  Bool_t fDirect;     //tell us if we use directly generalfunction
  Float_t fK3X;  
  Float_t fPadDistance; 
  //charge type
  char  fType[5];
  ClassDef(AliTPCRF1D,2)
}; 




#endif /* AliTPCRF1D_H */
  
