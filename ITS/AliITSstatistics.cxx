//////////////////////////////////////////////////////////////////////////
//  Alice ITS class to help keep statistical information                //
//                                                                      //
// version: 0.0.0 Draft.                                                //
// Date: April 18 1999                                                  //
// By: Bjorn S. Nilsen                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <stdio.h>
#include <math.h>
#include <TMath.h>

#include "AliITSstatistics.h"

ClassImp(AliITSstatistics)

//
AliITSstatistics::AliITSstatistics() : TObject(){
//
// default contructor
//
    fx = 0;
    fw = 0;
    fN = 0;
    fOrder = 0;
    return;
}


AliITSstatistics::AliITSstatistics(Int_t order) : TObject(){
//
// contructor to a specific order in the moments
//
    fOrder = order;
    fx = new Double_t[order];
    fw = new Double_t[order];
    Int_t i;
    for(i=0;i<order;i++) {fx[i] = 0.0; fw[i] = 0.0;}
    fN = 0;
    return;
}

AliITSstatistics::~AliITSstatistics(){
//
// default destructor
//
    if(fx!=0) delete[] fx;
    if(fw!=0) delete[] fw;
    fx = 0;
    fw = 0;
    fN = 0;
    fOrder = 0;
}
//_______________________________________________________________
AliITSstatistics& AliITSstatistics::operator=(AliITSstatistics &source){
// operator =

     Int_t i;
     if(this==&source) return *this;
     if(source.fOrder!=0){
       this->fOrder = source.fOrder;
       this->fN = source.fN;
       this->fx = new Double_t[this->fOrder];
       this->fw = new Double_t[this->fOrder];
       for(i=0;i<source.fOrder;i++){
	 this->fx[i] = source.fx[i];
	 this->fw[i] = source.fw[i];
       } // end for i
     }else{
       this->fx = 0;
       this->fw = 0;
       this->fN = 0;
       this->fOrder = 0;
     }// end if source.fOrder!=0
     return *this;
}
//_______________________________________________________________
AliITSstatistics::AliITSstatistics(AliITSstatistics &source){
// Copy constructor

  Int_t i;
  if(this==&source) return;
  if(source.fOrder!=0){
    this->fOrder = source.fOrder;
    this->fN = source.fN;
    this->fx = new Double_t[this->fOrder];
    this->fw = new Double_t[this->fOrder];
    for(i=0;i<source.fOrder;i++){
      this->fx[i] = source.fx[i];
      this->fw[i] = source.fw[i];
    } // end for i
  }else{
    this->fx = 0;
    this->fw = 0;
    this->fN = 0;
    this->fOrder = 0;
  }// end if source.fOrder!=0
}
//_______________________________________________________________
void AliITSstatistics::Reset(){
//
// reset all values to zero
//
    Int_t i;
    for(i=0;i<fOrder;i++) {fx[i] = 0.0; fw[i] = 0.0;}
    fN = 0;
    return;
}

void AliITSstatistics::AddValue(Double_t x,Double_t w=1.0){
//
// accumulate element x with weight w.
//
    Double_t y=1.0,z=1.0;
    Int_t i;


    if(isinf(y)!=0||isinf(x)!=0||isinf(w)!=0) return;
    if(isnan(y)!=0||isnan(x)!=0||isnan(w)!=0) return;
    fN++;
    for(i=0;i<fOrder;i++){
	y *= x;
	z *= w;
	fx[i] += y*w;
	fw[i] += z;
    } // end for i
}

Double_t AliITSstatistics::GetNth(Int_t order){
// This give the unbiased estimator for the RMS.
    Double_t s;

    if(fw[0]!=0.0&&order<=fOrder) s = fx[order-1]/fw[0];
    else {
	s = 0.0;
	printf("AliITSstatistics: error in GetNth: fOrder=%d fN=%d fw[0]=%f\n",
	       fOrder,fN,fw[0]);
    } // end else
    return s;
}

Double_t AliITSstatistics::GetRMS(){
// This give the unbiased estimator for the RMS.
    Double_t x,x2,w,ww,s;

    x  = GetMean(); // first order
    x2 = GetNth(2); // second order
    w  = fw[0];     // first order - 1.
    ww = fw[1];     // second order - 1.

    if(w*w==ww) return (-1.0);
    s = (x2-x*x)*w*w/(w*w-ww);
    return TMath::Sqrt(s);
}

Double_t AliITSstatistics::GetErrorMean(){
//This is the error in the mean or the square root of the variance of the mean.
    Double_t rms,w,ww,s;

    rms = GetRMS();
    w   = fw[0];
    ww  = fw[1];
    s   = rms*rms*ww/(w*w);
    return TMath::Sqrt(s);
}


Double_t AliITSstatistics::GetErrorRMS(){
//This is the error in the mean or the square root of the variance of the mean.
// at this moment this routine is only defined for weights=1.
    Double_t x,x2,x3,x4,w,ww,m2,m4,n,s;

    if(fw[0]!=(Double_t)fN||GetN()<4) return (-1.);
    x  = GetMean(); // first order
    x2 = GetNth(2); // second order
    w  = fw[0];     // first order - 1.
    ww = fw[1];     // second order - 1.
    if(w*w==ww) return (-1.0);
    s = (x2-x*x)*w*w/(w*w-ww);

    m2  = s;
    n   = (Double_t) GetN();
    x3  = GetNth(3);
    x4  = GetNth(4);
// This equation assumes that all of the weights are equal to 1.
    m4  = (n/(n-1.))*(x4-3.*x*x3+6.*x*x*x2-2.*x*x*x*x);
    s   = (m4-(n-3.)*m2*m2/(n-1.))/n;
    return TMath::Sqrt(s);
}
//_______________________________________________________________________
void AliITSstatistics::Streamer(TBuffer &R__b){
   // Stream an object of class AliITSstatistics.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> fN;
      R__b >> fOrder;
      R__b.ReadArray(fx);
      R__b.ReadArray(fw);
   } else {
      R__b.WriteVersion(AliITSstatistics::IsA());
      TObject::Streamer(R__b);
      R__b << fN;
      R__b << fOrder;
      R__b.WriteArray(fx,fOrder);
      R__b.WriteArray(fw,fOrder);
   }
}
