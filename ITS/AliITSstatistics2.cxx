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
#include "TMath.h"
#include "AliITSstatistics2.h"

ClassImp(AliITSstatistics2)

//
AliITSstatistics2::AliITSstatistics2() : TObject(){
//
// default constructor
//
    fX  = 0;
    fY  = 0;
    fYx = 0;
    fW  = 0;
    fN  = 0;
    fOrder = 0;
    return;
}


AliITSstatistics2::AliITSstatistics2(Int_t order) : TObject(){
//
// constructor to maximum moment/order order
//
    fOrder = order;
    fX     = new Double_t[order];
    fY     = new Double_t[order];
    fYx    = new Double_t[order];
    fW     = new Double_t[order];
    for(Int_t i=0;i<order;i++) {fX[i] = 0.0;fY[i] = 0.0;
                                fYx[i] = 0.0; fW[i] = 0.0;}
    fN = 0;
    return;
}

AliITSstatistics2::~AliITSstatistics2(){
//
// destructor
//
    if(fX!=0)  delete[] fX;
    if(fY!=0)  delete[] fY;
    if(fYx!=0) delete[] fYx;
    if(fW!=0)  delete[] fW;
    fX  = 0;
    fY  = 0;
    fYx = 0;
    fW  = 0;
    fN  = 0;
    fOrder = 0;
}

//_______________________________________________________________
AliITSstatistics2& AliITSstatistics2::operator=(AliITSstatistics2 &source){
// operator =

     if(this==&source) return *this;
	  if(source.fOrder!=0){
	       this->fOrder = source.fOrder;
			 this->fN = source.fN;
			 this->fX = new Double_t[this->fOrder];
			 this->fW = new Double_t[this->fOrder];
			 for(Int_t i=0;i<source.fOrder;i++){
			      this->fX[i] = source.fX[i];
					this->fW[i] = source.fW[i];
			 } // end for i
	  }else{
	       this->fX = 0;
			 this->fW = 0;
			 this->fN = 0;
			 this->fOrder = 0;
	  }// end if source.fOrder!=0
	  return *this;
}
//_______________________________________________________________
AliITSstatistics2::AliITSstatistics2(AliITSstatistics2 &source){
// Copy constructor

     if(this==&source) return;
	  if(source.fOrder!=0){
	       this->fOrder = source.fOrder;
			 this->fN = source.fN;
			 this->fX = new Double_t[this->fOrder];
			 this->fW = new Double_t[this->fOrder];
			 for(Int_t i=0;i<source.fOrder;i++){
			      this->fX[i] = source.fX[i];
					this->fW[i] = source.fW[i];
			 } // end for i
	  }else{
	       this->fX = 0;
			 this->fW = 0;
			 this->fN = 0;
			 this->fOrder = 0;
	  }// end if source.fOrder!=0
}

void AliITSstatistics2::Reset(){
//
// Reset/zero statistics
//
    for(Int_t i=0;i<fOrder;i++) {fX[i] = 0.0;fY[i] = 0.0;
                                fYx[i] = 0.0; fW[i] = 0.0;}
    fN = 0;
    return;
}

void AliITSstatistics2::AddValue(Double_t y,Double_t x,Double_t w=1.0){
//
// add next x,y pair to statistics
//
    Double_t xs=1.0,ys=1.0,yxs=1.0,ws=1.0;
    Int_t i;

    const Double_t kBig=1.0e+38;

    if(y>kBig || x>kBig || w>kBig) return;


    fN++;
    for(i=0;i<fOrder;i++){
	xs  *= x;
	ys  *= y;
	yxs *= x*y;
	ws  *= w;
	fX[i]  += xs*w;
	fY[i]  += ys*w;
	fYx[i] += yxs*w;
	fW[i]  += ws;
    } // end for i
}

Double_t AliITSstatistics2::GetXNth(Int_t order){
//
// This give the unbiased estimator for the RMS.
//

    Double_t s;

    if(fW[0]!=0.0&&order<=fOrder) s = fX[order-1]/fW[0];
    else {
	s = 0.0;
	printf("AliITSstatistics2: error in GetNth: fOrder=%d fN=%d fW[0]=%f\n",
	       fOrder,fN,fW[0]);
    } // end else
    return s;
}
Double_t AliITSstatistics2::GetYNth(Int_t order){
//
// This give the unbiased estimator for the RMS.
//
    Double_t s;

    if(fW[0]!=0.0&&order<=fOrder) s = fY[order-1]/fW[0];
    else {
	s = 0.0;
	printf("AliITSstatistics2: error in GetNth: fOrder=%d fN=%d fW[0]=%f\n",
	       fOrder,fN,fW[0]);
    } // end else
    return s;
}
Double_t AliITSstatistics2::GetYXNth(Int_t order){
// This give the unbiased estimator for the RMS.
    Double_t s;

    if(fW[0]!=0.0&&order<=fOrder) s = fYx[order-1]/fW[0];
    else {
	s = 0.0;
	printf("AliITSstatistics2: error in GetNth: fOrder=%d fN=%d fW[0]=%f\n",
	       fOrder,fN,fW[0]);
    } // end else
    return s;
}
Double_t AliITSstatistics2::GetRMSX(){
// This give the unbiased estimator for the RMS.
    Double_t x,x2,w,ww,s;

    x  = GetMeanX(); // first order
    x2 = GetXNth(2); // second order
    w  = fW[0];     // first order - 1.
    ww = fW[1];     // second order - 1.

    if(w*w==ww) return (-1.0);
    s = (x2-x*x)*w*w/(w*w-ww);
    return TMath::Sqrt(s);
}
Double_t AliITSstatistics2::GetRMSY(){
// This give the unbiased estimator for the RMS.
    Double_t x,x2,w,ww,s;

    x  = GetMeanY(); // first order
    x2 = GetYNth(2); // second order
    w  = fW[0];     // first order - 1.
    ww = fW[1];     // second order - 1.

    if(w*w==ww) return (-1.0);
    s = (x2-x*x)*w*w/(w*w-ww);
    return TMath::Sqrt(s);
}
Double_t AliITSstatistics2::GetRMSYX(){
// This give the unbiased estimator for the RMS.
    Double_t x,x2,w,ww,s;

    x  = GetMeanYX(); // first order
    x2 = GetYXNth(2); // second order
    w  = fW[0];     // first order - 1.
    ww = fW[1];     // second order - 1.

    if(w*w==ww) return (-1.0);
    s = (x2-x*x)*w*w/(w*w-ww);
    return TMath::Sqrt(s);
}
Double_t AliITSstatistics2::GetErrorMeanY(){
//This is the error in the mean or the square root of the variance of the mean.
    Double_t rms,w,ww,s;

    rms = GetRMSY();
    w   = fW[0];
    ww  = fW[1];
    s   = rms*rms*ww/(w*w);
    return TMath::Sqrt(s);
}
Double_t AliITSstatistics2::GetErrorMeanX(){
//This is the error in the mean or the square root of the variance of the mean.
    Double_t rms,w,ww,s;

    rms = GetRMSX();
    w   = fW[0];
    ww  = fW[1];
    s   = rms*rms*ww/(w*w);
    return TMath::Sqrt(s);
}
Double_t AliITSstatistics2::GetErrorMeanYX(){
//This is the error in the mean or the square root of the variance of the mean.
    Double_t rms,w,ww,s;

    rms = GetRMSYX();
    w   = fW[0];
    ww  = fW[1];
    s   = rms*rms*ww/(w*w);
    return TMath::Sqrt(s);
}


Double_t AliITSstatistics2::GetErrorRMSY(){
//This is the error in the mean or the square root of the variance of the mean.
// at this moment this routine is only defined for weights=1.
    Double_t x,x2,x3,x4,w,ww,m2,m4,n,s;

    if(fW[0]!=(Double_t)fN||GetN()<4) return (-1.);
    x  = GetMeanY(); // first order
    x2 = GetYNth(2); // second order
    w  = fW[0];     // first order - 1.
    ww = fW[1];     // second order - 1.
    if(w*w==ww) return (-1.0);
    s = (x2-x*x)*w*w/(w*w-ww);

    m2  = s;
    n   = (Double_t) GetN();
    x3  = GetYNth(3);
    x4  = GetYNth(4);
// This equation assumes that all of the weights are equal to 1.
    m4  = (n/(n-1.))*(x4-3.*x*x3+6.*x*x*x2-2.*x*x*x*x);
    s   = (m4-(n-3.)*m2*m2/(n-1.))/n;
    return TMath::Sqrt(s);
}
Double_t AliITSstatistics2::GetErrorRMSX(){
//This is the error in the mean or the square root of the variance of the mean.
// at this moment this routine is only defined for weights=1.
    Double_t x,x2,x3,x4,w,ww,m2,m4,n,s;

    if(fW[0]!=(Double_t)fN||GetN()<4) return (-1.);
    x  = GetMeanX(); // first order
    x2 = GetXNth(2); // second order
    w  = fW[0];     // first order - 1.
    ww = fW[1];     // second order - 1.
    if(w*w==ww) return (-1.0);
    s = (x2-x*x)*w*w/(w*w-ww);

    m2  = s;
    n   = (Double_t) GetN();
    x3  = GetXNth(3);
    x4  = GetXNth(4);
// This equation assumes that all of the weights are equal to 1.
    m4  = (n/(n-1.))*(x4-3.*x*x3+6.*x*x*x2-2.*x*x*x*x);
    s   = (m4-(n-3.)*m2*m2/(n-1.))/n;
    return TMath::Sqrt(s);
}
Double_t AliITSstatistics2::GetErrorRMSYX(){
//This is the error in the mean or the square root of the variance of the mean.
// at this moment this routine is only defined for weights=1.
    Double_t x,x2,x3,x4,w,ww,m2,m4,n,s;

    if(fW[0]!=(Double_t)fN||GetN()<4) return (-1.);
    x  = GetMeanYX(); // first order
    x2 = GetYXNth(2); // second order
    w  = fW[0];     // first order - 1.
    ww = fW[1];     // second order - 1.
    if(w*w==ww) return (-1.0);
    s = (x2-x*x)*w*w/(w*w-ww);

    m2  = s;
    n   = (Double_t) GetN();
    x3  = GetYXNth(3);
    x4  = GetYXNth(4);
// This equation assumes that all of the weights are equal to 1.
    m4  = (n/(n-1.))*(x4-3.*x*x3+6.*x*x*x2-2.*x*x*x*x);
    s   = (m4-(n-3.)*m2*m2/(n-1.))/n;
    return TMath::Sqrt(s);
}
//_______________________________________________________________________
Double_t AliITSstatistics2::FitToLine(Double_t &a,Double_t &b){
// fit to y = a+bx returns Chi squared or -1.0 if an error
    Double_t c,d,e,f,g,h;

    a = b = 0.0;
    if(fOrder<2 || fN<3) return -1.0;
    c = GetWN(1);
    d = GetYN(1);
    e = GetXN(1);
    f = GetYXN(1);
    g = GetXN(2);
    h = c*g-e*e;
    a = d*g-f*e;
    b = c*f-d*e;
    if(h!=0.0){
	a = a/h;
	b = b/h;
    }else{
	printf("AliITSstatistics2: Error in FitToLine vertical line\n");
	return -1.0;
    } // end if h
    h = GetYN(2)+a*a*c+b*b*g-2.0*a*d-2.0*b*f+2.0*a*b*e;
    h /= (Double_t)fN - 2.0;
    return h;
}
/*
//_______________________________________________________________________
void AliITSstatistics2::Streamer(TBuffer &R__b){
   // Stream an object of class AliITSstatistics2.

   if (R__b.IsReading()) {
      Version_t R__v = R__b.ReadVersion(); if (R__v) { }
      TObject::Streamer(R__b);
      R__b >> fN;
      R__b >> fOrder;
      R__b.ReadArray(fY);
      R__b.ReadArray(fX);
      R__b.ReadArray(fYx);
      R__b.ReadArray(fW);
   } else {
      R__b.WriteVersion(AliITSstatistics2::IsA());
      TObject::Streamer(R__b);
      R__b << fN;
      R__b << fOrder;
      R__b.WriteArray(fY,fOrder);
      R__b.WriteArray(fX,fOrder);
      R__b.WriteArray(fYx,fOrder);
      R__b.WriteArray(fW,fOrder);
   }
}
*/
