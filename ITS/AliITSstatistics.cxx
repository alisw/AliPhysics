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
AliITSstatistics::AliITSstatistics() : TObject(),
fN(0),
fOrder(0),
fX(0),
fW(0){
//
// default contructor
//
}


AliITSstatistics::AliITSstatistics(Int_t order) : TObject(),
fN(0),
fOrder(order),
fX(0),
fW(0){
//
// contructor to a specific order in the moments
//
    fX = new Double_t[order];
    fW = new Double_t[order];
    for(Int_t i=0;i<order;i++) {fX[i] = 0.0; fW[i] = 0.0;}
    return;
}

AliITSstatistics::~AliITSstatistics(){
//
// default destructor
//
    if(fX!=0) delete[] fX;
    if(fW!=0) delete[] fW;
    fX = 0;
    fW = 0;
    fN = 0;
    fOrder = 0;
}
//_______________________________________________________________
AliITSstatistics& AliITSstatistics::operator=(const AliITSstatistics &source){
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
AliITSstatistics::AliITSstatistics(const AliITSstatistics &source) : TObject(source),
fN(0),
fOrder(0),
fX(0),
fW(0){
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
//_______________________________________________________________
void AliITSstatistics::Reset(){
//
// reset all values to zero
//
    for(Int_t i=0;i<fOrder;i++) {fX[i] = 0.0; fW[i] = 0.0;}
    fN = 0;
    return;
}

//_______________________________________________________________
void AliITSstatistics::AddValue(Double_t x,Double_t w){
//
// accumulate element x with weight w.
//

  //it was AddValue(Double_t x,Double_t w=1.0);

    Double_t y=1.0,z=1.0;

    Int_t i;
    const Double_t kBig=1.0e+38;

    if(y>kBig || x>kBig || w>kBig) return;

    fN++;
    for(i=0;i<fOrder;i++){
	y *= x;
	z *= w;
	fX[i] += y*w;
	fW[i] += z;
    } // end for i
}

Double_t AliITSstatistics::GetNth(Int_t order){
// This give the unbiased estimator for the RMS.
    Double_t s;

    if(fW[0]!=0.0&&order<=fOrder) s = fX[order-1]/fW[0];
    else {
	s = 0.0;
	printf("AliITSstatistics: error in GetNth: fOrder=%d fN=%d fW[0]=%f\n",
	       fOrder,fN,fW[0]);
    } // end else
    return s;
}

Double_t AliITSstatistics::GetRMS(){
// This give the unbiased estimator for the RMS.
    Double_t x,x2,w,ww,s;

    x  = GetMean(); // first order
    x2 = GetNth(2); // second order
    w  = fW[0];     // first order - 1.
    ww = fW[1];     // second order - 1.

    if(w*w==ww) return (-1.0);
    s = (x2-x*x)*w*w/(w*w-ww);
    return TMath::Sqrt(s);
}

Double_t AliITSstatistics::GetErrorMean(){
//This is the error in the mean or the square root of the variance of the mean.
    Double_t rms,w,ww,s;

    rms = GetRMS();
    w   = fW[0];
    ww  = fW[1];
    s   = rms*rms*ww/(w*w);
    return TMath::Sqrt(s);
}


Double_t AliITSstatistics::GetErrorRMS(){
//This is the error in the mean or the square root of the variance of the mean.
// at this moment this routine is only defined for weights=1.
    Double_t x,x2,x3,x4,w,ww,m2,m4,n,s;

    if(fW[0]!=(Double_t)fN||GetN()<4) return (-1.);
    x  = GetMean(); // first order
    x2 = GetNth(2); // second order
    w  = fW[0];     // first order - 1.
    ww = fW[1];     // second order - 1.
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
