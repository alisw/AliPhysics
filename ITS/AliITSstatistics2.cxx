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
#include "Riostream.h"

ClassImp(AliITSstatistics2)

//
AliITSstatistics2::AliITSstatistics2() : 
TObject(), // Base Class
fN(-1),    // number of enetries -1 => Uninitilized
fOrder(0), // maximum moment of distributions (^n)
fX(0),    //[fOrder] array of sums of x^n
fYx(0),   //[fOrder] array of sums of (xy)^n
fY(0),    //[fOrder] array of sums of y^n
fW(0){    //[fOrder] array of sums of w^n (weights)
    // default constructor
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   A default constructed AliITSstatistics class

    return;
}
//______________________________________________________________________
AliITSstatistics2::AliITSstatistics2(Int_t order) : 
TObject(),     // Base Class
fN(0),         // number of enetries -1 => Uninitilized
fOrder(order), // maximum moment of distributions (^n)
fX(new Double_t[order]),    //[fOrder] array of sums of x^n
fYx(new Double_t[order]),   //[fOrder] array of sums of (xy)^n
fY(new Double_t[order]),    //[fOrder] array of sums of y^n
fW(new Double_t[order]){    //[fOrder] array of sums of w^n (weights)
    // constructor to maximum moment/order order
    // Inputs:
    //   Int_t order   The maximum moment of distributions {for example x^n}
    Int_t i;

    for(i=0;i<order;i++) {
        fX[i]  = 0.0;
        fY[i]  = 0.0;
        fYx[i] = 0.0;
        fW[i]  = 0.0;
    } // end for i
    fN = 0;
    return;
}
//______________________________________________________________________
AliITSstatistics2::~AliITSstatistics2(){
    // destructor
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.

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
    // Inputs:
    //   AliITSstaticstics2 &source The source of this copy.
    // Outputs:
    //   none.
    // Return:
    //   A copy of the source class

    if(this==&source) return *this;
    TObject::operator=(source);
    Reset(source.GetOrder());
    fN = source.GetN();
    fOrder=source.GetOrder();
    for(Int_t i=0;i<source.fOrder;i++){
        this->fX[i] = source.fX[i];
        this->fYx[i] = source.fYx[i];
        this->fY[i] = source.fY[i];
        this->fW[i] = source.fW[i];
    } // end for i
    return *this;
}
//_______________________________________________________________
AliITSstatistics2::AliITSstatistics2(AliITSstatistics2 &source): 
TObject(source),          // Base Class
fN(source.GetN()),        // number of enetries -1 => Uninitilized
fOrder(source.GetOrder()),// maximum moment of distributions (^n)
fX(new Double_t[source.GetOrder()]),//[fOrder] array of sums of x^n
fYx(new Double_t[source.GetOrder()]),//[fOrder] array of sums of (xy)^n
fY(new Double_t[source.GetOrder()]),//[fOrder] array of sums of y^n
fW(new Double_t[source.GetOrder()]){//[fOrder] array of sums of w^n (weights)
    // Copy constructor
    // Inputs:
    //   AliITSstatistics2 & source the source of this copy
    // Outputs:
    //   none.
    // Return:
    //   A copy of the source.

    if(this==&source) return;
    for(Int_t i=0;i<source.fOrder;i++){
        this->fX[i] = source.fX[i];
        this->fYx[i] = source.fYx[i];
        this->fY[i] = source.fY[i];
        this->fW[i] = source.fW[i];
    } // end for i
}
//______________________________________________________________________
void AliITSstatistics2::Reset(Int_t order){
    // Reset/zero all statistics variables statistics
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //   none.
    Int_t i;

    for(i=0;i<fOrder;i++) {
        fX[i] = 0.0;
        fY[i] = 0.0;
        fYx[i] = 0.0;
        fW[i] = 0.0;
    } // end for i
    fN = 0;
    if(order<0) return; // just zero
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
    if(order==0) return;
    fOrder = order;
    fX =  new Double_t[fOrder];
    fY =  new Double_t[fOrder];
    fYx = new Double_t[fOrder];
    fW =  new Double_t[fOrder];
    return;
}
//______________________________________________________________________
void AliITSstatistics2::AddValue(Double_t y,Double_t x,Double_t w=1.0){
    // add next x,y pair to statistics
    // Inputs:
    //   Double_t  y    y value of pair
    //   Double_t  x    x value of pair
    //   Double_t  w    weight of pair
    // Outputs:
    //   none.
    // Return:
    //   none.
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
//______________________________________________________________________
Double_t AliITSstatistics2::GetXNth(Int_t order)const{
    // This give the unbiased estimator for the RMS.
    // Inputs:
    //   Int_t order   the order of x^n value to be returned
    // Output:
    //   none.
    // Return:
    //   The value sum{x^n}.
    Double_t s;

    if(fW[0]!=0.0&&order<=fOrder) s = fX[order-1]/fW[0];
    else {
	s = 0.0;
	Error("GetXNth","error fOrder=%d fN=%d fW[0]=%f\n",
              fOrder,fN,fW[0]);
    } // end else
    return s;
}
//______________________________________________________________________
Double_t AliITSstatistics2::GetYNth(Int_t order)const{
    // This give the unbiased estimator for the RMS.
    // Inputs:
    //   Int_t order   the order of y^n value to be returned
    // Outputs:
    //   none.
    // Return:
    //  The value sum{y^n}
    Double_t s;

    if(fW[0]!=0.0&&order<=fOrder) s = fY[order-1]/fW[0];
    else {
	s = 0.0;
	Error("GetYNth","fOrder=%d fN=%d fW[0]=%f\n",
              fOrder,fN,fW[0]);
    } // end else
    return s;
}
//______________________________________________________________________
Double_t AliITSstatistics2::GetYXNth(Int_t order)const{
    // This give the unbiased estimator for the RMS.
    // Inputs:
    //   Int_t order   the order of (xy)^n value to be returned
    // Outputs:
    //   none.
    // Return:
    //  The value sum{(xy)^n}
    Double_t s;

    if(fW[0]!=0.0&&order<=fOrder) s = fYx[order-1]/fW[0];
    else {
	s = 0.0;
	Error("GetYXNth","fOrder=%d fN=%d fW[0]=%f\n",
	       fOrder,fN,fW[0]);
    } // end else
    return s;
}
//______________________________________________________________________
Double_t AliITSstatistics2::GetRMSX()const{
    // This give the unbiased estimator for the RMS.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //  The rms value
    Double_t x,x2,w,ww,s;

    x  = GetMeanX(); // first order
    x2 = GetXNth(2); // second order
    w  = fW[0];     // first order - 1.
    ww = fW[1];     // second order - 1.

    if(w*w==ww) return (-1.0);
    s = (x2-x*x)*w*w/(w*w-ww);
    return TMath::Sqrt(s);
}
//______________________________________________________________________
Double_t AliITSstatistics2::GetRMSY()const{
    // This give the unbiased estimator for the RMS.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //  The rms value
    Double_t x,x2,w,ww,s;

    x  = GetMeanY(); // first order
    x2 = GetYNth(2); // second order
    w  = fW[0];     // first order - 1.
    ww = fW[1];     // second order - 1.

    if(w*w==ww) return (-1.0);
    s = (x2-x*x)*w*w/(w*w-ww);
    return TMath::Sqrt(s);
}
//______________________________________________________________________
Double_t AliITSstatistics2::GetRMSYX()const{
    // This give the unbiased estimator for the RMS.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //  The rms value
    Double_t x,x2,w,ww,s;

    x  = GetMeanYX(); // first order
    x2 = GetYXNth(2); // second order
    w  = fW[0];     // first order - 1.
    ww = fW[1];     // second order - 1.

    if(w*w==ww) return (-1.0);
    s = (x2-x*x)*w*w/(w*w-ww);
    return TMath::Sqrt(s);
}
//______________________________________________________________________
Double_t AliITSstatistics2::GetErrorMeanY()const{
    //This is the error in the mean or the square root of the 
    // variance of the mean.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //  The error on the mean
    Double_t rms,w,ww,s;

    rms = GetRMSY();
    w   = fW[0];
    ww  = fW[1];
    s   = rms*rms*ww/(w*w);
    return TMath::Sqrt(s);
}
//______________________________________________________________________
Double_t AliITSstatistics2::GetErrorMeanX()const{
    //This is the error in the mean or the square root of the 
    // variance of the mean.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //  The error on the mean
    Double_t rms,w,ww,s;

    rms = GetRMSX();
    w   = fW[0];
    ww  = fW[1];
    s   = rms*rms*ww/(w*w);
    return TMath::Sqrt(s);
}
//______________________________________________________________________
Double_t AliITSstatistics2::GetErrorMeanYX()const{
    //This is the error in the mean or the square root of the 
    // variance of the mean.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //  The error on the mean
    Double_t rms,w,ww,s;

    rms = GetRMSYX();
    w   = fW[0];
    ww  = fW[1];
    s   = rms*rms*ww/(w*w);
    return TMath::Sqrt(s);
}
//______________________________________________________________________
Double_t AliITSstatistics2::GetErrorRMSY()const{
    // This is the error in the mean or the square root of the variance 
    // of the mean. at this moment this routine is only defined for 
    // weights=1.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //  The error on the rms
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
//______________________________________________________________________
Double_t AliITSstatistics2::GetErrorRMSX()const{
    // This is the error in the mean or the square root of the variance 
    // of the mean. at this moment this routine is only defined for 
    // weights=1.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //  The error on the rms
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
//______________________________________________________________________
Double_t AliITSstatistics2::GetErrorRMSYX()const{
    // This is the error in the mean or the square root of the variance 
    // of the mean. at this moment this routine is only defined for 
    // weights=1.
    // Inputs:
    //   none.
    // Outputs:
    //   none.
    // Return:
    //  The error on the rms
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
Double_t AliITSstatistics2::FitToLine(Double_t &a,Double_t &b)const{
    // fit to y = ax+b returns Chi squared or -1.0 if an error
    // Inputs:
    //   none.
    // Outputs:
    //   Double_t  a   The slope parameter
    //   Double_t  b   The intercept paramter
    // Return:
    //  The Chi^2 of the fit
    Double_t c,d,e,f,g,h;

    a = b = 0.0;
    if(fOrder<2 || fN<3){
        Error("FitToLine","Order=%d<2 or N=%d<3",fOrder,fN);
        return -1.0;
    } // end if
    c = GetWN(1);
    d = GetYN(1);
    e = GetXN(1);
    f = GetYXN(1);
    g = GetXN(2);
    h = c*g-e*e;
    b = d*g-f*e;
    a = c*f-d*e;
    if(h!=0.0){
	a = a/h;
	b = b/h;
    }else{
	Error("FitToLine","vertical line: fOrder=%d fN=%d "
              "GetWN(1)=%g X GetXN(2)=%g - GetXN(1)=%g^2 = 0",fOrder,fN,c,g,e);
	return -1.0;
    } // end if h
    return GetChiSquared(a,b);
}
//_______________________________________________________________________
Double_t AliITSstatistics2::GetChiSquared(Double_t a,Double_t b)const{
    //  returns Chi^2 value of data to line y=ax+b with given a,b
    // Inputs:
    //   Double_t  a   The slope parameter
    //   Double_t  b   The intercept paramter
    // Outputs::
    //   none.
    // Return:
    //  The Chi^2 of the fit
    Double_t c2;

    c2 = GetYN(2)+b*b*GetWN(1)+
        a*a*GetXN(2)-2.0*b*GetYN(1)-2.0*a*GetYXN(1)+2.0*b*a*GetXN(1);
    c2 /= (Double_t)fN - 2.0;
    return c2;
}
//______________________________________________________________________
void AliITSstatistics2::PrintAscii(ostream *os)const{
    // Print out class data values in Ascii Form to output stream
    // Inputs:
    //   ostream *os   Output stream where Ascii data is to be writen
    // Outputs:
    //   none.
    // Return:
    //   none.
    Int_t i;
#if defined __GNUC__
#if __GNUC__ > 2
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#else
#if defined __ICC || defined __ECC || defined __xlC__
    ios::fmtflags fmt;
#else
    Int_t fmt;
#endif
#endif

    *os << fN <<" "<< fOrder;
    fmt = os->setf(ios::scientific); // set scientific floating point output
    for(i=0;i<fOrder;i++) *os <<" "<< fX[i];
    for(i=0;i<fOrder;i++) *os <<" "<< fYx[i];
    for(i=0;i<fOrder;i++) *os <<" "<< fY[i];
    for(i=0;i<fOrder;i++) *os <<" "<< fW[i];
    os->flags(fmt); // reset back to old Formating.
    return;
}
//______________________________________________________________________
void AliITSstatistics2::ReadAscii(istream *is){
    // Read in class data values in Ascii Form to output stream
    // Inputs:
    //   istream *is   Input stream where Ascii data is to be read in from
    // Outputs:
    //   none.
    // Return:
    //   none.
    Int_t i;

    *is >> i >> fOrder;
    Reset(fOrder);
    fN = i;
    for(i=0;i<fOrder;i++) *is >> fX[i];
    for(i=0;i<fOrder;i++) *is >> fYx[i];
    for(i=0;i<fOrder;i++) *is >> fY[i];
    for(i=0;i<fOrder;i++) *is >> fW[i];
}
//______________________________________________________________________
ostream &operator<<(ostream &os,const AliITSstatistics2 &s){
    // Standard output streaming function
    // Inputs:
    //   ostream            &os  output steam
    //   AliITSstatistics2 &s class to be streamed.
    // Output:
    //   none.
    // Return:
    //   ostream &os  The stream pointer

    s.PrintAscii(&os);
    return os;
}
//______________________________________________________________________
istream &operator>>(istream &is,AliITSstatistics2 &s){
    // Standard inputput streaming function
    // Inputs:
    //   istream            &is  input steam
    //   AliITSstatistics2 &s class to be streamed.
    // Output:
    //   none.
    // Return:
    //   ostream &os  The stream pointer

    s.ReadAscii(&is);
    return is;
}
