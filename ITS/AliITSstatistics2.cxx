/**************************************************************************
 * Copyright(c) 2007-2009, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id$ */

//////////////////////////////////////////////////////////////////////////
//  Alice ITS class to help keep statistical information. Can also be   //
// used to fit data to lines and other 2 dimentional sytistical         //
// operations.                                                          //
//                                                                      //
// version: 0.0.0 Draft.                                                //
// Date: April 18 1999                                                  //
// By: Bjorn S. Nilsen                                                  //
// Updated: 1.0.0, Date: September 6 2007, By: Bjorn S. Nilsen          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <stdio.h>     //  ios::fmtflags fmt used in PrintAscii
#include "Riostream.h" // IO functions.
#include "TMath.h"     // TMath::Sqrt() function used.
#include "AliITSstatistics2.h" // Also defined TObject {base class}

ClassImp(AliITSstatistics2)

//
AliITSstatistics2::AliITSstatistics2() : 
TObject(), // Base Class
fN(-1),    // number of enetries -1 => Uninitilized
fOrder(0), // maximum moment of distributions (^n)
fX(0),     //[fOrder] array of sums of x^n
fYx(0),    //[fOrder] array of sums of (xy)^n
fY(0),     //[fOrder] array of sums of y^n
fW(0)      //[fOrder] array of sums of w^n (weights)
//,fDig(5)   // The number of significant digits to keep
//,fOver(0)  //! In case of numerical precistion problems
{
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
fW(new Double_t[order])     //[fOrder] array of sums of w^n (weights)
//,fDig(5)   // The number of significant digits to keep
//,fOver(0)                   //! In case of numeerical precistion problems
{   // constructor to maximum moment/order order
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
    // if(fOver!=0) delete fOver; fOver=0;
    // fDig=0;
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
    // this->fDig = source.fDig;
    // if(fOver!=0) this->fOver = new AliITSstatistics2(*(source.fOver));
    // else fOver=0;
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
fW(new Double_t[source.GetOrder()]) //[fOrder] array of sums of w^n (weights)
//,fDig(source.fDig)                  // The number of significant digits to keep
//,fOver(0)             //! In case of numerical precistion problems
{
    // Copy constructor
    // Inputs:
    //   AliITSstatistics2 & source the source of this copy
    // Outputs:
    //   none.
    // Return:
    //   A copy of the source.

    for(Int_t i=0;i<source.fOrder;i++){
        this->fX[i] = source.fX[i];
        this->fYx[i] = source.fYx[i];
        this->fY[i] = source.fY[i];
        this->fW[i] = source.fW[i];
    } // end for i
    //if(fOver!=0) this->fOver = new AliITSstatistics2(*(source.fOver));
    return;
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
    //if(fOver!=0) delete fOver; fOver = 0;
    return;
}
//----------------------------------------------------------------------
/*
void SetSignificantDigits(Int_t d){
    // Sets the number of significant digits. If adding a value to
    // one of this class' arrays looses significance at the fDig
    // level, a new instance of this class is created to keep
    // signigicance at or better than fDig level. if fDig<0, then
    // this feature it disabled and significance can be lost.
    // Inputs:
    //    Int_t  d   The new significance level
    // Outputs:
    //    none.
    // Return:
    //    none.

    fDig = d;
}
 */
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
    /*  If problem with precision, then creat/fill fOver
	as a partical sum to be added to "this" later.
	if(????fDig){
	if(fOver==0){
	    fOver = new AliITSstatistics2(fOrder);
	} // end if fOver==0
	fOver->AddValue(y,x,w);
	return;
	} // end if(???)
     */
    fN++;
    for(i=0;i<GetOrder();i++){
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

    if(GetWN(1)!=0.0 && order<=GetOrder()) s = GetXN(order)/GetWN(1);
    else {
	s = 0.0;
	Error("GetXNth","error fOrder=%d fN=%d fW[0]=%f\n",
              GetOrder(),GetN(),GetWN(1));
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

    if(GetWN(1)!=0.0&&order<=GetOrder()) s = GetYN(order)/GetWN(1);
    else {
	s = 0.0;
	Error("GetYNth","fOrder=%d fN=%d fW[0]=%f\n",
              GetOrder(),GetN(),GetWN(1));
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

    if(GetWN(1)!=0.0&&order<=GetOrder()) s = GetYXN(order)/GetWN(1);
    else {
	s = 0.0;
	Error("GetYXNth","fOrder=%d fN=%d fW[0]=%f\n",
	      GetOrder(),GetN(),GetWN(1));
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
    w  = GetWN(1);   // first order
    ww = GetWN(2);   // second order

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
    w  = GetWN(1);   // first order
    ww = GetWN(2);   // second order

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
    w  = GetWN(1);   // first order
    ww = GetWN(2);     // second order

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
    w   = GetWN(1);
    ww  = GetWN(2);
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
    w   = GetWN(1);
    ww  = GetWN(2);
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
    w   = GetWN(1);
    ww  = GetWN(2);
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

    if(GetWN(1)!=(Double_t)GetN()||GetN()<4) return (-1.);
    x  = GetMeanY(); // first order
    x2 = GetYNth(2); // second order
    w  = GetWN(1);   // first order
    ww = GetWN(2);   // second order
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

    if(GetWN(1)!=(Double_t)GetN()||GetN()<4) return (-1.);
    x  = GetMeanX(); // first order
    x2 = GetXNth(2); // second order
    w  = GetWN(1);   // first order
    ww = GetWN(2);   // second order
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

    if(GetWN(1)!=(Double_t)GetN()||GetN()<4) return (-1.);
    x  = GetMeanYX(); // first order
    x2 = GetYXNth(2); // second order
    w  = GetWN(1);    // first order
    ww = GetWN(2);    // second order
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
Double_t AliITSstatistics2::FitToLine(Double_t &a,Double_t &ea,
				      Double_t &b,Double_t &eb)const{
    // fit to y = ax+b returns Chi squared or -1.0 if an error.
    // The fitting is done by analitically minimizing 
    /*
      Begin_Latex
      \begin{equation*}
      \Chi^{2}=\sum_{i} (y_{i}-a x_{i} -b)^{2} w_{i}
      \end{equation*}
      Where if the weight used in 
      AliITSstatistics2::AddValue(Double_t y,Double_t x,Double_t w=1.0)
      is of the form
      \begin{equation*}
      w_{i}=\frac{1}{\delta y^{2}}.
      \end{equation*}
      Then we get the typicall chi square minimization.
      End_Latex
     */
    // Inputs:
    //   none.
    // Outputs:
    //   Double_t  a   The slope parameter
    //   Double_t  ea  Error on fitted slope parameter
    //   Double_t  b   The intercept paramter
    //   Double_t  eb  Error on fitted intercept parameter
    // Return:
    //  The Chi^2 of the fit
    Double_t c,d,e,f,g,h;

    a = ea = b = eb = 0.0;
    if(GetOrder()<2 || GetN()<3){
        Error("FitToLine","Order=%d<2 or N=%d<3",GetOrder(),GetN());
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
    if(h==0.0){
	Error("FitToLine","vertical line: fOrder=%d fN=%d "
              "GetWN(1)=%g X GetXN(2)=%g - GetXN(1)=%g^2 = 0",
	      GetOrder(),GetN(),c,g,e);
	return -1.0;
    } // end if h
    a = a/h;
    b = b/h;
    // Now for the errors.
    ea = c*c*g+(a*a-1.0)*c*e*e;
    ea = ea/(h*h);
    if(ea<0.0){
      Error("FitToLine","ea=%g is less than zero",ea);
      return -2.0;
    } // end if ea<0
    ea = TMath::Sqrt(ea);
    eb = c*g*g-2.0*d*e*g-2.0*(1.0-b)*c*e*e*g+2.0*(1.0-b)*d*e*e*e+
          GetYN(2)*e*e+(1.0-b)*(1.0-b)*c*e*e*e*e;
    eb = eb/(h*h);
    if(eb<0.0){
      Error("FitToLine","eb=%g is less than zero",eb);
      return -2.0;
    } // end if ea<0
    eb = TMath::Sqrt(eb);
    c = GetChiSquared(a,b);
    if(c<=0.0){ // must be a numerical precision problem.
    } // end if
    return c;
}
//_______________________________________________________________________
Double_t AliITSstatistics2::GetChiSquared(Double_t a,Double_t b)const{
    //  Returns Chi^2 value of data to line y=ax+b with given a,b.
    /*
      Begin_Latex
      Note: The Chi^2 value is computed from the expression
      \begin{equation*}
      \chi^{2}=\sum_{i}{w_{i}y_{i}^{2}} + b^{2}\sum_{i}{w_{i}}
                -2b\sum_{i}{w_{i}y_{i}}-2a\sum_{i}{w_{i}y_{i}x_{i}}
                +2ab\sum_{i}{w_{i}x_{i}}
                +a^{2}\sum_{i}w_{i}x_{i}^{2}
      \end{equation*}
      and not form the expression
      \begin{equation*}
      \chi^{2}= \sum_{i}{(y_{i}-ax_{i}-b)^{2}w_{i}.
      \end{equation*}
      Consiquently, there are occations when numerically these
      two expressions will not agree. In fact the form code here
      can give negitive values. This happens when the numerical
      significance is larger than the $\chi^{2}$ value. This should
      not be confused the the error values which can be returned.
      At present there is no check on the numberical significance
      of any results.
      End_Latex
     */
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
    c2 /= (Double_t)GetN() - 2.0;
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
    ios::fmtflags fmt;  // Standard IO format object, required for output.
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

    *os << fN <<" "<< GetOrder();
    fmt = os->setf(ios::scientific); // set scientific floating point output
    for(i=0;i<GetOrder();i++) *os <<" "<< GetXN(i+1);
    for(i=0;i<GetOrder();i++) *os <<" "<< GetYXN(i+1);
    for(i=0;i<GetOrder();i++) *os <<" "<< GetYN(i+1);
    for(i=0;i<GetOrder();i++) *os <<" "<< GetWN(i+1);
    //if(fOver!=0) { *os << " " << fDig;
    //*os << " " << fOver;
    // } else *os << " " << fDig;
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
    //*is >> fDig;
    // if(fDig>0) *is >> fOver;
    // else fDig *= -1;
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

