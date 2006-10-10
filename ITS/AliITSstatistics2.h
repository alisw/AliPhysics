#ifndef ALIITSSTATISTICS2_H
#define ALIITSSTATISTICS2_H
//////////////////////////////////////////////////////////////////////////
//  Alice ITS first detector alignment program.                         //
//                                                                      //
// version: 0.0.0 Draft.                                                //
// Date: April 18 1999                                                  //
// By: Bjorn S. Nilsen                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "TObject.h"

class AliITSstatistics2 : public TObject {
//
 public:
    AliITSstatistics2(); // default constructor
    AliITSstatistics2(Int_t order); // Standard constructor
    AliITSstatistics2(AliITSstatistics2 &source); // copy  constructor
    AliITSstatistics2& operator=(AliITSstatistics2 &source); // operator=
    virtual ~AliITSstatistics2();// destructor
    void Reset(Int_t order=-1); // Resets statictisc order=-1 do not change order.
    void AddValue(Double_t y,Double_t x,Double_t w);
    Double_t GetXNth (Int_t order)const;
    Double_t GetYNth (Int_t order)const;
    Double_t GetYXNth(Int_t order)const;
    Double_t GetMeanY()const{// return mean y
        return GetYNth(1);};
    Double_t GetMeanX()const{// return mean x
        return GetXNth(1);};
    Double_t GetMeanYX()const{// return mean Y*X
        return GetYXNth(1);};
    Int_t GetN()const{// retrun the number of entries
        return fN;};
    Int_t GetOrder()const{// return the maximum moment order
        return fOrder;};
    Double_t GetXN (Int_t order)const{// returns x^n
        return fX[order-1];};
    Double_t GetYN (Int_t order)const{// returns y^n
        return fY[order-1];};
    Double_t GetYXN(Int_t order)const{// returns (yx)^n
        return fYx[order-1];};
    Double_t GetWN (Int_t order)const{// returns w^n (weight)
        return fW[order-1];};
    Double_t GetRMSY()const; // Returns Y rms value
    Double_t GetRMSX()const; // Returns X rms value
    Double_t GetRMSYX()const;// Return XY rms value
    Double_t GetErrorMeanY()const; // returns error on y mean
    Double_t GetErrorMeanX()const; // returns error on x mean
    Double_t GetErrorMeanYX()const;// returns error on xy mean
    Double_t GetErrorRMSY()const;  // returns error on y rms
    Double_t GetErrorRMSX()const;  // returns error on x rms
    Double_t GetErrorRMSYX()const; // returns error on xy rms
    // returns Chi^2 value of data to line y=ax+b with given a,b
    Double_t GetChiSquared(Double_t a,Double_t b)const; 
    // Fits data to a line of the form y=ax+b and returns values
    // for a and b.
    Double_t FitToLine(Double_t &a,Double_t &b)const;
    // Print class in ascii form to stream
    virtual void PrintAscii(ostream *os)const;
    // Read in class in ascii form from stream
    virtual void ReadAscii(istream *is);

 private:
    Int_t  fN;       // number of enetries
    Int_t  fOrder;   // maximum moment of distributions (^n)
    Double_t *fX;    //[fOrder] array of sums of x^n
    Double_t *fYx;   //[fOrder] array of sums of (xy)^n
    Double_t *fY;    //[fOrder] array of sums of y^n
    Double_t *fW;    //[fOrder] array of sums of w^n (weights)

    ClassDef(AliITSstatistics2,1)  //compute usueful statistics in 2D
};
// Input and output function for standard C++ input/output.
ostream &operator<<(ostream &os,const AliITSstatistics2 &s);
istream &operator>>(istream &is,AliITSstatistics2 &s);
#endif
