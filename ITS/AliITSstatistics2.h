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
    AliITSstatistics2();
    AliITSstatistics2(Int_t order);
	 AliITSstatistics2(AliITSstatistics2 &source); // copy  constructor
	 AliITSstatistics2& operator=(AliITSstatistics2 &source); // operator=	 
    virtual ~AliITSstatistics2();
    void Reset();
    void AddValue(Double_t y,Double_t x,Double_t w);
    Double_t GetXNth (Int_t order);
    Double_t GetYNth (Int_t order);
    Double_t GetYXNth(Int_t order);
    Double_t GetMeanY()  {
	                       // return mean y
	                       return GetYNth(1);
								 };
    Double_t GetMeanX()  {
	                       // return mean x
	                       return GetXNth(1);
								 };
    Double_t GetMeanYX() {
	                       // return mean Y*X
	                       return GetYXNth(1);
								 };
    Int_t GetN(){
	              // retrun the number of entries
	              return fN;
					 };
    Int_t GetOrder(){
	                  // return the maximum moment order
	                  return fOrder;
						  };
    Double_t GetXN (Int_t order){
	                              // returns x^n
	                              return fX[order-1];
										  };
    Double_t GetYN (Int_t order){
	                              // returns y^n
	                              return fY[order-1];
										  };
    Double_t GetYXN(Int_t order){
	                              // returns (yx)^n
	                              return fYx[order-1];
										  };
    Double_t GetWN (Int_t order){
	                              // returns w^n (weight)
	                              return fW[order-1];
										  };
    Double_t GetRMSY();
    Double_t GetRMSX();
    Double_t GetRMSYX();
    Double_t GetErrorMeanY();
    Double_t GetErrorMeanX();
    Double_t GetErrorMeanYX();
    Double_t GetErrorRMSY();
    Double_t GetErrorRMSX();
    Double_t GetErrorRMSYX();
    Double_t FitToLine(Double_t &a,Double_t &b);

 private:
    Int_t  fN;       // number of enetries
    Int_t  fOrder;   // maximum moment of distributions (^n)
    Double_t *fX;    //[fOrder] array of sums of x^n
    Double_t *fYx;   //[fOrder] array of sums of (xy)^n
    Double_t *fY;    //[fOrder] array of sums of y^n
    Double_t *fW;    //[fOrder] array of sums of w^n (weights)


    ClassDef(AliITSstatistics2,1)  //compute usueful statistics in 2D
};
#endif
