#ifndef ALIITSSTATISTICS_H
#define ALIITSSTATISTICS_H
//////////////////////////////////////////////////////////////////////////
//  Alice ITS first detector alignment program.                         //
//                                                                      //
// version: 0.0.0 Draft.                                                //
// Date: April 18 1999                                                  //
// By: Bjorn S. Nilsen                                                  //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "TObject.h"

class AliITSstatistics : public TObject {
//

 public:
  AliITSstatistics();
  AliITSstatistics(Int_t order);
  AliITSstatistics(AliITSstatistics &source); // copy  constructor
  AliITSstatistics& operator=(AliITSstatistics &source); // operator=
  virtual ~AliITSstatistics();
  void Reset();
  void AddValue(Double_t x,Double_t w);
  void AddValue(Double_t x){
                            // Default weight of 1
                            AddValue(x,1.0);
									 } 
  void AddValue(Float_t x,Float_t w){
                                     //float
                                     AddValue((Double_t)x,(Double_t)w);
												 } 
  void AddValue(Float_t x){
                           // floats default weight of 1
                           AddValue((Double_t)x,(Double_t)1.0);} 
  Double_t GetNth(Int_t order);
  Double_t GetMean() {
                      // returns the mean
                      return GetNth(1);
							};
  Int_t GetN(){
               // returns the number of entries
               return fN;
              };
  Int_t GetOrder(){
                   // returns the order of the moment of the distribution
                   return 
						 fOrder;
                  };
  Double_t GetXN(Int_t order){
                              // returns X^N
                              return fx[order-1];
                              };
  Double_t GetWN(Int_t order){
                              // returns W^N
                              return fw[order-1];
                             };
  Double_t GetRMS();
  Double_t GetErrorMean();
  Double_t GetErrorRMS();

 private:
  Double_t *fx;   // fx array of x moments
  Double_t *fw;   // fw array of weight by moment
  Int_t    fN;    // fN number of enetries
  Int_t    fOrder;// fOrder maximum allowed moment

  ClassDef(AliITSstatistics,1)// A class to do simple statistics calculations
};
#endif
