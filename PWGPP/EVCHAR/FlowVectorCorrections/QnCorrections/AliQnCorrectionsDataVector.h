#ifndef ALIQNCORRECTIONS_DATAVECTOR_H
#define ALIQNCORRECTIONS_DATAVECTOR_H
/***************************************************************************
 * Package:       FlowVectorCorrections                                    *
 * Authors:       Jaap Onderwaater, GSI, jacobus.onderwaater@cern.ch       *
 *                Ilya Selyuzhenkov, GSI, ilya.selyuzhenkov@gmail.com      *
 *                Contributors are mentioned in the code where appropriate.*
 * Development:   2012-2015                                                *
 * See cxx source for GPL licence et. al.                                  *
 ***************************************************************************/
 
 
 

#include "AliQnCorrectionsConstants.h"
#include <TObject.h>
//#include <Rtypes.h>


class AliQnCorrectionsQnVector;
class TClonesArray;


//_____________________________________________________________________
class AliQnCorrectionsDataVector : public TObject {


 public:
  AliQnCorrectionsDataVector();
  ~AliQnCorrectionsDataVector();
    

  // setters
  void SetPhi(Float_t phi) {fPhi= phi;}
  void SetX(Float_t x) {fX= x;}
  void SetY(Float_t y) {fY= y;}
  void SetWeight(Float_t weight) {fWeight = weight;}
  void SetAverageEqualizedWeight( Float_t weight) {fEqualizedWeight[0] = weight;}
  void SetWidthEqualizedWeight(Float_t weight) {fEqualizedWeight[1] = weight;}
  void SetId(Int_t id) {fId = id;}
  void SetBin(Int_t bin) {fBin = bin;}

  // getters
  Float_t Phi()	    const  	    {return fPhi;}
  Float_t X()	    const  	    {return fX;}
  Float_t Y()	    const 	    {return fY;}
  Float_t Weight()const     {return fWeight;}
  Float_t Weight(Int_t method)	    const     {return fEqualizedWeight[method];}   // method 0: average equalized, 1: width equalized
  Int_t Id()	    const      {return fId;}
  Int_t Bin()	    const      {return fBin;}

  static void FillQvector(TClonesArray* dataVectorArray, AliQnCorrectionsQnVector* q, Int_t weight=-1);
  
 private:

  Float_t fPhi;
  Float_t fX;
  Float_t fY;
  Float_t fWeight;
  Float_t fEqualizedWeight[2];
  Int_t   fId;
  Int_t   fBin;
  //Float_t  fMinimumSignal;  


  ClassDef(AliQnCorrectionsDataVector, 2);
 

};





#endif
