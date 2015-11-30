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
#include <iostream>
#include <TObject.h>
//#include <Rtypes.h>


class TClonesArray;
class AliQnCorrectionsQnVector;



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

  void SetEventPlaneDetector(Int_t det)  {fEventPlaneDetectorMask |= (1<<det);}
  void SetEventPlaneDetectorMask(ULong64_t mask) {fEventPlaneDetectorMask=mask;}

  // getters
  Float_t Phi()	    const  	    {return fPhi;}
  Float_t X()	    const  	    {return fX;}
  Float_t Y()	    const 	    {return fY;}
  Float_t Weight()const     {return fWeight;}
  Float_t Weight(Int_t method)	    const     {return fEqualizedWeight[method];}   // method 0: average equalized, 1: width equalized
  Int_t Id()	    const      {return fId;}
  Int_t Bin()	    const      {return fBin;}
  Bool_t  CheckEventPlaneDetector(Int_t flag) const;
  ULong64_t EventPlaneDetectorFlag() const {return fEventPlaneDetectorMask;}

  static void FillQvector(TClonesArray* det, Int_t ep, AliQnCorrectionsQnVector* q, Int_t weight=-1);
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
  ULong64_t fEventPlaneDetectorMask;  // Bit maps for the event plane subdetectors


  ClassDef(AliQnCorrectionsDataVector, 1);
 

};





#endif
