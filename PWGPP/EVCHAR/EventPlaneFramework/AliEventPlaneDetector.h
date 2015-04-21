/*
***********************************************************
    Detector class that contains detector information (angles and weights)
    Contact: Jaap Onderwaater, j.onderwaater@gsi.de, jacobus.onderwaater@cern.ch
***********************************************************
*/

#ifndef ALIEVENTPLANEDETECTOR_H
#define ALIEVENTPLANEDETECTOR_H

//#include <iostream>
#include <TObject.h>
#include <Rtypes.h>

//const Int_t fgkEPMaxHarmonics = 6;
const Int_t fgkEPMaxDetectors = 20;



//_____________________________________________________________________
class AliEventPlaneDetector : public TObject {


 public:
  AliEventPlaneDetector();
  ~AliEventPlaneDetector();
    

  // setters
  void SetPhi(Float_t phi) {fPhi= phi;}
  void SetX(Float_t x) {fX= x;}
  void SetY(Float_t y) {fY= y;}
  void SetWeight(Float_t weight) {fWeight = weight;}
  void SetAverageEqualizedWeight(Int_t epDet, Float_t weight) {fEqualizedWeight[epDet][0] = weight;}
  void SetWidthEqualizedWeight(Int_t epDet, Float_t weight) {fEqualizedWeight[epDet][1] = weight;}
  void SetId(Int_t id) {fId = id;}
  void SetBin(Int_t bin) {fBin = bin;}

  void SetEventPlaneDetector(Int_t det) { 
      fEventPlaneDetectorMask |= (1<<det);
  }

  // getters
  Float_t Phi()	    const  	    {return fPhi;}
  Float_t X()	    const  	    {return fX;}
  Float_t Y()	    const 	    {return fY;}
  Float_t Weight()const     {return fWeight;}
  Float_t Weight(Int_t epDet, Int_t method)	    const     {return fEqualizedWeight[epDet][method];}   // method 0: average equalized, 1: width equalized
  Int_t Id()	    const      {return fId;}
  Int_t Bin()	    const      {return fBin;}
  Bool_t  CheckEventPlaneDetector(Int_t flag) const;
  ULong64_t EventPlaneDetectorFlag() const {return fEventPlaneDetectorMask;}

  
 private:

  Float_t fPhi;
  Float_t fX;
  Float_t fY;
  Float_t fWeight;
  Float_t fEqualizedWeight[fgkEPMaxDetectors][2];
  Int_t   fId;
  Int_t   fBin;
  //Float_t  fMinimumSignal;  
  ULong64_t fEventPlaneDetectorMask;  // Bit maps for the event plane subdetectors


  ClassDef(AliEventPlaneDetector, 1);
 

};




////_______________________________________________________________________________
inline Bool_t AliEventPlaneDetector::CheckEventPlaneDetector(Int_t det) const {
  //
  // Check the status of the event plane for a given detector and harmonic
  //
  return (det<fgkEPMaxDetectors ? (fEventPlaneDetectorMask&(1<<(det))) : kFALSE);
}

#endif
