#ifndef ALITOFT0V2_H
#define ALITOFT0V2_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//----------------------------------------------------------------------------//
//                                                                            //
//   Description: class to performe an event time measurment with TOF.        //
//                                                                            //
//----------------------------------------------------------------------------//

#include "TObject.h"

class AliESDtrack;
/* class AliTOFcalibHisto; */
class AliESDEvent;

class AliTOFT0v2: public TObject {
public:
  
  AliTOFT0v2() ;
  AliTOFT0v2(const AliTOFT0v2 & tzero);
  AliTOFT0v2 & operator=(const AliTOFT0v2 & tzero) ;
  AliTOFT0v2(AliESDEvent *event);
  virtual ~AliTOFT0v2() ; // dtor
 
  //  void SetCalib(AliTOFcalibHisto * const calib){fCalib = calib;};

  Double_t* DefineT0(Option_t *option); 
  
  void      SetTimeResolution(Double_t timeresolution);// timeresolution in [s] e.g. for 120 ps -> 1.2e-10
  
  Double_t GetTimeResolution() const {return fTimeResolution;}
  
  void          SetMomBounds(Float_t pLow, Float_t pUp) { fLowerMomBound=pLow; fUpperMomBound=pUp;} // momenta are expressed in [GeV/c]
  void          SetTimeCorr(Float_t timecorr) {fTimeCorr=timecorr;} //in ns!!!
  Float_t       GetMomError(Int_t index, Float_t mom, Float_t texp) const;
/*   void  Print(Option_t* option) const ; */

 private:

  Bool_t AcceptTrack(AliESDtrack *track); /* accept track */
  Float_t GetSigmaToVertex(AliESDtrack *track) const; /* get sigma to vertex */


  Float_t fLowerMomBound;   // momentum lower bound for selected primary tracks   
  Float_t fUpperMomBound;   // momentum upper bound for selected primary tracks 
  Double_t fTimeResolution;  // global time resolution used to calculate T0
  Float_t fTimeCorr;  // global time resolution used to calculate T0
  AliESDEvent* fEvent;      //evento per il quale si vuole calcolare il T0
/*   AliTOFcalibHisto *fCalib; // pointer to the class with the TOF time corrections */

  Double_t fT0SigmaT0def[4]; // array with the event information ([0]=event time -- [1] = sigma -- [2] = tracks on the TOF -- [3] = tracks used for the event time)
  
  ClassDef(AliTOFT0v2,2);  // Calculate the time zero using TOF detector */
  
};

#endif 
