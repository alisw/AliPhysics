#ifndef ALITOFT0V1_H
#define ALITOFT0V1_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//----------------------------------------------------------------------------//
//                                                                            //
//   Description: class to performe an event time measurment with TOF.        //
//                                                                            //
//----------------------------------------------------------------------------//

#include "TObject.h"

class AliESDtrack;
class AliESDEvent;

class AliTOFT0v1: public TObject {
public:
  
  AliTOFT0v1() ; // default constructor
  AliTOFT0v1(AliESDEvent *event); // overloaded constructor
  virtual ~AliTOFT0v1() ; // dtor
 
  Double_t* DefineT0(Option_t *option); 
  Double_t* DefineT0(Option_t *option,Float_t pMinCut,Float_t pMaxCut=1.5); 
  
  void      SetTimeResolution(Double_t timeresolution);// timeresolution in [s] e.g. for 120 ps -> 1.2e-10
  
  Double_t GetTimeResolution() const {return fTimeResolution;}
  
  void          SetMomBounds(Float_t pLow, Float_t pUp) { fLowerMomBound=pLow; fUpperMomBound=pUp;} // momenta are expressed in [GeV/c]
  void          SetTimeCorr(Float_t timecorr) {fTimeCorr=timecorr;} //in ns!!!
  Float_t       GetMomError(Int_t index, Float_t mom, Float_t texp) const;
  Double_t GetResult(Int_t i){if(i < 4) return fT0SigmaT0def[i]; else return -1.;};
/*   void  Print(Option_t* option) const ; */

  void Init(AliESDEvent *event); // init

 private:

  AliTOFT0v1(const AliTOFT0v1 &);
  AliTOFT0v1 & operator=(const AliTOFT0v1 &) ;

  Bool_t AcceptTrack(AliESDtrack *track); /* accept track */
  Float_t GetSigmaToVertex(AliESDtrack *track) const; /* get sigma to vertex */


  Float_t fLowerMomBound;   // momentum lower bound for selected primary tracks   
  Float_t fUpperMomBound;   // momentum upper bound for selected primary tracks 
  Double_t fTimeResolution;  // global time resolution used to calculate T0
  Float_t fTimeCorr;  // global time resolution used to calculate T0
  AliESDEvent* fEvent;      //evento per il quale si vuole calcolare il T0
  Double_t fT0SigmaT0def[4]; // array with the event information ([0]=event time -- [1] = sigma -- [2] = tracks on the TOF -- [3] = tracks used for the event time)
  
  ClassDef(AliTOFT0v1,2);  // Calculate the time zero using TOF detector */
  
};

#endif 
