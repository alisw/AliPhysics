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

#include "TObjArray.h"

class AliESDtrack;
class AliESDEvent;
class AliESDpid;

class TObjArray;

class AliTOFT0v1: public TObject {
public:
  
  AliTOFT0v1(AliESDpid *extPID=NULL); // default constructor
  AliTOFT0v1(AliESDEvent *event,AliESDpid *extPID=NULL); // overloaded constructor
  virtual ~AliTOFT0v1() ; // dtor
 
  Double_t* DefineT0(Option_t *option,Float_t pMinCut=3,Float_t pMaxCut=5); 
    
  void          SetMomBounds(Float_t pLow, Float_t pUp) { fLowerMomBound=pLow; fUpperMomBound=pUp;} // momenta are expressed in [GeV/c]
  void          SetTimeCorr(Float_t timecorr) {fTimeCorr=timecorr;} //in ns!!!
  Float_t       GetMomError(Int_t index, Float_t mom, Float_t texp) const;
  Double_t GetResult(Int_t i){if(i < 6) return fT0SigmaT0def[i]; else return -1.;};
/*   void  Print(Option_t* option) const ; */

  void   SetTimeResolution(Float_t /* timeres */){}; // obsolete

  void Init(AliESDEvent *event); // init

 private:

  Float_t ToCalculatePower(Float_t base, Int_t exponent) const ;
  Int_t   ToCalculatePower(Int_t base, Int_t exponent) const ;

  AliTOFT0v1(const AliTOFT0v1 &);
  AliTOFT0v1 & operator=(const AliTOFT0v1 &) ;

  Bool_t AcceptTrack(AliESDtrack *track); /* accept track */
  Float_t GetSigmaToVertex(AliESDtrack *track) const; /* get sigma to vertex */
  Bool_t CheckTPCMatching(AliESDtrack *track,Int_t imass) const;

  Float_t fLowerMomBound;   // momentum lower bound for selected primary tracks   
  Float_t fUpperMomBound;   // momentum upper bound for selected primary tracks 
  Float_t fTimeCorr;  // global time resolution used to calculate T0
  AliESDEvent* fEvent;      //evento per il quale si vuole calcolare il T0
  Double_t fT0SigmaT0def[6]; // array with the event information ([0]=event time -- [1] = sigma -- [2] = tracks on the TOF -- [3] = tracks used for the event time)
  
  AliESDpid *fPIDesd; // class with the detector response

  TObjArray *fTracks;   //! array of tracks
  TObjArray *fGTracks;  //! array of good tracks
  TObjArray *fTracksT0; //! array of tracks usefull for T0 estimate

  ClassDef(AliTOFT0v1,4);  // Calculate the time zero using TOF detector */
  
};

#endif 
