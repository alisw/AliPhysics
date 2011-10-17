#ifndef AliT0CalibSeasonTimeShift_H
#define AliT0CalibSeasonTimeShift_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

////////////////////////////////////////////////
//  class for T0 calibration                 //
////////////////////////////////////////////////

#include "TNamed.h"
class TH1F;
class AliT0CalibSeasonTimeShift: public TNamed {

 public:
  AliT0CalibSeasonTimeShift();
  AliT0CalibSeasonTimeShift(const char* name);
  AliT0CalibSeasonTimeShift(const AliT0CalibSeasonTimeShift &calibda);
   AliT0CalibSeasonTimeShift & operator= (const AliT0CalibSeasonTimeShift &calibda);
  virtual ~AliT0CalibSeasonTimeShift();
   virtual void  Print(Option_t* option= "") const; 
  Float_t  MeanAC()        const {return fMeanPar[0];}
  Float_t  MeanA()        const {return fMeanPar[1];}
  Float_t  MeanC()        const {return fMeanPar[2];}
  Float_t  T0resolution()        const {return fMeanPar[3];}

  Float_t  SigmaAC()        const {return fSigmaPar[0];}
  Float_t  SigmaA()        const {return fSigmaPar[01];}
  Float_t  SigmaC()        const {return fSigmaPar[2];}
  Float_t  SigmaT0resolution()        const {return fSigmaPar[3];}
  Float_t *GetT0Means() { return fMeanPar;}
  Float_t *GetT0Sigmas() { return fSigmaPar;};
 
  Bool_t SetT0Par(Float_t par[4],Float_t spar[4] );
  Int_t SetT0Par(const char* filePhys , Float_t *cdbtime);
  
  void GetMeanAndSigma(TH1F* hist, Float_t &mean, Float_t &sigma); 

 protected:
  Float_t  fMeanPar[4];     
// [0] (T0A+T0C)/2; 
// [1] T0A corrected by primary vertex
// [2] T0A corrected by primary vertex
// [3] T0resolution
   Float_t  fSigmaPar[4];     
// [0] sigma (T0A+T0C)/2; 
// [1] sigma T0A corrected by primary vertex
// [2]sigma T0A corrected by primary vertex
// [3]sigma T0resolution


  ClassDef(AliT0CalibSeasonTimeShift,1)    // T0 Sensor Calibration data
};


#endif

