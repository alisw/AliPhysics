#ifndef ALITOFT0_H
#define ALITOFT0_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//_____________________________________________________________________________//
//                                                                             //
//  Task Class for calculating the time zero of interaction using TOF          //
//  The input file need the track length till TOF detector                     //
//  It can be done modifyng the AliTOFvj StepManager and the AliTOFHit class   //
//  as follow                                                                  //
//                                                                             //
//-- Author: F. Pierella                                                       //
//                                                                             //
//_____________________________________________________________________________//

#include "TTask.h"

class TString;

class AliTOFT0: public TTask {

public:
  AliTOFT0() ;          // ctor
  AliTOFT0(char* headerFile, Int_t nEvents=0) ; 
  AliTOFT0(const AliTOFT0 & tzero);
//////                  {( (AliTOFT0 &)tzero ).Copy(*this) ;} 
  virtual ~AliTOFT0() ; // dtor
  /*
  AliTOFT0 & operator = (const AliTOFT0 & rvalue)  {
    // assignement operator requested by coding convention but not needed
    abort() ;
    return *this ; 
  }
  */
  const char*   GetTZeroFile() const {return fT0File.Data();}  
  virtual void  Exec(Option_t *option); 
  void          SetNEvents(Int_t nEvents) {fNevents = nEvents;}
  void          SetTimeResolution(Float_t timeresolution) { fTimeResolution=timeresolution;}// timeresolution in [s] e.g. for 120 ps -> 1.2e-10
  Int_t         GetNEvents() const {return fNevents;}
  void          SetTZeroFile(char* file) ;
  void          SetMomBounds(Float_t pLow, Float_t pUp) { fLowerMomBound=pLow; fUpperMomBound=pUp;} // momenta are expressed in [GeV/c]
  virtual void  Print(Option_t* option) const ;
  Bool_t   operator == (const AliTOFT0 & tzero) const ;

 protected:

 private:
  Int_t   fNevents;         // Number of events for which calculate the T0
  Float_t fTimeResolution;  // global time resolution used to calculate T0
  Float_t fLowerMomBound;   // momentum lower bound for selected primary tracks 
  Float_t fUpperMomBound;   // momentum upper bound for selected primary tracks 
  TString fT0File ;         // output file; it contains for time being only 3 histos 
  TString fHeadersFile;     // input file

  ClassDef(AliTOFT0,1)  // Calculate the time zero using TOF detector

};

#endif // AliTOFT0_H
