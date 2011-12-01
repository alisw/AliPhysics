#ifndef ALIAODTZERO_H
#define ALIAODTZERO_H

//-------------------------------------------------------------------------
//     Container class for AOD TZERO data
//     Author: Filip Krizek
//     filip.krizek@cern.ch 18/11/2011
//-------------------------------------------------------------------------

#include <TObject.h>

class AliAODTZERO : public TObject 
{
public:
  AliAODTZERO();
  AliAODTZERO(const AliAODTZERO& source);
  AliAODTZERO &operator=(const AliAODTZERO& source);

  virtual ~AliAODTZERO() {};

  // Getters 
  //1st
  Double32_t GetT0TOF(Int_t i)  const {return fT0TOF[i];}
  const Double32_t * GetT0TOF() const {return fT0TOF;}
  //best
  Double32_t GetT0TOFbest(Int_t i)  const {return fT0TOFbest[i];}
  const Double32_t * GetT0TOFbest() const {return fT0TOFbest;}
 
  Bool_t GetBackgroundFlag() const {return fBackground;}
  Bool_t GetPileupFlag()     const {return fPileup;}
  Bool_t GetSatellite()      const {return fSattelite;}
  
  //Setters
  void SetT0TOF(Int_t icase, Double32_t time) { fT0TOF[icase] = time;}
  void SetT0TOFbest(Int_t icase, Double32_t time) { fT0TOFbest[icase] = time;}
   
  void SetBackgroundFlag(Bool_t back = false) {fBackground = back;}
  void SetPileupFlag(Bool_t back = false) {fPileup  = back;}
  void SetSatelliteFlag(Bool_t sat = false) { fSattelite = sat;}
       
  
  
protected:
  Double32_t   fT0TOF[3];    // interaction time in ps with 1st time( A&C, A, C)
  Bool_t       fPileup;      // pile-up flag
  Bool_t       fSattelite;   // sattelite flag
  Bool_t       fBackground;  // sattelite flag
  Double32_t   fT0TOFbest[3];// interaction time in ps ( A&C, A, C) with best time

  ClassDef(AliAODTZERO,1)
};

#endif
