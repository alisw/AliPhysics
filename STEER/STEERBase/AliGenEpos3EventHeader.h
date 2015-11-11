#ifndef ALIGENEPOS3EVENTHEADER_H
#define ALIGENEPOS3EVENTHEADER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

// Header for EPOS v3.111 generated event.
//
// Author: Natalia Zhigareva <Natalia.Zhigareva@cern.ch>

#include "AliGenEventHeader.h"
#include "AliCollisionGeometry.h" 

class TGenerator;

class AliGenEpos3EventHeader : public AliGenEventHeader, public AliCollisionGeometry
{
public:
    AliGenEpos3EventHeader(const char* name);
    AliGenEpos3EventHeader();
    virtual ~AliGenEpos3EventHeader(){}

    // Getters:
    Int_t GetIversn()  const {return fIversn;}
    Int_t GetLaproj()  const {return fLaproj;}
    Int_t GetMaproj()  const {return fMaproj;}
    Int_t GetLatarg()  const {return fLatarg;}
    Int_t GetMatarg()  const {return fMatarg;}
    Float_t GetEngy()  const {return fEngy;}
    Int_t GetNfull()   const {return fNfull;}
    Int_t GetNfreeze() const {return fNfreeze;}
    Float_t GetBim()   const {return fBim;}

    //Setters:
    void SetIversn (Int_t iversn) {fIversn = iversn;}
    void SetLaproj (Int_t laproj) {fLaproj = laproj;}
    void SetMaproj (Int_t maproj) {fMaproj = maproj;}
    void SetLatarg (Int_t latarg) {fLatarg = latarg;}
    void SetMatarg (Int_t matarg) {fMatarg = matarg;}
    void SetEngy   (Float_t engy) {fEngy = engy;}
    void SetNfull  (Int_t nfull)  {fNfull = nfull;}
    void SetNfreeze (Int_t nfreeze) {fNfreeze = nfreeze;}
    void SetBim    (Int_t bim)    {fBim = bim;}

protected:
    
private:
  
//  Parameters in EPOS Tree:
//  --------------teposhead---------------
  Int_t fIversn;  //EPOS version number
  Int_t fLaproj;  //atomic number projectile
  Int_t fMaproj;  //mass number projectile
  Int_t fLatarg;  //atomic number target
  Int_t fMatarg;  //mass number target
  Float_t fEngy;  //energy in the CMS in GeV
  Int_t fNfull;   //number of full events
  Int_t fNfreeze; //number of freeze outs per full event
//--------------teposevent---------------
  Float_t fBim;   //impact parameter

  ClassDef(AliGenEpos3EventHeader,1) //event header for EPOS event
};

#endif /*ALIGENEPOS3EVENTHEADER_H_*/

  
  
  
  
  
  
  
  
  
  
