////////////////////////////////////////////////
//                                            //
//  Hit class for CRT                         //
//  Interface                                 //
//  Getters, Setters and member variables     //
//  declared here                             //
//                                            //
////////////////////////////////////////////////

#ifndef ALICRTHIT_H
#define ALICRTHIT_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliHit.h"

 
class AliCRThit : public AliHit {
  
public:
  AliCRThit() {}
  AliCRThit(Int_t shunt, Int_t track, Int_t* vol, Float_t *hits);
  AliCRThit(const AliCRThit & hit) ;
  virtual ~AliCRThit() {}

       // getters for AliCRThit object
  Float_t   Getnmou() const {return fnmou;}
  Float_t   Getmtyp()  const {return fmtyp;}
  Float_t Getxpit()     const {return fxpit;}
  Float_t Getypit()     const {return fypit;}
  Float_t Getzpit()     const {return fzpit;}
  Float_t Getpxug()     const {return fpxug;}
  Float_t Getpyug()     const {return fpyug;}
  Float_t Getpzug()     const {return fpzug;}
  Float_t Getlay()     const {return flay;}
  Float_t Getxver()     const {return fxver;}
  Float_t Getyver()     const {return fyver;}
  Float_t Getzver()     const {return fzver;}

  Int_t fVolume;
  Int_t fCopy;
  Float_t      fnmou; 
  Float_t      fmtyp;
  Float_t    fxpit;      //
  Float_t    fypit;    // 
  Float_t    fzpit;    // 
  Float_t    fpxug;     //
  Float_t    fpyug;      //
  Float_t    fpzug;      //
  Float_t    flay;    //
  Float_t    fxver;      //
  Float_t    fyver;      //
  Float_t    fzver;      //
  
private:
  
  ClassDef(AliCRThit,1)  // Hit for CRT (ACORDE)

};

#endif /* ALICRTHIT_H */
