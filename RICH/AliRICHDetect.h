#ifndef ALIRICHDETECT_H
#define ALIRICHDETECT_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////////
//   Reconstruction classes for set:RICH version 0     //
/////////////////////////////////////////////////////////

#include "AliRICH.h"
#include "TCanvas.h"

class AliRICHDetect;

class AliRICHDetect : public TObject {
    
 public:
  AliRICHDetect();
  AliRICHDetect(const char *name, const char *title);
  virtual       ~AliRICHDetect();
  void   Detect(Int_t nev);
  float Area(float theta,float OMEGA);
  Int_t  ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
  void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,long ndl, long ndh);
  Float_t SnellAngle(Float_t iangle);
  Float_t InvSnellAngle(Float_t rangle);
  void CreatePoints(Float_t theta, Float_t phi, Float_t omega, Float_t h);
  
 public:

  TCanvas *fc1;                   //Online reconstruction data
  TCanvas *fc2;                   //Online SPOT reconstruction data 
  TCanvas *fc3;                   //Online digits' coordinates data

  ClassDef(AliRICHDetect,1)  //Reconstruction module for :RICH version 0
	};


	
	
#endif
