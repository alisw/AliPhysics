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
  virtual void   Detect(Int_t nev, Int_t type);
  float Area(float theta,float OMEGA);
  Int_t Fiducial(Float_t x, Float_t y, Float_t theta, Float_t phi, Float_t height, Float_t maxOmega, Float_t minOmega);

  virtual Int_t  ***i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
  virtual void free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,long ndl, long ndh);
  Float_t SnellAngle(Float_t iangle);
  Float_t InvSnellAngle(Float_t rangle);
  void CreatePoints(Float_t theta, Float_t phi, Float_t omega, Float_t h);
  
 public:

  TCanvas *fc1;                   //Online reconstruction data
  TCanvas *fc2;                   //Online SPOT reconstruction data 
  TCanvas *fc3;                   //Online digits' coordinates data
  TCanvas *fc4;                   //Online mesh activation data

  ClassDef(AliRICHDetect,1)  //Reconstruction module for :RICH version 0
	};


	
	
#endif
