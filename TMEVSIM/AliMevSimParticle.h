
#ifndef ALIMEVSIMPARTICLE_H
#define ALIMEVSIMPARTICLE_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include <TPDGCode.h>

#include "TMevSimConverter.h"
#include "TMevSimPartTypeParams.h"

class AliMevSimParticle :public TMevSimPartTypeParams {

 public:
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  AliMevSimParticle();
  AliMevSimParticle(PDG_t pdg, Int_t multmean, Int_t multvc, 
		    Float_t tempmean, Float_t tempstdev, Float_t sigmamean,
		    Float_t sigmastdev, Float_t expvelmean, Float_t expvelstdev);
  
  virtual ~AliMevSimParticle();
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
  virtual void        SetPDG(PDG_t pdg);
  virtual PDG_t       GetPDG();
  
  ///////////////////////////////////////////////////////////////////////////////////////
  
 protected:
  
  PDG_t fPdg;
  TMevSimConverter *fConv;

 private:

  AliMevSimParticle(const AliMevSimParticle&); // Not implemented
  AliMevSimParticle& operator=(const AliMevSimParticle&); // Not implemented

  ClassDef(AliMevSimParticle,1)  
    
  ///////////////////////////////////////////////////////////////////////////////////////
    
};

#endif

