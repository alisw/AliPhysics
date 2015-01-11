/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

#ifndef ALITPCLASER_H
#define ALITPCLASER_H

////////////////////////////////////////////////
//            Laser for TPCv2                 //
////////////////////////////////////////////////

 
#include "AliTPCv2.h"

class AliTPCLaser : public AliTPCv2 {

public:
  AliTPCLaser():AliTPCv2(),   
    fNelPerCollision(10),
    fLaserPID(13),
    fCollisionsPerCm(20)  {}
  AliTPCLaser(const char *name, const char *title);
  virtual      ~AliTPCLaser() {}
  
  virtual void  StepManager();

  virtual Int_t   GetNelPerCollision() const {return fNelPerCollision;}
  virtual Int_t   GetLaserPID() const {return fLaserPID;}
  virtual Float_t GetCollisionsPerCm() const {return fCollisionsPerCm;}

  virtual void SetNelPerCollision(Int_t nel) {fNelPerCollision = nel;}
  virtual void SetLaserPID(Int_t pid) {fLaserPID = pid;}
  virtual void SetCollisionsPerCm(Int_t ncol) {fCollisionsPerCm = ncol;}
  
 private:
  Int_t   fNelPerCollision;  // Fixed number of electrons per collision 
  Int_t   fLaserPID;         // PID of laser  
  Float_t fCollisionsPerCm;  // Number of primary interactions per cm
  ClassDef(AliTPCLaser,2)  // For Laser
};

#endif
