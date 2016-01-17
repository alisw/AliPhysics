#ifndef ALIITSUCHIP_H
#define ALIITSUCHIP_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliITSUChip.h 53509 2011-12-10 18:55:52Z masera $ */
///////////////////////////////////////////////////////////////////////
//                                                                   //
//  Class AliITSUChip                                                //
//  The main function of chips is to simulate DIGITS from            //
//  GEANT HITS and produce POINTS from DIGITS                        //
//  It also make fast simulation without use of DIGITS               //
//                                                                   //
///////////////////////////////////////////////////////////////////////

#include "AliITSMFTChip.h"
#include "AliITSUGeomTGeo.h"

class AliITSUChip: public AliITSMFTChip {

public:
  AliITSUChip() : AliITSMFTChip() {}
  AliITSUChip(Int_t idx, AliITSUGeomTGeo* tg): AliITSMFTChip(idx,tg) {}
  virtual ~AliITSUChip() {}

 protected:
    AliITSUChip(const AliITSUChip &source); 
    AliITSUChip& operator=(const AliITSUChip &source); 
    ClassDef(AliITSUChip,2) // Copy the hits into a more useful order
};

#endif



