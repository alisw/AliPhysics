#ifndef MUONSegResV05_H
#define MUONSegResV05_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////
//  Segmentation and Response classes version 01   //
/////////////////////////////////////////////////////
 
#include "AliMUON.h"
#include "AliMUONSegResV02.h"
#include "TArrayF.h"
#include "TArrayI.h"
class AliMUONsegmentationV05 :
public AliMUONsegmentationV02 {
 public:
    AliMUONsegmentationV05(){}
    virtual ~AliMUONsegmentationV05(){}
    // Initialisation
    virtual void Init(AliMUONchamber*);
    // Test points for auto calibration
    void GiveTestPoints(Int_t &n, Float_t *x, Float_t *y);
    ClassDef(AliMUONsegmentationV05,1)
};
#endif






