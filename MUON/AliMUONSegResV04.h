#ifndef MUONSegResV04_H
#define MUONSegResV04_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////
//  Segmentation and Response classes version 01   //
/////////////////////////////////////////////////////
#include "AliMUON.h"
#include "AliMUONSegResV01.h"
#include "TArrayF.h"
#include "TArrayI.h"
class AliMUONsegmentationV04 :
public AliMUONsegmentationV01 {
 public:
    AliMUONsegmentationV04(){}
    virtual ~AliMUONsegmentationV04(){}
    // Initialisation
    virtual void Init(AliMUONchamber*);
    // Test points for auto calibration
    void GiveTestPoints(Int_t &n, Float_t *x, Float_t *y);
    ClassDef(AliMUONsegmentationV04,1)
};
	

#endif






