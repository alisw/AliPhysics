#ifndef ALIMUONSEGMENTATIONV05_H
#define ALIMUONSEGMENTATIONV05_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////
//  Segmentation                      version 05   //
/////////////////////////////////////////////////////
 
#include "AliMUONSegmentationV02.h"

class AliMUONSegmentationV05 :
public AliMUONSegmentationV02 {
 public:
    AliMUONSegmentationV05(){}
    virtual ~AliMUONSegmentationV05(){}
    // Initialisation
    virtual void Init(Int_t chamber);
    // Test points for auto calibration
    void GiveTestPoints(Int_t &n, Float_t *x, Float_t *y);
    ClassDef(AliMUONSegmentationV05,1)// Segmentation zones are rectangular modules 
};
#endif






