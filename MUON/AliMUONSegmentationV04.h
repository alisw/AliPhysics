#ifndef ALIMUONSEGMENTATIONV04_H
#define ALIMUONSEGMENTATIONV04_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////////
//  Segmentation and Response classes version 04   //
/////////////////////////////////////////////////////
 
#include "AliMUONSegmentationV01.h"

class AliMUONSegmentationV04 :
public AliMUONSegmentationV01 {
 public:
    AliMUONSegmentationV04(){}
    virtual ~AliMUONSegmentationV04(){}
    // Initialisation
    virtual void Init(Int_t chamber);
    // Test points for auto calibration
    void GiveTestPoints(Int_t &n, Float_t *x, Float_t *y) const;
    ClassDef(AliMUONSegmentationV04,1) // Segmentation zones are rectangular modules
};
	

#endif







