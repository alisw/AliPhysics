#ifndef ALIMUONFACTORY_H
#define ALIMUONFACTORY_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

////////////////////////////////////////////////////////////
//  Factory for muon chambers, segmentations and response //
////////////////////////////////////////////////////////////
#include "AliDetector.h"
// #include "AliMUONTriggerCircuit.h" // cp

class AliMUONChamber;
class AliMUON;

class AliMUONFactory : public  TObject {
 public:
    static void Build(AliMUON* where, const char* what);
 protected:
    ClassDef(AliMUONFactory,0)  // MUON Factory for Chambers and Segmentation
};
#endif















