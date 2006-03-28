#ifndef ALIITSTRIGGER_H
#define ALIITSTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliTriggerDetector.h"

class AliITSgeom;

////////////////////////////////////////////////////////////////////////
//
// Version 1
// Modified by D. Elia, C. Jorgensen
// March 2006
//
// Version 0
// Written by J. Conrad, E. Lopez Torres
// October 2005
//
// AliITSTrigger: implementation of the SPD Fast-OR based triggers.
//
////////////////////////////////////////////////////////////////////////

class AliITSTrigger : public AliTriggerDetector
{
 public:
                   AliITSTrigger();   // constructor
        virtual   ~AliITSTrigger(){}  // destructor
   virtual void    CreateInputs();
   virtual void    Trigger();

private:

   Int_t fGlobalFOThreshold;         // minimum number of FOs to fire Global FO trigger
   Int_t fHighMultFOThreshold;       // minimum number of FOs to fire High Mult FO trigger

   void MultiplicityTriggers(TObjArray* digDet, TTree* treeD, AliITSgeom* geom);
//   void GeometryTriggers(TObjArray* digDet, TTree* treeD, AliITSgeom* geom);
   void GeometryTriggers();

  ClassDef( AliITSTrigger, 1 )  // ITS SPD Trigger Detector class
};

#endif
