#ifndef ALIITSTRIGGER_H
#define ALIITSTRIGGER_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

/////////////////////////////////////////////////
//
//  ITS SPD Trigger Detector Class
//
//
//
/////////////////////////////////////////////////

#include "AliTriggerDetector.h"

#include "AliITSLoader.h"
#include "AliITSgeom.h"
#include "AliITSdigitSPD.h"

#include "AliTriggerInput.h"


class AliITSTrigger : public AliTriggerDetector
{
 public:
                   AliITSTrigger();   // constructor
        virtual   ~AliITSTrigger(){}  // destructor
   virtual void    CreateInputs();
   virtual void    Trigger();

private:

   Int_t fFODigistThreshold;         // minimum number of digits to fire the FO trigger
   Int_t fHighMultFODigistThreshold; // minimum number of digits to fire the FO high mult trigger

   void MultiplicityTriggers(TObjArray* digDet, TTree* treeD, AliITSgeom* geom);
   void GeometryTriggers(TObjArray* digDet, TTree* treeD, AliITSgeom* geom);

  ClassDef( AliITSTrigger, 1 )  // ITS SPD Trigger Detector class
};

#endif
