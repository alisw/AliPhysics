#ifndef ALIMUONCHAMBERTRIGGER_H
#define ALIMUONCHAMBERTRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004

#include "AliMUONChamber.h"


class AliMUONClusterFinder;
class AliMUONSegmentationTrigger;
class AliMUONResponseTrigger;
class AliMUONResponseTriggerV1;
class AliMUONHit;

class AliMUONChamberTrigger : public AliMUONChamber 
{
  public:
    AliMUONChamberTrigger();
    AliMUONChamberTrigger(Int_t id);
    virtual ~AliMUONChamberTrigger(){}
    
    // Cluster formation method (charge disintegration)
    
    virtual void   DisIntegration(AliMUONHit* hit,
       				Int_t& nnew, Float_t newclust[6][500]);

  ClassDef(AliMUONChamberTrigger,1) // Muon trigger chamber class
};
#endif












