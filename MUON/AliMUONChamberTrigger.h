#ifndef ALIMUONCHAMBERTRIGGER_H
#define ALIMUONCHAMBERTRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

#include "AliMUONChamber.h"


class AliMUONClusterFinder;
class AliMUONSegmentationTrigger ;
class AliMUONResponseTrigger ;

class AliMUONChamberTrigger:
public AliMUONChamber {
 public:
    AliMUONChamberTrigger();
    AliMUONChamberTrigger(Int_t id);
    virtual ~AliMUONChamberTrigger(){}
// Cluster formation method (charge disintegration)
    
    virtual void   DisIntegration(Float_t eloss, Float_t tof, Float_t xhit, Float_t yhit, Float_t zhit,
       				Int_t& nnew, Float_t newclust[6][500]);

  ClassDef(AliMUONChamberTrigger,1) // Muon trigger chamber class
      };
#endif












