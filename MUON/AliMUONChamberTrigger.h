#ifndef ALIMUONCHAMBERTRIGGER_H
#define ALIMUONCHAMBERTRIGGER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
// Revision of includes 07/05/2004
//
/// \ingroup base
/// \class AliMUONChamberTrigger
/// \brief Muon trigger chamber class

#include "AliMUONChamber.h"

class AliMUONClusterFinder;
class AliMUONSegmentationTrigger;
class AliMUONResponseTrigger;
class AliMUONResponseTriggerV1;
class AliMUONGeometryTransformer;
class AliMUONHit;

class AliMUONChamberTrigger : public AliMUONChamber 
{
  public:
    AliMUONChamberTrigger();
    AliMUONChamberTrigger(Int_t id, const AliMUONGeometryTransformer* kGeometry);
    virtual ~AliMUONChamberTrigger(){}
    
    // Cluster formation method (charge disintegration)
    
    virtual void   DisIntegration(AliMUONHit* hit,
       				 Int_t& nnew, Float_t newclust[6][500]);

  protected:   
    AliMUONChamberTrigger(const AliMUONChamberTrigger& right);
    AliMUONChamberTrigger&  operator = (const AliMUONChamberTrigger& right);

    const AliMUONGeometryTransformer* fkGeomTransformer;///< geometry transformations

  ClassDef(AliMUONChamberTrigger,2) // Muon trigger chamber class
};
#endif












