#ifndef ALIMUONRESPONSETRIGGERV1_H
#define ALIMUONRESPONSETRIGGERV1_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                                */

/* $Id$ */
// Revision of includes 07/05/2004

/// \ingroup sim
/// \class AliMUONResponseTriggerV1
/// \brief Implementation of RPC response

#include "AliMUONResponseTrigger.h"

class AliMUONResponseTriggerV1 : public AliMUONResponseTrigger 
{
  public:
    // default constructor
    AliMUONResponseTriggerV1();
    AliMUONResponseTriggerV1(Float_t hv);
    virtual ~AliMUONResponseTriggerV1(){} 
    // Charge disintegration
    virtual Float_t  IntXY(Int_t idDE, AliMUONGeometrySegmentation * segmentation) const;

    // Set the GenerCluster parameter       
    virtual Int_t SetGenerCluster();

  protected:
    Float_t fGenerCluster;   // Random number  
    Float_t fA;              // first parameter  of the cluster-size param
    Float_t fB;              // second parameter of the cluster-size param
    Float_t fC;              // third parameter  of the cluster-size param

  private:
    // initialize parameters
    void SetParameters(Float_t hv);
    // parametrization of the cluster-size
    Float_t FireStripProb(Float_t x4, Float_t theta) const;

  ClassDef(AliMUONResponseTriggerV1,1) // Implementation of RPC response
    
};
#endif













