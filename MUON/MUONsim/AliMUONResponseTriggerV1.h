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
#include "AliMUONHit.h"
#include "TArrayF.h"

class AliMUONResponseTriggerV1 : public AliMUONResponseTrigger 
{
  public:
    // default constructor
    AliMUONResponseTriggerV1();
    AliMUONResponseTriggerV1(Int_t mode);
    virtual ~AliMUONResponseTriggerV1();
    
    // Set the GenerCluster parameter       
    virtual Int_t SetGenerCluster();
    virtual void DisIntegrate(const AliMUONHit& hit, TList& digits, Float_t timeDif);
    
  protected:
    Float_t fGenerCluster;   ///< Random number
    TArrayF fHVvalues;       ///< Array containing HV values
    TArrayF fBValues;        ///< Array containing b parameters
    Int_t fWorkCondition;    ///< 1=streamer - 2=avalanche

  private:
    // parametrization of the cluster-size
    void SetHV();
    void SetBValues();
    Float_t FireStripProb(Float_t x4, Float_t theta,Int_t rpc,Int_t plane,Int_t cath) const;
    void Neighbours(const Int_t cath, const Int_t iX, const Int_t iY, Int_t Xlist[30], Int_t Ylist[30]) const;
    
  ClassDef(AliMUONResponseTriggerV1,2) // Implementation of RPC response
};
#endif
