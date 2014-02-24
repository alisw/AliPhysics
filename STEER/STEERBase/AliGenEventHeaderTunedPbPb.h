#ifndef ALIGENEVENTHEADERTUNEDPBPB_H
#define ALIGENEVENTHEADERTUNEDPBPB_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliGenEventHeaderTunedPbPb.h 50615 2013-08-19 23:19:19Z fnoferin $ */

//---------------------------------------------------------------------
// Event header base class for generator. 
// Stores generated event information
// Author: andreas.morsch@cern.ch
//---------------------------------------------------------------------

#include "AliGenEventHeader.h"

class AliGenEventHeaderTunedPbPb : public AliGenEventHeader
{
 public:

  AliGenEventHeaderTunedPbPb(const char* name);
  AliGenEventHeaderTunedPbPb();
  virtual ~AliGenEventHeaderTunedPbPb() {}
  // Getters
  Float_t GetCentrality() const {return fCentrality;};
  Float_t GetPsi2() const {return fPsi2;};
  Float_t GetPsi3() const {return fPsi3;};
  Float_t GetPsi4() const {return fPsi4;};
  // Setters
  void SetCentrality(Float_t centrality){fCentrality = centrality;};
  void SetPsi2(Float_t psi){fPsi2 = psi;};
  void SetPsi3(Float_t psi){fPsi3 = psi;};
  void SetPsi4(Float_t psi){fPsi4 = psi;};
	  
protected:
  Float_t fCentrality;// centrality
  Float_t fPsi2;     // psi_2 EP
  Float_t fPsi3;     // psi_3 EP
  Float_t fPsi4;     // psi_4 EP
  ClassDef(AliGenEventHeaderTunedPbPb, 2)        // Event header for primary event
};

#endif

