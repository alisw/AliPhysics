#ifndef ALIGENEVENTHEADER_H
#define ALIGENEVENTHEADER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

//---------------------------------------------------------------------
// Event header base class for generator. 
// Stores generated event information
// Author: andreas.morsch@cern.ch
//---------------------------------------------------------------------

#include <TNamed.h>
#include <TArrayF.h>

class AliGenEventHeader : public TNamed
{
 public:

  AliGenEventHeader(const char* name);
  AliGenEventHeader();
  virtual ~AliGenEventHeader() {}
  // Getters
  virtual Int_t           NProduced()  const   {return fNProduced;}
  virtual void            PrimaryVertex(TArrayF &o) const;
  
  // Setters
  virtual void   SetNProduced(Int_t nprod)         {fNProduced=nprod;}
  virtual void   SetPrimaryVertex(const TArrayF &o);
  
protected:
  Int_t     fNProduced;                 // Number stable or undecayed particles
  TArrayF   fVertex;                    // Primary Vertex Position
  ClassDef(AliGenEventHeader,2)         // Event header for primary event
};

#endif
