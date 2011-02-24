#ifndef AliAODHFUtil_H
#define AliAODHFUtil_H

/* Copyright(c) 1998-2011, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

//***********************************************************
// Class AliAODHFUtil
// class for enabling access to data not available for the moment in AODs
// Author: Carlos Perez
//***********************************************************
#include "TNamed.h"

class AliAODHFUtil : public TNamed {

 public:

 AliAODHFUtil();
 AliAODHFUtil(const char *pName);
 AliAODHFUtil(const AliAODHFUtil& pCopy);
 AliAODHFUtil& operator=(const AliAODHFUtil& pAssign);
 virtual ~AliAODHFUtil();

 //VZERO
 void SetVZERO(Float_t *pVzero);
 Float_t GetVZEROChannel(Int_t pCh) const;

 private:
 Float_t fVZERO[64]; // VZERO signals

 ClassDef(AliAODHFUtil,1) // class to keep event info need in HF analysis

};

#endif
