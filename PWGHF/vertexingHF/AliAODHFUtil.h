#ifndef ALIAODHFUTIL_H
#define ALIAODHFUTIL_H

/* Copyright(c) 1998-2011, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */ 

///***********************************************************
/// \class Class AliAODHFUtil
/// \brief class for enabling access to data not available for the moment in AODs
/// Author: Carlos Perez
///***********************************************************
#include "TNamed.h"

class AliAODHFUtil : public TNamed {

 public:

 AliAODHFUtil();
 AliAODHFUtil(const char *pName);
 AliAODHFUtil(const AliAODHFUtil& pCopy);
 virtual ~AliAODHFUtil();

 //VZERO
 void SetVZERO(Float_t *pVzero);
 Float_t GetVZEROChannel(Int_t pCh) const;

 private:

 AliAODHFUtil& operator=(const AliAODHFUtil& pAssign);

 Float_t fVZERO[64]; ///< VZERO signals

 /// \cond CLASSIMP
 ClassDef(AliAODHFUtil,1); // class to keep event info need in HF analysis
 /// \endcond

};

#endif
