#ifndef ALITPCRECO_H
#define ALITPCRECO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                   TPC reconstruction name space
//
//       Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------
#include <Rtypes.h>

//namespace AliTPCreco {    
   const Int_t kMaxClusterPerRow=2500;
   const Int_t kRowsToSkip=10;
   const Int_t kMaxRow=159;

   const Double_t kMaxCHI2=12.;
   const Double_t kMaxROAD=30.;
//}

//using namespace AliTPCreco;

#endif
