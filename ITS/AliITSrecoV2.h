#ifndef ALIITSRECO_H
#define ALIITSRECO_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

//-------------------------------------------------------------------------
//                   ITS reconstruction name space
//
//       Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch 
//-------------------------------------------------------------------------
#include <Rtypes.h>
#include <iostream.h>

//namespace AliITSreco {    
   const Int_t kMaxClusterPerLayer=3500*10;
   const Int_t kMaxDetectorPerLayer=1000;
   const Int_t kLayersToSkip=0;

   const Int_t kMaxLayer=6;
   const Double_t kSigmaY2[kMaxLayer]={
      1.44e-6, 1.44e-6, 1.444e-5, 1.444e-5, 4.0e-6, 4.0e-6 
   };
   const Double_t kSigmaZ2[kMaxLayer]={
     //4.9e-5, 4.9e-5, 7.84e-6, 7.84e-6, 0.006889, 0.006889
     1.44e-4, 1.44e-4, 7.84e-6, 7.84e-6, 0.006889, 0.006889
   };

   const Double_t kChi2PerCluster=7.;//10.;//7
   const Double_t kMaxChi2=15.;//20.; //15.
   const Double_t kMaxRoad=13.;

   const Double_t kSigmaYV=0.005e+0;
   const Double_t kSigmaZV=0.010e+0;

   const Double_t kConvConst=100/0.299792458/0.2; 
//}

//using namespace AliITSreco;

#endif
