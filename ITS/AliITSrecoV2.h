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

//namespace AliITSreco {    
   const Int_t kMaxClusterPerLayer=7000*10;
   const Int_t kMaxDetectorPerLayer=1000;

   const Int_t kLayersNotToSkip[]={0,0,0,0,0,0};
   const Int_t kLastLayerToTrackTo=0;

   const Int_t kMaxLayer=6;
   const Double_t kSigmaY2[kMaxLayer]={
      1.44e-6, 1.44e-6, 1.444e-5, 1.444e-5, 4.0e-6, 4.0e-6 
   };
   const Double_t kSigmaZ2[kMaxLayer]={
     //4.9e-5, 4.9e-5, 7.84e-6, 7.84e-6, 0.006889, 0.006889
     1.44e-4, 1.44e-4, 7.84e-6, 7.84e-6, 0.006889, 0.006889
   };

   const Double_t kChi2PerCluster=7.;
   const Double_t kMaxChi2=25.;
   const Double_t kMaxRoad=6.0;

   const Double_t kXV=0.0e+0;
   const Double_t kYV=0.0e+0;
   const Double_t kZV=0.0e+0;
   const Double_t kSigmaXV=0.005e+0;
   const Double_t kSigmaYV=0.005e+0;
   const Double_t kSigmaZV=0.010e+0;
//}

//using namespace AliITSreco;

#endif
