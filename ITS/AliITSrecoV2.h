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
   const Int_t kMaxClusterPerLayer5=7000*10*2/5;
   const Int_t kMaxClusterPerLayer10=7000*10*2/10;
   const Int_t kMaxClusterPerLayer20=7000*10*2/20;
   const Int_t kMaxDetectorPerLayer=1000;

   const Int_t kLayersNotToSkip[]={0,0,0,0,0,0};
   const Int_t kLastLayerToTrackTo=0;

   const Int_t kMaxLayer  = 6;
const Double_t kMaxSnp = 3.;
   const Double_t kSigmaY2[kMaxLayer]={
      1.44e-6, 1.44e-6, 1.444e-5, 1.444e-5, 4.0e-6, 4.0e-6 
   };
   const Double_t kSigmaZ2[kMaxLayer]={
     //4.9e-5, 4.9e-5, 7.84e-6, 7.84e-6, 0.006889, 0.006889
     1.44e-4, 1.44e-4, 7.84e-6, 7.84e-6, 0.006889, 0.006889
   };

   const Double_t kChi2PerCluster=9.;
//   const Double_t kMaxChi2PerCluster[5]={7.,5.,8.,8.,6.5};
   const Double_t kMaxChi2PerCluster[5]={11,12,12,5,12};
   const Double_t kMaxNormChi2NonC[6]  = {7,8,8,11,14,25};  //max norm chi2 for non constrained tracks
   const Double_t kMaxNormChi2C[6]  = {11,13,15,18,30,35};     //max norm chi2 for constrained tracks

   const Double_t kMaxChi2=35.;
//   const Double_t kMaxChi2s[6]={40,40,40,40,40,40};   
   const Double_t kMaxChi2s[6]={25,25,25,25,40,50};   
   const Double_t kMaxChi2sR[6]={10,10,10,10,30,40};   
   const Double_t kMaxChi2In=16.;
   const Double_t kVertexCut=25;
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
