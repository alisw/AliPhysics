#ifndef ALIMUONRAWCLUSTER_H
#define ALIMUONRAWCLUSTER_H

/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */

class TArrayF;

#include <TObject.h>
#include <TMath.h> // because of inline funtion GetRadius

class AliMUONRawCluster : public TObject {
public:

   Int_t       fTracks[3];        //labels of overlapped tracks
   Int_t       fQ[2]  ;           // Q of cluster (in ADC counts)     
   Float_t     fX[2]  ;           // X of cluster
   Float_t     fY[2]  ;           // Y of cluster
   Float_t     fZ[2]  ;           // Z of cluster
   Int_t       fPeakSignal[2];    // Peak signal 
   Int_t       fIndexMap[50][2];  // indeces of digits
   Int_t       fOffsetMap[50][2]; // Emmanuel special
   Float_t     fContMap[50][2];   // Contribution from digit
   Int_t       fPhysicsMap[50];   // Distinguish signal and background contr.
   Int_t       fMultiplicity[2];  // Cluster multiplicity
   Int_t       fNcluster[2];      // Number of clusters
   Int_t       fClusterType;      // Cluster type
   Float_t     fChi2[2];          // Chi**2 of fit
   Int_t       fGhost;            // 0 if not a ghost or ghost problem solved
                                  // >0 if ghost problem remains because
                                  // 1 both (true and ghost) satify 
                                  //   charge chi2 compatibility
                                  // 2 none give satisfactory chi2
 public:
   AliMUONRawCluster();
   virtual ~AliMUONRawCluster() {}
   Float_t GetRadius(Int_t i) {return TMath::Sqrt(fX[i]*fX[i]+fY[i]*fY[i]);}
   Bool_t IsSortable() const {return kTRUE;}
   Int_t  Compare(const TObject *obj) const;
   Int_t PhysicsContribution();
   static Int_t BinarySearch(Float_t r, TArrayF ccord, Int_t from, Int_t upto);
   static void  SortMin(Int_t *idx,Float_t *xdarray, Float_t *xarray,
			Float_t *yarray, Float_t *qarray,Int_t ntr);
   void DumpIndex();

   ClassDef(AliMUONRawCluster,1)  //Cluster class for MUON
};
#endif






