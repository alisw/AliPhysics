#ifndef ALIRICHRAWCLUSTER_H
#define ALIRICHRAWCLUSTER_H


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id$ */
#include <TObject.h>
#include <TMath.h>
#include <TArrayF.h>


class AliRICHRawCluster : public TObject {
public:
    Int_t       fTracks[3];      //labels of overlapped tracks
    Int_t       fQ  ;            // Q of cluster (in ADC counts)     
    Float_t     fX  ;            // X of cluster
    Float_t     fY  ;            // Y of cluster
    Int_t       fPeakSignal;     // Charge in the peak
    Int_t       fIndexMap[50];   //indeces of digits
    Int_t       fOffsetMap[50];  // offset map
    Float_t     fContMap[50];    //Contribution from digit
    Int_t       fPhysicsMap[50]; // physics processes
    Int_t       fMultiplicity;   //cluster multiplicity
    Int_t       fNcluster[2];    //number of clusters
    Int_t       fClusterType;    //??
    Int_t       fCtype;          //CL0, CL1, etc...
 public:
    AliRICHRawCluster() {
	fTracks[0]=fTracks[1]=fTracks[2]=-1; 
	fQ=0; fX=fY=0; fMultiplicity=0;
	for (int k=0;k<50;k++) {
	    fIndexMap[k]=-1;
	    fOffsetMap[k]=0;
	    fContMap[k]=0;
	    fPhysicsMap[k]=-1;
	    fCtype=-1;
	}
	fNcluster[0]=fNcluster[1]=-1;
    }
    virtual ~AliRICHRawCluster() {}
    
    Float_t GetRadius() {return TMath::Sqrt(fX*fX+fY*fY);}
    
    Bool_t IsSortable() const {return kTRUE;}
    Int_t  Compare(TObject *obj);
    Int_t PhysicsContribution();
    static Int_t BinarySearch(Float_t r, TArrayF coord, Int_t from, Int_t upto);
    static void  SortMin(Int_t *idx,Float_t *xdarray,Float_t *xarray,Float_t *yarray,Float_t *qarray,Int_t ntr);
    
   ClassDef(AliRICHRawCluster,1)  //Cluster object for set:RICH
};
#endif
