#ifndef AliRICHRawCluster_h
#define AliRICHRawCluster_h


/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

#include <TObject.h>


class AliRICHRawCluster : public TObject 
{
public:    
           AliRICHRawCluster();
  virtual ~AliRICHRawCluster(){;}
  Bool_t IsSortable() const {return kTRUE;}  //virtual
  Int_t  Compare(const TObject *obj) const;  //virtual
  void   Print(const Option_t *option="")const; //virtual  
  Int_t  PhysicsContribution();
public:
  Int_t       fTracks[3];      //labels of overlapped tracks
  Int_t       fQ  ;            // Q of cluster (in ADC counts)     
  Float_t     fX  ;            // X of cluster
  Float_t     fY  ;            // Y of cluster
  Int_t       fPeakSignal;     // Charge in the peak
  Int_t       fIndexMap[50];   //indeces of digits
  Float_t     fContMap[50];    //Contribution from digit
  Int_t       fPhysicsMap[50]; // physics processes
  Int_t       fMultiplicity;   //cluster multiplicity
  Int_t       fNcluster[2];    //number of clusters
  Int_t       fClusterType;    //??
  Int_t       fCtype;          //CL0, CL1, etc...    
  ClassDef(AliRICHRawCluster,2)  //Cluster object for set:RICH
};
#endif
