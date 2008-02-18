#ifndef ALITPCCLUSTERKR_H
#define ALITPCCLUSTERKR_H
/* Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 * See cxx source for full Copyright notice                               */

/* $Id: AliTPCclusterKr.h,v 1.8 2008/01/22 16:07:15 matyja Exp $ */

//-------------------------------------------------------
//                    TPC Kr Cluster Class
//
//   Origin: Adam Matyja, INP PAN, adam.matyja@ifj.edu.pl
//-------------------------------------------------------

#include <vector>
#include "TObject.h"
#include "AliTPCvtpr.h"

//_____________________________________________________________________________
class AliTPCclusterKr: public TObject{
public:
  AliTPCclusterKr();
  AliTPCclusterKr(const AliTPCclusterKr & param);//copy constructor
  AliTPCclusterKr &operator = (const AliTPCclusterKr & param); 
  ~AliTPCclusterKr();


  void SetMax(AliTPCvtpr q){fMax=q;}//set values of max. in cluster
  void SetADCcluster(Int_t q){fADCcluster=q;}
  void SetSec(Short_t q){fNsec=q;}
  void SetNpads(Short_t q){fNpads=q;}
  void SetSize(){fSize=fCluster.size();}

  AliTPCvtpr GetMax(){return fMax;}
  Int_t GetADCcluster(){return  fADCcluster;}
  Short_t GetSec(){return fNsec;}
  Short_t GetNpads(){return fNpads;}
  Short_t GetSize(){return fSize;}

  std::vector< AliTPCvtpr*> fCluster;//cluster contents(adc,nt,np,nr)

private:
  AliTPCvtpr fMax;//max (ADC,timebin,pad,row) in cluster
  Int_t fADCcluster; //ADC of cluster
  Short_t fNsec;  //sector of the cluster
  Short_t fNpads; //number of pads in cluster
  Short_t fSize; //size of vector

  ClassDef(AliTPCclusterKr,2)  // Time Projection Chamber Kr clusters
};


#endif


