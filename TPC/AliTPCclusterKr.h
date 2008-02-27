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
#include <list>
#include "TObject.h"
#include "AliTPCvtpr.h"

//_____________________________________________________________________________
class AliTPCclusterKr: public TObject{
public:
  AliTPCclusterKr();
  AliTPCclusterKr(const AliTPCclusterKr & param);//copy constructor
  AliTPCclusterKr &operator = (const AliTPCclusterKr & param); 
  ~AliTPCclusterKr();

  void SetCenter();//set center of the cluster weighted by charge

  void SetMax(AliTPCvtpr q){fMax=q;}//set values of max. in cluster
  void SetADCcluster(Int_t q){fADCcluster=q;}
  void SetSec(Short_t q){fSec=q;}
  void SetNPads(Short_t q){fNPads=q;}
  void SetNRows(Short_t q){fNRows=q;}
  void SetSize(){fSize=fCluster.size();}
  void SetCenterX(Double_t q){fCenterX=q;}
  void SetCenterY(Double_t q){fCenterY=q;}
  void SetCenterT(Double_t q){fCenterT=q;}

  AliTPCvtpr GetMax(){return fMax;}
  Int_t GetADCcluster(){return  fADCcluster;}
  Short_t GetSec(){return fSec;}
  Short_t GetNPads(){return fNPads;}
  Short_t GetNRows(){return fNRows;}
  Short_t GetSize(){return fSize;}
  Double_t GetCenterX(){return fCenterX;}
  Double_t GetCenterY(){return fCenterY;}
  Double_t GetCenterT(){return fCenterT;}

  std::vector< AliTPCvtpr*> fCluster;//cluster contents(adc,nt,np,nr)

private:
  AliTPCvtpr fMax;//max (ADC,timebin,pad,row) in cluster
  Int_t fADCcluster; //ADC of cluster
  Short_t fSec;  //sector of the cluster
  Short_t fNPads; //number of pads in cluster
  Short_t fNRows; //number of rows in cluster
  Short_t fSize; //size of vector
  Double_t fCenterX;// X coordinate of the cluster center in cm
  Double_t fCenterY;// Y coordinate of the cluster center in cm
  Double_t fCenterT;// time coordinate of the cluster center in timebins

  ClassDef(AliTPCclusterKr,2)  // Time Projection Chamber Kr clusters
};


#endif


