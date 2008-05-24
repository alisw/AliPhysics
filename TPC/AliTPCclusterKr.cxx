/**************************************************************************
 * Copyright(c) 1998-1999, ALICE Experiment at CERN, All rights reserved. *
 *                                                                        *
 * Author: The ALICE Off-line Project.                                    *
 * Contributors are mentioned in the code where appropriate.              *
 *                                                                        *
 * Permission to use, copy, modify and distribute this software and its   *
 * documentation strictly for non-commercial purposes is hereby granted   *
 * without fee, provided that the above copyright notice appears in all   *
 * copies and that both the copyright notice and this permission notice   *
 * appear in the supporting documentation. The authors make no claims     *
 * about the suitability of this software for any purpose. It is          *
 * provided "as is" without express or implied warranty.                  *
 **************************************************************************/

/* $Id: AliTPCclusterKr.cxx,v 1.7 2008/01/22 17:24:53 matyja Exp $ */

//-----------------------------------------------------------------
//           Implementation of the TPC Kr cluster class
//
// Origin: Adam Matyja, INP PAN, adam.matyja@ifj.edu.pl
//-----------------------------------------------------------------

#include "AliTPCclusterKr.h"
#include "AliCluster.h"
#include "AliTPCvtpr.h"
#include "TObjArray.h"

ClassImp(AliTPCclusterKr)


AliTPCclusterKr::AliTPCclusterKr()
:AliCluster(),
 fMax(),
 fADCcluster(0),
 fSec(0),
 fNPads(0),
 fNRows(0),
 fSize(0),
 fCenterX(0),
 fCenterY(0),
 fCenterT(0),
 fCluster(0)
{
//
// default constructor
//
  fCluster=new TObjArray();
}

AliTPCclusterKr::AliTPCclusterKr(const AliTPCclusterKr &param)
:AliCluster(param),
 fMax(),
 fADCcluster(0),
 fSec(0),
 fNPads(0),
 fNRows(0),
 fSize(0),
 fCenterX(0),
 fCenterY(0),
 fCenterT(0),
 fCluster(0)
{
//
// copy constructor
//
  fADCcluster = param.fADCcluster;
  fSec  = param.fSec ;
  fNPads = param.fNPads;
  fNRows = param.fNRows;
  fMax = param.fMax;
  //  fCluster = param.fCluster;
  fCenterX = param.fCenterX;
  fCenterY = param.fCenterY;
  fCenterT = param.fCenterT;
  fCluster=new TObjArray(*(param.fCluster));
  fSize = param.fSize;
} 

AliTPCclusterKr &AliTPCclusterKr::operator = (const AliTPCclusterKr & param)
{
  //
  // assignment operator
  // 
  (AliCluster&)(*this) = (AliCluster&)param;
  fADCcluster = param.fADCcluster;
  fSec  = param.fSec ;
  fNPads = param.fNPads;
  fNRows = param.fNRows;
  fMax = param.fMax;
  //  fCluster=param.fCluster;
  fCenterX = param.fCenterX;
  fCenterY = param.fCenterY;
  fCenterT = param.fCenterT;
  delete fCluster;
  fCluster=new TObjArray(*(param.fCluster));
  fSize=param.fSize;
  return (*this);
}

AliTPCclusterKr::~AliTPCclusterKr()
{
  //
  // destructor
  //
  if(fCluster) {
    fCluster->SetOwner(kTRUE);
    fCluster->Delete();
    delete fCluster;
  }
  fCluster=0;
}

////____________________________________________________________________________
void AliTPCclusterKr::SetCenter(){
  //
  // calculate geometrical center of the cluster
  //
  Double_t rX=0;
  Double_t rY=0;
  Double_t rT=0;

  Short_t adc;
  fADCcluster=0;
  for(Int_t iter = 0; iter < fCluster->GetEntriesFast(); ++iter) {
    AliTPCvtpr *iclus=(AliTPCvtpr *)fCluster->At(iter);

    //for( std::vector<AliTPCvtpr*>::iterator iclus  = fCluster.begin();
    //iclus != fCluster.end(); ++iclus ) {
    adc = (iclus)->GetAdc();
    fADCcluster+=adc;
    rX += ((iclus)->GetX() * adc);
    rY += ((iclus)->GetY() * adc);
    rT += ((iclus)->GetT() * adc);
  }
  fCenterX=rX/fADCcluster;
  fCenterY=rY/fADCcluster;
  fCenterT=rT/fADCcluster;

  return;
}
