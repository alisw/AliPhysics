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

//-------------------------------------------------------------------------
//                Implementation of the HLT ITS track class
//
//          Origin: Cvetan Cheshkov, CERN, Cvetan.Cheshkov@cern.ch
//-------------------------------------------------------------------------

#include <TMath.h>

#include "AliL3StandardIncludes.h"

#include "AliESDHLTtrack.h"
#include "AliL3ITStrack.h"

ClassImp(AliL3ITStrack)

//____________________________________________________________________________
AliL3ITStrack::AliL3ITStrack():AliITStrackV2(),
  fESDHLTtrack(0)
{
  //------------------------------------------------------------------
  //Constructor
  //------------------------------------------------------------------
}

//____________________________________________________________________________
AliL3ITStrack::AliL3ITStrack(const AliL3ITStrack& t) : AliITStrackV2(t) {
  //------------------------------------------------------------------
  //Copy constructor
  //------------------------------------------------------------------
  fESDHLTtrack=t.fESDHLTtrack;
}

//____________________________________________________________________________
AliL3ITStrack::AliL3ITStrack(AliESDHLTtrack& t, Double_t zvertex) throw (const Char_t *) {
  // The method constructs an AliL3ITStrack object from an ESD HLT track
 
  SetChi2(0.);
  if(t.GetNHits()==1)
    SetNumberOfClusters(0);
  else
    SetNumberOfClusters(t.GetNHits());
  SetLabel(t.GetMCid());
  SetMass(0.13957);

  fdEdx=0;
  fAlpha = fmod((t.GetSector()+0.5)*(2*TMath::Pi()/18),2*TMath::Pi());
  if      (fAlpha < -TMath::Pi()) fAlpha += 2*TMath::Pi();
  else if (fAlpha >= TMath::Pi()) fAlpha -= 2*TMath::Pi();

  //First the emiision angle
  Double_t psi = t.GetPsi()-(t.GetSector()+0.5)*(2*TMath::Pi()/18);

  //Then local x,y coordinates
  Double_t radius = t.GetPt()*GetConvConst();
  Double_t xhit = 82.97; //Position at first TPC padrow
  Double_t trackphi0 = psi + (-t.GetCharge())*TMath::Pi()/2;
  Double_t x0 = t.GetFirstPointX()*TMath::Cos(fAlpha) + t.GetFirstPointY()*TMath::Sin(fAlpha);
  Double_t y0 = t.GetFirstPointY()*TMath::Cos(fAlpha) - t.GetFirstPointX()*TMath::Sin(fAlpha);
  Double_t centerx = radius *  cos(trackphi0) + x0;
  Double_t centery = radius *  sin(trackphi0) + y0;
  Double_t aa = (xhit - centerx)*(xhit - centerx);
  Double_t r2 = radius*radius;
  if(aa > r2) throw "AliITStrackV2: conversion failed !\n";
  Double_t aa2 = sqrt(r2 - aa);
  Double_t y1 = centery + aa2;
  Double_t y2 = centery - aa2;
  Double_t yhit = y1;
  if(fabs(y2) < fabs(y1)) yhit = y2;

  //Local z coordinate
  Double_t angle1 = atan2((yhit - centery),(xhit - centerx));
  if(angle1 < 0) angle1 += 2.*TMath::Pi();
  Double_t angle2 = atan2((x0-centery),(y0-centerx));
  if(angle2 < 0) angle2 += 2.*TMath::Pi();
  Double_t diffangle = angle1 - angle2;
  diffangle = fmod(diffangle,2.*TMath::Pi());
  if(((-t.GetCharge())*diffangle) < 0) diffangle = diffangle - (-t.GetCharge())*2.*TMath::Pi();
  Double_t stot = fabs(diffangle)*radius;
  Double_t zhit;
  if(t.GetNHits()==1)
    zhit = zvertex + stot*t.GetTgl();
  else
    zhit = t.GetFirstPointZ() + stot*t.GetTgl();

  //Local sine of track azimuthal angle
  if((-t.GetCharge())<0)
    radius = -radius;
  Double_t sinbeta = -1.*(centerx - xhit)/radius;

  //Filling of the track paramaters
  fX=xhit;
  fP0=yhit;
  fP1=zhit;
  fP2=sinbeta;
  fP3=t.GetTgl();
  fP4=1./radius;

  //and covariance matrix
  fC22=0.005*0.005;
  fC33=0.005*0.005;
  fC00=fC22*82.97*82.97;
  fC11=fC33*82.97*82.97;
  //  fC44=(0.005+0.025*t.GetPt())*(0.005+0.025*t.GetPt())*fP4*fP4;
  //  fC44=(0.01+0.01*t.GetPt())*(0.01+0.01*t.GetPt())*fP4*fP4;
  fC44=0.01*0.01*fP4*fP4;

  fC10=0;
  fC20=0;   fC21=0;
  fC30=0;   fC31=0;   fC32=0;
  fC40=0;   fC41=0;   fC42=0;   fC43=0;

  fESDHLTtrack=&t;
  fESDtrack=0;

  SetFakeRatio(0.);

}

//_____________________________________________________________________________
Int_t AliL3ITStrack::Compare(const TObject *o) const {
  //-----------------------------------------------------------------
  // This function compares tracks according to the their curvature
  //-----------------------------------------------------------------
  AliL3ITStrack *t=(AliL3ITStrack*)o;
  Double_t co=TMath::Abs(t->Get1Pt());
  Double_t c =TMath::Abs(Get1Pt());
  //  Double_t co=t->GetSigmaY2()*t->GetSigmaZ2();
  //  Double_t c =GetSigmaY2()*GetSigmaZ2();
  if (c>co) return 1;
  else if (c<co) return -1;
  return 0;
}
