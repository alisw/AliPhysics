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

/* $Id$ */


// Generates n particles with in the same phi angle, varies theta
// in equidistant intervals
// This class is intended to use for studies of TPC response
// via merging with background event.
// Note that for a given theta pt and p are not independent 
// Range for only one variable (pt or p) should be given.
// Based on the AliGenBox class written by andreas.morsch@cern.ch
//
// Comments and suggestions: Jiri.Chudoba@cern.ch


#include <TPDGCode.h>

#include "AliConst.h"
#include "AliGenThetaSlice.h"
#include "AliRun.h"

ClassImp(AliGenThetaSlice)

//_____________________________________________________________________________
AliGenThetaSlice::AliGenThetaSlice()
    :AliGenerator(),
     fIpart(0)
{
  //
  // Default constructor
  //
}

//_____________________________________________________________________________
AliGenThetaSlice::AliGenThetaSlice(Int_t npart)
    :AliGenerator(npart),
     fIpart(kProton)
{
  //
  // Standard constructor
  //
  fName  = "ThetaSlice";
  fTitle = "Particle generator - const. phi, slices in theta";
}

//_____________________________________________________________________________

void AliGenThetaSlice::Generate()
{
  //
  // Generate one trigger
  //
  
    Float_t polar[3]= {0,0,0};
    Float_t origin[3];
    Float_t time;
    Float_t p[3];
    Int_t i, j, nt;
    Double_t pmom, theta, phi, pt;
    Float_t random[6];

    if (fNpart == 0) return;

    for (j=0;j<3;j++) origin[j]=fOrigin[j];
    time = fTimeOrigin;
    if(fVertexSmear==kPerEvent) {
	Rndm(random,6);
	for (j=0;j<3;j++) {
	    origin[j]+=fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
		TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
	}
	Rndm(random,2);
	time += fOsigma[2]/TMath::Ccgs()*
	  TMath::Cos(2*random[0]*TMath::Pi())*
	  TMath::Sqrt(-2*TMath::Log(random[1]));
    }
    Float_t thetaInterval = 0.;
    if (fNpart > 1) {
      thetaInterval = (fThetaMax-fThetaMin)/(fNpart-1);
    }
    Rndm(random,1);
    phi=fPhiMin+random[0]*(fPhiMax-fPhiMin);
    for(i=0;i<fNpart;i++) {
	Rndm(random,1);
	theta=fThetaMin+i*thetaInterval;
	if(TestBit(kMomentumRange)) {
	    pmom=fPMin+random[0]*(fPMax-fPMin);
	    pt=pmom*TMath::Sin(theta);
	} else {
	    pt=fPtMin+random[0]*(fPtMax-fPtMin);
	    pmom=pt/TMath::Sin(theta);
	}
	p[0] = pt*TMath::Cos(phi);
	p[1] = pt*TMath::Sin(phi);
	p[2] = pmom*TMath::Cos(theta);

	if(fVertexSmear==kPerTrack) {
	    Rndm(random,6);
	    for (j=0;j<3;j++) {
		origin[j]=fOrigin[j]+fOsigma[j]*TMath::Cos(2*random[2*j]*TMath::Pi())*
		    TMath::Sqrt(-2*TMath::Log(random[2*j+1]));
	    }
	    Rndm(random,2);
	    time = fTimeOrigin + fOsigma[2]/TMath::Ccgs()*
	      TMath::Cos(2*random[0]*TMath::Pi())*
	      TMath::Sqrt(-2*TMath::Log(random[1]));
	}
	PushTrack(fTrackIt,-1,fIpart,p,origin,polar,time,kPPrimary,nt);
    }
}

//_____________________________________________________________________________

void AliGenThetaSlice::Init()
{
// Initialisation, check consistency of selected ranges
  if(TestBit(kPtRange)&&TestBit(kMomentumRange)) 
    Fatal("Init","You should not set the momentum range and the pt range!\n");
  if((!TestBit(kPtRange))&&(!TestBit(kMomentumRange))) 
    Fatal("Init","You should set either the momentum or the pt range!\n");
}

