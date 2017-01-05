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


// Uncorrelated muon pairs generator 
// Author: alessandro.de.falco@ca.infn.it


#include "TPDGCode.h"

#include "AliConst.h"
#include "AliGenMuonUncorr.h"
#include "AliRun.h"
#include "AliGenEventHeader.h"
#include "TDatabasePDG.h"
#include "AliPDG.h"

ClassImp(AliGenMuonUncorr)

//_____________________________________________________________________________
AliGenMuonUncorr::AliGenMuonUncorr()
    :AliGenerator(), 
     fIpart(13),
  fY(0),
  fPt(0) 
{
  fY = new TF1("fY","gaus(0)+gaus(3)",-4.5,-2);
  fY->SetParameters(3.44458e+02,-3.84918e+00,5.92541e-01,4.57350e+02,-2.61815e+00,5.95914e-01);
  fPt = new TF1("fPt","[0]*x*exp(-x/[1])+[2]*exp(-x/[3])",0.5,5); 
  fPt->SetParameters(5.82192e+03,1.53286e-01,7.57423e+01,7.37488e-01); 
  fName  = "Uncorr";
  fTitle = "Uncorr particle generator";

}


//_____________________________________________________________________________

void AliGenMuonUncorr::Generate() {
    Float_t polar[3]= {0,0,0};
  //
    Float_t origin[3];
    Float_t time;
    Float_t p[3];
    Int_t i, j, nt;
    Double_t pmom, theta, phi, pt;
    Double_t y, mt;
    //
    Float_t random[6];
  //
    for (j=0;j<3;j++) origin[j]=fOrigin[j];
    time = fTimeOrigin;
    if(fVertexSmear==kPerEvent) {
	Vertex();
	for (j=0;j<3;j++) origin[j]=fVertex[j];
	time = fTime;
    }

    Double_t m = TDatabasePDG::Instance()->GetParticle(fIpart)->Mass();

    Int_t charge = 1; 

    TLorentzVector p4mu[2]; 
    for(i=0;i<2;i++) {
      if (gRandom->Rndm()>0.5) charge = 1;
      else charge = -1;
      
      y = fY->GetRandom();
      pt = fPt->GetRandom();
      phi = gRandom->Rndm()*2*TMath::Pi(); 
      mt = TMath::Sqrt(m*m + pt*pt); 
      
      p[0] = pt * TMath::Cos(phi); 
      p[1] = pt * TMath::Sin(phi);
      p[2] = mt * TMath::SinH(y); 
      
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
      PushTrack(fTrackIt,-1,charge*fIpart,p,origin,polar,time,kPPrimary,nt, 1., 1);
    }
    
    AliGenEventHeader* header = new AliGenEventHeader("UNCORR");
    header->SetPrimaryVertex(fVertex);
    header->SetNProduced(2);
    header->SetInteractionTime(fTime);
    
 // Passes header either to the container or to gAlice
    if (fContainer) {
        header->SetName(fName);
	fContainer->AddHeader(header);
    } else {
	gAlice->SetGenEventHeader(header);	
    }
}

//_____________________________________________________________________________

void AliGenMuonUncorr::Init()
{;}

