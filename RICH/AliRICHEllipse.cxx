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

/*
  $Log$
*/

#include "AliRICHEllipse.h"
#include "AliRICH.h"
#include "AliRun.h"

#include <TRandom.h>

ClassImp(AliRICHEllipse)

//________________________________________________________________________________
AliRICHEllipse::AliRICHEllipse()
{ 

//  Default Constructor for a RICH ellipse

    fCx = 0;
    fCy = 0;
    fOmega = 0;
    fTheta = 0;
    fPhi = 0;
    fh= 0;
}

//________________________________________________________________________________
AliRICHEllipse::~AliRICHEllipse()
{ 

// Destructor

    fCx = 0;
    fCy = 0;
    fOmega = 0;
    fTheta = 0;
    fPhi = 0;
    fh= 0;
}


//________________________________________________________________________________
AliRICHEllipse::AliRICHEllipse(Float_t cx, Float_t cy, Float_t omega, Float_t theta, Float_t phi)
{ 

//  Constructor for a RICH ellipse

    fCx = cx;
    fCy = cy;
    fOmega = omega;
    fTheta = theta;
    fPhi = phi;
    fh=11.25;
}

//________________________________________________________________________________
void AliRICHEllipse::CreatePoints(Int_t chamber)
{

// Create points along the ellipse equation

  Int_t s1,s2;
  Float_t fiducial=fh*TMath::Tan(fOmega+fTheta), l=fh/TMath::Cos(fTheta), xtrial, y=0, c0, c1, c2;
  //TRandom *random=new TRandom();

  AliRICH *pRICH  = (AliRICH*)gAlice->GetModule("RICH");
  AliRICHChamber*       iChamber;
  
  iChamber = &(pRICH->Chamber(chamber));
  //cout<<"fiducial="<<fiducial<<endl;
  
  for(Float_t i=0;i<1000;i++)
    {
      
      Float_t counter=0;
      
      c0=0;c1=0;c2=0;
      while((c1*c1-4*c2*c0)<=0 && counter<1000)
	{
	  //Choose which side to go...
	  if(i>250 && i<750) s1=1; 
	  //if (gRandom->Rndm(1)>.5) s1=1;
	  else s1=-1;
	  //printf("s1:%d\n",s1);
	  //Trial a y
	  y=s1*i*gRandom->Rndm(Int_t(fiducial/50));
	  //printf("Fiducial %f  for omega:%f theta:%f phi:%f\n",fiducial,fOmega,fTheta,fPhi);
	  Float_t alfa1=fTheta;
	  Float_t theta1=fPhi;
	  Float_t omega1=fOmega;
	  
	  //Solve the eq for a trial x
	  c0=-TMath::Power(y*TMath::Cos(alfa1)*TMath::Cos(theta1),2)-TMath::Power(y*TMath::Sin(alfa1),2)+TMath::Power(l*TMath::Tan(omega1),2)+2*l*y*TMath::Cos(alfa1)*TMath::Sin(theta1)*TMath::Power(TMath::Tan(omega1),2)+TMath::Power(y*TMath::Cos(alfa1)*TMath::Sin(theta1)*TMath::Tan(omega1),2);
	  c1=2*y*TMath::Cos(alfa1)*TMath::Sin(alfa1)-2*y*TMath::Cos(alfa1)*TMath::Power(TMath::Cos(theta1),2)*TMath::Sin(alfa1)+2*l*TMath::Sin(alfa1)*TMath::Sin(theta1)*TMath::Power(TMath::Tan(omega1),2)+2*y*TMath::Cos(alfa1)*TMath::Sin(alfa1)*TMath::Power(TMath::Sin(theta1),2)*TMath::Power(TMath::Tan(omega1),2);
	  c2=-TMath::Power(TMath::Cos(alfa1),2)-TMath::Power(TMath::Cos(theta1)*TMath::Sin(alfa1),2)+TMath::Power(TMath::Sin(alfa1)*TMath::Sin(theta1)*TMath::Tan(omega1),2);
	  //cout<<"Trial: y="<<y<<"c0="<<c0<<" c1="<<c1<<" c2="<<c2<<endl;
	  //printf("Result:%f\n\n",c1*c1-4*c2*c0);
	  //i+=.01;
	  counter +=1;
	}
      
      if (counter>=1000)
	y=0; 

      //Choose which side to go...
      //if(gRandom->Rndm(1)>.5) s=1; 
      //else s=-1;
      if(i>500) s2=1;
      //if (gRandom->Rndm(1)>.5) s2=1;
      else s2=-1;
      xtrial=fCx+(-c1+s2*TMath::Sqrt(c1*c1-4*c2*c0))/(2*c2);
      //cout<<"x="<<xtrial<<" y="<<cy+y<<endl;
      //printf("Coordinates: %f %f\n",xtrial,fCy+y);

      Float_t vectorLoc[3]={xtrial,6.276,(fCy+y)};
      Float_t  vectorGlob[3];
      iChamber->LocaltoGlobal(vectorLoc,vectorGlob);
      SetPoint(i,vectorGlob[0],vectorGlob[1],vectorGlob[2]);
      //printf("Coordinates: %f %f %f\n",vectorGlob[0],vectorGlob[1],vectorGlob[2]);
    }
}














