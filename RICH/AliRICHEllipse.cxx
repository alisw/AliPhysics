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
  Revision 1.4  2001/05/10 12:34:43  jbarbosa
  Changed drwaing routines.

  Revision 1.3  2000/11/01 15:37:44  jbarbosa
  Removed verbose output.

  Revision 1.2  2000/06/30 16:31:51  dibari
  New drawing routine from Nico and Daniela.

  Revision 1.1  2000/06/12 15:21:57  jbarbosa
  Cleaned up version.

*/

#include "TMath.h"
#include "AliRICHEllipse.h"
#include "AliRICH.h"
#include "AliRun.h"
#include "AliRICHPatRec.h"
#include "AliRICHGeometry.h"

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
AliRICHEllipse::AliRICHEllipse(Float_t cx, Float_t cy, Float_t omega, Float_t theta, Float_t phi, Float_t emiss)
{ 

//  Constructor for a RICH ellipse

    fCx = cx;
    fCy = cy;
    fOmega = omega;
    fTheta = theta;
    fPhi = phi;
    fEmissPoint = emiss;
    fh=9.25;
}

//________________________________________________________________________________
void AliRICHEllipse::CreatePoints(Int_t chamber)
{

// Create points along the ellipse equation
  
  Float_t x, y, rotx, roty, h, cx, cy, phi, omega, theta, omega1, theta1, phiinc;
  Float_t a,b,c,r,e, offset;
  
  Float_t kPi=TMath::Pi();

  AliRICH *pRICH  = (AliRICH*)gAlice->GetModule("RICH");
  AliRICHChamber*       iChamber;
  AliRICHGeometry*  geometry;
  
  iChamber = &(pRICH->Chamber(chamber));
  geometry=iChamber->GetGeometryModel();
  
  //h = 2.3 * geometry->GetRadiatorToPads();
  h = geometry->GetRadiatorToPads();
  //printf("h: %f",h);

  cx = fCx;
  cy = fCy;
  theta = fTheta;
  omega = fOmega;
  phiinc = fPhi+kPi/2;
  
  printf("Omega: %f, Theta: %f, Phi: %f\n", omega, theta, phiinc); 


  for(Float_t i=0;i<1000;i++)
    {
      phi=((Float_t)(i)/1000)*2*kPi;
      //printf("Phi: %f\n",phi);

      //theta1=SnellAngle(theta1);
      
      //if(phi<=TMath::Pi())
      omega1=SnellAngle(omega);
      theta1=SnellAngle(theta);
      //omega1=SnellAngle(omega+cos(phi)*theta);
      //if(phi>TMath::Pi())
	//omega1=SnellAngle(omega+(1-2*(2*TMath::Pi()-phi)/(TMath::Pi()))*theta);

   
      //Omega1->Fill(omega1,(float) 1);

      a = h*(tan(omega1+theta1)+tan(omega1-theta1))/2;
      b = h*tan(omega1);
      e = sqrt(1 - (b*b)/(a*a));
      c = a*e;
      r = (a*(1-e*e))/(1+e*cos(e));
      offset = h*(tan(omega1+theta1)-tan(omega1-theta1))/2;
	
      x = b* sin(phi);
      y = a* cos(phi) + offset;
		
      rotx = x*cos(phiinc)-y*sin(phiinc);
      roty = x*sin(phiinc)+y*cos(phiinc);
   
      //x = h * 1/(tan(omega1)) * sin(phi+phiinc);
      //y = x * 1/(tan(phi+phiinc));

      

      //Rings->Fill(x,y,(float) 1);

      rotx += cx;
      roty += cy;

      //printf("x:%f, y: %f\n",x,y);

      Float_t vectorLoc[3]={rotx,6.276,roty};
      Float_t  vectorGlob[3];
      iChamber->LocaltoGlobal(vectorLoc,vectorGlob);
      SetPoint((Int_t) i,vectorGlob[0],vectorGlob[1],vectorGlob[2]);
      //printf("Coordinates: %f %f %f\n",vectorGlob[0],vectorGlob[1],vectorGlob[2]);
      
    }

}

void AliRICHEllipse::CerenkovRingDrawing(Int_t chamber,Int_t track)
{

//to draw Cherenkov ring by known Cherenkov angle

  Int_t nmaxdegrees;
  Int_t Nphpad;
  Float_t phpad;
  Float_t nfreonave, nquartzave;
  Float_t aveEnerg;
  Float_t energy[2];
  Float_t e1, e2, f1, f2;
  Float_t pointsOnCathode[3];

  //printf("Drawing ring in chamber:%d\n",chamber);


  AliRICH *pRICH  = (AliRICH*)gAlice->GetModule("RICH");
  AliRICHChamber*       iChamber;
  
  iChamber = &(pRICH->Chamber(chamber));

  AliRICHPatRec *PatRec = new AliRICHPatRec;
  PatRec->TrackParam(track,chamber,fTheta,fOmega);

  //printf("Just created PateRec\n");

//parameters to calculate freon window refractive index vs. energy

    Float_t a = 1.177;
    Float_t b = 0.0172;
    
//parameters to calculate quartz window refractive index vs. energy
/*
   Energ[0]  = 5.6;
   Energ[1]  = 7.7;
*/	
    energy[0]  = 5.0;
    energy[1]  = 8.0;
    e1  = 10.666;
    e2  = 18.125;
    f1  = 46.411;
    f2  = 228.71;
  

    /*Float_t nquartz = 1.585;
      Float_t ngas    = 1.;
      Float_t nfreon  = 1.295;
      Float_t value;
    */



   nmaxdegrees = 360;
   
   for (Nphpad=0; Nphpad<nmaxdegrees;Nphpad++) { 

       phpad = (360./(Float_t)nmaxdegrees)*(Float_t)Nphpad;
      
       aveEnerg =  (energy[0]+energy[1])/2.;
       //aveEnerg = 6.5;
       
       
       nfreonave  = a+b*aveEnerg;
       nquartzave = sqrt(1+(f1/(TMath::Power(e1,2)-TMath::Power(aveEnerg,2)))+
			 (f2/(TMath::Power(e2,2)-TMath::Power(aveEnerg,2))));

       //nfreonave = 1.295;
       //nquartzave = 1.585;
       
       ///printf("Calling DistancefromMip %f %f \n",fEmissPoint,fOmega);
       
       //Float_t dummy = 
	 PatRec->DistanceFromMip(nfreonave, nquartzave,fEmissPoint,fOmega, phpad, pointsOnCathode,fTheta,fPhi);

       //Float_t points[3];

       //points = pointsOnCathode;


       //printf(" values %f %f %f\n",points[0],points[1],points[2]);
       
       Float_t vectorLoc[3]={pointsOnCathode[0],pointsOnCathode[2],pointsOnCathode[1]};
       Float_t  vectorGlob[3];
       iChamber->LocaltoGlobal(vectorLoc,vectorGlob);
       SetPoint(Nphpad,vectorGlob[0],vectorGlob[1],vectorGlob[2]);
      //fCoordEllipse[0][Nphpad] = pointsOnCathode[0];
      //fCoordEllipse[1][Nphpad] = pointsOnCathode[1];
       
       //printf(" values %f %f \n",pointsOnCathode[0],pointsOnCathode[1]);
       
   }

}



Float_t AliRICHEllipse:: SnellAngle(Float_t iangle)
{ 

// Compute the Snell angle

  Float_t nfreon  = 1.295;
  Float_t nquartz = 1.585;
  Float_t ngas = 1;

  Float_t sinrangle;
  Float_t rangle;
  Float_t a1, a2;

  a1=nfreon/nquartz;
  a2=nquartz/ngas;

  sinrangle = a1*a2*sin(iangle);

  if(sinrangle>1.) {
     rangle = 999.;
     return rangle;
  }
  
  rangle = asin(sinrangle);  
  //printf("iangle %f, a1*a2, %f, sinranlge, %f, rangle, %f\n", iangle, a1*a2, sinrangle, rangle);
  return rangle;

}








