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

#include <stdlib.h> 

#include "AliRICHDetectV1.h"
#include "AliRICH.h"
#include "AliRICHPoints.h"
#include "AliRICHDetect.h"
#include "AliRICHHit.h"
#include "AliRICHDigit.h"
#include "AliRICHRawCluster.h"
#include "AliRICHCerenkov.h"
#include "AliRICHSegmentationV0.h"
#include "AliRun.h"
#include "TParticle.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom.h"
#include "TH3.h"
#include "TH2.h"
#include "TCanvas.h"
#include <TStyle.h>



ClassImp(AliRICHDetectV1)


//___________________________________________
AliRICHDetectV1::AliRICHDetectV1() : AliRICHDetect()
{

// Default constructor 

  fc1 = 0;
  fc2 = 0;
  fc3 = 0;

}

//___________________________________________
AliRICHDetectV1::AliRICHDetectV1(const char *name, const char *title)
    : AliRICHDetect()
{

  TStyle *mystyle=new TStyle("Plain","mystyle");
  mystyle->SetPalette(1,0);
  mystyle->cd();


  fc1= new TCanvas("c1","Reconstructed points",50,50,300,350);
  fc1->Divide(2,2);
  fc2= new TCanvas("c2","Reconstructed points after SPOT",370,50,300,350);
  fc2->Divide(2,2); 
  fc3= new TCanvas("c3","Used Digits",690,50,300,350);
  fc4= new TCanvas("c4","Mesh activation data",50,430,600,350);
  fc4->Divide(2,1);
 

}

//___________________________________________
AliRICHDetectV1::~AliRICHDetectV1()
{
    
// Destructor

}


void AliRICHDetectV1::Detect(Int_t nev, Int_t type)
{  	
    
//
// Detection algorithm


  //printf("Detection started!\n");
  Float_t omega,omega1,theta1,x,y,q=0,z,cx,cy,max,radius=0,meanradius=0;
  Int_t maxi,maxj,maxk;
  Float_t originalOmega, originalPhi, originalTheta;
  Float_t steptheta,stepphi,stepomega;
  Float_t binomega, bintheta, binphi;
  Int_t intomega, inttheta, intphi;
  Float_t maxRadius,minRadius,eccentricity,angularadius,offset,phi_relative;
  Int_t i,j,k;

  AliRICH *pRICH  = (AliRICH*)gAlice->GetDetector("RICH");
  AliRICHSegmentationV0*  segmentation;
  AliRICHChamber*       iChamber;
  AliRICHGeometry*  geometry;
  
  iChamber = &(pRICH->Chamber(0));
  segmentation=(AliRICHSegmentationV0*) iChamber->GetSegmentationModel(0);
  geometry=iChamber->GetGeometryModel();
 
  
  //const Float_t Noise_Level=0;                       //Noise Level in percentage of mesh points
  //const Float_t t=0.6;			       //Softening of Noise Correction (factor)
  
  const Float_t kPi=TMath::Pi();		
  
  const Float_t kHeight=geometry->GetRadiatorToPads(); //Distance from Radiator to Pads in centimeters
  //printf("Distance to Pads:%f\n",kHeight);
 
  const Int_t kSpot=2;                                 //number of passes with spot algorithm
  const Int_t activ_tresh=0;                           //activation treshold to identify a track
  
  const Int_t kDimensionTheta=2;		       //Matrix dimension for angle Detection
  const Int_t kDimensionPhi=2;
  const Int_t kDimensionOmega=50;
  
  const Float_t SPOTp=.25;		              //Percentage of spot action
  const Float_t kMinOmega=.6;
  const Float_t kMaxOmega=.7;		      //Maximum Cherenkov angle to identify
  const Float_t kMinTheta=0;
  const Float_t kMaxTheta=0.5*kPi/180;	
  const Float_t kMinPhi=0;
  const Float_t kMaxPhi=20*kPi/180;

  const Float_t sigma=0.5;                          //half thickness of fiducial band in cm

  Float_t rechit[6];                                 //Reconstructed point data

  Int_t ***point = i3tensor(0,kDimensionTheta,0,kDimensionPhi,0,kDimensionOmega);
  Int_t ***point1 = i3tensor(0,kDimensionTheta,0,kDimensionPhi,0,kDimensionOmega);
  
  steptheta=(kMaxTheta-kMinTheta)/kDimensionTheta;
  stepphi=(kMaxPhi-kMinPhi)/kDimensionPhi;
  stepomega=(kMaxOmega-kMinOmega)/kDimensionOmega;

  static TH3F *Points = new TH3F("Points","Reconstructed points 3D",kDimensionTheta,0,kDimensionTheta,kDimensionPhi,0,kDimensionPhi,kDimensionOmega,0,kDimensionOmega);
  static TH2F *ThetaPhi = new TH2F("ThetaPhi","Theta-Phi projection",kDimensionTheta,0,kDimensionTheta,kDimensionPhi,0,kDimensionPhi);
  static TH2F *OmegaTheta = new TH2F("OmegaTheta","Omega-Theta projection",kDimensionTheta,0,kDimensionTheta,kDimensionOmega,0,kDimensionOmega);
  static TH2F *OmegaPhi = new TH2F("OmegaPhi","Omega-Phi projection",kDimensionPhi,0,kDimensionPhi,kDimensionOmega,0,kDimensionOmega);
  static TH3F *SpotPoints = new TH3F("Points","Reconstructed points 3D, spot",kDimensionTheta,0,kDimensionTheta,kDimensionPhi,0,kDimensionPhi,kDimensionOmega,0,kDimensionOmega);
  static TH2F *SpotThetaPhi = new TH2F("ThetaPhi","Theta-Phi projection, spot",kDimensionTheta,0,kDimensionTheta,kDimensionPhi,0,kDimensionPhi);
  static TH2F *SpotOmegaTheta = new TH2F("OmegaTheta","Omega-Theta projection, spot",kDimensionTheta,0,kDimensionTheta,kDimensionOmega,0,kDimensionOmega);
  static TH2F *SpotOmegaPhi = new TH2F("OmegaPhi","Omega-Phi projection, spot",kDimensionPhi,0,kDimensionPhi,kDimensionOmega,0,kDimensionOmega);
  static TH2F *DigitsXY = new TH2F("DigitsXY","Pads used for reconstruction",150,-25,25,150,-25,25);
  static TH1F *AngleAct = new TH1F("AngleAct","Activation per angle",100,.45,1);
  static TH1F *Activation = new TH1F("Activation","Activation per ring",100,0,25);
  Points->SetXTitle("theta");
  Points->SetYTitle("phi");
  Points->SetZTitle("omega");
  ThetaPhi->SetXTitle("theta");
  ThetaPhi->SetYTitle("phi");
  OmegaTheta->SetXTitle("theta");
  OmegaTheta->SetYTitle("omega");
  OmegaPhi->SetXTitle("phi");
  OmegaPhi->SetYTitle("omega");
  SpotPoints->SetXTitle("theta");
  SpotPoints->SetYTitle("phi");
  SpotPoints->SetZTitle("omega");
  SpotThetaPhi->SetXTitle("theta");
  SpotThetaPhi->SetYTitle("phi");
  SpotOmegaTheta->SetXTitle("theta");
  SpotOmegaTheta->SetYTitle("omega");
  SpotOmegaPhi->SetXTitle("phi");
  SpotOmegaPhi->SetYTitle("omega");
  AngleAct->SetFillColor(5);
  AngleAct->SetXTitle("rad");
  AngleAct->SetYTitle("activation");
  Activation->SetFillColor(5);
  Activation->SetXTitle("activation");

  Int_t ntracks = (Int_t)pRICH->TreeH()->GetEntries();
   
  Float_t trackglob[3];
  Float_t trackloc[3];

  //printf("Area de uma elipse com teta 0 e Omega 45:%f",Area(0,45));
   
  Int_t track;
	
  for (track=0; track<ntracks;track++) {
    gAlice->ResetHits();
    pRICH->TreeH()->GetEvent(track);
    TClonesArray *pHits  = pRICH->Hits();
    if (pHits == 0) return;
    Int_t nhits = pHits->GetEntriesFast();
    if (nhits == 0) continue;
    //Int_t nent=(Int_t)gAlice->TreeD()->GetEntries();
    AliRICHHit *mHit = 0;
    //Int_t npoints=0;
    
    Int_t counter=0, counter1=0;
    //Initialization
    for(i=0;i<kDimensionTheta;i++)
      {
	for(j=0;j<kDimensionPhi;j++)
	  {
	    for(k=0;k<kDimensionOmega;k++)
	      {
		counter++;
		point[i][j][k]=0;
		//printf("Dimensions theta:%d, phi:%d, omega:%d",kDimensionTheta,kDimensionPhi,kDimensionOmega);
		//printf("Resetting %d %d %d, time %d\n",i,j,k,counter);
		//-Noise_Level*(Area(i*kPi/(18*dimension),k*kMaxOmega/dimension)-Area((i-1)*kPi/(18*dimension),(k-1)*kMaxOmega/dimension));
		//printf("n-%f",-Noise_Level*(Area(i*kPi/(18*dimension),k*kMaxOmega/dimension)-Area((i-1)*kPi/(18*dimension),(k-1)*kMaxOmega/dimension)));
	      }
	  }
      }

    Int_t ncerenkovs = pRICH->Cerenkovs()->GetEntriesFast();
    
    
    originalOmega = 0;
    counter = 0;
    
    for (Int_t hit=0;hit<ncerenkovs;hit++) {
      AliRICHCerenkov* cHit = (AliRICHCerenkov*) pRICH->Cerenkovs()->UncheckedAt(hit);
      Float_t loss = cHit->fLoss;                 //did it hit the CsI?
      Float_t production = cHit->fProduction;     //was it produced in freon?  
      Float_t cherenkov = cHit->fCerenkovAngle;   //production cerenkov angle
      if (loss == 4 && production == 1)
	{
	  counter +=1;
	  originalOmega += cherenkov;
	  //printf("%f\n",cherenkov);
	}
    }

    originalOmega = originalOmega/counter;
    
    //printf("Cerenkovs       : %d\n",counter);
    
    mHit = (AliRICHHit*) pHits->UncheckedAt(0);
    Int_t nch  = mHit->Chamber();
    originalTheta = mHit->Theta();
    originalPhi = mHit->Phi();
    trackglob[0] = mHit->X();
    trackglob[1] = mHit->Y();
    trackglob[2] = mHit->Z();
    
    
    printf("\n--------------------------------------\n");
    printf("Chamber %d, track %d\n", nch, track);
    printf("Original omega: %f\n",originalOmega);
    
    iChamber = &(pRICH->Chamber(nch-1));
    
    printf("Nch:%d x:%f y:%f\n",nch,trackglob[0],trackglob[2]);
    
    iChamber->GlobaltoLocal(trackglob,trackloc);
    
    iChamber->LocaltoGlobal(trackloc,trackglob);
    
    
    cx=trackloc[0];
    cy=trackloc[2];

    printf("cy:%f  ", cy);
    
    if(counter != 0)    //if there are cerenkovs
      {
	
	AliRICHDigit *points = 0;
	TClonesArray *pDigits = pRICH->DigitsAddress(nch-1);   
	
	AliRICHRawCluster *cluster =0;
	TClonesArray *pClusters = pRICH->RawClustAddress(nch-1);

	Int_t maxcycle=0;

	//digitize from digits
	
	if(type==0)
	  {
	    gAlice->TreeD()->GetEvent(0);
	    Int_t ndigits = pDigits->GetEntriesFast();
	    maxcycle=ndigits;
	    //printf("Got %d digits\n",ndigits);
	  }
	
	//digitize from clusters
	
	if(type==1)
	  {
	    Int_t nent=(Int_t)gAlice->TreeR()->GetEntries();
	    gAlice->TreeR()->GetEvent(nent-1);
	    Int_t nclusters = pClusters->GetEntriesFast();
	    maxcycle=nclusters;
	    //printf("Got %d clusters\n",nclusters);
	  }



	
	counter=0;
	printf("Starting calculations\n");
	printf("           Start                                                                                              Finish\n");
	printf("Progress:     ");


	for(Float_t theta=0;theta<kMaxTheta;theta+=steptheta)
	  {		
	    printf(".");
	    for(Float_t phi=0;phi<=kMaxPhi;phi+=stepphi)
	      {		
		for(omega=kMinOmega;omega<=kMaxOmega;omega+=stepomega)
		  {		
		    //printf("Entering angle cycle\n");
		    omega1=SnellAngle(omega);
		    theta1=SnellAngle(theta);
		    
		    maxRadius = kHeight*(tan(omega1+theta1)+tan(omega1-theta1))/2;
		    minRadius = kHeight*tan(omega1);
		    eccentricity = sqrt(1-(minRadius*minRadius)/(maxRadius*maxRadius));
		


		    offset = kHeight*(tan(omega1+theta1)-tan(omega1-theta1))/2;
		    
		    //printf("phi:%f theta:%f omega:%f \n", phi,theta,omega);

		    //printf("offset:%f cx:%f cy:%f \n", offset,cx,cy);
		     
		    Float_t cxn = cx + offset * sin(phi);
		    Float_t cyn = cy + offset * cos(phi);

		    //printf("cxn:%f cyn:%f\n", cxn, cyn);

		    for (Int_t cycle=0;cycle<maxcycle;cycle++)
		      {	
			//printf("Entering point cycle");
			if(type==0)
			  {
			    points=(AliRICHDigit*) pDigits->UncheckedAt(cycle);
			    segmentation->GetPadC(points->PadX(), points->PadY(),x, y, z);
			  }
			
			if(type==1)
			  {
			    cluster=(AliRICHRawCluster*) pClusters->UncheckedAt(cycle);
			    x=cluster->fX;
			    y=cluster->fY;
			    q=cluster->fQ;
			  }
			
			if(type ==0 || q > 100)
			  
			  {

			    x=x-cxn;
			    y=y-cyn;
			    radius=TMath::Sqrt(TMath::Power(x,2)+TMath::Power(y,2));

			    phi_relative = asin(y/radius);
			    phi_relative = TMath::Abs(phi_relative - phi);
			    
			    angularadius = maxRadius*sqrt((1-eccentricity*eccentricity)/(1-eccentricity*eccentricity*cos(phi_relative)*cos(phi_relative)));
			
			    //printf("omega:%f min:%f rad:%f max:%f\n",omega, angularadius-sigma,radius,angularadius+sigma);


			    if((angularadius-sigma)<radius && (angularadius+sigma)>radius)
			      {
				printf("omega:%f min:%f rad:%f max:%f\n",omega, angularadius-sigma,radius,angularadius+sigma);

				bintheta=theta*kDimensionTheta/kMaxTheta;
				binphi=phi*kDimensionPhi/kMaxPhi;
				binomega=omega*kDimensionOmega/(kMaxOmega-kMinOmega);
				
				if(Int_t(bintheta+0.5)==Int_t(bintheta))
				  inttheta=Int_t(bintheta);
				else
				  inttheta=Int_t(bintheta+0.5);
				
				if(Int_t(binomega+0.5)==Int_t(binomega))
				  intomega=Int_t(binomega);
				else
			      intomega=Int_t(binomega+0.5);
				
				if(Int_t(binphi+0.5)==Int_t(binphi))
				  intphi=Int_t(binphi);
				else
				  intphi=Int_t(binphi+0.5);
				
				//printf("Point added at %d %d %d\n",inttheta,intphi,intomega);
				
				//if(type==0)
				point[inttheta][intphi][intomega]+=1;
				//if(type==1)
				//point[inttheta][intphi][intomega]+=(Int_t)(q);
				
				//printf("Omega stored:%d\n",intomega);
				Points->Fill(inttheta,intphi,intomega,(float) 1);
				ThetaPhi->Fill(inttheta,intphi,(float) 1);
				OmegaTheta->Fill(inttheta,intomega,(float) 1);
				OmegaPhi->Fill(intphi,intomega,(float) 1);
				//printf("Filling at %d %d %d\n",Int_t(theta*kDimensionTheta/kMaxTheta),Int_t(phi*kDimensionPhi/kMaxPhi),Int_t(omega*kDimensionOmega/kMaxOmega));
			      }
				//if(omega<kMaxOmega)point[Int_t(theta)][Int_t(phi)][Int_t(omega)]+=1;
			  }
		      }
		  }
	      }
	    
	  }
		//printf("Used %d digits for theta %3.1f\n",counter1, theta*180/kPi);
	
	printf("\n");

	meanradius=meanradius/counter;
	//printf("Mean radius:%f, counter:%d\n",meanradius,counter);
	rechit[5]=meanradius;
	//printf("Used %d digits\n",counter1);
	//printf("\n");

	if(nev<2)
	  {
	    if(nev==0)
	      {
		fc1->cd(1);
		Points->Draw("colz");
		fc1->cd(2);
		ThetaPhi->Draw("colz");
		fc1->cd(3);
		OmegaTheta->Draw("colz");
		fc1->cd(4);
		OmegaPhi->Draw("colz");
		fc3->cd();
		DigitsXY->Draw("colz");
	      }
	    else
	      {
		//fc1->cd(1);
		//Points->Draw("same");
		//fc1->cd(2);
		//ThetaPhi->Draw("same");
		//fc1->cd(3);
		//OmegaTheta->Draw("same");
		//fc1->cd(4);
		//OmegaPhi->Draw("same");	
	      }
	  }
	
    
	//SPOT execute twice
	for(Int_t s=0;s<kSpot;s++)
	  {
	    printf("     Applying Spot algorithm, pass %d\n", s);
	
	//buffer copy
	for(i=0;i<=kDimensionTheta;i++)
	  {
	    for(j=0;j<=kDimensionPhi;j++)
	      {
		for(k=0;k<=kDimensionOmega;k++)
		  {
		    point1[i][j][k]=point[i][j][k];	
		  }
	      }
	  }

	//SPOT algorithm			
	for(i=1;i<kDimensionTheta-1;i++)
	  {
	    for(j=1;j<kDimensionPhi-1;j++)
	      {
		for(k=1;k<kDimensionOmega-1;k++)
		  {
		    if((point[i][k][j]>point[i-1][k][j])&&(point[i][k][j]>point[i+1][k][j])&&
		       (point[i][k][j]>point[i][k-1][j])&&(point[i][k][j]>point[i][k+1][j])&&
		       (point[i][k][j]>point[i][k][j-1])&&(point[i][k][j]>point[i][k][j+1]))
		      {
			//cout<<"SPOT"<<endl;
			//Execute SPOT on point											   	
			point1[i][j][k]+=Int_t(SPOTp*(point[i-1][k][j]+point[i+1][k][j]+point[i][k-1][j]+point[i][k+1][j]+point[i][k][j-1]+point[i][k][j+1]));    
			point1[i-1][k][j]=Int_t(SPOTp*point[i-1][k][j]);
			point1[i+1][k][j]=Int_t(SPOTp*point[i+1][k][j]);
			point1[i][k-1][j]=Int_t(SPOTp*point[i][k-1][j]);
			point1[i][k+1][j]=Int_t(SPOTp*point[i][k+1][j]);
			point1[i][k][j-1]=Int_t(SPOTp*point[i][k][j-1]);
			point1[i][k][j+1]=Int_t(SPOTp*point[i][k][j+1]);
		      }
		  }
	      }
	  }
	
	//copy from buffer copy
	counter1=0;
	for(i=1;i<kDimensionTheta;i++)
	  {
	    for(j=1;j<kDimensionPhi;j++)
	      {
		for(k=1;k<kDimensionOmega;k++)
		  {
		    point[i][j][k]=point1[i][j][k];
		    if(nev<20)
		      {
			if(s==kSpot-1)
			  {
			    if(point1[i][j][k] != 0)
			      {
				SpotPoints->Fill(i,j,k,(float) point1[i][j][k]);
				//printf("Random number %f\n",random->Rndm(2));
				//if(random->Rndm() < .2)
				  //{
				SpotThetaPhi->Fill(i,j,(float) point1[i][j][k]);
				SpotOmegaTheta->Fill(i,k,(float) point1[i][j][k]);
				SpotOmegaPhi->Fill(j,k,(float) point1[i][j][k]);
				    counter1++;
				  //}
				//printf("Filling at %d %d %d value %f\n",i,j,k,(float) point1[i][j][k]);
			      }
			  }
		      }
		    //if(point1[i][j][k] != 0)
		      //printf("Last transfer point: %d, point1, %d\n",point[i][j][k],point1[i][j][k]);
		  }
	      }
	  }
      }
    
    //printf("Filled %d cells\n",counter1);

	if(nev<2)
	  {
	    if(nev==0)
	      {
		fc2->cd(1);
		SpotPoints->Draw("colz");
		fc2->cd(2);
		SpotThetaPhi->Draw("colz");
		fc2->cd(3);
		SpotOmegaTheta->Draw("colz");
		fc2->cd(4);
		SpotOmegaPhi->Draw("colz");
	      }
	    else
	      {
		//fc2->cd(1);
		//SpotPoints->Draw("same");
		//fc2->cd(2);
		//SpotThetaPhi->Draw("same");
		//fc2->cd(3);
		//SpotOmegaTheta->Draw("same");
		//fc2->cd(4);
		//SpotOmegaPhi->Draw("same");	
	      }
	  }
    
    
	//Identification is equivalent to maximum determination
	max=0;maxi=0;maxj=0;maxk=0;
	
	printf("     Proceeding to identification");
	
	for(i=0;i<kDimensionTheta;i++)
	  for(j=0;j<kDimensionPhi;j++)
	    for(k=0;k<kDimensionOmega;k++)
	      if(point[i][j][k]>max)
		{
		  //cout<<"maxi="<<i*90/dimension<<" maxj="<<j*90/dimension<<" maxk="<<k*kMaxOmega/dimension*180/kPi<<" max="<<max<<endl;
		maxi=i;maxj=j;maxk=k;
		max=point[i][j][k];
		printf(".");
		//printf("Max Omega %d, Max Theta %d, Max Phi %d (%d counts)\n",maxk,maxi,maxj,max);
		}
	printf("\n");
 
	Float_t FinalOmega; 
	Float_t FinalTheta;
	Float_t FinalPhi;

   
	if(max>activ_tresh)
	  {
	    FinalOmega = maxk*(kMaxOmega-kMinOmega)/kDimensionOmega; 
	    FinalTheta = maxi*kMaxTheta/kDimensionTheta;
	    FinalPhi = maxj*kMaxPhi/kDimensionPhi;
	    
	    FinalOmega += kMinOmega;
	  }
	else
	  {
	    FinalOmega = 0; 
	    FinalTheta = 0;
	    FinalPhi = 0;
	    
	    printf("     Ambiguous data!\n");
	  }

	    

	//printf("Detected angle for height %3.1f and for center %3.1f %3.1f:%f\n",h,cx,cy,maxk*kPi/(kDimensionTheta*4));
	printf("     Indentified angles: cerenkov - %f, theta - %3.1f, phi - %3.1f (%f activation)\n", FinalOmega, FinalTheta*180/kPi, FinalPhi*180/kPi, max);
	//printf("Detected angle for height %3.1f and for center %3.1f %3.1f:%f\n",kHeight,cx,cy,maxk);

	AngleAct->Fill(FinalOmega, (float) max);
	Activation->Fill(max, (float) 1);

	//fscanf(omegas,"%f",&realomega);
	//fscanf(thetas,"%f",&realtheta);
	//printf("Real Omega: %f",realomega);			
	//cout<<"Detected:theta="<<maxi*90/kDimensionTheta<<"phi="<<maxj*90/kDimensionPhi<<"omega="<<maxk*kMaxOmega/kDimensionOmega*180/kPi<<" OmegaError="<<fabs(maxk*kMaxOmega/kDimensionOmega*180/kPi-realomega)<<" ThetaError="<<fabs(maxi*90/kDimensionTheta-realtheta)<<endl<<endl;		
    
	//fprintf(results,"Center Coordinates, cx=%6.2f cy=%6.2f, Real Omega=%6.2f, Detected Omega=%6.2f, Omega Error=%6.2f Theta Error=%6.2f\n",cx,cy,realomega,maxk*kMaxOmega/kDimensionOmega*180/kPi,fabs(maxk*kMaxOmega/kDimensionOmega*180/kPi-realomega),fabs(maxi*90/kDimensionTheta-realtheta));
    
    /*for(j=0;j<np;j++)
      pointpp(maxj*90/kDimensionTheta,maxi*90/kDimensionPhi,maxk*kMaxOmega/kDimensionOmega*180/kPi,cx,cy);//Generates a point on the elipse*/		    


    //Start filling rec. hits
    
    rechit[0] = FinalTheta;
    rechit[1] = 90*kPi/180 + FinalPhi;
    rechit[2] = FinalOmega;
    rechit[3] = cx;
    rechit[4] = cy;

    //CreatePoints(FinalTheta, 270*kPi/180 + FinalPhi, FinalOmega, kHeight);
       
    //printf ("track %d, theta %f, phi %f, omega %f\n\n\n",track,rechit[0],rechit[1],rechit[2]);
    //printf("rechit %f %f %f %f %f\n",rechit[0],rechit[1],rechit[2],rechit[3],rechit[4]);
    //printf("Chamber:%d",nch);


      }

    else   //if no cerenkovs
      
      {
	
	rechit[0] = 0;
	rechit[1] = 0;
	rechit[2] = 0;
	rechit[3] = 0;
	rechit[4] = 0;
	originalOmega = 0; 
	originalTheta = 0;
	originalPhi =0;
      }

    
    // fill rechits
    pRICH->AddRecHit3D(nch-1,rechit,originalOmega, originalTheta, originalPhi);
    printf("track %d, theta r:%f o:%f, phi r:%f o:%f, omega r:%f o:%f cx:%f cy%f\n\n\n", track, rechit[0], originalTheta, rechit[1], originalPhi, rechit[2], originalOmega, cx, cy);

  }

  if(type==1)  //reco from clusters
    {
      pRICH->ResetRawClusters();
      //Int_t nent=(Int_t)gAlice->TreeR()->GetEntries();
      //gAlice->TreeR()->GetEvent(track);
      //printf("Going to branch %d\n",track);
      //gAlice->GetEvent(nev);
    }

	
  //printf("\n\n\n\n");
  gAlice->TreeR()->Fill();
  TClonesArray *fRec;
  for (i=0;i<kNCH;i++) {
    fRec=pRICH->RecHitsAddress3D(i);
    int ndig=fRec->GetEntriesFast();
    printf ("Chamber %d, rings %d\n",i+1,ndig);
  }

  fc4->cd(1);
  AngleAct->Draw();
  fc4->cd(2);
  Activation->Draw();

  pRICH->ResetRecHits3D();

  free_i3tensor(point,0,kDimensionTheta,0,kDimensionPhi,0,kDimensionOmega);
  free_i3tensor(point1,0,kDimensionTheta,0,kDimensionPhi,0,kDimensionOmega);
  

}

Int_t ***AliRICHDetectV1::i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
// allocate a Int_t 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] 
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
  Int_t ***t;
  
  int NR_END=1; 

  // allocate pointers to pointers to rows 
  t=(Int_t ***) malloc((size_t)((nrow+NR_END)*sizeof(Int_t**)));
  if (!t) printf("allocation failure 1 in f3tensor()");
  t += NR_END;
  t -= nrl;
  
  // allocate pointers to rows and set pointers to them 
  t[nrl]=(Int_t **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(Int_t*)));
  if (!t[nrl]) printf("allocation failure 2 in f3tensor()");
  t[nrl] += NR_END;
  t[nrl] -= ncl;
  
  // allocate rows and set pointers to them 
  t[nrl][ncl]=(Int_t *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(Int_t)));
  if (!t[nrl][ncl]) printf("allocation failure 3 in f3tensor()");
  t[nrl][ncl] += NR_END;
  t[nrl][ncl] -= ndl;
  
  for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
  for(i=nrl+1;i<=nrh;i++) {
    t[i]=t[i-1]+ncol;
    t[i][ncl]=t[i-1][ncl]+ncol*ndep;
    for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
  }
  
  // return pointer to array of pointers to rows 
  return t;
}


void AliRICHDetectV1::free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,long ndl, long ndh)
// free a Int_t f3tensor allocated by i3tensor()
{
  int NR_END=1; 

  free((char*) (t[nrl][ncl]+ndl-NR_END));
  free((char*) (t[nrl]+ncl-NR_END));
  free((char*) (t+nrl-NR_END));
}
