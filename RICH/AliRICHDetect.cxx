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


#include "AliRICH.h"
#include "AliRICHPoints.h"
#include "AliRICHDetect.h"
#include "AliRICHDigit.h"
#include "AliRICHRawCluster.h"
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



ClassImp(AliRICHDetect)
//___________________________________________
AliRICHDetect::AliRICHDetect() : TObject()
{

// Default constructor 

  fc1 = 0;
  fc2 = 0;
  fc3 = 0;

}

//___________________________________________
AliRICHDetect::AliRICHDetect(const char *name, const char *title)
    : TObject()
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
AliRICHDetect::~AliRICHDetect()
{
    
// Destructor

}


void AliRICHDetect::Detect(Int_t nev, Int_t type)
{  	
    
//
// Detection algorithm


  //printf("Detection started!\n");
  Float_t omega,omega1,theta1,phi_relative,steptheta,stepphi,x,y,q=0,z,cx,cy,l,aux1,aux2,aux3,max,radius=0,meanradius=0;
  Int_t maxi,maxj,maxk;
  Float_t originalOmega, originalPhi, originalTheta;
  Float_t binomega, bintheta, binphi;
  Int_t intomega, inttheta, intphi;
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
 
  const Int_t kSpot=0;                                 //number of passes with spot algorithm
  
  const Int_t kDimensionTheta=100;		       //Matrix dimension for angle Detection
  const Int_t kDimensionPhi=100;
  const Int_t kDimensionOmega=100;
  
  const Float_t SPOTp=.25;		              //Percentage of spot action
  const Float_t kMinOmega=.4;
  const Float_t kMaxOmega=.7;		      //Maximum Cherenkov angle to identify
  const Float_t kMinTheta=0;
  const Float_t kMaxTheta=10*kPi/180;	
  const Float_t kMinPhi=0;
  const Float_t kMaxPhi=360*kPi/180;

  Float_t rechit[6];                                 //Reconstructed point data

  Int_t ***point = i3tensor(0,kDimensionTheta,0,kDimensionPhi,0,kDimensionOmega);
  Int_t ***point1 = i3tensor(0,kDimensionTheta,0,kDimensionPhi,0,kDimensionOmega);
  
  steptheta=(kMaxTheta-kMinTheta)/kDimensionTheta;
  stepphi=(kMaxPhi-kMinPhi)/kDimensionPhi;

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
    AliRICHhit *mHit = 0;
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
    
    printf("Cerenkovs       : %d\n",counter);
    
    if(counter != 0)    //if there are cerenkovs
      {
	originalOmega = originalOmega/counter;
	printf("Original omega: %f\n",originalOmega);
     


	mHit = (AliRICHhit*) pHits->UncheckedAt(0);
	Int_t nch  = mHit->Chamber();
	originalTheta = mHit->Theta();
	originalPhi = mHit->Phi();
	trackglob[0] = mHit->X();
	trackglob[1] = mHit->Y();
	trackglob[2] = mHit->Z();
	
	printf("\n--------------------------------------\n");
	printf("Chamber %d, track %d\n", nch, track);

	
	iChamber = &(pRICH->Chamber(nch-1));
    
	//printf("Nch:%d\n",nch);

	iChamber->GlobaltoLocal(trackglob,trackloc);
    
	iChamber->LocaltoGlobal(trackloc,trackglob);
       
       
	cx=trackloc[0];
	cy=trackloc[2];
     
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
	    printf("Got %d digits\n",ndigits);
	  }
	
	//digitize from clusters
	
	if(type==1)
	  {
	    Int_t nent=(Int_t)gAlice->TreeR()->GetEntries();
	    gAlice->TreeR()->GetEvent(nent-1);
	    Int_t nclusters = pClusters->GetEntriesFast();
	    maxcycle=nclusters;
	    printf("Got %d clusters\n",nclusters);
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
		//printf("Phi:%3.1f\n", phi*180/kPi);
		counter1=0;
		for (Int_t cycle=0;cycle<maxcycle;cycle++)
		  {	
		    
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
			
			x=x-cx;
			y=y-cy;
			radius=TMath::Sqrt(TMath::Power(x,2)+TMath::Power(y,2));
			
			//calculation of relative phi angle of digit
			
			phi_relative = acos(y/radius);
			phi_relative = TMath::Abs(phi_relative - phi);
			
			
			if(radius>4)
			  {
			    meanradius+=radius;
			    counter++;
			    
			    
			    //if (radius < (2*kHeight*tan(theta+kMaxOmega)))
			    if (radius < (2*kHeight*tan(kMaxOmega)))
			      //if(Fiducial(x,y,theta,phi,kHeight,kMaxOmega,kMinOmega))
			      {
				
				if(phi==0)
				  {
				//printf("Radius: %f, Max Radius: %f\n",radius,2*kHeight*tan(theta+kMaxOmega)*3/4);
				//printf("Loaded digit %d with coordinates x:%f, y%f\n",dig,x,y);
				//printf("Using digit %d, for theta %f\n",dig,theta);
				  }
				
				counter1++;
				
				l=kHeight/cos(theta);
				
				//main calculation
				
				DigitsXY->Fill(x,y,(float) 1);
				
				theta1=SnellAngle(theta);
				
				aux1=-y*sin(phi)+x*cos(phi);
				aux2=y*cos(phi)+x*sin(phi);
				aux3=( TMath::Power(aux1,2)+TMath::Power(cos(theta1)*aux2 ,2))/TMath::Power(sin(theta1)*aux2+l,2);
				omega=atan(sqrt(aux3));
				
				//omega is distorted, theta1 is distorted
				
				if(InvSnellAngle(omega+TMath::Abs(cos(phi_relative))*theta1)<999)
				  {
				    omega1=InvSnellAngle(omega+TMath::Abs(cos(phi_relative))*theta);
				    theta1=InvSnellAngle(theta1);
				    
				  }
				else
				  {
				    omega1=0;
				    theta1=0;
				  }
				
				if(omega1<kMaxOmega && omega1>kMinOmega)
				  {
				//printf("Omega found:%f\n",omega);
				    omega1=omega1-kMinOmega;
				    
				//printf("Omega: %f Theta: %3.1f Phi:%3.1f\n",omega, theta*180/kPi, phi*180/kPi);
				    
				    bintheta=theta1*kDimensionTheta/kMaxTheta;
				    binphi=phi*kDimensionPhi/kMaxPhi;
				    binomega=omega1*kDimensionOmega/(kMaxOmega-kMinOmega);
				    
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
				    
				    if(type==0)
				      point[inttheta][intphi][intomega]+=1;
				    if(type==1)
				      point[inttheta][intphi][intomega]+=(Int_t)(q);
				    
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
		//printf("Used %d digits for theta %3.1f\n",counter1, theta*180/kPi);
	      }
	    
	  }

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
    
	Float_t FinalOmega = maxk*(kMaxOmega-kMinOmega)/kDimensionOmega; 
	Float_t FinalTheta = maxi*kMaxTheta/kDimensionTheta;
	Float_t FinalPhi = maxj*kMaxPhi/kDimensionPhi;

	FinalOmega += kMinOmega;
    
	//printf("Detected angle for height %3.1f and for center %3.1f %3.1f:%f\n",h,cx,cy,maxk*kPi/(kDimensionTheta*4));
	printf("     Indentified angles: cerenkov - %f, theta - %3.1f, phi - %3.1f (%f activation)\n", FinalOmega, FinalTheta*180/kPi, FinalPhi*180/kPi, max);
	//printf("Detected angle for height %3.1f and for center %3.1f %3.1f:%f\n",kHeight,cx,cy,maxk);

	AngleAct->Fill(FinalOmega, (float) max);
	Activation->Fill(max, (float) 1);

	if(nev==0)
	      {
		fc4->cd(1);
		AngleAct->Draw();
		fc4->cd(2);
		Activation->Draw();
	      }
	    else
	      {
		fc4->cd(1);
		AngleAct->Draw("same");
		fc4->cd(2);
		Activation->Draw("same");
	      }
	  

	//fscanf(omegas,"%f",&realomega);
	//fscanf(thetas,"%f",&realtheta);
	//printf("Real Omega: %f",realomega);			
	//cout<<"Detected:theta="<<maxi*90/kDimensionTheta<<"phi="<<maxj*90/kDimensionPhi<<"omega="<<maxk*kMaxOmega/kDimensionOmega*180/kPi<<" OmegaError="<<fabs(maxk*kMaxOmega/kDimensionOmega*180/kPi-realomega)<<" ThetaError="<<fabs(maxi*90/kDimensionTheta-realtheta)<<endl<<endl;		
    
	//fprintf(results,"Center Coordinates, cx=%6.2f cy=%6.2f, Real Omega=%6.2f, Detected Omega=%6.2f, Omega Error=%6.2f Theta Error=%6.2f\n",cx,cy,realomega,maxk*kMaxOmega/kDimensionOmega*180/kPi,fabs(maxk*kMaxOmega/kDimensionOmega*180/kPi-realomega),fabs(maxi*90/kDimensionTheta-realtheta));
    
    /*for(j=0;j<np;j++)
      pointpp(maxj*90/kDimensionTheta,maxi*90/kDimensionPhi,maxk*kMaxOmega/kDimensionOmega*180/kPi,cx,cy);//Generates a point on the elipse*/		    


    //Start filling rec. hits
    
    rechit[0] = FinalTheta;
    rechit[1] = FinalPhi - 90*kPi/180;
    rechit[2] = FinalOmega;
    rechit[3] = cx;
    rechit[4] = cy;

    //CreatePoints(FinalTheta, 270*kPi/180 + FinalPhi, FinalOmega, kHeight);
       
    printf ("track %d, theta %f, phi %f, omega %f\n\n\n",track,rechit[0],rechit[1]*180/kPi,rechit[2]);
    
    // fill rechits
    pRICH->AddRecHit3D(nch-1,rechit,originalOmega, originalTheta, originalPhi);
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

      }

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
  pRICH->ResetRecHits3D();

  free_i3tensor(point,0,kDimensionTheta,0,kDimensionPhi,0,kDimensionOmega);
  free_i3tensor(point1,0,kDimensionTheta,0,kDimensionPhi,0,kDimensionOmega);
}



Float_t AliRICHDetect:: Area(Float_t theta,Float_t omega)
{

//
// Calculates area of an ellipse for given incidence angles    


    Float_t area;
    const Float_t kHeight=9.25;                       //Distance from Radiator to Pads in pads
    
    area=TMath::Pi()*TMath::Power(kHeight*tan(omega),2)/TMath::Power(TMath::Power(cos(theta),2)-TMath::Power(tan(omega)*sin(theta),2),3/2);
    
    return (area);
}


Int_t ***AliRICHDetect::i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
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

void AliRICHDetect::free_i3tensor(int ***t, long nrl, long nrh, long ncl, long nch,long ndl, long ndh)
// free a Int_t f3tensor allocated by i3tensor()
{
  int NR_END=1; 

  free((char*) (t[nrl][ncl]+ndl-NR_END));
  free((char*) (t[nrl]+ncl-NR_END));
  free((char*) (t+nrl-NR_END));
}


Float_t AliRICHDetect:: SnellAngle(Float_t iangle)
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
  return rangle;
}

Float_t AliRICHDetect:: InvSnellAngle(Float_t rangle)
{ 

// Compute the inverse Snell angle

  Float_t nfreon  = 1.295;
  Float_t nquartz = 1.585;
  Float_t ngas = 1;

  Float_t siniangle;
  Float_t iangle;
  Float_t a1,a2;

  a1=nfreon/nquartz;
  a2=nquartz/ngas;

  siniangle = sin(rangle)/(a1*a2);
  iangle = asin(siniangle);

  if(siniangle>1.) {
     iangle = 999.;
     return iangle;
  }
  
  iangle = asin(siniangle);
  return iangle;
}



//________________________________________________________________________________
void AliRICHDetect::CreatePoints(Float_t theta, Float_t phi, Float_t omega, Float_t h)
{
  
  // Create points along the ellipse equation
  
  Int_t s1,s2;
  Float_t fiducial=h*TMath::Tan(omega+theta), l=h/TMath::Cos(theta), xtrial, y=0, c0, c1, c2;
  //TRandom *random=new TRandom();
  
  static TH2F *REllipse = new TH2F("REllipse","Reconstructed ellipses",150,-25,25,150,-25,25);

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
	  //printf("Fiducial %f  for omega:%f theta:%f phi:%f\n",fiducial,omega,theta,fphi);
	  Float_t alfa1=theta;
	  Float_t theta1=phi;
	  Float_t omega1=omega;
	  
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
      xtrial=(-c1+s2*TMath::Sqrt(c1*c1-4*c2*c0))/(2*c2);
      //cout<<"x="<<xtrial<<" y="<<cy+y<<endl;
      //printf("Coordinates: %f %f\n",xtrial,fCy+y);

      REllipse->Fill(xtrial,y);

      //printf("Coordinates: %f %f %f\n",vectorGlob[0],vectorGlob[1],vectorGlob[2]);
    }
  
  fc3->cd(2);
  REllipse->Draw();
}

Int_t AliRICHDetect::Fiducial(Float_t rotx, Float_t roty, Float_t theta, Float_t phi, Float_t height, Float_t maxOmega, Float_t minOmega)
{

  Float_t x,y,y1,h,omega1,omega2;

  Float_t a,b, offset;
  Float_t a1,b1, offset1;

  h = height;

  //refraction calculation
  
  theta = SnellAngle(theta);
  phi = phi - TMath::Pi();
  omega1 = SnellAngle(maxOmega);
  omega2 = SnellAngle(minOmega);

  //maximum zone
  a = h*(tan(omega1+theta)+tan(omega1-theta))/2;
  b = h*tan(omega1);
 
  offset = h*(tan(omega1+theta)-tan(omega1-theta))/2;
  
 
  //minimum zone
  a1 = h*(tan(omega2+theta)+tan(omega2-theta))/2;
  b1 = h*tan(omega2);
 
  offset1 = h*(tan(omega2+theta)-tan(omega2-theta))/2;
  

  //rotation to phi=0
  x = rotx*cos(phi)+roty*sin(phi);
  y = -rotx*sin(phi)+roty*cos(phi) - offset;
  y1 = -rotx*sin(phi)+roty*cos(phi) - offset1;

  
  if(x*x/a+y*y/b<1 && x*x/a1+y1*y1/b1>1)
    return 1;
  else
    return 0;

}
