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
  Revision 1.10  2001/02/13 20:39:06  jbarbosa
  Changes to make it work with new IO.

  Revision 1.9  2001/01/22 21:39:11  jbarbosa
  Several tune-ups

  Revision 1.8  2000/11/15 15:52:53  jbarbosa
  Turned on spot algorithm.

  Revision 1.7  2000/11/01 15:37:05  jbarbosa
  Updated to use its own rec. point object.

  Revision 1.6  2000/10/02 21:28:12  fca
  Removal of useless dependecies via forward declarations

  Revision 1.5  2000/06/30 16:30:28  dibari
  Disabled writing to rechits.

  Revision 1.4  2000/06/15 15:46:59  jbarbosa
  Corrected compilation errors on HP-UX (replaced pow with TMath::Power)

  Revision 1.3  2000/06/13 13:15:41  jbarbosa
  Still some code cleanup done (variable names)

  Revision 1.2  2000/06/12 15:19:30  jbarbosa
  Cleaned up version.

  Revision 1.1  2000/04/19 13:05:14  morsch
  J. Barbosa's spot reconstruction algorithm.

*/


#include "AliRICH.h"
#include "AliRICHPoints.h"
#include "AliRICHDetect.h"
#include "AliRICHHit.h"
#include "AliRICHDigit.h"
#include "AliRun.h"
#include "TParticle.h"
#include "TTree.h"
#include "TMath.h"
#include "TRandom.h"



ClassImp(AliRICHDetect)
//___________________________________________
AliRICHDetect::AliRICHDetect() : TObject()
{

// Default constructor 

    //fChambers = 0;
}

//___________________________________________
AliRICHDetect::AliRICHDetect(const char *name, const char *title)
    : TObject()
{
    
// Constructor

    /*fChambers = new TObjArray(7);
    for (Int_t i=0; i<7; i++) {
    
	(*fChambers)[i] = new AliRICHchamber();  
	
    } */     
}


void AliRICHDetect::Detect()
{  	
    
//
// Detection algorithm


  //printf("Detection started!\n");
  Float_t omega,steptheta,stepphi,x,y,cx,cy,l,aux1,aux2,aux3,maxi,maxj,maxk,max;
  //Float_t theta,phi,realomega,realtheta;
  Int_t i,j,k;
 
  
  //const Float_t Noise_Level=0;                //Noise Level in percentage of mesh points
  //const Float_t t=0.6;			//Softening of Noise Correction (factor)
  
  const Float_t kPi=3.1415927;		
  
  const Float_t kHeight=10;                     //Distance from Radiator to Pads in pads
 
  const Int_t kSpot=0;                          //number of passes with spot algorithm
  
  const Int_t kDimensionTheta=50;		//Matrix dimension for angle Detection
  const Int_t kDimensionPhi=50;
  const Int_t kDimensionOmega=50;
  
  const Float_t SPOTp=.2;		        //Percentage of spot action
  //const Int_t np=500;		                //Number of points to reconstruct elipse 
  const Float_t kMinOmega=30*kPi/180;
  const Float_t kMaxOmega=65*kPi/180;		//Maximum Cherenkov angle to identify
 
  const Float_t kCorr=.5;                       //Correction factor, accounting for aberration, refractive index, etc.
   
  Int_t point[kDimensionTheta][kDimensionPhi][kDimensionOmega];
  Int_t point1[kDimensionTheta][kDimensionPhi][kDimensionOmega];
  
  steptheta=kPi/kDimensionTheta;
  stepphi=kPi/kDimensionPhi;

  AliRICHChamber*       iChamber;
  
  AliRICH *pRICH  = (AliRICH*)gAlice->GetDetector("RICH");
  Int_t ntracks = (Int_t)gAlice->TreeH()->GetEntries();
  //Int_t ntrks = gAlice->GetNtrack();
  
  Float_t trackglob[3];
  Float_t trackloc[3];

  //printf("Got ntracks:%d\n",ntracks);
  /*TVector *xp = new TVector(1000);
  TVector *yp = new TVector(1000);
  TVector *zp = new TVector(1000);
  TVector *ptrk = new TVector(1000);
  TVector *phit = new TVector(1000);*/
  
  //printf("Area de uma elipse com teta 0 e Omega 45:%f",Area(0,45));
    
  Int_t track;
	
  for (track=0; track<ntracks;track++) {
    gAlice->ResetHits();
    gAlice->TreeH()->GetEvent(track);
    TClonesArray *pHits  = pRICH->Hits();
    if (pHits == 0) return;
    Int_t nhits = pHits->GetEntriesFast();
    if (nhits == 0) continue;
    Int_t nent=(Int_t)gAlice->TreeD()->GetEntries();
    gAlice->TreeD()->GetEvent(1);
    AliRICHHit *mHit = 0;
    AliRICHDigit *points = 0;
    //Int_t npoints=0;
    
    Int_t counter=0;
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
    mHit = (AliRICHHit*) pHits->UncheckedAt(0);
    //printf("Aqui vou eu\n");
    Int_t nch  = mHit->fChamber;
    //printf("Aqui fui eu\n");
    trackglob[0] = mHit->X();
    trackglob[1] = mHit->Y();
    trackglob[2] = mHit->Z();

    cx=trackglob[0];
    cy=trackglob[2];
    
    
    //printf("Chamber processed:%d\n",nch);

    printf("\nChamber %d, particle at: %3.1f %3.1f,\n",nch,trackglob[0],trackglob[2]);

    iChamber = &(pRICH->Chamber(nch-1));
    
    //printf("Nch:%d\n",nch);

    iChamber->GlobaltoLocal(trackglob,trackloc);
    
    //printf("Transformation 1: %3.1f %3.1f %3.1f\n",trackloc[0],trackloc[1],trackloc[2]);


    iChamber->LocaltoGlobal(trackloc,trackglob);
       
    //printf("Transformation 2: %3.1f %3.1f %3.1f\n",trackglob[0],trackglob[1],trackglob[2]);
    
    
     

    TClonesArray *pDigits = pRICH->DigitsAddress(nch-1);   
    Int_t ndigits = pDigits->GetEntriesFast();
    
    //printf("Got %d digits\n",ndigits);

    //printf("Starting calculations\n");
    
    for(Float_t theta=0;theta<kPi/18;theta+=steptheta)
      {			
	for(Float_t phi=0;phi<=kPi/3;phi+=stepphi)
	  {		       
	    for (Int_t dig=0;dig<ndigits;dig++)
	      {	
		points=(AliRICHDigit*) pDigits->UncheckedAt(dig);
		
		x=points->fPadX-cx;
		y=points->fPadY-cy;
		//printf("Loaded digit %d with coordinates x:%f, y%f\n",dig,x,y);
		//cout<<"x="<<x<<" y="<<y<<endl;
			
		if (sqrt(TMath::Power(x,2)+TMath::Power(y,2))<kHeight*tan(theta+kMaxOmega)*3/4)
		  {
		    
	
		    l=kHeight/cos(theta);
		    
		    aux1=-y*sin(phi)+x*cos(phi);
		    aux2=y*cos(phi)+x*sin(phi);
		    aux3=( TMath::Power(aux1,2)+TMath::Power(cos(theta)*aux2 ,2))/TMath::Power(sin(theta)*aux2+l,2);
			//cout<<"aux1="<<aux1<<" aux2="<<aux2<<" aux3="<<aux3;
			    
		    omega=atan(sqrt(aux3));
		    //printf("Omega: %f\n",omega);
		    
		    //cout<<"\ni="<<i<<" theta="<<Int_t(2*theta*dimension/kPi)<<" phi="<<Int_t(2*phi*dimension/kPi)<<" omega="<<Int_t(2*omega*dimension/kPi)<<endl<<endl;
		    //{Int_t lixo;cin>>lixo;}
		    if(omega<kMaxOmega && omega>kMinOmega)
		      {
			omega=omega-kMinOmega;
			//point[Int_t(2*theta*kDimensionTheta/kPi)][Int_t(2*phi*kDimensionPhi/kPi)][Int_t(kCorr*2*omega*kDimensionOmega/kMaxOmega)]+=1;	
			point[Int_t(2*theta*kDimensionTheta/kPi)][Int_t(2*phi*kDimensionPhi/kPi)][Int_t(kCorr*(omega/(kMaxOmega-kMinOmega)*kDimensionOmega))]+=1;		
		      }
		    //if(omega<kMaxOmega)point[Int_t(theta)][Int_t(phi)][Int_t(omega)]+=1;
		  }
		}
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
	for(i=1;i<kDimensionTheta;i++)
	  {
	    for(j=1;j<kDimensionPhi;j++)
	      {
		for(k=1;k<kDimensionOmega;k++)
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
	for(i=1;i<kDimensionTheta;i++)
	  {
	    for(j=1;j<kDimensionPhi;j++)
	      {
		for(k=1;k<kDimensionOmega;k++)
		  {
		    point[i][j][k]=point1[i][j][k];
		    //if(point1[i][j][k] != 0)
		      //printf("Last transfer point: %d, point1, %d\n",point[i][j][k],point1[i][j][k]);
		  }
	      }
	  }
      }
    
    
    //Identification is equivalent to maximum determination
    max=0;maxi=0;maxj=0;maxk=0;
    
    printf("     Proceeding to identification");
    
    for(i=1;i<kDimensionTheta-3;i++)
      for(j=1;j<=kDimensionPhi-3;j++)
	for(k=0;k<=kDimensionOmega;k++)
	    if(point[i][j][k]>max)
	      {
		//cout<<"maxi="<<i*90/dimension<<" maxj="<<j*90/dimension<<" maxk="<<k*kMaxOmega/dimension*180/kPi<<" max="<<max<<endl;
		maxi=i;maxj=j;maxk=k;
		max=point[i][j][k];
		printf(".");
		//printf("Max Omega %f, Max Theta %f, Max Phi %f\n",maxk,maxi,maxj);
	      }
    printf("\n");
    
    maxk=maxk*(kMaxOmega-kMinOmega)/kDimensionOmega + kMinOmega;

    
    //printf("Detected angle for height %3.1f and for center %3.1f %3.1f:%f\n",h,cx,cy,maxk*kPi/(kDimensionTheta*4));
    printf("     Indentified cerenkov angle: %f\n", maxk);
    //printf("Detected angle for height %3.1f and for center %3.1f %3.1f:%f\n",kHeight,cx,cy,maxk);


    //fscanf(omegas,"%f",&realomega);
    //fscanf(thetas,"%f",&realtheta);
    //printf("Real Omega: %f",realomega);			
    //cout<<"Detected:theta="<<maxi*90/kDimensionTheta<<"phi="<<maxj*90/kDimensionPhi<<"omega="<<maxk*kMaxOmega/kDimensionOmega*180/kPi<<" OmegaError="<<fabs(maxk*kMaxOmega/kDimensionOmega*180/kPi-realomega)<<" ThetaError="<<fabs(maxi*90/kDimensionTheta-realtheta)<<endl<<endl;		
    
    //fprintf(results,"Center Coordinates, cx=%6.2f cy=%6.2f, Real Omega=%6.2f, Detected Omega=%6.2f, Omega Error=%6.2f Theta Error=%6.2f\n",cx,cy,realomega,maxk*kMaxOmega/kDimensionOmega*180/kPi,fabs(maxk*kMaxOmega/kDimensionOmega*180/kPi-realomega),fabs(maxi*90/kDimensionTheta-realtheta));
    
    /*for(j=0;j<np;j++)
      pointpp(maxj*90/kDimensionTheta,maxi*90/kDimensionPhi,maxk*kMaxOmega/kDimensionOmega*180/kPi,cx,cy);//Generates a point on the elipse*/		    


    //Start filling rec. hits
    
    Float_t rechit[6];
    
    rechit[0] = (Float_t)( maxi*kPi/(kDimensionTheta*4));
    rechit[1]   = (Float_t)( maxj*kPi/(kDimensionPhi*4));
    rechit[2] = (Float_t)( maxk);
    //rechit[0] = (Float_t)( maxi);
    //rechit[1]   = (Float_t)( maxj);
    //rechit[2] = (Float_t)( maxk);
    rechit[3] = cx;
    rechit[4] = cy;
    rechit[5] = 0.5;
    
    //printf ("track %d, theta %f, phi %f, omega %f\n\n\n",track,rechit[0],rechit[1],rechit[2]);
    
    // fill rechits
    pRICH->AddRecHit3D(nch-1,rechit);
    //printf("Chamber:%d",nch);
  }			
  //printf("\n\n\n\n");
  gAlice->TreeR()->Fill();
  //TTree *TR=gAlice->TreeR();
  //Stat_t ndig=TR->GetEntries();
  TClonesArray *fRec;
  for (i=0;i<kNCH;i++) {
    fRec=pRICH->RecHitsAddress3D(i);
    int ndig=fRec->GetEntriesFast();
    printf ("Chamber %d, rings %d\n",i,ndig);
  }
  //printf("Number of rec. hits: %d",ndig);
  pRICH->ResetRecHits3D();
  //char hname[30];
  //sprintf(hname,"TreeR%d",track);
  //gAlice->TreeR()->Write(hname);
	
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

/*Int_t ***AliRICHDetect::i3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
// allocate a Float_t 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] 
{
long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
Int_t ***t;

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
}*/

/*void pointpp(Float_t alfa,Float_t theta,Float_t omega,Float_t cx,Float_t cy)
  {
  Int_t s;
  Float_t fiducial=h*tan((omega+theta)*kPi/180),l=h/cos(theta*kPi/180),xtrial,y,c0,c1,c2;
  
  //cout<<"fiducial="<<fiducial<<endl;
  
  c0=0;c1=0;c2=0;
  while((c1*c1-4*c2*c0)<=0)
  {	
  //Choose which side to go...
  if(aleat(1)>.5) s=1; else s=-1;
  //Trial a y
  y=s*aleat(fiducial);		
  Float_t alfa1=alfa*kPi/180;
  Float_t theta1=theta*kPi/180;
  Float_t omega1=omega*kPi/180;
  //Solve the eq for a trial x
  c0=-TMath::Power(y*cos(alfa1)*cos(theta1),2)-TMath::Power(y*sin(alfa1),2)+TMath::Power(l*tan(omega1),2)+2*l*y*cos(alfa1)*sin(theta1)*TMath::Power(tan(omega1),2)+TMath::Power(y*cos(alfa1)*sin(theta1)*tan(omega1),2);
  c1=2*y*cos(alfa1)*sin(alfa1)-2*y*cos(alfa1)*TMath::Power(cos(theta1),2)*sin(alfa1)+2*l*sin(alfa1)*sin(theta1)*TMath::Power(tan(omega1),2)+2*y*cos(alfa1)*sin(alfa1)*TMath::Power(sin(theta1),2)*TMath::Power(tan(omega1),2);
  c2=-TMath::Power(cos(alfa1),2)-TMath::Power(cos(theta1)*sin(alfa1),2)+TMath::Power(sin(alfa1)*sin(theta1)*tan(omega1),2);
  //cout<<"Trial: y="<<y<<"c0="<<c0<<" c1="<<c1<<" c2="<<c2<<endl;
  }
  //Choose which side to go...
  if(aleat(1)>.5) s=1; else s=-1;
  xtrial=cx+(-c1+s*sqrt(c1*c1-4*c2*c0))/(2*c2);
  //cout<<"x="<<xtrial<<" y="<<cy+y<<endl;
  fprintf(final,"%f %f\n",xtrial,cy+y);
  }*/






