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


#include "AliRICH.h"
#include "AliRICHPoints.h"
#include "AliRICHDetect.h"
#include "DataStructures.h"
#include "AliRun.h"
#include "TParticle.h"
#include "TMath.h"
#include "TRandom.h"



ClassImp(AliRICHDetect)
//___________________________________________
AliRICHDetect::AliRICHDetect() : TObject()
{
    //fChambers = 0;
}

//___________________________________________
AliRICHDetect::AliRICHDetect(const char *name, const char *title)
    : TObject()
{
    
    /*fChambers = new TObjArray(7);
    for (Int_t i=0; i<7; i++) {
    
	(*fChambers)[i] = new AliRICHchamber();  
	
    } */     
}


void AliRICHDetect::Detect()
{  	
    
  //printf("Detection started!\n");
  Float_t OMEGA,steptheta,stepphi,x,y,cx,cy,l,aux1,aux2,aux3,maxi,maxj,maxk,max;
  //Float_t theta,phi,realomega,realtheta;
  Int_t i,j,k;
  
  //const Float_t Noise_Level=0;          //Noise Level in percentage of mesh points
  //const Float_t t=0.6;			//Softening of Noise Correction (factor)
  
  const Float_t Pii=3.1415927;		
  
  const Float_t h=10;                       //Distance from Radiator to Pads in pads
  
  
  const Int_t dimensiontheta=100;		//Matrix dimension for angle Detection
  const Int_t dimensionphi=100;
  const Int_t dimensionOMEGA=100;
  
  //const Float_t SPOTp=.2;		//Percentage of spot action
  //const Int_t np=500;		//Number of points to reconstruct elipse 
  const Float_t maxOMEGA=65*Pii/180;		//Maximum Cherenkov angle to identify
  
  Int_t Point[dimensiontheta][dimensionphi][dimensionOMEGA];
  //Int_t Point1[dimensiontheta][dimensionphi][dimensionOMEGA];
  
  steptheta=Pii/dimensiontheta;
  stepphi=Pii/dimensionphi;

  AliRICHChamber*       iChamber;
  
  AliRICH *RICH  = (AliRICH*)gAlice->GetDetector("RICH");
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
    
	
  for (Int_t track=0; track<ntracks;track++) {
    gAlice->ResetHits();
    gAlice->TreeH()->GetEvent(track);
    TClonesArray *Hits  = RICH->Hits();
    if (Hits == 0) return;
    Int_t nhits = Hits->GetEntriesFast();
    if (nhits == 0) continue;
    Int_t nent=(Int_t)gAlice->TreeD()->GetEntries();
    gAlice->TreeD()->GetEvent(nent-1);
    AliRICHHit *mHit = 0;
    AliRICHDigit *points = 0;
    //Int_t npoints=0;
    
    Int_t counter=0;
    //Initialization
    for(i=0;i<dimensiontheta;i++)
      {
	for(j=0;j<dimensionphi;j++)
	  {
	    for(k=0;k<dimensionOMEGA;k++)
	      {
		counter++;
		Point[i][j][k]=0;
		//printf("Dimensions theta:%d, phi:%d, omega:%d",dimensiontheta,dimensionphi,dimensionOMEGA);
		//printf("Resetting %d %d %d, time %d\n",i,j,k,counter);
		//-Noise_Level*(Area(i*Pii/(18*dimension),k*maxOMEGA/dimension)-Area((i-1)*Pii/(18*dimension),(k-1)*maxOMEGA/dimension));
		//printf("n-%f",-Noise_Level*(Area(i*Pii/(18*dimension),k*maxOMEGA/dimension)-Area((i-1)*Pii/(18*dimension),(k-1)*maxOMEGA/dimension)));
	      }
	  }
      }
    mHit = (AliRICHHit*) Hits->UncheckedAt(0);
    //printf("Aqui vou eu\n");
    Int_t nch  = mHit->fChamber;
    //printf("Aqui fui eu\n");
    trackglob[0] = mHit->fX;
    trackglob[1] = mHit->fY;
    trackglob[2] = mHit->fZ;

    cx=trackglob[0];
    cy=trackglob[2];
    
    
    //printf("Chamber processed:%d\n",nch);
    printf("Center processed: %3.1f %3.1f %3.1f\n",trackglob[0],trackglob[1],trackglob[2]);

    iChamber = &(RICH->Chamber(nch-1));
    
    //printf("Nch:%d\n",nch);

    iChamber->GlobaltoLocal(trackglob,trackloc);
    
    printf("Transformation 1: %3.1f %3.1f %3.1f\n",trackloc[0],trackloc[1],trackloc[2]);


    iChamber->LocaltoGlobal(trackloc,trackglob);
       
    printf("Transformation 2: %3.1f %3.1f %3.1f\n",trackglob[0],trackglob[1],trackglob[2]);
    
    
     

    TClonesArray *Digits = RICH->DigitsAddress(nch-1);   
    Int_t ndigits = Digits->GetEntriesFast();
    
    //printf("Got %d digits\n",ndigits);

    //printf("Starting calculations\n");
    
    for(Float_t theta=0;theta<Pii/18;theta+=steptheta)
      {			
	for(Float_t phi=0;phi<=Pii/3;phi+=stepphi)
	  {		       
	    for (Int_t dig=0;dig<ndigits;dig++)
	      {	
		points=(AliRICHDigit*) Digits->UncheckedAt(dig);
		
		x=points->fPadX-cx;
		y=points->fPadY-cy;
		//printf("Loaded digit %d with coordinates x:%f, y%f\n",dig,x,y);
		//cout<<"x="<<x<<" y="<<y<<endl;
		
		if (sqrt(pow(x,2)+pow(y,2))<h*tan(theta+maxOMEGA)*3/4)
		  {
		    
		    l=h/cos(theta);
		    
		    aux1=-y*sin(phi)+x*cos(phi);
		    aux2=y*cos(phi)+x*sin(phi);
		    aux3=( pow(aux1,2)+pow(cos(theta)*aux2 ,2))/pow(sin(theta)*aux2+l,2);
		    //cout<<"aux1="<<aux1<<" aux2="<<aux2<<" aux3="<<aux3;
		    
		    OMEGA=atan(sqrt(aux3));
		    //printf("Omega: %f\n",OMEGA);
		    
		    //cout<<"\ni="<<i<<" theta="<<Int_t(2*theta*dimension/Pii)<<" phi="<<Int_t(2*phi*dimension/Pii)<<" OMEGA="<<Int_t(2*OMEGA*dimension/Pii)<<endl<<endl;
		    //{Int_t lixo;cin>>lixo;}
		    if(OMEGA<maxOMEGA)Point[Int_t(2*theta*dimensiontheta/Pii)][Int_t(2*phi*dimensionphi/Pii)][Int_t(OMEGA*dimensionOMEGA/maxOMEGA)]+=1;	
		    //if(OMEGA<maxOMEGA)Point[Int_t(theta)][Int_t(phi)][Int_t(OMEGA)]+=1;
		  }
		}
	  }
      }	
    
    
    
    //SPOT execute twice
    /*for(s=1;i<=2;s++)
      {
	//buffer copy
	for(i=0;i<=dimensiontheta;i++)
	  for(j=0;j<=dimensionphi;j++)
	    for(k=0;k<=dimensionOMEGA;k++)
	      Point1[i][j][k]=Point[i][j][k];	
	
	cout<<"COM SPOT!"<<endl;{Int_t lixo;cin>>lixo;}					
	//SPOT algorithm			
	for(i=1;i<dimensiontheta;i++)
	  for(j=1;j<dimensionphi;j++)
	    for(k=1;k<dimensionOMEGA;k++)
	      {
		if((Point[i][k][j]>Point[i-1][k][j])&&(Point[i][k][j]>Point[i+1][k][j])&&
		   (Point[i][k][j]>Point[i][k-1][j])&&(Point[i][k][j]>Point[i][k+1][j])&&
		   (Point[i][k][j]>Point[i][k][j-1])&&(Point[i][k][j]>Point[i][k][j+1]))
		  {
		    //cout<<"SPOT"<<endl;
		    //Execute SPOT on point											   	
		    Point1[i][j][k]+=int(SPOTp*(Point[i-1][k][j]+Point[i+1][k][j]+Point[i][k-1][j]+Point[i][k+1][j]+Point[i][k][j-1]+Point[i][k][j+1]));    
		    Point1[i-1][k][j]=int(SPOTp*Point[i-1][k][j]);
		    Point1[i+1][k][j]=Int_t(SPOTp*Point[i+1][k][j]);
		    Point1[i][k-1][j]=Int_t(SPOTp*Point[i][k-1][j]);
		    Point1[i][k+1][j]=Int_t(SPOTp*Point[i][k+1][j]);
		    Point1[i][k][j-1]=Int_t(SPOTp*Point[i][k][j-1]);
		    Point1[i][k][j+1]=Int_t(SPOTp*Point[i][k][j+1]);
		  }
	      }
	//copy from buffer copy
	for(i=1;i<dimensiontheta;i++)
	  for(j=1;j<dimensionphi;j++)
	    for(k=1;k<dimensionOMEGA;k++)
	      Point[i][j][k]=Point1[i][j][k];										
	  
	  }*/
    
    
    //Identification is equivalent to maximum determination
    max=0;maxi=0;maxj=0;maxk=0;
    
    //cout<<"Proceeding to Identification"<<endl;
    
    for(i=1;i<dimensiontheta-3;i++)
      for(j=1;j<=dimensionphi-3;j++)
	for(k=0;k<=dimensionOMEGA;k++)
	  if(Point[i][j][k]>max)
	    {
	      //cout<<"maxi="<<i*90/dimension<<" maxj="<<j*90/dimension<<" maxk="<<k*maxOMEGA/dimension*180/Pii<<" max="<<max<<endl;
	      maxi=i;maxj=j;maxk=k;
	      max=Point[i][j][k];
	      //printf("Max Omega %f, Max Theta %f, Max Phi %f\n",maxk,maxi,maxj);
	    }
    
    //printf("Detected angle for height %3.1f and for center %3.1f %3.1f:%f\n",h,cx,cy,maxk*Pii/(dimensiontheta*4));
    //printf("Detected angle for height %3.1f and for center %3.1f %3.1f:%f\n",h,cx,cy,maxk);


    //fscanf(omegas,"%f",&realomega);
    //fscanf(thetas,"%f",&realtheta);
    //printf("Real Omega: %f",realomega);			
    //cout<<"Detected:theta="<<maxi*90/dimensiontheta<<"phi="<<maxj*90/dimensionphi<<"OMEGA="<<maxk*maxOMEGA/dimensionOMEGA*180/Pii<<" OmegaError="<<fabs(maxk*maxOMEGA/dimensionOMEGA*180/Pii-realomega)<<" ThetaError="<<fabs(maxi*90/dimensiontheta-realtheta)<<endl<<endl;		
    
    //fprintf(results,"Center Coordinates, cx=%6.2f cy=%6.2f, Real Omega=%6.2f, Detected Omega=%6.2f, Omega Error=%6.2f Theta Error=%6.2f\n",cx,cy,realomega,maxk*maxOMEGA/dimensionOMEGA*180/Pii,fabs(maxk*maxOMEGA/dimensionOMEGA*180/Pii-realomega),fabs(maxi*90/dimensiontheta-realtheta));
    
    /*for(j=0;j<np;j++)
      Pointpp(maxj*90/dimensiontheta,maxi*90/dimensionphi,maxk*maxOMEGA/dimensionOMEGA*180/Pii,cx,cy);//Generates a point on the elipse*/		    


    //Start filling rec. hits
    
    Float_t rechit[5];
    
    rechit[0] = (Float_t)( maxi*Pii/(dimensiontheta*4));
    rechit[1]   = (Float_t)( maxj*Pii/(dimensionphi*4));
    rechit[2] = (Float_t)( maxk*Pii/(dimensionOMEGA*4));
    //rechit[0] = (Float_t)( maxi);
    //rechit[1]   = (Float_t)( maxj);
    //rechit[2] = (Float_t)( maxk);
    rechit[3] = cx;
    rechit[4] = cy;
    
    //printf ("track %d, theta %f, phi %f, omega %f\n\n\n",track,rechit[0],rechit[1],rechit[2]);
    
    // fill rechits
    RICH->AddRecHit(nch-1,rechit);
  }			
  //printf("\n\n\n\n");
  gAlice->TreeR()->Fill();
  //TTree *TR=gAlice->TreeR();
  //Stat_t ndig=TR->GetEntries();
  TClonesArray *fRec;
  for (i=0;i<7;i++) {
    fRec=RICH->RecHitsAddress(i);
    int ndig=fRec->GetEntriesFast();
    printf ("Chamber %d, rings %d\n",i,ndig);
  }
  //printf("Number of rec. hits: %d",ndig);
  RICH->ResetRecHits();
  //char hname[30];
  //sprintf(hname,"TreeR%d",track);
  //gAlice->TreeR()->Write(hname);
	
}

Float_t AliRICHDetect:: Area(Float_t theta,Float_t OMEGA)
{
    
    Float_t area;
    const Float_t h=9.25;                       //Distance from Radiator to Pads in pads
    
    area=TMath::Pi()*pow(h*tan(OMEGA),2)/pow(pow(cos(theta),2)-pow(tan(OMEGA)*sin(theta),2),3/2);
    
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

/*void Pointpp(Float_t alfa,Float_t theta,Float_t OMEGA,Float_t cx,Float_t cy)
  {
  Int_t s;
  Float_t fiducial=h*tan((OMEGA+theta)*Pii/180),l=h/cos(theta*Pii/180),xtrial,y,c0,c1,c2;
  
  //cout<<"fiducial="<<fiducial<<endl;
  
  c0=0;c1=0;c2=0;
  while((c1*c1-4*c2*c0)<=0)
  {	
  //Choose which side to go...
  if(aleat(1)>.5) s=1; else s=-1;
  //Trial a y
  y=s*aleat(fiducial);		
  Float_t alfa1=alfa*Pii/180;
  Float_t theta1=theta*Pii/180;
  Float_t OMEGA1=OMEGA*Pii/180;
  //Solve the eq for a trial x
  c0=-pow(y*cos(alfa1)*cos(theta1),2)-pow(y*sin(alfa1),2)+pow(l*tan(OMEGA1),2)+2*l*y*cos(alfa1)*sin(theta1)*pow(tan(OMEGA1),2)+pow(y*cos(alfa1)*sin(theta1)*tan(OMEGA1),2);
  c1=2*y*cos(alfa1)*sin(alfa1)-2*y*cos(alfa1)*pow(cos(theta1),2)*sin(alfa1)+2*l*sin(alfa1)*sin(theta1)*pow(tan(OMEGA1),2)+2*y*cos(alfa1)*sin(alfa1)*pow(sin(theta1),2)*pow(tan(OMEGA1),2);
  c2=-pow(cos(alfa1),2)-pow(cos(theta1)*sin(alfa1),2)+pow(sin(alfa1)*sin(theta1)*tan(OMEGA1),2);
  //cout<<"Trial: y="<<y<<"c0="<<c0<<" c1="<<c1<<" c2="<<c2<<endl;
  }
  //Choose which side to go...
  if(aleat(1)>.5) s=1; else s=-1;
  xtrial=cx+(-c1+s*sqrt(c1*c1-4*c2*c0))/(2*c2);
  //cout<<"x="<<xtrial<<" y="<<cy+y<<endl;
  fprintf(final,"%f %f\n",xtrial,cy+y);
  }*/






