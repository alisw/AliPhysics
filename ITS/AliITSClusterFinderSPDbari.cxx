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


#include "AliITSClusterFinderSPDbari.h"
#include "AliITS.h"
#include "AliITSgeom.h"
#include "AliITSdigit.h"
#include "AliITSRawCluster.h"
#include "AliITSRecPoint.h"
#include "AliITSsegmentation.h"
#include "AliITSresponse.h"
#include "AliRun.h"




ClassImp(AliITSClusterFinderSPDbari)

//----------------------------------------------------------
AliITSClusterFinderSPDbari::AliITSClusterFinderSPDbari
(AliITSsegmentation *seg, TClonesArray *digits, TClonesArray *recp)   
{
  // constructor
    fSegmentation=seg;
    fDigits=digits;
    fClusters=recp;
    fNclusters= fClusters->GetEntriesFast();
    SetDx();
    SetDz();
}

//_____________________________________________________________________________
AliITSClusterFinderSPDbari::AliITSClusterFinderSPDbari()
{
  // constructor
  fSegmentation=0;
  fDigits=0;
  fClusters=0;
  fNclusters=0;
  SetDx();
  SetDz();
  
}

//__________________________________________________________________________
AliITSClusterFinderSPDbari::AliITSClusterFinderSPDbari(const
AliITSClusterFinderSPDbari &source){
  //     Copy Constructor 
  if(&source == this) return;
  this->fClusters = source.fClusters ;
  this->fNclusters = source.fNclusters ;
  this->fDz = source.fDz ;
  this->fDx = source.fDx ;
  return;
}

//_________________________________________________________________________
AliITSClusterFinderSPDbari& 
  AliITSClusterFinderSPDbari::operator=(const AliITSClusterFinderSPDbari &source) {
  //    Assignment operator
  if(&source == this) return *this;
  this->fClusters = source.fClusters ;
  this->fNclusters = source.fNclusters ;
  this->fDz = source.fDz ;
  this->fDx = source.fDx ;
  return *this;
}

//_____________________________________________________________________________
void AliITSClusterFinderSPDbari::FindRawClusters(){
   
    // input of Cluster Finder  
    Int_t digitcount=0;
    Int_t numberd=10000;
    Int_t   *digx       = new Int_t[numberd];
    Int_t   *digz       = new Int_t[numberd];
    Int_t   *digtr1     = new Int_t[numberd];
    Int_t   *digtr2     = new Int_t[numberd];
    Int_t   *digtr3     = new Int_t[numberd];
    
    //  output of Cluster Finder    
    Int_t numberc=1000;
    Float_t *xcenterl   = new Float_t[numberc];
    Float_t *zcenterl   = new Float_t[numberc];
    Float_t *errxcenter = new Float_t[numberc];
    Float_t *errzcenter = new Float_t[numberc];
    Int_t   *tr1clus    = new Int_t[numberc];
    Int_t   *tr2clus    = new Int_t[numberc];
    Int_t   *tr3clus    = new Int_t[numberc];

    Int_t nclus;

    digitcount=0;
    Int_t ndigits = fDigits->GetEntriesFast();  
    if (!ndigits) return;


    AliITSdigit *dig;
    AliITSdigitSPD *dig1;
    Int_t ndig;
    for(ndig=0; ndig<ndigits; ndig++) {
         dig= (AliITSdigit*)fDigits->UncheckedAt(ndig);
     
         digx[digitcount] = dig->fCoord2+1;  //starts at 1
         digz[digitcount] = dig->fCoord1+1;  //starts at 1

         dig1= (AliITSdigitSPD*)fDigits->UncheckedAt(ndig);

         digtr1[digitcount] = dig1->fTracks[0];
         digtr2[digitcount] = dig1->fTracks[1];
         digtr3[digitcount] = dig1->fTracks[2];

         digitcount++;
    }


        ClusterFinder(digitcount,digx,digz,digtr1,digtr2,digtr3,
              nclus,xcenterl,zcenterl,errxcenter,errzcenter,
              tr1clus, tr2clus, tr3clus);
 
        DigitToPoint(nclus,xcenterl,zcenterl,errxcenter,errzcenter,
              tr1clus, tr2clus, tr3clus);


  delete[] digx       ;
  delete[] digz       ;
  delete[] digtr1     ;
  delete[] digtr2     ;
  delete[] digtr3     ;
  delete[] xcenterl   ;
  delete[] zcenterl   ;
  delete[] errxcenter ;
  delete[] errzcenter ;
  delete[] tr1clus    ;
  delete[] tr2clus    ;
  delete[] tr3clus    ;
  
}
//-----------------------------------------------------------------
void AliITSClusterFinderSPDbari::ClusterFinder(Int_t ndigits,
                    Int_t digx[],Int_t digz[],
                    Int_t digtr1[],Int_t digtr2[],Int_t digtr3[],
                    Int_t &nclus, Float_t xcenter[],Float_t zcenter[],
				    Float_t errxcenter[],Float_t errzcenter[],
                    Int_t tr1clus[], Int_t tr2clus[], Int_t tr3clus[]) {
//
// Search for clusters of fired pixels (digits). Two digits are linked
// inside a cluster if they are countiguous both in row or column
// direction.  Diagonal digits are not linked.
// xcenter, ycenter, zcenter are the coordinates of the center
// of each found cluster, calculated from the averaging the corresponding
// coordinate of the center of the linked digits. The coordinates are
// given in the local reference sistem. 
// errxcenter, errycenter, errzcenter are the errors associated to
// the corresponding average.
//
//

  Int_t if1, min, max, nd;
  Int_t x1, z1, t1, t2, t3;
  Int_t ndx, ndz, ndxmin, ndxmax, ndzmin, ndzmax;
  Float_t dx, dz; 
  Int_t i,k,ipos=0;
  Float_t xdum, zdum;      

  Int_t numberd=10000;
  Int_t *ifpad = new Int_t[numberd];
  Int_t *xpad  = new Int_t[numberd];
  Int_t *zpad  = new Int_t[numberd];
  Int_t *tr1pad  = new Int_t[numberd];
  Int_t *tr2pad  = new Int_t[numberd];
  Int_t *tr3pad  = new Int_t[numberd];
  Int_t *iclus   = new Int_t[numberd];

 
  nclus=1;
  for (i=0; i < ndigits ; i++) ifpad[i] = -1;

  ifpad[0]=0;
  for (i=0; i < ndigits-1 ; i++) 
  {
    if ( ifpad[i] == -1 ) 
    { 
	   	nclus++;
		ipos++;
     	ifpad[i]=nclus-1;
    }
    for (Int_t j=i+1 ; j < ndigits ; j++)  
    {  
      if (ifpad[j]== -1 )
      {
   	     dx = TMath::Abs(digx[i]-digx[j]);
   	     dz = TMath::Abs(digz[i]-digz[j]);

// 	     if ( ( dx+dz )==1 )  //  clusters are not diagonal
   	     if ( ( dx+dz )==1 || (dx==1 && dz==1) )  //  diagonal clusters allowed
   	     {
   		    ipos++;
		    ifpad[j]=ifpad[i];

   		    x1=digx[j];
   		    z1=digz[j];
   		    digx[j]=digx[ipos];
   		    digz[j]=digz[ipos];
   		    digx[ipos]=x1;
   		    digz[ipos]=z1;

   		    t1=digtr1[j];
   		    t2=digtr2[j];
   		    t3=digtr3[j];
   		    digtr1[j]=digtr1[ipos];
   		    digtr2[j]=digtr2[ipos];
   		    digtr3[j]=digtr3[ipos];
   		    digtr1[ipos]=t1;
   		    digtr2[ipos]=t2;
   		    digtr3[ipos]=t3;

   		    if1=ifpad[j];
   		    ifpad[j]=ifpad[ipos];
   		    ifpad[ipos]=if1;
	     }
      }
    }
   }   
   if ( ifpad[ndigits-1] == -1 )
   {
	  nclus++;
	  ifpad[ndigits-1]=nclus-1;
   }
   for (i=0 ; i < ndigits ; i++) iclus[ifpad[i]]++;

   min=0;
   max=0;

   // loop on found clusters 

   for (i=0 ; i < nclus ; i++)  
   {
      min = max;
      max += iclus[i];
      Float_t deltax = fSegmentation->Dpx(0);
      if (iclus[i]!=1) 
      {
        //cluster with more than one digit
        nd=iclus[i];
	    Int_t count=0;
        for (k=min;k<min+nd;k++)
        {
	       xpad[count] = digx[k];	   
	       zpad[count] = digz[k];

	       tr1pad[count] = digtr1[k];	   
	       tr2pad[count] = digtr2[k];	   
	       tr3pad[count] = digtr3[k];	   

	       count++; 
        }
        ndxmin = xpad[TMath::LocMin(nd,xpad)];
        ndxmax = xpad[TMath::LocMax(nd,xpad)];
        ndzmin = zpad[TMath::LocMin(nd,zpad)];
        ndzmax = zpad[TMath::LocMax(nd,zpad)];
        ndx = ndxmax - ndxmin+1;
        ndz = ndzmax - ndzmin+1;

        // calculate x and z coordinates of the center of the cluster
        fSegmentation->GetPadCxz(digx[min],digz[min]-1,xdum, zdum);


	if (ndx == 1) {	    
	    xcenter[i] = xdum;
	}    
	else{ 
           xcenter[i] = 0.;
	   for (k=0;k<nd;k++) {
             fSegmentation->GetPadCxz(xpad[k],zpad[k]-1,xdum,zdum);
	     xcenter[i] += (xdum / nd);
	   }	               
	}

	if (ndz == 1) {
	    zcenter[i] = zdum;
	}   
	else {
	   zcenter[i] = 0.;
	   for (k=0;k<nd;k++) {	      
             fSegmentation->GetPadCxz(xpad[k],zpad[k]-1,xdum,zdum);
	     zcenter[i] += (zdum / nd);
	   }
	}

        // error on points in x and z directions

        if (ndx == 1) {
             errxcenter[i] = deltax / TMath::Sqrt(12.);
        }
        else {
             errxcenter[i] = 0.;	 		
             for (k=0;k<nd;k++){ 
               fSegmentation->GetPadCxz(xpad[k],zpad[k]-1,xdum,zdum);
               errxcenter[i] += ((xdum-xcenter[i])*(xdum-xcenter[i]))/(nd*(nd-1)); 
             }   
	     errxcenter[i] = TMath::Sqrt(errxcenter[i]);
        }
	
	if (ndz == 1) {
            Float_t deltaz = fSegmentation->Dpz(digz[min]);	              
	    errzcenter[i] = deltaz / TMath::Sqrt(12.);
        }
	else {
	     errzcenter[i] = 0.;
	     for (k=0;k<nd;k++){ 
               fSegmentation->GetPadCxz(xpad[k],zpad[k]-1,xdum,zdum);
	       errzcenter[i] += ((zdum-zcenter[i])*(zdum-zcenter[i]))/(nd*(nd-1));
	     }
	     errzcenter[i] = TMath::Sqrt(errzcenter[i]);
	}    

        // take three track numbers for the cluster
        for (k=0;k<nd;k++){
          if(tr1pad[k] != -2) tr1clus[i]=tr1pad[k];
          if(tr2pad[k] != -2) tr2clus[i]=tr2pad[k];
          if(tr3pad[k] != -2) tr3clus[i]=tr3pad[k];
        }
        if(tr1clus[i] == 0) tr1clus[i]= -2;
        if(tr2clus[i] == 0) tr2clus[i]= -2;
        if(tr3clus[i] == 0) tr3clus[i]= -2;
      }
      else  {
      
        // cluster with single digit
	ndx = 1;
        ndz = 1;
        fSegmentation->GetPadCxz(digx[min],digz[min]-1,xdum,zdum);
        xcenter[i] = xdum;
        zcenter[i] = zdum;
        tr1clus[i]=digtr1[min];
        tr2clus[i]=digtr2[min];
        tr3clus[i]=digtr3[min];
	    Float_t deltaz = fSegmentation->Dpz(digz[min]);
	    errxcenter[i] = deltax / TMath::Sqrt(12.);
  	    errzcenter[i] = deltaz / TMath::Sqrt(12.);
     }

     // store the cluster information to the AliITSRawCLusterSPD object
     AliITS *iTS=(AliITS*)gAlice->GetModule("ITS");

     //put the cluster center in local reference frame of the detector
     // and in microns
     xcenter[i] = xcenter[i] - fSegmentation->Dx()/2.; 
     zcenter[i] = zcenter[i] - fSegmentation->Dz()/2.;

     AliITSRawClusterSPD *clust = new AliITSRawClusterSPD(zcenter[i],xcenter[i],1.,ndz,ndx,0.,0.,0.,0.,0.,0.,0.);
     iTS->AddCluster(0,clust);
   }	  
   delete[] ifpad;
   delete[] xpad ;
   delete[] zpad ;
   delete[] iclus;
}
//______________________________________________________
void AliITSClusterFinderSPDbari::DigitToPoint(Int_t nclus,
                   Float_t *xcenter,   Float_t *zcenter,
		   Float_t *errxcenter,Float_t *errzcenter, 
		   Int_t *tr1clus, Int_t *tr2clus, Int_t *tr3clus){
 //
 // A point is associated to each cluster of SPD digits. The points
 // and their associated errors are stored in the file galiceSP.root.
 //
 //
 
     Float_t l[3],xg,zg;
     const Float_t kconv = 1.0e-4; // micron -> cm

     // get rec points
     AliITS *iTS=(AliITS*)gAlice->GetModule("ITS");

     for (Int_t i=0; i<nclus; i++)
     {
        l[0] = kconv*xcenter[i];
        l[1] = kconv*fSegmentation->Dy()/2.;
        l[2] = kconv*zcenter[i];

        xg = l[0]; 
        zg = l[2]; 

	    Float_t sigma2x = (kconv*errxcenter[i]) * (kconv*errxcenter[i]);
	    Float_t sigma2z = (kconv*errzcenter[i]) * (kconv*errzcenter[i]);
        AliITSRecPoint rnew;
        rnew.SetX(xg);
        rnew.SetZ(zg);
        rnew.SetQ(1.);
        rnew.SetdEdX(0.);
        rnew.SetSigmaX2(sigma2x);
        rnew.SetSigmaZ2(sigma2z);
        rnew.fTracks[0]=tr1clus[i];
        rnew.fTracks[1]=tr2clus[i];
        rnew.fTracks[2]=tr3clus[i];
        iTS->AddRecPoint(rnew); 
     }

}
