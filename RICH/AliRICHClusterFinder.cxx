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
  Revision 1.9  2001/01/26 20:00:27  hristov
  Major upgrade of AliRoot code

  Revision 1.8  2000/11/02 09:11:12  jbarbosa
  Removed AliRICHRecHit.h from include.

  Revision 1.7  2000/10/03 21:44:09  morsch
  Use AliSegmentation and AliHit abstract base classes.

  Revision 1.6  2000/10/02 21:28:12  fca
  Removal of useless dependecies via forward declarations

  Revision 1.5  2000/10/02 15:45:58  jbarbosa
  Fixed forward declarations.

  Revision 1.4  2000/06/12 19:01:29  morsch
  Clean-up bug in Centered() corrected.

  Revision 1.3  2000/06/12 15:49:44  jbarbosa
  Removed verbose output.

  Revision 1.2  2000/06/12 15:18:19  jbarbosa
  Cleaned up version.

  Revision 1.1  2000/04/19 13:01:48  morsch
  A cluster finder and hit reconstruction class for RICH (adapted from MUON).
  Cluster Finders for MUON and RICH should derive from the same class in the
  future (JB, AM).

*/


#include "AliRICHClusterFinder.h"
#include "AliRun.h"
#include "AliRICH.h"
#include "AliRICHHit.h"
#include "AliRICHHitMapA1.h"
#include "AliRICHCerenkov.h"
#include "AliRICHSDigit.h"
#include "AliRICHDigit.h"
#include "AliRICHRawCluster.h"

#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TPad.h>
#include <TGraph.h> 
#include <TPostScript.h> 
#include <TMinuit.h> 

//----------------------------------------------------------
static AliSegmentation* gSegmentation;
static AliRICHResponse*     gResponse;
static Int_t                gix[500];
static Int_t                giy[500];
static Float_t              gCharge[500];
static Int_t                gNbins;
static Int_t                gFirst=kTRUE;
static TMinuit *gMyMinuit ;
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
static Int_t                gChargeTot;

ClassImp(AliRICHClusterFinder)

AliRICHClusterFinder::AliRICHClusterFinder
(AliSegmentation *segmentation, AliRICHResponse *response, 
 TClonesArray *digits, Int_t chamber)   
{

// Constructor for Cluster Finder object

    fSegmentation=segmentation;
    fResponse=response;
    
    fDigits=digits;
    fNdigits = fDigits->GetEntriesFast();
    fChamber=chamber;
    fRawClusters=new TClonesArray("AliRICHRawCluster",10000);
    fNRawClusters=0;
    fCogCorr = 0;
    SetNperMax();
    SetClusterSize();
    SetDeclusterFlag();
    fNPeaks=-1;
}

AliRICHClusterFinder::AliRICHClusterFinder()
{

// Default constructor

    fSegmentation=0;
    fResponse=0;
    
    fDigits=0;
    fNdigits = 0;
    fChamber=-1;
    fRawClusters=new TClonesArray("AliRICHRawCluster",10000);
    fNRawClusters=0;
    fHitMap = 0;
    fCogCorr = 0;
    SetNperMax();
    SetClusterSize();
    SetDeclusterFlag();
    fNPeaks=-1;
}

AliRICHClusterFinder::AliRICHClusterFinder(const AliRICHClusterFinder& ClusterFinder)
{
// Copy Constructor
}

AliRICHClusterFinder::~AliRICHClusterFinder()
{

// Destructor

delete fRawClusters;
}

void AliRICHClusterFinder::AddRawCluster(const AliRICHRawCluster c)
{
  //
  // Add a raw cluster copy to the list
  //
    AliRICH *pRICH=(AliRICH*)gAlice->GetModule("RICH");
    pRICH->AddRawCluster(fChamber,c); 
    fNRawClusters++;
}



void AliRICHClusterFinder::Decluster(AliRICHRawCluster *cluster)
{

//
// Decluster algorithm
    
    Int_t mul = cluster->fMultiplicity;
//    printf("Decluster - multiplicity   %d \n",mul);

    if (mul == 1 || mul ==2) {
//
// Nothing special for 1- and 2-clusters
	if (fNPeaks != 0) {
            cluster->fNcluster[0]=fNPeaks;
            cluster->fNcluster[1]=0;
        } 
	AddRawCluster(*cluster); 
        fNPeaks++;
    } else if (mul ==3) {
//
// 3-cluster, check topology
//	printf("\n 3-cluster, check topology \n");
	if (fDeclusterFlag) {
	    if (Centered(cluster)) {
		// ok, cluster is centered 
	    } else {
		// cluster is not centered, split into 2+1
	    }
	} else {
	    if (fNPeaks != 0) {
		cluster->fNcluster[0]=fNPeaks;
		cluster->fNcluster[1]=0;
	    } 
	    AddRawCluster(*cluster); 
	    fNPeaks++;
	}	    
    } else {
// 
// 4-and more-pad clusters
//
	if (mul <= fClusterSize) {
	    if (fDeclusterFlag) {
		SplitByLocalMaxima(cluster);
	    } else {
		if (fNPeaks != 0) {
		    cluster->fNcluster[0]=fNPeaks;
		    cluster->fNcluster[1]=0;
		} 
		AddRawCluster(*cluster);
		fNPeaks++;
	    }	
	}
    } // multiplicity 
}


Bool_t AliRICHClusterFinder::Centered(AliRICHRawCluster *cluster)
{

// Is the cluster centered?

    AliRICHDigit* dig;
    dig= (AliRICHDigit*)fDigits->UncheckedAt(cluster->fIndexMap[0]);
    Int_t ix=dig->fPadX;
    Int_t iy=dig->fPadY;
    Int_t nn;
    Int_t x[kMaxNeighbours], y[kMaxNeighbours], xN[kMaxNeighbours], yN[kMaxNeighbours];
    
    fSegmentation->Neighbours(ix,iy,&nn,x,y);
    Int_t nd=0;
    for (Int_t i=0; i<nn; i++) {
	if (fHitMap->TestHit(x[i],y[i]) == kUsed) {
	    xN[nd]=x[i];
	    yN[nd]=y[i];
	    nd++;
	    
	    //printf("Getting: %d %d %d\n",i,x[i],y[i]);
	}
    }
    if (nd==2) {
//
// cluster is centered !
	if (fNPeaks != 0) {
            cluster->fNcluster[0]=fNPeaks;
            cluster->fNcluster[1]=0;
        }  
	cluster->fCtype=0;
	AddRawCluster(*cluster);
	fNPeaks++;
	return kTRUE;
    } else if (nd ==1) {
//
// Highest signal on an edge, split cluster into 2+1
//
// who is the neighbour ?
      
      //printf("Calling GetIndex with x:%d y:%d\n",xN[0], yN[0]);
      
	Int_t nind=fHitMap->GetHitIndex(xN[0], yN[0]);
	Int_t i1= (nind==cluster->fIndexMap[1]) ? 1:2;
	Int_t i2= (nind==cluster->fIndexMap[1]) ? 2:1;    
//
// 2-cluster
	AliRICHRawCluster cnew;
	if (fNPeaks == 0) {
            cnew.fNcluster[0]=-1;
            cnew.fNcluster[1]=fNRawClusters;
        } else {
            cnew.fNcluster[0]=fNPeaks;
            cnew.fNcluster[1]=0;
        }
	cnew.fMultiplicity=2;
	cnew.fIndexMap[0]=cluster->fIndexMap[0];
	cnew.fIndexMap[1]=cluster->fIndexMap[i1];
	FillCluster(&cnew);
	cnew.fClusterType=cnew.PhysicsContribution();
	AddRawCluster(cnew);
        fNPeaks++;
//
// 1-cluster
	cluster->fMultiplicity=1;
	cluster->fIndexMap[0]=cluster->fIndexMap[i2];
	cluster->fIndexMap[1]=0;
	cluster->fIndexMap[2]=0;	
	FillCluster(cluster);
        if (fNPeaks != 0) {
            cluster->fNcluster[0]=fNPeaks;
            cluster->fNcluster[1]=0;
        }  
	cluster->fClusterType=cluster->PhysicsContribution();
	AddRawCluster(*cluster);
	fNPeaks++;
	return kFALSE;
    } else {
	printf("\n Completely screwed up %d !! \n",nd);
	
    }
    
	return kFALSE;
}
void AliRICHClusterFinder::SplitByLocalMaxima(AliRICHRawCluster *c)
{

//
// Split the cluster according to the number of maxima inside


    AliRICHDigit* dig[100], *digt;
    Int_t ix[100], iy[100], q[100];
    Float_t x[100], y[100], zdum;
    Int_t i; // loops over digits
    Int_t j; // loops over local maxima
    //    Float_t xPeak[2];
    //    Float_t yPeak[2];
    //    Int_t threshold=500;
    Int_t mul=c->fMultiplicity;
//
//  dump digit information into arrays
//
    for (i=0; i<mul; i++)
    {
	dig[i]= (AliRICHDigit*)fDigits->UncheckedAt(c->fIndexMap[i]);
	ix[i]= dig[i]->fPadX;
	iy[i]= dig[i]->fPadY;
	q[i] = dig[i]->fSignal;
	fSegmentation->GetPadC(ix[i], iy[i], x[i], y[i], zdum);
    }
//
//  Find local maxima
//
    Bool_t isLocal[100];
    Int_t nLocal=0;
    Int_t associatePeak[100];
    Int_t indLocal[100];
    Int_t nn;
    Int_t xNei[kMaxNeighbours], yNei[kMaxNeighbours];
    for (i=0; i<mul; i++) {
	fSegmentation->Neighbours(ix[i], iy[i], &nn, xNei, yNei);
	isLocal[i]=kTRUE;
	for (j=0; j<nn; j++) {
	    if (fHitMap->TestHit(xNei[j], yNei[j])==kEmpty) continue;
	    digt=(AliRICHDigit*) fHitMap->GetHit(xNei[j], yNei[j]);
	    if (digt->fSignal > q[i]) {
		isLocal[i]=kFALSE;
		break;
//
// handle special case of neighbouring pads with equal signal
	    } else if (digt->fSignal == q[i]) {
		if (nLocal >0) {
		    for (Int_t k=0; k<nLocal; k++) {
			if (xNei[j]==ix[indLocal[k]] && yNei[j]==iy[indLocal[k]]){
			    isLocal[i]=kFALSE;
			}
		    }
		}
	    } 
	} // loop over next neighbours
	// Maxima should not be on the edge
	if (isLocal[i]) {
	    indLocal[nLocal]=i;
	    nLocal++;
	} 
    } // loop over all digits
//    printf("Found %d local Maxima",nLocal);
//
// If only one local maximum found but multiplicity is high 
// take global maximum from the list of digits.    
    if (nLocal==1 && mul>5) {
	Int_t nnew=0;
	for (i=0; i<mul; i++) {
	    if (!isLocal[i]) {
		indLocal[nLocal]=i;
		isLocal[i]=kTRUE;
		nLocal++;
		nnew++;
	    }
	    if (nnew==1) break;
	}
    }
    
// If number of local maxima is 2 try to fit a double gaussian
    if (nLocal==-100) {
//
//  Initialise global variables for fit
	gFirst=1;
	gSegmentation=fSegmentation;
	gResponse    =fResponse;
	gNbins=mul;
	
	for (i=0; i<mul; i++) {
	    gix[i]=ix[i];
	    giy[i]=iy[i];
	    gCharge[i]=Float_t(q[i]);
	}
//
	if (gFirst) {
	    gFirst=kFALSE;
	    gMyMinuit = new TMinuit(5);
	}
	gMyMinuit->SetFCN(fcn);
	gMyMinuit->mninit(5,10,7);
	Double_t arglist[20];
	Int_t ierflag=0;
	arglist[0]=1;
//	gMyMinuit->mnexcm("SET ERR",arglist,1,ierflag);
// Set starting values 
	static Double_t vstart[5];
	vstart[0]=x[indLocal[0]];
	vstart[1]=y[indLocal[0]];	
	vstart[2]=x[indLocal[1]];
	vstart[3]=y[indLocal[1]];	
	vstart[4]=Float_t(q[indLocal[0]])/
	    Float_t(q[indLocal[0]]+q[indLocal[1]]);
// lower and upper limits
	static Double_t lower[5], upper[5];
	Int_t isec=fSegmentation->Sector(ix[indLocal[0]], iy[indLocal[0]]);
	lower[0]=vstart[0]-fSegmentation->Dpx(isec)/2;
	lower[1]=vstart[1]-fSegmentation->Dpy(isec)/2;
//	lower[1]=vstart[1];
	
	upper[0]=lower[0]+fSegmentation->Dpx(isec);
	upper[1]=lower[1]+fSegmentation->Dpy(isec);
//	upper[1]=vstart[1];
	
	isec=fSegmentation->Sector(ix[indLocal[1]], iy[indLocal[1]]);
	lower[2]=vstart[2]-fSegmentation->Dpx(isec)/2;
	lower[3]=vstart[3]-fSegmentation->Dpy(isec)/2;
//	lower[3]=vstart[3];
	
	upper[2]=lower[2]+fSegmentation->Dpx(isec);
	upper[3]=lower[3]+fSegmentation->Dpy(isec);
//	upper[3]=vstart[3];
	
	lower[4]=0.;
	upper[4]=1.;
// step sizes
	static Double_t step[5]={0.005, 0.03, 0.005, 0.03, 0.01};
	
	gMyMinuit->mnparm(0,"x1",vstart[0],step[0],lower[0],upper[0],ierflag);
	gMyMinuit->mnparm(1,"y1",vstart[1],step[1],lower[1],upper[1],ierflag);
	gMyMinuit->mnparm(2,"x2",vstart[2],step[2],lower[2],upper[2],ierflag);
	gMyMinuit->mnparm(3,"y2",vstart[3],step[3],lower[3],upper[3],ierflag);
	gMyMinuit->mnparm(4,"a0",vstart[4],step[4],lower[4],upper[4],ierflag);
// ready for minimisation	
	gMyMinuit->SetPrintLevel(-1);
	gMyMinuit->mnexcm("SET OUT", arglist, 0, ierflag);
	arglist[0]= -1;
	arglist[1]= 0;
	
	gMyMinuit->mnexcm("SET NOGR", arglist, 0, ierflag);
	gMyMinuit->mnexcm("SCAN", arglist, 0, ierflag);
	gMyMinuit->mnexcm("EXIT" , arglist, 0, ierflag);
// Print results
//	Double_t amin,edm,errdef;
//	Int_t nvpar,nparx,icstat;
//	gMyMinuit->mnstat(amin,edm,errdef,nvpar,nparx,icstat);
//	gMyMinuit->mnprin(3,amin);
// Get fitted parameters

	Double_t xrec[2], yrec[2], qfrac;
	TString chname;
	Double_t epxz, b1, b2;
	Int_t ierflg;
	gMyMinuit->mnpout(0, chname, xrec[0], epxz, b1, b2, ierflg);	
	gMyMinuit->mnpout(1, chname, yrec[0], epxz, b1, b2, ierflg);	
	gMyMinuit->mnpout(2, chname, xrec[1], epxz, b1, b2, ierflg);	
	gMyMinuit->mnpout(3, chname, yrec[1], epxz, b1, b2, ierflg);	
	gMyMinuit->mnpout(4, chname, qfrac,   epxz, b1, b2, ierflg);	
	//printf("\n %f %f %f %f %f\n", xrec[0], yrec[0], xrec[1], yrec[1],qfrac);
//	delete gMyMinuit;
	
	
 //
 // One cluster for each maximum
 //
	for (j=0; j<2; j++) {
	    AliRICHRawCluster cnew;
	    if (fNPeaks == 0) {
		cnew.fNcluster[0]=-1;
		cnew.fNcluster[1]=fNRawClusters;
	    } else {
		cnew.fNcluster[0]=fNPeaks;
		cnew.fNcluster[1]=0;
	    }
	    cnew.fMultiplicity=0;
	    cnew.fX=Float_t(xrec[j]);
	    cnew.fY=Float_t(yrec[j]);
	    if (j==0) {
		cnew.fQ=Int_t(gChargeTot*qfrac);
	    } else {
		cnew.fQ=Int_t(gChargeTot*(1-qfrac));
	    }
	    gSegmentation->SetHit(xrec[j],yrec[j],0);
	    for (i=0; i<mul; i++) {
		cnew.fIndexMap[cnew.fMultiplicity]=c->fIndexMap[i];
		gSegmentation->SetPad(gix[i], giy[i]);
		Float_t q1=gResponse->IntXY(gSegmentation);
		cnew.fContMap[cnew.fMultiplicity]=Float_t(q[i])/(q1*cnew.fQ);
		cnew.fMultiplicity++;
	    }
	    FillCluster(&cnew,0);
	    //printf("\n x,y %f %f ", cnew.fX, cnew.fY);
	    cnew.fClusterType=cnew.PhysicsContribution();
	    AddRawCluster(cnew);
	    fNPeaks++;
	}
    }

    Bool_t fitted=kTRUE;

    if (nLocal !=-100 || !fitted) {
	// Check if enough local clusters have been found,
	// if not add global maxima to the list 
	//
	Int_t nPerMax;
	if (nLocal!=0) {
	    nPerMax=mul/nLocal;
	} else {
	    printf("\n Warning, no local maximum found \n");
	    nPerMax=fNperMax+1;
	}
	
	if (nPerMax > fNperMax) {
	    Int_t nGlob=mul/fNperMax-nLocal+1;
	    if (nGlob > 0) {
		Int_t nnew=0;
		for (i=0; i<mul; i++) {
		    if (!isLocal[i]) {
			indLocal[nLocal]=i;
			isLocal[i]=kTRUE;
			nLocal++;
			nnew++;
		    }
		    if (nnew==nGlob) break;
		}
	    }
	}
	//
	// Associate hits to peaks
	//
	for (i=0; i<mul; i++) {
	    Float_t dmin=1.E10;
	    Float_t qmax=0;
	    if (isLocal[i]) continue;
	    for (j=0; j<nLocal; j++) {
		Int_t il=indLocal[j];
		Float_t d=TMath::Sqrt((x[i]-x[il])*(x[i]-x[il])
				      +(y[i]-y[il])*(y[i]-y[il]));
		Float_t ql=q[il];
		//
		// Select nearest peak
		//
		if (d<dmin) {
		    dmin=d;
		    qmax=ql;
		    associatePeak[i]=j;
		} else if (d==dmin) {
		    //
		    // If more than one take highest peak
		    //
		    if (ql>qmax) {
			dmin=d;
			qmax=ql;
			associatePeak[i]=j;
		    }
		}
	    }
	}
	

 //
 // One cluster for each maximum
 //
	for (j=0; j<nLocal; j++) {
	    AliRICHRawCluster cnew;
	    if (fNPeaks == 0) {
		cnew.fNcluster[0]=-1;
		cnew.fNcluster[1]=fNRawClusters;
	    } else {
		cnew.fNcluster[0]=fNPeaks;
		cnew.fNcluster[1]=0;
	    }
	    cnew.fIndexMap[0]=c->fIndexMap[indLocal[j]];
	    cnew.fMultiplicity=1;
	    for (i=0; i<mul; i++) {
		if (isLocal[i]) continue;
		if (associatePeak[i]==j) {
		    cnew.fIndexMap[cnew.fMultiplicity]=c->fIndexMap[i];
		    cnew.fMultiplicity++;
		}
	    }
	    FillCluster(&cnew);
	    cnew.fClusterType=cnew.PhysicsContribution();
	    AddRawCluster(cnew);
	    fNPeaks++;
	}
    }
}


void  AliRICHClusterFinder::FillCluster(AliRICHRawCluster* c, Int_t flag) 
{
//
//  Completes cluster information starting from list of digits
//
    AliRICHDigit* dig;
    Float_t x, y, z;
    Int_t  ix, iy;
    Float_t frac=0;
    
    c->fPeakSignal=0;
    if (flag) {
	c->fX=0;
	c->fY=0;
	c->fQ=0;
    }
    //c->fQ=0;
 

    for (Int_t i=0; i<c->fMultiplicity; i++)
    {
	dig= (AliRICHDigit*)fDigits->UncheckedAt(c->fIndexMap[i]);
	ix=dig->fPadX+c->fOffsetMap[i];
	iy=dig->fPadY;
	Int_t q=dig->fSignal;
	if (dig->fPhysics >= dig->fSignal) {
	  c->fPhysicsMap[i]=2;
	} else if (dig->fPhysics == 0) {
	  c->fPhysicsMap[i]=0;
	} else  c->fPhysicsMap[i]=1;
//
//
// peak signal and track list
	if (flag) {
	   if (q>c->fPeakSignal) {
	      c->fPeakSignal=q;
/*
	    c->fTracks[0]=dig->fTracks[0];
	    c->fTracks[1]=dig->fTracks[1];
	    c->fTracks[2]=dig->fTracks[2];
*/
	      //c->fTracks[0]=dig->fTrack;
	    c->fTracks[0]=dig->fHit;
	    c->fTracks[1]=dig->fTracks[0];
	    c->fTracks[2]=dig->fTracks[1];
	   }
	} else {
	   if (c->fContMap[i] > frac) {
              frac=c->fContMap[i];
	      c->fPeakSignal=q;
/*
	    c->fTracks[0]=dig->fTracks[0];
	    c->fTracks[1]=dig->fTracks[1];
	    c->fTracks[2]=dig->fTracks[2];
*/
	      //c->fTracks[0]=dig->fTrack;
	    c->fTracks[0]=dig->fHit;
	    c->fTracks[1]=dig->fTracks[0];
	    c->fTracks[2]=dig->fTracks[1];
	   }
	}
//
	if (flag) {
	    fSegmentation->GetPadC(ix, iy, x, y, z);
	    c->fX += q*x;
	    c->fY += q*y;
	    c->fQ += q;
	}

    } // loop over digits

 if (flag) {
     
     c->fX/=c->fQ;
     c->fX=fSegmentation->GetAnod(c->fX);
     c->fY/=c->fQ; 
//
//  apply correction to the coordinate along the anode wire
//
     x=c->fX;   
     y=c->fY;
     fSegmentation->GetPadI(x, y, 0, ix, iy);
     fSegmentation->GetPadC(ix, iy, x, y, z);
     Int_t isec=fSegmentation->Sector(ix,iy);
     TF1* cogCorr = fSegmentation->CorrFunc(isec-1);
     
     if (cogCorr) {
	 Float_t yOnPad=(c->fY-y)/fSegmentation->Dpy(isec);
	 c->fY=c->fY-cogCorr->Eval(yOnPad, 0, 0);
     }
 }
}


void  AliRICHClusterFinder::FindCluster(Int_t i, Int_t j, AliRICHRawCluster &c){
//
//  Find clusters
//
//
//  Add i,j as element of the cluster
//    
    
    Int_t idx = fHitMap->GetHitIndex(i,j);
    AliRICHDigit* dig = (AliRICHDigit*) fHitMap->GetHit(i,j);
    Int_t q=dig->fSignal;
    if (q > TMath::Abs(c.fPeakSignal)) {
	c.fPeakSignal=q;
/*
	c.fTracks[0]=dig->fTracks[0];
	c.fTracks[1]=dig->fTracks[1];
	c.fTracks[2]=dig->fTracks[2];
*/
	//c.fTracks[0]=dig->fTrack;
	c.fTracks[0]=dig->fHit;
	c.fTracks[1]=dig->fTracks[0];
	c.fTracks[2]=dig->fTracks[1];
    }
//
//  Make sure that list of digits is ordered 
// 
    Int_t mu=c.fMultiplicity;
    c.fIndexMap[mu]=idx;

    if (dig->fPhysics >= dig->fSignal) {
        c.fPhysicsMap[mu]=2;
    } else if (dig->fPhysics == 0) {
        c.fPhysicsMap[mu]=0;
    } else  c.fPhysicsMap[mu]=1;

    if (mu > 0) {
	for (Int_t ind=mu-1; ind>=0; ind--) {
	    Int_t ist=(c.fIndexMap)[ind];
	    Int_t ql=((AliRICHDigit*)fDigits
		      ->UncheckedAt(ist))->fSignal;
	    if (q>ql) {
		c.fIndexMap[ind]=idx;
		c.fIndexMap[ind+1]=ist;
	    } else {
		break;
	    }
	}
    }
    
    c.fMultiplicity++;
    
    if (c.fMultiplicity >= 50 ) {
	printf("FindCluster - multiplicity >50  %d \n",c.fMultiplicity);
	c.fMultiplicity=49;
    }

// Prepare center of gravity calculation
    Float_t x, y, z;
    fSegmentation->GetPadC(i, j, x, y, z);
    c.fX += q*x;
    c.fY += q*y;
    c.fQ += q;
// Flag hit as taken  
    fHitMap->FlagHit(i,j);
//
//  Now look recursively for all neighbours
//  
    Int_t nn;
    Int_t xList[kMaxNeighbours], yList[kMaxNeighbours];
    fSegmentation->Neighbours(i,j,&nn,xList,yList);
    for (Int_t in=0; in<nn; in++) {
	Int_t ix=xList[in];
	Int_t iy=yList[in];
	if (fHitMap->TestHit(ix,iy)==kUnused) FindCluster(ix, iy, c);
    }
}

//_____________________________________________________________________________

void AliRICHClusterFinder::FindRawClusters()
{
  //
  // simple RICH cluster finder from digits -- finds neighbours and 
  // fill the tree with raw clusters
  //
    if (!fNdigits) return;

    fHitMap = new AliRICHHitMapA1(fSegmentation, fDigits);

    AliRICHDigit *dig;

    //printf ("Now I'm here");

    Int_t ndig;
    Int_t nskip=0;
    Int_t ncls=0;
    fHitMap->FillHits();
    for (ndig=0; ndig<fNdigits; ndig++) {
	dig = (AliRICHDigit*)fDigits->UncheckedAt(ndig);
	Int_t i=dig->fPadX;
	Int_t j=dig->fPadY;
	if (fHitMap->TestHit(i,j)==kUsed ||fHitMap->TestHit(i,j)==kEmpty) {
	    nskip++;
	    continue;
	}
	AliRICHRawCluster c;
	c.fMultiplicity=0;
	c.fPeakSignal=dig->fSignal;
/*
	c.fTracks[0]=dig->fTracks[0];
	c.fTracks[1]=dig->fTracks[1];
	c.fTracks[2]=dig->fTracks[2];
*/
	//c.fTracks[0]=dig->fTrack;
	c.fTracks[0]=dig->fHit;
	c.fTracks[1]=dig->fTracks[0];
	c.fTracks[2]=dig->fTracks[1];
        // tag the beginning of cluster list in a raw cluster
        c.fNcluster[0]=-1;
	FindCluster(i,j, c);
	// center of gravity
	c.fX /= c.fQ;
	c.fX=fSegmentation->GetAnod(c.fX);
	c.fY /= c.fQ;
//
//  apply correction to the coordinate along the anode wire
//
	Int_t ix,iy;
	Float_t x=c.fX;   
	Float_t y=c.fY;
	Float_t z;
	
	fSegmentation->GetPadI(x, y, 0, ix, iy);
	fSegmentation->GetPadC(ix, iy, x, y, z);
	Int_t isec=fSegmentation->Sector(ix,iy);
	TF1* cogCorr=fSegmentation->CorrFunc(isec-1);
	if (cogCorr) {
	    Float_t yOnPad=(c.fY-y)/fSegmentation->Dpy(isec);
	    c.fY=c.fY-cogCorr->Eval(yOnPad,0,0);
	}

//
//      Analyse cluster and decluster if necessary
//	
    ncls++;
    c.fNcluster[1]=fNRawClusters;
    c.fClusterType=c.PhysicsContribution();
    Decluster(&c);
    fNPeaks=0;
//
//
//
//      reset Cluster object
	for (int k=0;k<c.fMultiplicity;k++) {
	    c.fIndexMap[k]=0;
	}
	c.fMultiplicity=0;
    } // end loop ndig    
    delete fHitMap;
}

void AliRICHClusterFinder::
CalibrateCOG()
{

// Calibration

    Float_t x[5];
    Float_t y[5];
    Int_t n, i;
    TF1 func;
    if (fSegmentation) {
	fSegmentation->GiveTestPoints(n, x, y);
	for (i=0; i<n; i++) {
	    Float_t xtest=x[i];
	    Float_t ytest=y[i];	    
	    SinoidalFit(xtest, ytest, func);
	    fSegmentation->SetCorrFunc(i, new TF1(func));
	}
    }
}


void AliRICHClusterFinder::
SinoidalFit(Float_t x, Float_t y, TF1 &func)
{
// Sinoidal fit


    static Int_t count=0;
    char canvasname[3];
    Float_t z;
    
    count++;
    sprintf(canvasname,"c%d",count);

    const Int_t kNs=101;
    Float_t xg[kNs], yg[kNs], xrg[kNs], yrg[kNs];
    Float_t xsig[kNs], ysig[kNs];
   
    AliSegmentation *segmentation=fSegmentation;

    Int_t ix,iy;
    segmentation->GetPadI(x,y,0,ix,iy);   
    segmentation->GetPadC(ix,iy,x,y,z);   
    Int_t isec=segmentation->Sector(ix,iy);
// Pad Limits    
    Float_t xmin = x-segmentation->Dpx(isec)/2;
    Float_t ymin = y-segmentation->Dpy(isec)/2;
//      	
//      Integration Limits
    Float_t dxI=fResponse->SigmaIntegration()*fResponse->ChargeSpreadX();
    Float_t dyI=fResponse->SigmaIntegration()*fResponse->ChargeSpreadY();

//
//  Scanning
//
    Int_t i;
    Float_t qp;
//
//  y-position
    Float_t yscan=ymin;
    Float_t dy=segmentation->Dpy(isec)/(kNs-1);

    for (i=0; i<kNs; i++) {
//
//      Pad Loop
//      
	Float_t sum=0;
	Float_t qcheck=0;
	segmentation->SigGenInit(x, yscan, 0);
	
	for (segmentation->FirstPad(x, yscan,0, dxI, dyI); 
	     segmentation->MorePads(); 
	     segmentation->NextPad()) 
	{
	    qp=fResponse->IntXY(segmentation);
	    qp=TMath::Abs(qp);
//
//
	    if (qp > 1.e-4) {
		qcheck+=qp;
		Int_t ixs=segmentation->Ix();
		Int_t iys=segmentation->Iy();
		Float_t xs,ys,zs;
		segmentation->GetPadC(ixs,iys,xs,ys,zs);
		sum+=qp*ys;
	    }
	} // Pad loop
	Float_t ycog=sum/qcheck;
	yg[i]=(yscan-y)/segmentation->Dpy(isec);
	yrg[i]=(ycog-y)/segmentation->Dpy(isec);
	ysig[i]=ycog-yscan;
	yscan+=dy;
    } // scan loop
//
//  x-position
    Float_t xscan=xmin;
    Float_t dx=segmentation->Dpx(isec)/(kNs-1);

    for (i=0; i<kNs; i++) {
//
//      Pad Loop
//      
	Float_t sum=0;
	Float_t qcheck=0;
	segmentation->SigGenInit(xscan, y, 0);
	
	for (segmentation->FirstPad(xscan, y, 0, dxI, dyI); 
	     segmentation->MorePads(); 
	     segmentation->NextPad()) 
	{
	    qp=fResponse->IntXY(segmentation);
	    qp=TMath::Abs(qp);
//
//
	    if (qp > 1.e-2) {
		qcheck+=qp;
		Int_t ixs=segmentation->Ix();
		Int_t iys=segmentation->Iy();
		Float_t xs,ys,zs;
		segmentation->GetPadC(ixs,iys,xs,ys,zs);
		sum+=qp*xs;
	    }
	} // Pad loop
	Float_t xcog=sum/qcheck;
	xcog=segmentation->GetAnod(xcog);
	
	xg[i]=(xscan-x)/segmentation->Dpx(isec);
	xrg[i]=(xcog-x)/segmentation->Dpx(isec);
	xsig[i]=xcog-xscan;
	xscan+=dx;
    }
//
// Creates a Root function based on function sinoid above
// and perform the fit
//
    //    TGraph *graphx = new TGraph(kNs,xg ,xsig);
    //    TGraph *graphxr= new TGraph(kNs,xrg,xsig);   
    //    TGraph *graphy = new TGraph(kNs,yg ,ysig);
    TGraph *graphyr= new TGraph(kNs,yrg,ysig);

    Double_t sinoid(Double_t *x, Double_t *par);
    new TF1("sinoidf",sinoid,0.5,0.5,5);
    graphyr->Fit("sinoidf","Q");
    func = *((TF1*)((graphyr->GetListOfFunctions())->At(0)));
/*
    
    TCanvas *c1=new TCanvas(canvasname,canvasname,400,10,600,700);
    TPad* pad11 = new TPad("pad11"," ",0.01,0.51,0.49,0.99);
    TPad* pad12 = new TPad("pad12"," ",0.51,0.51,0.99,0.99);
    TPad* pad13 = new TPad("pad13"," ",0.01,0.01,0.49,0.49);
    TPad* pad14 = new TPad("pad14"," ",0.51,0.01,0.99,0.49);
    pad11->SetFillColor(11);
    pad12->SetFillColor(11);
    pad13->SetFillColor(11);
    pad14->SetFillColor(11);
    pad11->Draw();
    pad12->Draw();
    pad13->Draw();
    pad14->Draw();
    
//
    pad11->cd();
    graphx->SetFillColor(42);
    graphx->SetMarkerColor(4);
    graphx->SetMarkerStyle(21);
    graphx->Draw("AC");
    graphx->GetHistogram()->SetXTitle("x on pad");
    graphx->GetHistogram()->SetYTitle("xcog-x");   


    pad12->cd();
    graphxr->SetFillColor(42);
    graphxr->SetMarkerColor(4);
    graphxr->SetMarkerStyle(21);
    graphxr->Draw("AP");
    graphxr->GetHistogram()->SetXTitle("xcog on pad");
    graphxr->GetHistogram()->SetYTitle("xcog-x");   
    

    pad13->cd();
    graphy->SetFillColor(42);
    graphy->SetMarkerColor(4);
    graphy->SetMarkerStyle(21);
    graphy->Draw("AF");
    graphy->GetHistogram()->SetXTitle("y on pad");
    graphy->GetHistogram()->SetYTitle("ycog-y");   
    


    pad14->cd();
    graphyr->SetFillColor(42);
    graphyr->SetMarkerColor(4);
    graphyr->SetMarkerStyle(21);
    graphyr->Draw("AF");
    graphyr->GetHistogram()->SetXTitle("ycog on pad");
    graphyr->GetHistogram()->SetYTitle("ycog-y");   
    
    c1->Update();
*/
}

Double_t sinoid(Double_t *x, Double_t *par)
{

// Sinoid function

    Double_t arg = -2*TMath::Pi()*x[0];
    Double_t fitval= par[0]*TMath::Sin(arg)+
	par[1]*TMath::Sin(2*arg)+
	par[2]*TMath::Sin(3*arg)+
	par[3]*TMath::Sin(4*arg)+
	par[4]*TMath::Sin(5*arg);
    return fitval;
 }


Double_t DoubleGauss(Double_t *x, Double_t *par)
{

// Doublr gaussian function

    Double_t arg1 = (x[0]-par[1])/0.18;
    Double_t arg2 = (x[0]-par[3])/0.18;
    Double_t fitval= par[0]*TMath::Exp(-arg1*arg1/2)
	+par[2]*TMath::Exp(-arg2*arg2/2);
    return fitval;
 }

Float_t DiscrCharge(Int_t i,Double_t *par) 
{
// par[0]    x-position of first  cluster
// par[1]    y-position of first  cluster
// par[2]    x-position of second cluster
// par[3]    y-position of second cluster
// par[4]    charge fraction of first  cluster
// 1-par[4]  charge fraction of second cluster

    static Float_t qtot;
    if (gFirst) {
	qtot=0;
	for (Int_t jbin=0; jbin<gNbins; jbin++) {
	    qtot+=gCharge[jbin];
	}
	gFirst=0;
	//printf("\n sum of charge from DiscrCharge %f\n", qtot);
	gChargeTot=Int_t(qtot);
	
    }
    gSegmentation->SetPad(gix[i], giy[i]);
//  First Cluster
    gSegmentation->SetHit(par[0],par[1],0);
    Float_t q1=gResponse->IntXY(gSegmentation);
    
//  Second Cluster
    gSegmentation->SetHit(par[2],par[3],0);
    Float_t q2=gResponse->IntXY(gSegmentation);
    
    Float_t value = qtot*(par[4]*q1+(1.-par[4])*q2);
    return value;
}

//
// Minimisation function
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    Int_t i;
    Float_t delta;
    Float_t chisq=0;
    Float_t qcont=0;
    Float_t qtot=0;
    
    for (i=0; i<gNbins; i++) {
	Float_t q0=gCharge[i];
	Float_t q1=DiscrCharge(i,par);
	delta=(q0-q1)/TMath::Sqrt(q0);
	chisq+=delta*delta;
	qcont+=q1;
	qtot+=q0;
    }
    chisq=chisq+=(qtot-qcont)*(qtot-qcont)*0.5;
    f=chisq;
}


void AliRICHClusterFinder::SetDigits(TClonesArray *RICHdigits)
{

// Get all the digits

    fDigits=RICHdigits;
    fNdigits = fDigits->GetEntriesFast();
}

AliRICHClusterFinder& AliRICHClusterFinder::operator=(const AliRICHClusterFinder& rhs)
{
// Assignment operator
    return *this;
    
}
