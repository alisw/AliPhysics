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
Revision 1.5  2000/06/15 07:58:48  morsch
Code from MUON-dev joined

Revision 1.4.4.2  2000/06/09 21:58:15  morsch
Most coding rule violations corrected.

*/

#include "AliMUONClusterFinder.h"
#include "AliMUON.h"
#include "AliMUONHitMap.h"
#include "AliMUONHitMapA1.h"
#include "AliMUONSegmentation.h"
#include "AliMUONResponse.h"
#include "AliMUONDigit.h"
#include "AliMUONRawCluster.h"
#include "AliRun.h"


#include <TTree.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TPad.h>
#include <TGraph.h> 
#include <TPostScript.h> 
#include <TMinuit.h> 
#include <TClonesArray.h> 

//_____________________________________________________________________
static AliMUONSegmentation* fgSegmentation;
static AliMUONResponse*     fgResponse;
static Int_t                fgix[500];
static Int_t                fgiy[500];
static Float_t              fgCharge[500];
static Int_t                fgNbins;
static Int_t                fgFirst=kTRUE;
static Int_t                fgChargeTot;
static Float_t              fgQtot;
static TMinuit*             fgMyMinuit ;
// This function is minimized in the double-Mathieson fit
void fcn2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
void fcn1(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);



ClassImp(AliMUONClusterFinder)

    AliMUONClusterFinder::AliMUONClusterFinder
	(AliMUONSegmentation *segmentation, 
	AliMUONResponse *response, 
 	TClonesArray *digits, Int_t chamber)   
{
// Constructor    
    fSegmentation[0]=segmentation;
    fResponse=response;
    
    fDigits=digits;
    fNdigits = fDigits->GetEntriesFast();
    fChamber=chamber;
    fRawClusters=new TClonesArray("AliMUONRawCluster",10000);
    fNRawClusters=0;
    fCogCorr = 0;
    SetNperMax();
    SetClusterSize();
    SetDeclusterFlag();
    fNPeaks=-1;
}

    AliMUONClusterFinder::AliMUONClusterFinder()
{
// Default constructor
    fSegmentation[0]=0;
    fResponse=0;
    fDigits=0;
    fNdigits = 0;
    fChamber=-1;
    fRawClusters=new TClonesArray("AliMUONRawCluster",10000);
    fNRawClusters=0;
    fHitMap = 0;
    fCogCorr = 0;
    SetNperMax();
    SetClusterSize();
    SetDeclusterFlag();
    fNPeaks=-1;
}

AliMUONClusterFinder::AliMUONClusterFinder(
    const AliMUONClusterFinder & clusterFinder)
{
// Dummy copy Constructor
    ;
}

AliMUONClusterFinder::~AliMUONClusterFinder()
{
// Destructor
    delete fRawClusters;
}

void AliMUONClusterFinder::SetDigits(TClonesArray *MUONdigits) 
{
// Set pointer to digits
    fDigits=MUONdigits;
    fNdigits = fDigits->GetEntriesFast();
}

void AliMUONClusterFinder::AddRawCluster(const AliMUONRawCluster c)
{
  //
  // Add a raw cluster copy to the list
  //
    AliMUON *pMUON=(AliMUON*)gAlice->GetModule("MUON");
    pMUON->AddRawCluster(fChamber,c); 
    fNRawClusters++;
}



void AliMUONClusterFinder::Decluster(AliMUONRawCluster *cluster)
{
// Decluster composite clusters
    Bool_t fitted;
    
    Int_t mul = cluster->fMultiplicity[0];
    if (mul == 1 || mul ==2) {
	if (mul==2) {
	    fitted=SingleMathiesonFit(cluster);
	}
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
// 
// This part could be activated again in the future
	if (fDeclusterFlag) {
	    if (Centered(cluster)) {
		// ok, cluster is centered
		fitted = SingleMathiesonFit(cluster); 	    
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
		fitted=SingleMathiesonFit(cluster);
		AddRawCluster(*cluster);
		fNPeaks++;
	    } // if Declustering selected
	} // if < maximum clustersize for deconvolution 
    } // multiplicity 
}


Bool_t AliMUONClusterFinder::Centered(AliMUONRawCluster *cluster)
{
// True if cluster is centered
    AliMUONDigit* dig;
    dig= (AliMUONDigit*)fDigits->UncheckedAt(cluster->fIndexMap[0][0]);
    Int_t ix=dig->fPadX;
    Int_t iy=dig->fPadY;
    Int_t nn;
    Int_t x[kMaxNeighbours], y[kMaxNeighbours], xN[kMaxNeighbours], yN[kMaxNeighbours];
    
    fSegmentation[0]->Neighbours(ix,iy,&nn,x,y);
    Int_t nd=0;
    for (Int_t i=0; i<nn; i++) {
	if (fHitMap->TestHit(x[i],y[i]) == kUsed) {
	    xN[nd]=x[i];
	    yN[nd]=y[i];
	    nd++;
	}
    }
    if (nd==2) {
//
// cluster is centered !
	if (fNPeaks != 0) {
            cluster->fNcluster[0]=fNPeaks;
            cluster->fNcluster[1]=0;
	}  
	AddRawCluster(*cluster);
	fNPeaks++;
	return kTRUE;
    } else if (nd == 1) {
//
// Highest signal on an edge, split cluster into 2+1
//
// who is the neighbour ?
	Int_t nind=fHitMap->GetHitIndex(xN[0], yN[0]);
	Int_t i1= (nind==cluster->fIndexMap[1][0]) ? 1:2;
	Int_t i2= (nind==cluster->fIndexMap[1][0]) ? 2:1;    
//
// 2-cluster
	AliMUONRawCluster cnew;
	if (fNPeaks == 0) {
	    cnew.fNcluster[0]=-1;
            cnew.fNcluster[1]=fNRawClusters;
        } else {
            cnew.fNcluster[0]=fNPeaks;
            cnew.fNcluster[1]=0;
        }
	cnew.fMultiplicity[0]=2;
	cnew.fIndexMap[0][0]=cluster->fIndexMap[0][0];
	cnew.fIndexMap[1][0]=cluster->fIndexMap[i1][0];
	FillCluster(&cnew);
	cnew.fClusterType=cnew.PhysicsContribution();
	AddRawCluster(cnew);
        fNPeaks++;
//
// 1-cluster
	cluster->fMultiplicity[0]=1;
	cluster->fIndexMap[0][0]=cluster->fIndexMap[i2][0];
	cluster->fIndexMap[1][0]=0;
	cluster->fIndexMap[2][0]=0;	
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

void AliMUONClusterFinder::SplitByLocalMaxima(AliMUONRawCluster *c)
{
// Search for local maxima and split cluster accordingly
    Bool_t fitted;
    
    AliMUONDigit* digt;
     Int_t i; // loops over digits
     Int_t j; // loops over local maxima
     fMul=c->fMultiplicity[0];
//
//  dump digit information into arrays
//
    for (i=0; i<fMul; i++)
    {
	fDig[i]= (AliMUONDigit*)fDigits->UncheckedAt(c->fIndexMap[i][0]);
	fIx[i]= fDig[i]->fPadX;
	fIy[i]= fDig[i]->fPadY;
	fQ[i] = fDig[i]->fSignal;
	fSegmentation[0]->GetPadCxy(fIx[i], fIy[i], fX[i], fY[i]);
    }
//
//  Find local maxima
//
    fNLocal=0;
    Bool_t isLocal[100];
    Int_t assocPeak[100];
    Int_t nn;
    Int_t x[kMaxNeighbours], y[kMaxNeighbours];
    for (i=0; i<fMul; i++) {
	fSegmentation[0]->Neighbours(fIx[i], fIy[i], &nn, x, y);
	isLocal[i]=kTRUE;
	for (j=0; j<nn; j++) {
	    if (fHitMap->TestHit(x[j], y[j])==kEmpty) continue;
	    digt=(AliMUONDigit*) fHitMap->GetHit(x[j], y[j]);
	    if (digt->fSignal > fQ[i]) {
		isLocal[i]=kFALSE;
		break;
//
// handle special case of neighbouring pads with equal signal
	    } else if (digt->fSignal == fQ[i]) {
		if (fNLocal >0) {
		    for (Int_t k=0; k<fNLocal; k++) {
			if (x[j]==fIx[fIndLocal[k]] && y[j]==fIy[fIndLocal[k]])
			{
			    isLocal[i]=kFALSE;
			} 
		    } // loop over local maxima
		} // are there are already local maxima
	    } 
	} // loop over next neighbours
	// Maxima should not be on the edge
	if (isLocal[i]) {
	    fIndLocal[fNLocal]=i;
	    fNLocal++;
	} 
    } // loop over all digits
//    printf("Found %d local Maxima",fNLocal);
//
// If only one local maximum found but multiplicity is high 
// take global maximum from the list of digits.    
// 12 should not be hard wired but a parameter 
// 
    if (fNLocal==1 && fMul>=1) {
	Int_t nnew=0;
	for (i=0; i<fMul; i++) {
	    if (!isLocal[i]) {
		fIndLocal[fNLocal]=i;
		isLocal[i]=kTRUE;
		fNLocal++;
		nnew++;
	    }
	    if (nnew==1) break;
	}
    }
    
// If number of local maxima is 2 try to fit a double mathieson
    if (fNLocal==2) {
	fitted = DoubleMathiesonFit(c);
    }
    
    if (fNLocal !=2 || !fitted) {
	// Check if enough local clusters have been found,
	// if not add global maxima to the list 
	//
	Int_t nPerMax;
	if (fNLocal!=0) {
	    nPerMax=fMul/fNLocal;
	} else {
	    printf("\n Warning, no local maximum found \n");
	    nPerMax=fNperMax+1;
	}

	if (nPerMax > fNperMax) {
	    Int_t nGlob=fMul/fNperMax-fNLocal+1;
	    if (nGlob > 0) {
		Int_t nnew=0;
		for (i=0; i<fMul; i++) {
		    if (!isLocal[i]) {
			fIndLocal[fNLocal]=i;
			isLocal[i]=kTRUE;
			fNLocal++;
			nnew++;
		    }
		    if (nnew==nGlob) break;
		}
	    }
	}
	//
	// Associate hits to peaks
	//
	for (i=0; i<fMul; i++) {
	    Float_t dmin=1.E10;
	    Float_t qmax=0;
	    if (isLocal[i]) continue;
	    for (j=0; j<fNLocal; j++) {
		Int_t il=fIndLocal[j];
		Float_t d=TMath::Sqrt((fX[i]-fX[il])*(fX[i]-fX[il])
				      +(fY[i]-fY[il])*(fY[i]-fY[il]));
		Float_t ql=fQ[il];
 //
 // Select nearest peak
 //
		if (d<dmin) {
		    dmin=d;
		    qmax=ql;
		    assocPeak[i]=j;
		} else if (d==dmin) {
 //
 // If more than one take highest peak
 //
		    if (ql>qmax) {
			dmin=d;
			qmax=ql;
			assocPeak[i]=j;
		    }
		}
	    }
	}
	//
	// One cluster for each maximum
	//
	for (j=0; j<fNLocal; j++) {
	    AliMUONRawCluster cnew;
	    if (fNPeaks == 0) {
		cnew.fNcluster[0]=-1;
		cnew.fNcluster[1]=fNRawClusters;
	    } else {
		cnew.fNcluster[0]=fNPeaks;
		cnew.fNcluster[1]=0;
	    }
	    cnew.fIndexMap[0][0]=c->fIndexMap[fIndLocal[j]][0];
	    cnew.fMultiplicity[0]=1;
	    for (i=0; i<fMul; i++) {
		if (isLocal[i]) continue;
		if (assocPeak[i]==j) {
		    cnew.fIndexMap[cnew.fMultiplicity[0]][0]=c->fIndexMap[i][0];
		    cnew.fMultiplicity[0]++;
		}
	    }
	    FillCluster(&cnew);
	    cnew.fClusterType=cnew.PhysicsContribution();
	    AddRawCluster(cnew);
	    fNPeaks++;
	}
    }
}


void  AliMUONClusterFinder::FillCluster(AliMUONRawCluster* c, Int_t flag) 
{
//
//  Completes cluster information starting from list of digits
//
    AliMUONDigit* dig;
    Float_t x, y;
    Int_t  ix, iy;
    Float_t frac=0;
    
    c->fPeakSignal[0]=0;
    if (flag) {
	c->fX[0]=0;
	c->fY[0]=0;
	c->fQ[0]=0;
    }
 

    for (Int_t i=0; i<c->fMultiplicity[0]; i++)
    {
	dig= (AliMUONDigit*)fDigits->UncheckedAt(c->fIndexMap[i][0]);
	ix=dig->fPadX+c->fOffsetMap[i][0];
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
	    if (q>c->fPeakSignal[0]) {
		c->fPeakSignal[0]=q;
		c->fTracks[0]=dig->fHit;
		c->fTracks[1]=dig->fTracks[0];
		c->fTracks[2]=dig->fTracks[1];
	    }
	} else {
	    if (c->fContMap[i][0] > frac) {
		frac=c->fContMap[i][0];
		c->fPeakSignal[0]=q;
		c->fTracks[0]=dig->fHit;
		c->fTracks[1]=dig->fTracks[0];
		c->fTracks[2]=dig->fTracks[1];
	    }
	}
//
	if (flag) {
	    fSegmentation[0]->GetPadCxy(ix, iy, x, y);
	    c->fX[0] += q*x;
	    c->fY[0] += q*y;
	    c->fQ[0] += q;
	}
    } // loop over digits
    
    if (flag) {
    	c->fX[0]/=c->fQ[0];
     	c->fX[0]=fSegmentation[0]->GetAnod(c->fX[0]);
     	c->fY[0]/=c->fQ[0]; 
//
//  apply correction to the coordinate along the anode wire
//
     	x=c->fX[0];   
     	y=c->fY[0];
     	fSegmentation[0]->GetPadIxy(x, y, ix, iy);
     	fSegmentation[0]->GetPadCxy(ix, iy, x, y);
     	Int_t isec=fSegmentation[0]->Sector(ix,iy);
    	TF1* cogCorr = fSegmentation[0]->CorrFunc(isec-1);
	
     	if (cogCorr) {
	    Float_t yOnPad=(c->fY[0]-y)/fSegmentation[0]->Dpy(isec);
	    c->fY[0]=c->fY[0]-cogCorr->Eval(yOnPad, 0, 0);
     	}
    }
}


void  AliMUONClusterFinder::FindCluster(Int_t i, Int_t j, AliMUONRawCluster &c){
//
//  Find clusters
//
//
//  Add i,j as element of the cluster
//    
    Int_t idx = fHitMap->GetHitIndex(i,j);
    AliMUONDigit* dig = (AliMUONDigit*) fHitMap->GetHit(i,j);
    Int_t q=dig->fSignal;
    if (q > TMath::Abs(c.fPeakSignal[0])) {
	c.fPeakSignal[0]=q;
	c.fTracks[0]=dig->fHit;
	c.fTracks[1]=dig->fTracks[0];
	c.fTracks[2]=dig->fTracks[1];
    }
//
//  Make sure that list of digits is ordered 
// 
    Int_t mu=c.fMultiplicity[0];
    c.fIndexMap[mu][0]=idx;
    
    if (dig->fPhysics >= dig->fSignal) {
        c.fPhysicsMap[mu]=2;
    } else if (dig->fPhysics == 0) {
        c.fPhysicsMap[mu]=0;
    } else  c.fPhysicsMap[mu]=1;
    
    if (mu > 0) {
	for (Int_t ind=mu-1; ind>=0; ind--) {
	    Int_t ist=(c.fIndexMap)[ind][0];
	    Int_t ql=((AliMUONDigit*)fDigits
		      ->UncheckedAt(ist))->fSignal;
	    if (q>ql) {
		c.fIndexMap[ind][0]=idx;
		c.fIndexMap[ind+1][0]=ist;
	    } else {
		break;
	    }
	}
    }
    
    c.fMultiplicity[0]++;
    
    if (c.fMultiplicity[0] >= 50 ) {
	printf("FindCluster - multiplicity >50  %d \n",c.fMultiplicity[0]);
	c.fMultiplicity[0]=49;
    }

// Prepare center of gravity calculation
    Float_t x, y;
    fSegmentation[0]->GetPadCxy(i, j, x, y);
    c.fX[0] += q*x;
    c.fY[0] += q*y;
    c.fQ[0] += q;
// Flag hit as taken  
    fHitMap->FlagHit(i,j);
//
//  Now look recursively for all neighbours
//  
    Int_t nn;
    Int_t xList[kMaxNeighbours], yList[kMaxNeighbours];
    fSegmentation[0]->Neighbours(i,j,&nn,xList,yList);
    for (Int_t in=0; in<nn; in++) {
	Int_t ix=xList[in];
	Int_t iy=yList[in];
	if (fHitMap->TestHit(ix,iy)==kUnused) FindCluster(ix, iy, c);
    }
}

//_____________________________________________________________________________

void AliMUONClusterFinder::FindRawClusters()
{
  //
  // simple MUON cluster finder from digits -- finds neighbours and 
  // fill the tree with raw clusters
  //
    if (!fNdigits) return;

    fHitMap = new AliMUONHitMapA1(fSegmentation[0], fDigits);

    AliMUONDigit *dig;

    Int_t ndig;
    Int_t nskip=0;
    Int_t ncls=0;
    fHitMap->FillHits();
    for (ndig=0; ndig<fNdigits; ndig++) {
	dig = (AliMUONDigit*)fDigits->UncheckedAt(ndig);
	Int_t i=dig->fPadX;
	Int_t j=dig->fPadY;
	if (fHitMap->TestHit(i,j)==kUsed ||fHitMap->TestHit(i,j)==kEmpty) {
	    nskip++;
	    continue;
	}
	AliMUONRawCluster c;
	c.fMultiplicity[0]=0;
	c.fPeakSignal[0]=dig->fSignal;
	c.fTracks[0]=dig->fHit;
	c.fTracks[1]=dig->fTracks[0];
	c.fTracks[2]=dig->fTracks[1];
	// tag the beginning of cluster list in a raw cluster
	c.fNcluster[0]=-1;
	FindCluster(i,j, c);
	// center of gravity
	c.fX[0] /= c.fQ[0];
	c.fX[0]=fSegmentation[0]->GetAnod(c.fX[0]);
	c.fY[0] /= c.fQ[0];
//
//  apply correction to the coordinate along the anode wire
//
	Int_t ix,iy;
	Float_t x=c.fX[0];   
	Float_t y=c.fY[0];
	fSegmentation[0]->GetPadIxy(x, y, ix, iy);
	fSegmentation[0]->GetPadCxy(ix, iy, x, y);
	Int_t isec=fSegmentation[0]->Sector(ix,iy);
	TF1* cogCorr=fSegmentation[0]->CorrFunc(isec-1);
	if (cogCorr) {
	    Float_t yOnPad=(c.fY[0]-y)/fSegmentation[0]->Dpy(isec);
	    c.fY[0]=c.fY[0]-cogCorr->Eval(yOnPad,0,0);
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
	for (int k=0;k<c.fMultiplicity[0];k++) {
	    c.fIndexMap[k][0]=0;
	}
	c.fMultiplicity[0]=0;
    } // end loop ndig    
    delete fHitMap;
}

void AliMUONClusterFinder::
CalibrateCOG()
{
// Calibrate the cog method
    Float_t x[5];
    Float_t y[5];
    Int_t n, i;
    TF1 func;
    if (fSegmentation[0]) {
	fSegmentation[0]->GiveTestPoints(n, x, y);
	for (i=0; i<n; i++) {
	    Float_t xtest=x[i];
	    Float_t ytest=y[i];	    
	    SinoidalFit(xtest, ytest, func);
	    fSegmentation[0]->SetCorrFunc(i, new TF1(func));
	}
    }
}


void AliMUONClusterFinder::
SinoidalFit(Float_t x, Float_t y, TF1 &func)
{
// Perform a senoidal fit to the residuals of the cog method
    static Int_t count=0;
    char canvasname[3];
    count++;
    sprintf(canvasname,"c%d",count);
    
    const Int_t kns=101;
    Float_t xg[kns], yg[kns], xrg[kns], yrg[kns];
    Float_t xsig[kns], ysig[kns];
    
    AliMUONSegmentation *segmentation=fSegmentation[0];
    
    Int_t ix,iy;
    segmentation->GetPadIxy(x,y,ix,iy);   
    segmentation->GetPadCxy(ix,iy,x,y);   
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
    Float_t dy=segmentation->Dpy(isec)/(kns-1);

    for (i=0; i<kns; i++) {
//
//      Pad Loop
//      
	Float_t sum=0;
	Float_t qcheck=0;
	segmentation->SigGenInit(x, yscan, 0);
	
	for (segmentation->FirstPad(x, yscan, dxI, dyI); 
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
		Float_t xs,ys;
		segmentation->GetPadCxy(ixs,iys,xs,ys);
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
    Float_t dx=segmentation->Dpx(isec)/(kns-1);
    
    for (i=0; i<kns; i++) {
//
//      Pad Loop
//      
	Float_t sum=0;
	Float_t qcheck=0;
	segmentation->SigGenInit(xscan, y, 0);
	
	for (segmentation->FirstPad(xscan, y, dxI, dyI); 
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
		Float_t xs,ys;
		segmentation->GetPadCxy(ixs,iys,xs,ys);
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
/*
        TGraph *graphx = new TGraph(kns,xg ,xsig);
        TGraph *graphxr= new TGraph(kns,xrg,xsig);   
        TGraph *graphy = new TGraph(kns,yg ,ysig);
*/
    TGraph *graphyr= new TGraph(kns,yrg,ysig);

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

Bool_t AliMUONClusterFinder::SingleMathiesonFit(AliMUONRawCluster *c)
{
//
//  Initialise global variables for fit
    Int_t i;
    fMul=c->fMultiplicity[0];
    fgSegmentation=fSegmentation[0];
    fgResponse    =fResponse;
    fgNbins=fMul;
    Float_t qtot=0;
//
//  dump digit information into arrays
//
    for (i=0; i<fMul; i++)
    {
	fDig[i]= (AliMUONDigit*)fDigits->UncheckedAt(c->fIndexMap[i][0]);
	fIx[i]= fDig[i]->fPadX;
	fIy[i]= fDig[i]->fPadY;
	fQ[i] = fDig[i]->fSignal;
	fSegmentation[0]->GetPadCxy(fIx[i], fIy[i], fX[i], fY[i]);
	fgix[i]=fIx[i];
	fgiy[i]=fIy[i];
	fgCharge[i]=Float_t(fQ[i]);
	qtot+=fgCharge[i];
    }

    fgQtot=qtot;
    fgChargeTot=Int_t(qtot);
    
//
    if (fgFirst) {
	fgFirst=kFALSE;
	fgMyMinuit = new TMinuit(5);
    }

    fgMyMinuit->SetFCN(fcn1);
    fgMyMinuit->mninit(2,10,7);
    Double_t arglist[20];
    Int_t ierflag=0;
    arglist[0]=1;
//	fgMyMinuit->mnexcm("SET ERR",arglist,1,ierflag);
// Set starting values 
    static Double_t vstart[2];
    vstart[0]=c->fX[0];
    vstart[1]=c->fY[0];	
// lower and upper limits
    static Double_t lower[2], upper[2];
    Int_t ix,iy;
    fSegmentation[0]->GetPadIxy(c->fX[0], c->fY[0], ix, iy);
    Int_t isec=fSegmentation[0]->Sector(ix, iy);
    lower[0]=vstart[0]-fSegmentation[0]->Dpx(isec)/2;
    lower[1]=vstart[1]-fSegmentation[0]->Dpy(isec)/2;
    
    upper[0]=lower[0]+fSegmentation[0]->Dpx(isec);
    upper[1]=lower[1]+fSegmentation[0]->Dpy(isec);
    
// step sizes
    static Double_t step[2]={0.0005, 0.0005};
    
    fgMyMinuit->mnparm(0,"x1",vstart[0],step[0],lower[0],upper[0],ierflag);
    fgMyMinuit->mnparm(1,"y1",vstart[1],step[1],lower[1],upper[1],ierflag);
// ready for minimisation	
    fgMyMinuit->SetPrintLevel(1);
    fgMyMinuit->mnexcm("SET OUT", arglist, 0, ierflag);
    arglist[0]= -1;
    arglist[1]= 0;
    
    fgMyMinuit->mnexcm("SET NOGR", arglist, 0, ierflag);
    fgMyMinuit->mnexcm("MIGRAD", arglist, 0, ierflag);
    fgMyMinuit->mnexcm("EXIT" , arglist, 0, ierflag);
// Print results
// Get fitted parameters
    Double_t xrec, yrec;
    TString chname;
    Double_t epxz, b1, b2;
    Int_t ierflg;
    fgMyMinuit->mnpout(0, chname, xrec, epxz, b1, b2, ierflg);	
    fgMyMinuit->mnpout(1, chname, yrec, epxz, b1, b2, ierflg);	
    c->fX[0]=xrec;
    c->fY[0]=yrec;
    return kTRUE;
}

Bool_t AliMUONClusterFinder::DoubleMathiesonFit(AliMUONRawCluster *c)
{
//
//  Initialise global variables for fit
    Int_t i,j;
    
    fgSegmentation=fSegmentation[0];
    fgResponse    =fResponse;
    fgNbins=fMul;
    Float_t qtot=0;
	
    for (i=0; i<fMul; i++) {
	fgix[i]=fIx[i];
	fgiy[i]=fIy[i];
	fgCharge[i]=Float_t(fQ[i]);
	qtot+=fgCharge[i];
    }
    fgQtot=qtot;
    fgChargeTot=Int_t(qtot);
    
//
    if (fgFirst) {
	fgFirst=kFALSE;
	fgMyMinuit = new TMinuit(5);
    }
    fgMyMinuit->SetFCN(fcn2);
    fgMyMinuit->mninit(5,10,7);
    Double_t arglist[20];
    Int_t ierflag=0;
    arglist[0]=1;
//	fgMyMinuit->mnexcm("SET ERR",arglist,1,ierflag);
// Set starting values 
    static Double_t vstart[5];
    vstart[0]=fX[fIndLocal[0]];
    vstart[1]=fY[fIndLocal[0]];	
    vstart[2]=fX[fIndLocal[1]];
    vstart[3]=fY[fIndLocal[1]];	
    vstart[4]=Float_t(fQ[fIndLocal[0]])/
	Float_t(fQ[fIndLocal[0]]+fQ[fIndLocal[1]]);
// lower and upper limits
    static Double_t lower[5], upper[5];
    Int_t isec=fSegmentation[0]->Sector(fIx[fIndLocal[0]], fIy[fIndLocal[0]]);
    lower[0]=vstart[0]-fSegmentation[0]->Dpx(isec);
    lower[1]=vstart[1]-fSegmentation[0]->Dpy(isec);
    
    upper[0]=lower[0]+2.*fSegmentation[0]->Dpx(isec);
    upper[1]=lower[1]+2.*fSegmentation[0]->Dpy(isec);
    
    isec=fSegmentation[0]->Sector(fIx[fIndLocal[1]], fIy[fIndLocal[1]]);
    lower[2]=vstart[2]-fSegmentation[0]->Dpx(isec)/2;
    lower[3]=vstart[3]-fSegmentation[0]->Dpy(isec)/2;
    
    upper[2]=lower[2]+fSegmentation[0]->Dpx(isec);
    upper[3]=lower[3]+fSegmentation[0]->Dpy(isec);
    
    lower[4]=0.;
    upper[4]=1.;
// step sizes
    static Double_t step[5]={0.0005, 0.0005, 0.0005, 0.0005, 0.01};
    
    fgMyMinuit->mnparm(0,"x1",vstart[0],step[0],lower[0],upper[0],ierflag);
    fgMyMinuit->mnparm(1,"y1",vstart[1],step[1],lower[1],upper[1],ierflag);
    fgMyMinuit->mnparm(2,"x2",vstart[2],step[2],lower[2],upper[2],ierflag);
    fgMyMinuit->mnparm(3,"y2",vstart[3],step[3],lower[3],upper[3],ierflag);
    fgMyMinuit->mnparm(4,"a0",vstart[4],step[4],lower[4],upper[4],ierflag);
// ready for minimisation	
    fgMyMinuit->SetPrintLevel(-1);
    fgMyMinuit->mnexcm("SET OUT", arglist, 0, ierflag);
    arglist[0]= -1;
    arglist[1]= 0;
    
    fgMyMinuit->mnexcm("SET NOGR", arglist, 0, ierflag);
    fgMyMinuit->mnexcm("MIGRAD", arglist, 0, ierflag);
    fgMyMinuit->mnexcm("EXIT" , arglist, 0, ierflag);
// Get fitted parameters
    Double_t xrec[2], yrec[2], qfrac;
    TString chname;
    Double_t epxz, b1, b2;
    Int_t ierflg;
    fgMyMinuit->mnpout(0, chname, xrec[0], epxz, b1, b2, ierflg);	
    fgMyMinuit->mnpout(1, chname, yrec[0], epxz, b1, b2, ierflg);	
    fgMyMinuit->mnpout(2, chname, xrec[1], epxz, b1, b2, ierflg);	
    fgMyMinuit->mnpout(3, chname, yrec[1], epxz, b1, b2, ierflg);	
    fgMyMinuit->mnpout(4, chname, qfrac,   epxz, b1, b2, ierflg);	

    Double_t fmin, fedm, errdef;
    Int_t   npari, nparx, istat;
      
    fgMyMinuit->mnstat(fmin, fedm, errdef, npari, nparx, istat);  
    
    printf("\n fmin %f \n", fmin);
    
//
// One cluster for each maximum
//
    for (j=0; j<2; j++) {
	AliMUONRawCluster cnew;
	cnew.fChi2[0]=Float_t(fmin);
	
	if (fNPeaks == 0) {
	    cnew.fNcluster[0]=-1;
	    cnew.fNcluster[1]=fNRawClusters;
	} else {
	    cnew.fNcluster[0]=fNPeaks;
	    cnew.fNcluster[1]=0;
	}
	cnew.fMultiplicity[0]=0;
	cnew.fX[0]=Float_t(xrec[j]);
	cnew.fY[0]=Float_t(yrec[j]);
	if (j==0) {
	    cnew.fQ[0]=Int_t(fgChargeTot*qfrac);
	} else {
	    cnew.fQ[0]=Int_t(fgChargeTot*(1-qfrac));
	}
	fgSegmentation->SetHit(xrec[j],yrec[j]);
	for (i=0; i<fMul; i++) {
	    cnew.fIndexMap[cnew.fMultiplicity[0]][0]=c->fIndexMap[i][0];
	    fgSegmentation->SetPad(fgix[i], fgiy[i]);
	    Float_t q1=fgResponse->IntXY(fgSegmentation);
	    cnew.fContMap[cnew.fMultiplicity[0]][0]=(q1*cnew.fQ[0])/Float_t(fQ[i]);
	    cnew.fMultiplicity[0]++;
	}
	FillCluster(&cnew,0);
	cnew.fClusterType=cnew.PhysicsContribution();
	AddRawCluster(cnew);
	fNPeaks++;
    }
    return kTRUE;
}

AliMUONClusterFinder& AliMUONClusterFinder
::operator = (const AliMUONClusterFinder& rhs)
{
// Dummy assignment operator
    return *this;
}


Double_t sinoid(Double_t *x, Double_t *par)
{
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
    Double_t arg1 = (x[0]-par[1])/0.18;
    Double_t arg2 = (x[0]-par[3])/0.18;
    Double_t fitval= par[0]*TMath::Exp(-arg1*arg1/2)
	+par[2]*TMath::Exp(-arg2*arg2/2);
    return fitval;
 }

Float_t DiscrCharge1(Int_t i,Double_t *par) 
{
// par[0]    x-position of cluster
// par[1]    y-position of cluster

   fgSegmentation->SetPad(fgix[i], fgiy[i]);
//  First Cluster
   fgSegmentation->SetHit(par[0],par[1]);
   Float_t q1=fgResponse->IntXY(fgSegmentation);
    
   Float_t value = fgQtot*q1;
   return value;
}


Float_t DiscrCharge2(Int_t i,Double_t *par) 
{
// par[0]    x-position of first  cluster
// par[1]    y-position of first  cluster
// par[2]    x-position of second cluster
// par[3]    y-position of second cluster
// par[4]    charge fraction of first  cluster
// 1-par[4]  charge fraction of second cluster

   fgSegmentation->SetPad(fgix[i], fgiy[i]);
//  First Cluster
   fgSegmentation->SetHit(par[0],par[1]);
   Float_t q1=fgResponse->IntXY(fgSegmentation);
    
//  Second Cluster
   fgSegmentation->SetHit(par[2],par[3]);
   Float_t q2=fgResponse->IntXY(fgSegmentation);
    
   Float_t value = fgQtot*(par[4]*q1+(1.-par[4])*q2);
   return value;
}

//
// Minimisation functions
// Single Mathieson
void fcn1(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    Int_t i;
    Float_t delta;
    Float_t chisq=0;
    Float_t qcont=0;
    Float_t qtot=0;
    
    for (i=0; i<fgNbins; i++) {
	Float_t q0=fgCharge[i];
	Float_t q1=DiscrCharge1(i,par);
	delta=(q0-q1)/TMath::Sqrt(q0);
	chisq+=delta*delta;
	qcont+=q1;
	qtot+=q0;
    }
    f=chisq;
}

// Double Mathieson
void fcn2(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
    Int_t i;
    Float_t delta;
    Float_t chisq=0;
    Float_t qcont=0;
    Float_t qtot=0;
    
    for (i=0; i<fgNbins; i++) {

	Float_t q0=fgCharge[i];
	Float_t q1=DiscrCharge2(i,par);
	delta=(q0-q1)/TMath::Sqrt(q0);
	chisq+=delta*delta;
	qcont+=q1;
	qtot+=q0;
    }
//    chisq=chisq+=(qtot-qcont)*(qtot-qcont)*0.5;
    f=chisq;
}











