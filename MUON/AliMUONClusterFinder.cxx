#include "AliMUONClusterFinder.h"
#include "TTree.h"
#include "AliRun.h"
#include <TCanvas.h>
#include <TH1.h>
#include <TPad.h>
#include <TGraph.h> 
#include <TPostScript.h> 
#include <TMinuit.h> 

ClassImp(AliMUONRecCluster)
//_____________________________________________________________________
static AliMUONsegmentation* gSegmentation;
static AliMUONresponse*     gResponse;
static Int_t                gix[500];
static Int_t                giy[500];
static Float_t              gCharge[500];
static Int_t                gNbins;
static Int_t                gFirst=kTRUE;
static TMinuit *gMyMinuit ;
void fcn(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag);
static Int_t                gChargeTot;


AliMUONRecCluster::AliMUONRecCluster()

{
    fDigits=0;
    fNdigit=-1;
}

AliMUONRecCluster::
AliMUONRecCluster(Int_t FirstDigit,Int_t Ichamber, Int_t Icathod)
{
    fX = 0.;
    fY = 0.;
    fDigits = new TArrayI(10);
    fNdigit=0;
    AddDigit(FirstDigit);
    fChamber=Ichamber;
    fCathod=Icathod;
}

void AliMUONRecCluster::AddDigit(Int_t Digit)
{
    if (fNdigit==fDigits->GetSize()) {
	//enlarge the list by hand!
	Int_t *array= new Int_t[fNdigit*2];
	for(Int_t i=0;i<fNdigit;i++)
	    array[i] = fDigits->At(i);
	fDigits->Adopt(fNdigit*2,array);
    }
    fDigits->AddAt(Digit,fNdigit);
    fNdigit++;
}


AliMUONRecCluster::~AliMUONRecCluster()
{
    if (fDigits)
	delete fDigits;
}

Int_t AliMUONRecCluster::FirstDigitIndex()
{
    fCurrentDigit=0;
    return fDigits->At(fCurrentDigit);
}

Int_t AliMUONRecCluster::NextDigitIndex()
{
    fCurrentDigit++;
    if (fCurrentDigit<fNdigit)
	return fDigits->At(fCurrentDigit);
    else 
	return InvalidDigitIndex();
}

Int_t AliMUONRecCluster::NDigits()
{
    return fNdigit;
}
void AliMUONRecCluster::Finish()
{
    // In order to reconstruct coordinates, one has to
    // get back to the digits which is not trivial here,
    // because we don't know where digits are stored!
    // Center Of Gravity, or other method should be
    // a property of AliMUON class!
}


//----------------------------------------------------------
ClassImp(AliMUONClusterFinder)

    AliMUONClusterFinder::AliMUONClusterFinder
(AliMUONsegmentation *segmentation, AliMUONresponse *response, 
 TClonesArray *digits, Int_t chamber)   
{
    fSegmentation=segmentation;
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
    fSegmentation=0;
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

void AliMUONClusterFinder::AddRawCluster(const AliMUONRawCluster c)
{
  //
  // Add a raw cluster copy to the list
  //
    AliMUON *MUON=(AliMUON*)gAlice->GetModule("MUON");
    MUON->AddRawCluster(fChamber,c); 
    fNRawClusters++;
}



void AliMUONClusterFinder::Decluster(AliMUONRawCluster *cluster)
{
//    AliMUONdigit *dig;
//    Int_t q;

    
    Int_t mul = cluster->fMultiplicity;
//    printf("Decluster - multiplicity   %d \n",mul);

    if (mul == 1 || mul ==2) {
//	printf("\n Nothing special for 1- and 2-clusters \n");
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
	    //	    printf("\n ok, cluster is centered \n");
	  } else {
	    // cluster is not centered, split into 2+1
	    //	    printf("\n cluster is not centered, split into 2+1 \n");
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
      //if (mul < 12) {
	  //	  printf("Decluster - multiplicity > 45   %d \n",mul);
	  //printf("Decluster - multiplicity < 25   %d \n",mul);
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
	//}
    } // multiplicity 
}


Bool_t AliMUONClusterFinder::Centered(AliMUONRawCluster *cluster)
{
    AliMUONdigit* dig;
    dig= (AliMUONdigit*)fDigits->UncheckedAt(cluster->fIndexMap[0]);
    Int_t ix=dig->fPadX;
    Int_t iy=dig->fPadY;
    Int_t nn;
    Int_t X[kMaxNeighbours], Y[kMaxNeighbours], XN[kMaxNeighbours], YN[kMaxNeighbours];
    
    fSegmentation->Neighbours(ix,iy,&nn,X,Y);
    Int_t nd=0;
    for (Int_t i=0; i<nn; i++) {
	if (fHitMap->TestHit(X[i],Y[i]) == used) {
	    XN[nd]=X[i];
	    YN[nd]=Y[i];
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
    } else if (nd ==1) {
//
// Highest signal on an edge, split cluster into 2+1
//
// who is the neighbour ?
	Int_t nind=fHitMap->GetHitIndex(XN[0], YN[0]);
	Int_t i1= (nind==cluster->fIndexMap[1]) ? 1:2;
	Int_t i2= (nind==cluster->fIndexMap[1]) ? 2:1;    
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
void AliMUONClusterFinder::SplitByLocalMaxima(AliMUONRawCluster *c)
{
    AliMUONdigit* dig[100], *digt;
    Int_t ix[100], iy[100], q[100];
    Float_t x[100], y[100];
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
	dig[i]= (AliMUONdigit*)fDigits->UncheckedAt(c->fIndexMap[i]);
	ix[i]= dig[i]->fPadX;
	iy[i]= dig[i]->fPadY;
	q[i] = dig[i]->fSignal;
	fSegmentation->GetPadCxy(ix[i], iy[i], x[i], y[i]);
    }
//
//  Find local maxima
//
    Bool_t IsLocal[100];
    Int_t NLocal=0;
    Int_t AssocPeak[100];
    Int_t IndLocal[100];
    Int_t nn;
    Int_t X[kMaxNeighbours], Y[kMaxNeighbours];
    for (i=0; i<mul; i++) {
	fSegmentation->Neighbours(ix[i], iy[i], &nn, X, Y);
	IsLocal[i]=kTRUE;
	for (j=0; j<nn; j++) {
	    if (fHitMap->TestHit(X[j], Y[j])==empty) continue;
	    digt=(AliMUONdigit*) fHitMap->GetHit(X[j], Y[j]);
	    if (digt->fSignal > q[i]) {
		IsLocal[i]=kFALSE;
		break;
//
// handle special case of neighbouring pads with equal signal
	    } else if (digt->fSignal == q[i]) {
		if (NLocal >0) {
		    for (Int_t k=0; k<NLocal; k++) {
			if (X[j]==ix[IndLocal[k]] && Y[j]==iy[IndLocal[k]]){
			    IsLocal[i]=kFALSE;
			}
		    }
		}
	    } 
	} // loop over next neighbours
	// Maxima should not be on the edge
	if (IsLocal[i]) {
	    IndLocal[NLocal]=i;
	    NLocal++;
	} 
    } // loop over all digits
//    printf("Found %d local Maxima",NLocal);
//
// If only one local maximum found but multiplicity is high 
// take global maximum from the list of digits.    
    if (NLocal==1 && mul>12) {
	Int_t nnew=0;
	for (i=0; i<mul; i++) {
	    if (!IsLocal[i]) {
		IndLocal[NLocal]=i;
		IsLocal[i]=kTRUE;
		NLocal++;
		nnew++;
	    }
	    if (nnew==1) break;
	}
    }
    
// If number of local maxima is 2 try to fit a double gaussian
    if (NLocal==2) {
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
	vstart[0]=x[IndLocal[0]];
	vstart[1]=y[IndLocal[0]];	
	vstart[2]=x[IndLocal[1]];
	vstart[3]=y[IndLocal[1]];	
	vstart[4]=Float_t(q[IndLocal[0]])/
	    Float_t(q[IndLocal[0]]+q[IndLocal[1]]);
// lower and upper limits
	static Double_t lower[5], upper[5];
	Int_t isec=fSegmentation->Sector(ix[IndLocal[0]], iy[IndLocal[0]]);
	lower[0]=vstart[0]-fSegmentation->Dpx(isec)/2;
	lower[1]=vstart[1]-fSegmentation->Dpy(isec)/2;
//	lower[1]=vstart[1];
	
	upper[0]=lower[0]+fSegmentation->Dpx(isec);
	upper[1]=lower[1]+fSegmentation->Dpy(isec);
//	upper[1]=vstart[1];
	
	isec=fSegmentation->Sector(ix[IndLocal[1]], iy[IndLocal[1]]);
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
	printf("\n %f %f %f %f %f\n", xrec[0], yrec[0], xrec[1], yrec[1],qfrac);
//	delete gMyMinuit;
	
	
 //
 // One cluster for each maximum
 //
	 for (j=0; j<2; j++) {
	     AliMUONRawCluster cnew;
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
	     gSegmentation->SetHit(xrec[j],yrec[j]);
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

     if (NLocal !=2 || !fitted) {
 // Check if enough local clusters have been found,
 // if not add global maxima to the list 
 //
	 Int_t nPerMax;
	 if (NLocal!=0) {
	     nPerMax=mul/NLocal;
	 } else {
	     printf("\n Warning, no local maximum found \n");
	     nPerMax=fNperMax+1;
	 }

	 if (nPerMax > fNperMax) {
	     Int_t nGlob=mul/fNperMax-NLocal+1;
	     if (nGlob > 0) {
		 Int_t nnew=0;
		 for (i=0; i<mul; i++) {
		     if (!IsLocal[i]) {
			 IndLocal[NLocal]=i;
			 IsLocal[i]=kTRUE;
			 NLocal++;
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
	     if (IsLocal[i]) continue;
	     for (j=0; j<NLocal; j++) {
		 Int_t il=IndLocal[j];
		 Float_t d=TMath::Sqrt((x[i]-x[il])*(x[i]-x[il])
				   +(y[i]-y[il])*(y[i]-y[il]));
		 Float_t ql=q[il];
 //
 // Select nearest peak
 //
		 if (d<dmin) {
		     dmin=d;
		     qmax=ql;
		     AssocPeak[i]=j;
		 } else if (d==dmin) {
 //
 // If more than one take highest peak
 //
		     if (ql>qmax) {
			 dmin=d;
			 qmax=ql;
			 AssocPeak[i]=j;
		     }
		 }
	     }
	 }


 //
 // One cluster for each maximum
 //
	 for (j=0; j<NLocal; j++) {
	     AliMUONRawCluster cnew;
	     if (fNPeaks == 0) {
		 cnew.fNcluster[0]=-1;
		 cnew.fNcluster[1]=fNRawClusters;
	     } else {
		 cnew.fNcluster[0]=fNPeaks;
		 cnew.fNcluster[1]=0;
	     }
	     cnew.fIndexMap[0]=c->fIndexMap[IndLocal[j]];
	     cnew.fMultiplicity=1;
	     for (i=0; i<mul; i++) {
		 if (IsLocal[i]) continue;
		 if (AssocPeak[i]==j) {
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


void  AliMUONClusterFinder::FillCluster(AliMUONRawCluster* c, Int_t flag) 
{
//
//  Completes cluster information starting from list of digits
//
    AliMUONdigit* dig;
    Float_t x, y;
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
	dig= (AliMUONdigit*)fDigits->UncheckedAt(c->fIndexMap[i]);
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
	    fSegmentation->GetPadCxy(ix, iy, x, y);
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
     fSegmentation->GetPadIxy(x, y, ix, iy);
     fSegmentation->GetPadCxy(ix, iy, x, y);
     Int_t isec=fSegmentation->Sector(ix,iy);
     TF1* CogCorr = fSegmentation->CorrFunc(isec-1);
     
     if (CogCorr) {
	 Float_t YonPad=(c->fY-y)/fSegmentation->Dpy(isec);
	 c->fY=c->fY-CogCorr->Eval(YonPad, 0, 0);
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
    AliMUONdigit* dig = (AliMUONdigit*) fHitMap->GetHit(i,j);
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
	    Int_t ql=((AliMUONdigit*)fDigits
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
    Float_t x, y;
    fSegmentation->GetPadCxy(i, j, x, y);
    c.fX += q*x;
    c.fY += q*y;
    c.fQ += q;
// Flag hit as taken  
    fHitMap->FlagHit(i,j);
//
//  Now look recursively for all neighbours
//  
    Int_t nn;
    Int_t Xlist[kMaxNeighbours], Ylist[kMaxNeighbours];
    fSegmentation->Neighbours(i,j,&nn,Xlist,Ylist);
    for (Int_t in=0; in<nn; in++) {
	Int_t ix=Xlist[in];
	Int_t iy=Ylist[in];
	if (fHitMap->TestHit(ix,iy)==unused) FindCluster(ix, iy, c);
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

    fHitMap = new AliMUONHitMapA1(fSegmentation, fDigits);

    AliMUONdigit *dig;

    Int_t ndig;
    Int_t nskip=0;
    Int_t ncls=0;
    fHitMap->FillHits();
    for (ndig=0; ndig<fNdigits; ndig++) {
	dig = (AliMUONdigit*)fDigits->UncheckedAt(ndig);
	Int_t i=dig->fPadX;
	Int_t j=dig->fPadY;
	if (fHitMap->TestHit(i,j)==used ||fHitMap->TestHit(i,j)==empty) {
	    nskip++;
	    continue;
	}
	AliMUONRawCluster c;
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
	fSegmentation->GetPadIxy(x, y, ix, iy);
	fSegmentation->GetPadCxy(ix, iy, x, y);
	Int_t isec=fSegmentation->Sector(ix,iy);
	TF1* CogCorr=fSegmentation->CorrFunc(isec-1);
	if (CogCorr) {
	    Float_t YonPad=(c.fY-y)/fSegmentation->Dpy(isec);
	    c.fY=c.fY-CogCorr->Eval(YonPad,0,0);
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

void AliMUONClusterFinder::
CalibrateCOG()
{
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


void AliMUONClusterFinder::
SinoidalFit(Float_t x, Float_t y, TF1 &func)
{
//
    static Int_t count=0;
    char canvasname[3];
    count++;
    sprintf(canvasname,"c%d",count);

    const Int_t ns=101;
    Float_t xg[ns], yg[ns], xrg[ns], yrg[ns];
    Float_t xsig[ns], ysig[ns];
   
    AliMUONsegmentation *segmentation=fSegmentation;

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
    Float_t dy=segmentation->Dpy(isec)/(ns-1);

    for (i=0; i<ns; i++) {
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
    Float_t dx=segmentation->Dpx(isec)/(ns-1);

    for (i=0; i<ns; i++) {
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
    //    TGraph *graphx = new TGraph(ns,xg ,xsig);
    //    TGraph *graphxr= new TGraph(ns,xrg,xsig);   
    //    TGraph *graphy = new TGraph(ns,yg ,ysig);
    TGraph *graphyr= new TGraph(ns,yrg,ysig);

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
    gSegmentation->SetHit(par[0],par[1]);
    Float_t q1=gResponse->IntXY(gSegmentation);
    
//  Second Cluster
    gSegmentation->SetHit(par[2],par[3]);
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











