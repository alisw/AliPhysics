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

#include "AliMUONClusterFinderV0.h"
#include "AliMUONSegResV1.h"
//#include "TTree.h"
//#include "AliRun.h"
//#include <TCanvas.h>
//#include <TH1.h>
//#include <TPad.h>
//#include <TGraph.h> 

//----------------------------------------------------------
ClassImp(AliMUONClusterFinderV0)

    AliMUONClusterFinderV0::AliMUONClusterFinderV0
(AliMUONSegmentation *segmentation, AliMUONResponse *response, 
 TClonesArray *digits, Int_t chamber) : AliMUONClusterFinder(segmentation,response,digits,chamber)
{;}

    AliMUONClusterFinderV0::AliMUONClusterFinderV0():AliMUONClusterFinder()
{;}

/*
void AliMUONClusterFinder::AddRawCluster(const AliMUONRawCluster c)
{
  //
  // Add a raw cluster copy to the list
  //
    AliMUON *MUON=(AliMUON*)gAlice->GetModule("MUON");
    MUON->AddRawCluster(fChamber,c); 
    fNRawClusters++;
}
*/



void AliMUONClusterFinderV0::Decluster(AliMUONRawCluster *cluster)
{
//    AliMUONDigit *dig;
//    Int_t q;
    static int done=0;
    if (!done) {
	printf("Calling decluster\n");
	done=1;
    }
    

    
    Int_t mul = cluster->fMultiplicity;
//    printf("Decluster - multiplicity   %d \n",mul);

    if (mul == 1) {
//	printf("\n Nothing special for 1-clusters \n");
//
// Nothing special for 1-clusters
//
	AddRawCluster(*cluster); 
    } else if (mul ==2) {
//
// 2-cluster, compute offset
//
        SetOffset(cluster);
	FillCluster(cluster);
        AddRawCluster(*cluster); 
    } else if (mul ==3) {
//
// 3-cluster, check topology
//	printf("\n 3-cluster, check topology \n");
//
	if (Centered(cluster)) {
//
// ok, cluster is centered 
//	    printf("\n ok, cluster is centered \n");
	} else {
//
// cluster is not centered, split into 2+1
//	    printf("\n cluster is not centered, split into 2+1 \n");
	}
	    
    } else {
	if (mul >(50-5)) printf("Decluster - multiplicity %d approaching 50\n",mul);
// 
// 4-and more-pad clusters
//
	SplitByLocalMaxima(cluster);
    } // multiplicity 
}

Int_t AliMUONClusterFinderV0::PeakOffsetAndCoordinates(Int_t DigitIndex, Float_t *X, Float_t *Y)
//
// Computes for which allowed offsets the digit has the highest neighbouring charge
// Returns the value of the offset, and sets the pyisical coordinates of that pad
// Loop on physical neighbours is specific to AliMUONSegmentationV1
{
Int_t nPara, offset, returnOffset=0 ;
AliMUONDigit* dig= (AliMUONDigit*)fDigits->UncheckedAt(DigitIndex);
AliMUONSegmentationV1* seg = (AliMUONSegmentationV1*) fSegmentation;
seg->GetNParallelAndOffset(dig->fPadX,dig->fPadY,&nPara,&offset);
if (nPara>1)
  {
  Float_t qMax=0;
  for (Int_t i=0;i<nPara; i++)
    {
    // Compute the charge on the 9 neighbouring pads
    // We assume that there are no pads connected in parallel in the neighbourhood
    //
    Float_t q=0;
    for (Int_t dx=-1;dx<2;dx++)
      for (Int_t dy=-1;dy<2;dy++)
        {
        if (dx==dy && dy==0)
          continue; 
        Int_t padY=dig->fPadY+dy;
        Int_t padX=seg->Ix((Int_t) (dig->fPadX+dx+i*offset) , padY);
	if (fHitMap->TestHit(padX, padY)==empty) 
	   continue;
        AliMUONDigit* digt = (AliMUONDigit*) fHitMap->GetHit(padX,padY);
        q += digt->fSignal;
        }
    if (q>qMax)
      {
      returnOffset=i*offset;
      qMax=q;
      }
    }
  }
fSegmentation->GetPadCxy(dig->fPadX+returnOffset,dig->fPadY,*X,*Y);
return returnOffset;
}


void AliMUONClusterFinderV0::SetOffset(AliMUONRawCluster *cluster)
// compute the offsets assuming that there is only one peak !
{
//DumpCluster(cluster);
Float_t X,Y;
cluster->fOffsetMap[0]=PeakOffsetAndCoordinates(cluster->fIndexMap[0],&X,&Y);
for (Int_t i=1;i<cluster->fMultiplicity;i++) {
  AliMUONDigit* dig= (AliMUONDigit*)fDigits->UncheckedAt(cluster->fIndexMap[i]);
  fSegmentation->Distance2AndOffset(dig->fPadX,dig->fPadY,X,Y,&(cluster->fOffsetMap[i]));
  }
}

void AliMUONClusterFinderV0::DumpCluster(AliMUONRawCluster *cluster)
{    
printf ("other cluster\n");
for (Int_t i=0; i<cluster->fMultiplicity; i++)
    {
	AliMUONDigit* dig= (AliMUONDigit*)fDigits->UncheckedAt(cluster->fIndexMap[i]);
	Int_t nPara, offset;
	fSegmentation->GetNParallelAndOffset(dig->fPadX,dig->fPadY,&nPara,&offset);

	printf("X %d Y %d Q %d NPara %d \n",dig->fPadX, dig->fPadY,dig->fSignal, nPara);
    }
}

Bool_t AliMUONClusterFinderV0::Centered(AliMUONRawCluster *cluster)
{
    AliMUONDigit* dig;
    dig= (AliMUONDigit*)fDigits->UncheckedAt(cluster->fIndexMap[0]);
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
	SetOffset(cluster);
	FillCluster(cluster);
	AddRawCluster(*cluster);
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
	cnew.fMultiplicity=2;
	cnew.fIndexMap[0]=cluster->fIndexMap[0];
	cnew.fIndexMap[1]=cluster->fIndexMap[i1];
	SetOffset(&cnew);
	FillCluster(&cnew);
	AddRawCluster(cnew);
//
// 1-cluster
	cluster->fMultiplicity=1;
	cluster->fIndexMap[0]=cluster->fIndexMap[i2];
	cluster->fIndexMap[1]=0;
	cluster->fIndexMap[2]=0;	
	FillCluster(cluster);
	AddRawCluster(*cluster);
	return kFALSE;
    } else {
	printf("\n Completely screwed up %d !! \n",nd);
	
    }
    
	return kFALSE;
}


void AliMUONClusterFinderV0::SplitByLocalMaxima(AliMUONRawCluster *c)
{
    AliMUONDigit* dig[50], *digt;
    Int_t ix[50], iy[50], q[50];
    Float_t x[50], y[50];
    Int_t i; // loops over digits
    Int_t j; // loops over local maxima
    
    Int_t mul=c->fMultiplicity;
//
//  dump digit information into arrays
//
    for (i=0; i<mul; i++)
    {
	dig[i]= (AliMUONDigit*)fDigits->UncheckedAt(c->fIndexMap[i]);
	ix[i]= dig[i]->fPadX;
	iy[i]= dig[i]->fPadY;
	q[i] = dig[i]->fSignal;
	fSegmentation->GetPadCxy(ix[i], iy[i], x[i], y[i]);
    }
//
//  Find local maxima
//
    Bool_t IsLocal[50];
    Int_t NLocal=0;
    Int_t AssocPeak[50];
    Int_t IndLocal[50];
    Int_t nn;
    Int_t X[kMaxNeighbours], Y[kMaxNeighbours];
    for (i=0; i<mul; i++) {
	fSegmentation->Neighbours(ix[i], iy[i], &nn, X, Y);
	IsLocal[i]=kTRUE;
	for (j=0; j<nn; j++) {
	    if (fHitMap->TestHit(X[j], Y[j])==empty) continue;
	    digt=(AliMUONDigit*) fHitMap->GetHit(X[j], Y[j]);
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
	if (IsLocal[i]) {
	    IndLocal[NLocal]=i;
            // New for LYON : we guess which is the actual position of the pad hit
	    // But this would run like that for normal chamber !
            c->fOffsetMap[i]=PeakOffsetAndCoordinates(c->fIndexMap[i], &(x[i]), &(y[i]));
	    NLocal++;
	} 
    } // loop over all digits
//    printf("Found %d local Maxima",NLocal);
//
// Associate hits to peaks
//
    for (i=0; i<mul; i++) {
        //
        // loop on digits
        //
	Float_t dmin=1.E10;
	Float_t qmax=0;
	Int_t offset;
	if (IsLocal[i]) continue;
	for (j=0; j<NLocal; j++) {
            //
            // Loop on peaks
            //
	    Int_t il=IndLocal[j];
//	    Float_t d=TMath::Sqrt((x[i]-x[il])*(x[i]-x[il])
//				  +(y[i]-y[il])*(y[i]-y[il]));
            // Can run like that for non-Lyon chambers
	    Float_t d = fSegmentation->Distance2AndOffset(ix[i],iy[i],x[il],y[il], &offset);
	    Float_t ql=q[il];
//
// Select nearest peak
//
	    if (d<dmin) {
		dmin=d;
		qmax=ql;
		AssocPeak[i]=j;
		c->fOffsetMap[i]=offset;
	    } else if (d==dmin) {
//
// If more than one take highest peak
//
		if (ql>qmax) {
		    dmin=d;
		    qmax=ql;
		    AssocPeak[i]=j;
		    c->fOffsetMap[i]=offset;
		}
	    } // end if
	} // End loop on peaks
    } // end loop on digits
//
// One cluster for each maximum
//
    for (j=0; j<NLocal; j++) {
	AliMUONRawCluster cnew;
	cnew.fIndexMap[0]=c->fIndexMap[IndLocal[j]];
	cnew.fOffsetMap[0]=c->fOffsetMap[IndLocal[j]];
	cnew.fMultiplicity=1;
	for (i=0; i<mul; i++) {
	    if (IsLocal[i]) continue;
	    if (AssocPeak[i]==j) {
		cnew.fIndexMap[cnew.fMultiplicity]=c->fIndexMap[i];
		cnew.fOffsetMap[cnew.fMultiplicity]=c->fOffsetMap[i];
		cnew.fMultiplicity++;
	    }
	}
	FillCluster(&cnew);
	AddRawCluster(cnew);
    }
}

/*
void  AliMUONClusterFinderV0::FillCluster(AliMUONRawCluster* c)
{
//
//  Completes cluster information starting from list of digits
//
    AliMUONDigit* dig;
    Float_t x, y;
    Int_t  ix, iy;
    
    c->fPeakSignal=0;
    c->fX=0;
    c->fY=0;
    c->fQ=0;
    for (Int_t i=0; i<c->fMultiplicity; i++)
    {
	dig= (AliMUONDigit*)fDigits->UncheckedAt(c->fIndexMap[i]);
	ix=dig->fPadX + c.fOffsetMap[i]; // should be 0 for non-LYON
	iy=dig->fPadY;
	Int_t q=dig->fSignal;
//
//
// peak signal and track list
	if (q>c->fPeakSignal) {
	    c->fPeakSignal=0;
	    c->fTracks[0]=dig->fTracks[0];
	    c->fTracks[1]=dig->fTracks[1];
	    c->fTracks[2]=dig->fTracks[2];
	}
//
// centre of gravity
	fSegmentation->GetPadCxy(ix, iy, x, y);
	c->fX += q*x;
	c->fY += q*y;
	c->fQ += q;
    }
    c->fX/=c->fQ;
// Not valid for inclined tracks in X !!! (Manu)
//    c->fX=fSegmentation->GetAnod(c->fX);
    c->fY/=c->fQ; 
//
//  apply correction to the coordinate along the anode wire
//
    if (fCogCorr) {
	x=c->fX;   
	y=c->fY;
	fSegmentation->GetPadIxy(x, y, ix, iy);
	fSegmentation->GetPadCxy(ix, iy, x, y);
	Float_t YonPad=(c->fY-y)/fSegmentation->Dpy();
	c->fY=y-fCogCorr->Eval(YonPad, 0, 0);
    }

}
*/
/*
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
    if (q > TMath::Abs(c.fPeakSignal)) {
	c.fPeakSignal=q;
	c.fTracks[0]=dig->fTracks[0];
	c.fTracks[1]=dig->fTracks[1];
	c.fTracks[2]=dig->fTracks[2];
    }
//
//  Make sure that list of digits is ordered 
// 
    Int_t mu=c.fMultiplicity;
    c.fIndexMap[mu]=idx;

    if (mu > 0) {
	for (Int_t ind=mu-1; ind>=0; ind--) {
	    Int_t ist=(c.fIndexMap)[ind];
	    Int_t ql=((AliMUONDigit*)fDigits
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
	c.fMultiplicity=50-1;
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
*/

//_____________________________________________________________________________

void AliMUONClusterFinderV0::FindRawClusters()
{
  //
  // simple MUON cluster finder from digits -- finds neighbours and 
  // fill the tree with raw clusters
  //
    if (!fNdigits) return;

    fHitMap = new AliMUONHitMapA1(fSegmentation, fDigits);

    AliMUONDigit *dig;

    int ndig;
    int nskip=0;

    fHitMap->FillHits();
    for (ndig=0; ndig<fNdigits; ndig++) {
	dig = (AliMUONDigit*)fDigits->UncheckedAt(ndig);
	Int_t i=dig->fPadX;
	Int_t j=dig->fPadY;
	if (fHitMap->TestHit(i,j)==used ||fHitMap->TestHit(i,j)==empty) {
	    nskip++;
	    continue;
	}
	AliMUONRawCluster c;
	c.fMultiplicity=0;
//	c.fPeakSignal=dig->fSignal;
//	c.fTracks[0]=dig->fTracks[0];
//	c.fTracks[1]=dig->fTracks[1];
//	c.fTracks[2]=dig->fTracks[2];
	c.fPeakSignal=0;
	FindCluster(i,j, c);
	// center of gravity
	c.fX /= c.fQ;
	c.fX=fSegmentation->GetAnod(c.fX);
	c.fY /= c.fQ;
//
//  apply correction to the coordinate along the anode wire
//

	
    if (fCogCorr) {
	Int_t ix,iy;
	Float_t x=c.fX;   
	Float_t y=c.fY;
	fSegmentation->GetPadIxy(x, y, ix, iy);
	fSegmentation->GetPadCxy(ix, iy, x, y);
	Float_t YonPad=(c.fY-y)/fSegmentation->Dpy();
	c.fY=y-fCogCorr->Eval(YonPad,0,0);
    }
//
//      Analyse cluster and decluster if necessary
//	
	Decluster(&c);
//
//
//
//      reset Cluster object
	for (int k=0;k<c.fMultiplicity;k++) {
	    c.fIndexMap[k]=0;
 	    c.fOffsetMap[k]=0;
	}
	c.fMultiplicity=0;
    } // end loop ndig    
    delete fHitMap;
}

/*

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
	}
	fCogCorr = new TF1(func);
    }
}
*/
/*

void AliMUONClusterFinder::
SinoidalFit(Float_t x, Float_t y, TF1 &func)
{
//
    static Int_t count=0;
    char canvasname[3];
    count++;
    sprintf(canvasname,"c%d",count);

// MANU : without const, error on HP
    const Int_t ns=101;
    Float_t xg[ns], yg[ns], xrg[ns], yrg[ns];
    Float_t xsig[ns], ysig[ns];
   
    AliMUONSegmentation *segmentation=fSegmentation;

    Int_t ix,iy;
    segmentation->GetPadIxy(x,y,ix,iy);   
    segmentation->GetPadCxy(ix,iy,x,y);   
    Int_t isec=segmentation->Sector(ix,iy);
// Pad Limits    
    Float_t xmin = x-segmentation->GetRealDpx(isec)/2;
    Float_t ymin = y-segmentation->Dpy()/2;
//      	
//      Integration Limits
    Float_t dxI=fResponse->Nsigma()*fResponse->ChwX();
    Float_t dyI=fResponse->Nsigma()*fResponse->ChwY();

//
//  Scanning
//
    Int_t i;
    Float_t qp;
//
//  y-position
    Float_t yscan=ymin;
    Float_t dy=segmentation->Dpy()/(ns-1);

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
	yg[i]=(yscan-y)/segmentation->Dpy();
	yrg[i]=(ycog-y)/segmentation->Dpy();
	ysig[i]=ycog-yscan;
	yscan+=dy;
    } // scan loop
//
//  x-position
    Float_t xscan=xmin;
    Float_t dx=segmentation->GetRealDpx(isec)/(ns-1);

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
	
	xg[i]=(xscan-x)/segmentation->GetRealDpx(isec);
	xrg[i]=(xcog-x)/segmentation->GetRealDpx(isec);
	xsig[i]=xcog-xscan;
	xscan+=dx;
    }

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
   TGraph *graphx = new TGraph(ns,xg ,xsig);
   TGraph *graphxr= new TGraph(ns,xrg,xsig);   
   TGraph *graphy = new TGraph(ns,yg ,ysig);
   TGraph *graphyr= new TGraph(ns,yrg,ysig);
//
// Creates a Root function based on function sinoid above
// and perform the fit
//
   Double_t sinoid(Double_t *x, Double_t *par);
   TF1 *sinoidf = new TF1("sinoidf",sinoid,0.5,0.5,5);
   graphyr->Fit("sinoidf","V");
   sinoidf->Copy(func);
   func.Eval(0,0,0);
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
   
}
*/
/*
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
*/







