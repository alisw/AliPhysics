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
Revision 1.2  2000/06/30 12:07:49  kowal2
Updated from the TPC-PreRelease branch

Revision 1.1.2.1  2000/06/25 08:52:51  kowal2
replacing AliClusterFinder

*/

//-----------------------------------------------------------------------------
//
//  Implementation of class ALICLUSTERFINDER
// 
//Class for cluster finding in two dimension.
//In the present there are implemented two algorithm
//primitive recursion algorithm. (FindPeaks) 
//Algorithm is not working in case of overlaping clusters
//Maximum - minimum in direction algoritm (Find clusters)
//In this algoritm we suppose that each cluster has local 
//maximum. From this local maximum I mus see each point 
//of cluster.
//From maximum i can accept every point in radial 
//direction which is before border in direction
//Border in direction occur if we have next in
//direction nder threshold or response begin
//to increase in given radial direction
//-----------------------------------------------------------------------------

#include "TMinuit.h"
#include "AliArrayI.h"
#include "TClonesArray.h"
#include "AliTPC.h"
#include "TRandom.h"
#include "AliH2F.h"
#include "TMarker.h"
#include "AliCluster.h"
#include "AliTPCClusterFinder.h"
#include <fstream.h>

//direction constants possible direction in 8 different sectors
//


const Int_t kClStackSize =1000;




static AliTPCClusterFinder * gClusterFinder; //for fitting routine

void gauss(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag)
{
  AliArrayI * points = gClusterFinder->GetStack();
  const Int_t nbins = gClusterFinder->GetStackIndex();
  Int_t i;
  //calculate chisquare
   Double_t chisq = 0;
   Double_t delta;
   for (i=0;i<nbins; i++) {
     Float_t x = points->At(i*3);
     Float_t y = points->At(i*3+1);
     Float_t z = points->At(i*3+2);
     Float_t deltax2 = (x-par[1]);
     deltax2*=deltax2;
     deltax2*=par[3];
     Float_t deltay2 = (y-par[2]);
     deltay2*=deltay2;
     deltay2*=par[4];
     
     delta  = z-par[0]*TMath::Exp(-deltax2-deltay2);
     chisq += delta*delta;
   }
   f = chisq;   
}


ClassImp(AliTPCClusterFinder)
  //ClassImp(AliCell)

AliTPCClusterFinder::AliTPCClusterFinder()
{
  fDigits =0;
  fDimX = 0;
  fDimY = 0;
  fNoiseTh = 3;
  fMulSigma2 = 16; //4 sigma
  fDirSigmaFac = 1.4;
  fDirAmpFac =1.3;
  fNType=8;
  fThreshold = 2;
  fStack = new AliArrayI;
  fStack->Set(kClStackSize);  
  fClustersArray =0;
  SetSigmaX(1,0,0);
  SetSigmaY(1,0,0);


  fDetectorParam = 0;
  ResetStatus(); 
  fBFit = kFALSE;
  fMinuit= new TMinuit(5);
  fMinuit->SetFCN(gauss);
  gClusterFinder = this;
  
}


AliTPCClusterFinder::~AliTPCClusterFinder()
{
 if (fDigits  != 0) delete fDigits;
}

void AliTPCClusterFinder::SetSigmaX(Float_t s0, Float_t s1x, Float_t s1y)
{
  fSigmaX[0]=s0;
  fSigmaX[1]=s1x;
  fSigmaX[2]=s1y;

}
void AliTPCClusterFinder::SetSigmaY(Float_t s0, Float_t s1x, Float_t s1y)
{
  fSigmaY[0]=s0;
  fSigmaY[1]=s1x;
  fSigmaY[2]=s1y;
}



Bool_t AliTPCClusterFinder::SetSigma2(Int_t i, Int_t j, Float_t & sigmax2, Float_t &sigmay2)
{
  //
  //set sigmax2 and sigma y2  accordig i and j position of cell 
  //
  
  //  Float_t x[3] = {ItoX(i),JtoY(j),0};
  Float_t x= ItoX(i);
  Float_t y= JtoY(j);

  sigmax2= fSigmaX[0]+fSigmaX[1]*x+fSigmaX[2]*y;
  sigmay2= fSigmaY[0]+fSigmaY[1]*x+fSigmaY[2]*y;
  return kTRUE;  
}

/*
Bool_t AliTPCClusterFinder::SetSigma2(Int_t i, Int_t j, Float_t & sigmax2, Float_t &sigmay2)
{
  //
  //set sigmax2 and sigma y2  accordig i and j position of cell 
  //
  if (fDetectorParam==0) {
    sigmax2=1;
    sigmay2=1;
    return kFALSE;
  }
  Float_t x[3] = {ItoX(i),JtoY(j),0};
  Float_t sigma[2];
  fDetectorParam->GetClusterSize(x,fDetectorIndex,0,0,sigma);
  sigmax2=sigma[0]*(fX2-fX1)*(fX2-fX1)/(fDimX*fDimX);
  sigmay2=sigma[1]*(fY2-fY1)*(fY2-fY1)/(fDimY*fDimY);
  return kTRUE;
}
*/


void AliTPCClusterFinder::GetHisto(TH2F * his2)
{
  
  UInt_t idim =his2->GetNbinsX();
  UInt_t jdim =his2->GetNbinsY();
  fX1 = his2->GetXaxis()->GetXmin();
  fX2 = his2->GetXaxis()->GetXmax();
  fY1 = his2->GetYaxis()->GetXmin();
  fY2 = his2->GetYaxis()->GetXmax();
 
  if ( (idim>0) && (jdim>0))
    {
      rOK = kTRUE;
      fDimX = idim;
      fDimY = jdim;
      Int_t size =idim*jdim;       
      if (fDigits !=0) delete fDigits;
      fDigits  = (Int_t*) new Int_t[size];
      fCells  = (AliCell*) new AliCell[size];

    }  else 
      rOK=kFALSE;
  for (Int_t i = 0; i<(Int_t)idim;i++)    
    for (Int_t j = 0; j<(Int_t)jdim;j++)
      {
	Int_t index = his2->GetBin(i+1,j+1);
	//AliCell * cell = GetCell(i,j);
	//if (cell!=0) cell->SetSignal(his2->GetBinContent(index));
	SetSignal(his2->GetBinContent(index),i,j);
      }
   
}




void AliTPCClusterFinder::FindMaxima()
{
  for (Int_t i=0; i<fDimX; i++) 
    for (Int_t j=0;j<fDimY; j++)      
      if (IsMaximum(i,j)) cout<<i<<"   "<<j<<"\n"; 		     
}

 
void  AliTPCClusterFinder::Transform(AliDigitCluster * c)
{
  //transform coordinata from bin coordinata to "normal coordinata"
  //for example if we initialize finder with histogram
  //it transform values from bin coordinata to the histogram coordinata
  c->fX=ItoX(c->fX);
  c->fY=JtoY(c->fY);
  c->fMaxX=ItoX(c->fMaxX);
  c->fMaxY=JtoY(c->fMaxY);

  c->fSigmaX2=c->fSigmaX2*(fX2-fX1)*(fX2-fX1)/(fDimX*fDimX);
  c->fSigmaY2=c->fSigmaY2*(fY2-fY1)*(fY2-fY1)/(fDimY*fDimY);  
  c->fArea   =c->fArea*(fX2-fX1)*(fY2-fY1)/(fDimX*fDimY); 
}
void  AliTPCClusterFinder::AddToStack(Int_t i, Int_t j, Int_t signal)
{
  //
  //add digit to stack
  //
  if ( ((fStackIndex+2)>=kClStackSize) || (fStackIndex<0) ) return; 
  fStack->AddAt(i,fStackIndex);
  fStack->AddAt(j,fStackIndex+1);
  fStack->AddAt(signal,fStackIndex+2);
  fStackIndex+=3;  
}

void AliTPCClusterFinder::GetClusterStatistic(AliDigitCluster & cluster)
{
  //
  //calculate statistic of cluster 
  //
  Double_t sumxw,sumyw,sumx2w,sumy2w,sumxyw,sumw;
  Int_t minx,maxx,miny,maxy;
  sumxw=sumyw=sumx2w=sumy2w=sumxyw=sumw=0;
  minx=fDimX;
  maxx=-fDimX;
  miny=fDimY;
  maxy=-fDimY;
  Int_t x0=fStack->At(0);
  Int_t y0=fStack->At(1);
  Int_t maxQx =x0;
  Int_t maxQy =y0;  
  Int_t maxQ=fStack->At(2);
  

  for (Int_t i = 0; i<fStackIndex;i+=3){
    Int_t x = fStack->At(i);
    Int_t y = fStack->At(i+1);
    Int_t dx=x-x0;
    Int_t dy=y-y0;
    Int_t w = fStack->At(i+2);
    if (w>maxQ){
      maxQ = w;
      maxQx = x;
      maxQy=y;
    }	
    if (x<minx) minx=x;
    if (y<miny) miny=y;
    if (x>maxx) maxx=x;
    if (y>maxy) maxy=y;   
    sumxw+=dx*w;
    sumyw+=dy*w;
    sumx2w+=dx*dx*w;
    sumy2w+=dy*dy*w;
    sumxyw+=dx*dy*w;
    sumw+=w;    
  }
  cluster.fQ = sumw;
  if (sumw>0){
    cluster.fX = sumxw/sumw;
    cluster.fY = sumyw/sumw;
    cluster.fQ = sumw;
    cluster.fSigmaX2 = sumx2w/sumw-cluster.fX*cluster.fX;
    cluster.fSigmaY2 = sumy2w/sumw-cluster.fY*cluster.fY;
    cluster.fSigmaXY = sumxyw/sumw-cluster.fX*cluster.fY;
    cluster.fMaxX = maxQx;
    cluster.fMaxY = maxQy; 
    cluster.fMax = maxQ;
    cluster.fArea = fStackIndex/3;  
    cluster.fNx = maxx-minx+1;
    cluster.fNy = maxy-miny+1;
    cluster.fX +=x0; 
    cluster.fY +=y0;
  }
}
void AliTPCClusterFinder::GetClusterFit(AliDigitCluster & cluster)
{
  //
  //calculate statistic of cluster 
  //
  Double_t arglist[10];
  Int_t ierflg = 0;
  
  arglist[0] = 1;
  fMinuit->mnexcm("SET ERR", arglist ,1,ierflg);
  
  //fistly find starting parameters
  Int_t minx,maxx,miny,maxy,maxQ,maxQx,maxQy;
  Int_t over =0;
  Float_t sumxw,sumyw,sumw;
  sumxw=sumyw=sumw=0;
  minx=fDimX;
  maxx=-fDimX;
  miny=fDimY;
  maxy=-fDimY;
  maxQx=fStack->At(0);
  maxQy=fStack->At(1);
  maxQ=fStack->At(2);
 
  for (Int_t i = 0; i<fStackIndex;i+=3){
    Int_t x = fStack->At(i);
    Int_t y = fStack->At(i+1);
    Int_t w = fStack->At(i+2);
    if (w>fThreshold) {
      over++;
      sumw+=w;    
      sumxw+=x*w;
      sumyw+=y*w;
      if (x<minx) minx=x;
      if (y<miny) miny=y;
      if (x>maxx) maxx=x;
      if (y>maxy) maxy=y;
      if (w>maxQ) {
	maxQ=w;   
	maxQx=x;
	maxQy=y;
      }    
    }
  }
  Int_t nx = maxx-minx+1;
  Int_t ny = maxy-miny+1;
  
  SetSigma2(maxQx,maxQy,fCurrentSigmaX2,fCurrentSigmaY2);
  Double_t vstart[5]={maxQ,sumxw/sumw,sumyw/sumw,1/(2*fCurrentSigmaX2),1/(2*fCurrentSigmaY2)};
  Double_t step[5]={1.,0.01,0.01,0.01,0.01};
  fMinuit->mnparm(0, "amp", vstart[0], step[0], 0,0,ierflg);
  fMinuit->mnparm(1, "x0", vstart[1], step[1], 0,0,ierflg);
  fMinuit->mnparm(2, "y0", vstart[2], step[2], 0,0,ierflg);
  fMinuit->mnparm(3, "sx2", vstart[3], step[3], 0,0,ierflg);
  fMinuit->mnparm(4, "sy2", vstart[4], step[4], 0,0,ierflg);
  arglist[0] = 500;
  arglist[1] = 1.;

  fMinuit->mnfree(0);  //set unfixed all parameters
  //if we have area less then
  if (over<=21) {  //if we dont't have more then  7  points
    fMinuit->FixParameter(3); 
    fMinuit->FixParameter(4);
  }
  else {
    if (nx<3)  fMinuit->FixParameter(3); //fix sigma x if no data in x direction
    if (ny<3)  fMinuit->FixParameter(4);  //fix sigma y if no data in y direction
  }
  fMinuit->mnexcm("MIGRAD", arglist ,2,ierflg);
 
  if (sumw>0){
    Double_t  x[5];
    Double_t  error[5];
    fMinuit->GetParameter(0,x[0],error[0]);
    fMinuit->GetParameter(1,x[1],error[1]);
    fMinuit->GetParameter(2,x[2],error[2]);
    fMinuit->GetParameter(3,x[3],error[3]);
    fMinuit->GetParameter(4,x[4],error[4]);

    cluster.fX = x[1];
    cluster.fY = x[2];
    cluster.fMaxX = maxQx;
    cluster.fMaxY = maxQy;
    
    cluster.fQ = sumw;
    cluster.fSigmaX2 = 1/TMath::Sqrt(2*x[3]);
    cluster.fSigmaY2 = 1/TMath::Sqrt(2*x[4]);
    cluster.fSigmaXY = 0;
    cluster.fMax = x[0];
    cluster.fArea = over;  
    cluster.fNx = nx;
    cluster.fNy = ny;
  }
}

Bool_t   AliTPCClusterFinder::CheckIfDirBorder(Float_t x, Float_t y, 
					  Int_t i,Int_t j)
{
  //
  //function which control if given cell with index i, j is the 
  //minimum in direction
  // x and y are estimate of local maximum 
  //direction is given by the 
  Float_t virtualcell;
  AliCell * cellor= GetCell(i,j);
  Int_t     sigor = GetSignal(i,j);

  //control derivation in direction
  //if function grows up in direction then there is border
  Float_t dx = i-x;
  Float_t dy = j-y; 
  Float_t dd = TMath::Sqrt(dx*dx+dy*dy);
  Float_t ddx = TMath::Abs(dx);
  ddx =   (ddx>0.5) ? ddx-0.5: 0;
  ddx*=ddx;
  Float_t ddy = TMath::Abs(dy);
  ddy = (ddy>0.5) ? ddy-0.5: 0;
  ddy*=ddy;
  Float_t d2 = ddx/(2*fDirSigmaFac*fCurrentSigmaX2)+ddy/(2*fDirSigmaFac*fCurrentSigmaY2); //safety factor 
  //I accept sigmax and sigma y bigge by factor sqrt(fDirsigmaFac)
  Float_t amp = TMath::Exp(-d2)*fCurrentMaxAmp*fDirAmpFac; //safety factor fDirFac>1

  if (sigor>amp) return kTRUE; 
  if (dd==0) return kFALSE;

  dx/=dd;
  dy/=dd;  
  virtualcell = GetVirtualSignal(i+dx,j+dy);
  if (virtualcell <=fThreshold) return kFALSE;
  if (virtualcell>sigor)
    if (virtualcell>(sigor+fNoiseTh))
      {cellor->SetDirBorder(fIndex); return kTRUE;}
    else
      {
	virtualcell = GetVirtualSignal(i+2*dx,j+2*dy);
	if (virtualcell>sigor)
	  { cellor->SetDirBorder(fIndex); return kTRUE;}       
      };
  return kFALSE;  
}





Bool_t  AliTPCClusterFinder::IsMaximum(Int_t i, Int_t  j)
{
  //there is maximum if given digits is 1 sigma over all adjacent
  //in 8 neighborow 
  //or ther exist virual maximum
  //is maximum on 24 points neighboring
  //  Bool_t res = kFALSE;
  Int_t over =0;
  Int_t overth=0;
  Int_t oversigma =0;
  AliCell * cell = GetCell(i,j); 
  Int_t signal = GetSignal(i,j);
  if (cell == 0) return kFALSE;
  for ( Int_t di=-1;di<=1;di++)
    for ( Int_t dj=-1;dj<=1;dj++){      
      if ( (di!=0) || (dj!=0))
	{
	  AliCell * cell2=GetCell(i+di,j+dj);
	  Int_t signal2 = GetSignal(i+di,j+dj);
	  if (cell2 == 0) {
	    over+=1;
	    oversigma+=1;
	  }
	  else
	    {
	      if (signal2>signal) return kFALSE;
	      if (signal2>fThreshold) overth++;
	      if (signal2==signal) {
		if (di<0) return kFALSE; 
		if ( (di+dj)<0) return kFALSE;
	      }
	      //	      if (signal>=signal2){
	      over+=1;
	      if (signal>fNoiseTh+signal2)	     
		oversigma+=1;		
	      //}
	    }
	}
    }
  //if I have only one neighborough over threshold 
  if (overth<2) return kFALSE;
  if (over<8)   return kFALSE;
  if (oversigma==8) {
    fCurrentMaxX = i;
    fCurrentMaxY = j;
    fCurrentMaxAmp =signal;
    SetMaximum(fIndex,i,j);
    return kTRUE;
  }
  //check if there exist virtual maximum
  for (Float_t ddi=0.;(ddi<1.);ddi+=0.5)
    for (Float_t ddj=0.;(ddj<1.);ddj+=0.5)	   	  
      if (IsVirtualMaximum(Float_t(i)+ddi,Float_t(j)+ddj)){
	fCurrentMaxX = i+ddi;
	fCurrentMaxY = j+ddj;
	fCurrentMaxAmp =signal; 
	SetMaximum(fIndex,i,j);
	return kTRUE;	
      }
  return kFALSE;
}

Bool_t  AliTPCClusterFinder::IsVirtualMaximum(Float_t x, Float_t  y)
{
  //there is maximum if given digits is 1 sigma over all adjacent
  //in 8 neighborow or 
  //is maximum on 24 points neighboring
  Bool_t res = kFALSE;
  Int_t over =0;
  Int_t overth=0;
  Int_t oversigma =0;
  Float_t virtualcell = GetVirtualSignal(x,y); 
  if (virtualcell < 0) return kFALSE;
  for ( Int_t di=-1;di<=1;di++)
    for ( Int_t dj=-1;dj<=1;dj++)
      if ( (di!=0) || (dj!=0))
	{
	  Float_t virtualcell2=GetVirtualSignal(x+di,y+dj);
	  if (virtualcell2 < 0) {
	    over+=1;
	    oversigma+=1;
	  }
	  else
	    {
	      if (virtualcell2>fThreshold) overth++;
	      if (virtualcell>=virtualcell2){
		over+=1;
		if (virtualcell>fNoiseTh+virtualcell2)	     
		  oversigma+=1;
	      }
	    }
	}
  if (overth<2) return kFALSE;
  //if there exist only one or less neighboring above threshold
  if (oversigma==8)  res = kTRUE;
  else if ((over==8)&&(GetNType()==8)) res=kTRUE;
  else if (over ==8 )
    for ( Int_t di=-2;di<=2;di++)
      for ( Int_t dj=-2;dj<=2;dj++)
	if ( (di==2)||(di==-2) || (dj==2)|| (dj==-2) )
	  {
	    Float_t virtualcell2=GetVirtualSignal(x+di,y+dj);
	    if (virtualcell2 < 0) {
	      over+=1;
	      oversigma+=1;
	    }
	    else
	      {
		if (virtualcell>=virtualcell2) over+=1;
	      }
	  }	
  if (over == 24) res=kTRUE;
  return res;
  
}


void AliTPCClusterFinder::ResetSignal()
{
   //reset dignals to 0
  Int_t size = fDimX*fDimY;
  AliCell *dig=fCells;
  if (rOK==kTRUE) for (Int_t i=0 ; i<size;i++) dig[i] = 0; 
}



void AliTPCClusterFinder::ResetStatus()
{
   //reset status of signals to not used
  Int_t size = fDimX*fDimY;
  AliCell *dig=fCells;
  if (rOK==kTRUE) for (Int_t i=0 ; i<size;i++) 
      dig[i].SetStatus(0);     
} 


AliCell  *  AliTPCClusterFinder::GetCell(Int_t i, Int_t j)
{
  //return reference to the cell with index i,j 
  if (rOK == kTRUE)
    if ( (i>=0) && (i<fDimX) && (j>=0) && (j<fDimY) )
      return &fCells[i+j*fDimX];
  return 0; 
}

Float_t   AliTPCClusterFinder::GetVirtualSignal(Float_t ri, Float_t rj)
{
  //it generate virtual cell as mean value from different cels
  //after using it must be destructed !!!  
  Int_t i=(Int_t)ri;
  Int_t j=(Int_t)rj;
  Int_t ddi = (ri>i)? 1:0;
  Int_t ddj = (rj>j)? 1:0;
  Float_t sum = 0;
  Float_t sumw= 0;
  for (Int_t di=0;di<=ddi;di++)   
    for (Int_t dj=0;dj<=ddj;dj++)
      {             
	Float_t w = (ri-i-di)*(ri-i-di)+(rj-j-dj)*(rj-j-dj);	
	if (w>0) w=1/TMath::Sqrt(w);
	else w=9999999;
	AliCell * cel2 =GetCell(i+di,j+dj);
	Int_t signal2 = GetSignal(i+di,j+dj);
        if (cel2!=0) {
	  sumw+=w;
	  sum+= signal2*w;
	}
      }
  if (sumw>0)  return (sum/sumw);
  else 
    return -1;
}



void AliTPCClusterFinder::Streamer(TBuffer & R__b)
{
  if (R__b.IsReading()) {
    //      Version_t R__v = R__b.ReadVersion();
   } else {
      R__b.WriteVersion(AliTPCClusterFinder::IsA());    
   } 
}



void AliTPCClusterFinder::SetBlockIndex(Int_t * index)
{
  //
  //calculate which indexes we must check for border
  //
  if (TMath::Abs(index[0])<2) index[2] = 0;
  else {
    index[2] = TMath::Abs(index[0])-1;
    if (index[0]<0) index[2]*=-1;   //first x block
  } 
  if (TMath::Abs(index[1])<2) index[3] = 0;
  else {
    index[3] = TMath::Abs(index[1])-1;
    if (index[1]<0) index[3]*=-1;   //first y block
  } 
  if (TMath::Abs(index[0])<TMath::Abs(index[1])){
    index[4]=index[0];
    index[5]=index[3];
  }
  else
    if (index[0]==index[1]) {
      index[4]=0;
      index[5]=0;
    }
    else{
      index[4]=index[2];
      index[5]=index[1]; 
    }
  return;  
}

//***********************************************************************
//***********************************************************************

TClonesArray * AliTPCClusterFinder::FindPeaks1(TClonesArray *arr)
{
  //find peaks and write it in form of AliTPCcluster to array
  if (arr==0){
    fClustersArray=new TClonesArray("AliDigitCluster",300);
    fIndex=1;
  }
  else {
    fClustersArray = arr;
    fIndex = fClustersArray->GetEntriesFast();
  }
 
  AliDigitCluster c;           
  ResetStatus();  
   for (Int_t i=0; i<fDimX; i++) 
     for (Int_t j=0;j<fDimY; j++) 
       {	
	 fStackIndex=0;          
	 fBDistType = kFALSE;
	 AliCell * cell = GetCell(i,j);
	 if (!(cell->IsChecked()))  Adjacent(i,j);
         //if there exists more then 2 digits cluster 
	 if (fStackIndex >2 ){	   
	   if (fBFit==kFALSE) GetClusterStatistic(c);
	     else GetClusterFit(c);
	   //write some important chracteristic area of cluster
	   //	   
	   Transform(&c);
	   //write cluster information to array
	   TClonesArray &lclusters = *fClustersArray;
	   new (lclusters[fIndex++])  AliDigitCluster(c);
	   //             cout<<"fx="<<c.fX<<"   fy"<<c.fY<<"\n";	   
	 } 
       }
   return fClustersArray;
}


TClonesArray * AliTPCClusterFinder::FindPeaks2(TClonesArray *arr)
{
  //find peaks and write it in form of AliTPCcluster to array
  if (arr==0){
    fClustersArray=new TClonesArray("AliDigitCluster",300);
    fIndex=1;
  }
  else {
    fClustersArray = arr;
    fIndex = fClustersArray->GetEntriesFast();
  }

  AliDigitCluster c;           
  ResetStatus();  
  
   for (Int_t i=0; i<fDimX; i++) 
     for (Int_t j=0;j<fDimY; j++) 
       {
	 fStackIndex=0;  
	 if (IsMaximum(i,j) == kTRUE){
	   SetSigma2(i,j,fCurrentSigmaX2,fCurrentSigmaY2);
	   fBDistType = kTRUE;
	   Adjacent(i,j);
	   //if there exists more then 2 digits cluster 
	   if (fStackIndex >2 ){
	     if (fBFit==kFALSE) GetClusterStatistic(c);
	     else GetClusterFit(c);
	     //write some important chracteristic area of cluster
	     //	   
	     Transform(&c);
	     //write cluster information to array
	     TClonesArray &lclusters = *fClustersArray;
	     new(lclusters[fIndex++]) AliDigitCluster(c);
	     //             cout<<"fx="<<c.fX<<"   fy"<<c.fY<<"\n";	   
	   } 
	 }
       }
   return fClustersArray;
}


TClonesArray * AliTPCClusterFinder::FindPeaks3(TClonesArray *arr)
{
  //find peaks and write it in form of AliTPCcluster to array
  if (arr==0){
    fClustersArray=new TClonesArray("AliDigitCluster",300);
    fIndex=1;
  }
  else {
    fClustersArray = arr;
    fIndex = fClustersArray->GetEntriesFast();
  }
  
  AliDigitCluster c;    
  ResetStatus();  
  
  Int_t dmax=5;
  Int_t naccepted =1;
   for (Int_t i=0; i<fDimX; i++) 
     for (Int_t j=0;j<fDimY; j++) 
       {
	 fStackIndex=0;  
	 if (IsMaximum(i,j) == kTRUE){
	   SetSigma2(i,j,fCurrentSigmaX2,fCurrentSigmaY2);
           AddToStack(i,j,GetSignal(i,j));
	   
	   //loop over different distance 
	   naccepted =1;
	   for ( Int_t dd =1;((dd<=dmax) && (naccepted>0));dd++){
             naccepted=0; 
	     for (Int_t di = -dd;di<=dd;di++){
	       Int_t ddj = dd-abs(di);
	       Int_t sigstart = (ddj>0) ?  -1 : 0;
	       for (Int_t sig = sigstart;sig<=1;sig+=2){
		 Int_t dj= sig*ddj; 
		 AliCell *cell= GetCell(i+di,j+dj);
		 Int_t signal = GetSignal(i+di,j+dj);
		 if (cell==0) continue;
		 Int_t index[6];
		 index[0]=di;
		 index[1]=dj;
		 if (dd>2) {
		   SetBlockIndex(index);  //adjust index to control	       
		   if ( IsBorder(fIndex,i+index[2],j+index[3]) || 
			IsBorder(fIndex,i+index[4],j+index[5])) {
		     cell->SetBorder(fIndex);   
		     continue;
		   }
		 }
		 if ( signal<=fThreshold ){
		   //if under threshold
		   cell->SetThBorder(fIndex);
		   if (fBFit==kTRUE)  AddToStack(i+di,j+dj,signal);
		   continue;
		 }
		 naccepted++;	       
		 if (CheckIfDirBorder(fCurrentMaxX,fCurrentMaxY,i+di,j+dj) == kTRUE) {
		   if (fBFit==kFALSE) AddToStack(i+di,j+dj,signal/2);
		   continue; 
		 }
		 AddToStack(i+di,j+dj,signal);

	       } //loop over sig dj 
	     } //loop over di
	     
	   }//loop over dd
	 } //if there is maximum
	 //if there exists more then 2 digits cluster 
	 if (fStackIndex >2 ){
	   if (fBFit==kFALSE) GetClusterStatistic(c);
	   else GetClusterFit(c);
	   //write some important chracteristic area of cluster
	   //	   
	   Transform(&c);
	   //write cluster information to array
	   TClonesArray &lclusters = *fClustersArray;
	   new(lclusters[fIndex++]) AliDigitCluster(c);
	   //             cout<<"fx="<<c.fX<<"   fy"<<c.fY<<"\n";	   
	 }
       } //lopp over all digits

   return fClustersArray;
}






void AliTPCClusterFinder::Adjacent(Int_t i,Int_t j)
{
  //
  //recursive agorithm program
  //
  if (fBDistType==kTRUE) {
    Float_t delta = (i-fCurrentMaxX)*(i-fCurrentMaxX)/fCurrentSigmaX2;
    delta+=(j-fCurrentMaxY)*(j-fCurrentMaxY)/fCurrentSigmaY2;
    if (delta > fMulSigma2) {
       SetDirBorder(fIndex,i,j);
      return;
    }
  }
  AliCell *cell = GetCell(i,j);
  Int_t signal = GetSignal(i,j);
  Int_t q=signal;  
  cell->SetChecked(fIndex);  
  if ( (q>fThreshold) || (fBFit==kTRUE))   AddToStack(i,j,q);
  if ( q >fThreshold )
    {
      
      AliCell * newcel;      
      newcel = GetCell(i-1,j);
      if (newcel !=0) if (!newcel->IsChecked(fIndex) ) Adjacent(i-1,j);
      newcel = GetCell(i,j-1);
      if (newcel !=0) if (!newcel->IsChecked(fIndex) ) Adjacent(i,j-1);
      newcel = GetCell(i+1,j);
      if (newcel !=0) if (!newcel->IsChecked(fIndex) ) Adjacent(i+1,j);
      newcel = GetCell(i,j+1);
      if (newcel !=0) if (!newcel->IsChecked(fIndex) ) Adjacent(i,j+1);
    }      
  else cell->SetThBorder(fIndex);
}



AliH2F *  AliTPCClusterFinder::DrawHisto( const char *option=0, 
			      Float_t x1, Float_t x2, Float_t y1, Float_t y2)
{
  //
  //draw digits in given array
  //  
  //make digits histo 
  char ch[30];
  sprintf(ch,"Cluster finder digits ");
  if ( (fDimX<1)|| (fDimY<1)) {
    return 0;
  }
  AliH2F * his  = new AliH2F(ch,ch,fDimX,fX1,fX2,fDimY,fY1,fY2);
  //set histogram  values
  for (Int_t i = 0; i<fDimX;i++)    
    for (Int_t j = 0; j<fDimY;j++){
      Float_t x = ItoX(i);
      Float_t y= JtoY(j);
      his->Fill(x,y,GetSignal(i,j));
    }
  if (x1>=0) {
      AliH2F *h2fsub = his->GetSubrange2d(x1,x2,y1,y2);
      delete his;
      his=h2fsub;
  }
  if (his==0) return 0;
  if (option!=0) his->Draw(option);
  else his->Draw();
  return his;  
}


void AliTPCClusterFinder::DrawCluster(
				  Int_t color, Int_t size, Int_t style)
{

  if (fClustersArray==0) return;  
  //draw marker for each of cluster
  Int_t ncl=fClustersArray->GetEntriesFast();
  for (Int_t i=0;i<ncl;i++){
    AliCluster *cl = (AliCluster*)fClustersArray->UncheckedAt(i);
    TMarker * marker = new TMarker;
    marker->SetX(cl->fX);
    marker->SetY(cl->fY);
    marker->SetMarkerSize(size);
    marker->SetMarkerStyle(style);
    marker->SetMarkerColor(color);
    marker->Draw();    
  }
}



AliH2F *  AliTPCClusterFinder::DrawBorders( const char *option,  AliH2F *h, Int_t type ,
			      Float_t x1, Float_t x2, Float_t y1, Float_t y2)
{
  //
  //draw digits in given array
  //  
  //make digits histo 
  char ch[30];
  sprintf(ch,"Cluster finder digits borders");
  if ( (fDimX<1)|| (fDimY<1)) {
    return 0;
  }
  AliH2F * his;
  if (h!=0) his =h;
  else his  = new AliH2F(ch,ch,fDimX,fX1,fX2,fDimY,fY1,fY2);
  //set histogram  values
  for (Int_t i = 0; i<fDimX;i++)    
    for (Int_t j = 0; j<fDimY;j++){      
      Float_t x = ItoX(i);
      Float_t y= JtoY(j);
      if (((type==1)||(type==0))&&IsMaximum(0,i,j)) his->Fill(x,y,16);   
      if (((type==3)||(type==0))&&(IsDirBorder(0,i,j))) his->Fill(x,y,8);
      if (((type==4)||(type==0))&&(IsThBorder(0,i,j))) his->Fill(x,y,4);       
      if (((type==2)||(type==0))&&IsBorder(0,i,j)) his->Fill(x,y,1);

    }
	 
  if (x1>=0) {
      AliH2F *h2fsub = his->GetSubrange2d(x1,x2,y1,y2);
      delete his;
      his=h2fsub;
  }
  if (his==0) return 0;
  if (option!=0) his->Draw(option);
  else his->Draw();
  return his;  
}
