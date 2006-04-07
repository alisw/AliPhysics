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

//-------------------------------------------------------
//          Implementation of the TPC clusterer
//
//   Origin: Marian Ivanov 
//-------------------------------------------------------

#include "AliTPCReconstructor.h"
#include "AliTPCclustererMI.h"
#include "AliTPCclusterMI.h"
#include <TObjArray.h>
#include <TFile.h>
#include "AliTPCClustersArray.h"
#include "AliTPCClustersRow.h"
#include "AliTPCRawStream.h"
#include "AliDigits.h"
#include "AliSimDigits.h"
#include "AliTPCParam.h"
#include "AliRawReader.h"
#include "AliTPCRawStream.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "Riostream.h"
#include <TTree.h>

#include "AliTPCcalibDB.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"


ClassImp(AliTPCclustererMI)



AliTPCclustererMI::AliTPCclustererMI(const AliTPCParam* par)
{
  fInput =0;
  fOutput=0;
  fParam = par;
}
void AliTPCclustererMI::SetInput(TTree * tree)
{
  //
  // set input tree with digits
  //
  fInput = tree;  
  if  (!fInput->GetBranch("Segment")){
    cerr<<"AliTPC::Digits2Clusters(): no porper input tree !\n";
    fInput=0;
    return;
  }
}

void AliTPCclustererMI::SetOutput(TTree * tree) 
{
  //
  //
  fOutput= tree;  
  AliTPCClustersRow clrow;
  AliTPCClustersRow *pclrow=&clrow;  
  clrow.SetClass("AliTPCclusterMI");
  clrow.SetArray(1); // to make Clones array
  fOutput->Branch("Segment","AliTPCClustersRow",&pclrow,32000,200);    
}


Float_t  AliTPCclustererMI::GetSigmaY2(Int_t iz){
  // sigma y2 = in digits  - we don't know the angle
  Float_t z = iz*fParam->GetZWidth()+fParam->GetNTBinsL1()*fParam->GetZWidth();
  Float_t sd2 = (z*fParam->GetDiffL()*fParam->GetDiffL())/
    (fPadWidth*fPadWidth);
  Float_t sres = 0.25;
  Float_t res = sd2+sres;
  return res;
}


Float_t  AliTPCclustererMI::GetSigmaZ2(Int_t iz){
  //sigma z2 = in digits - angle estimated supposing vertex constraint
  Float_t z = iz*fZWidth+fParam->GetNTBinsL1()*fParam->GetZWidth();
  Float_t sd2 = (z*fParam->GetDiffL()*fParam->GetDiffL())/(fZWidth*fZWidth);
  Float_t angular = fPadLength*(fParam->GetZLength()-z)/(fRx*fZWidth);
  angular*=angular;
  angular/=12.;
  Float_t sres = fParam->GetZSigma()/fZWidth;
  sres *=sres;
  Float_t res = angular +sd2+sres;
  return res;
}

void AliTPCclustererMI::MakeCluster(Int_t k,Int_t max,Float_t *bins, UInt_t /*m*/,
AliTPCclusterMI &c) 
{
  Int_t i0=k/max;  //central pad
  Int_t j0=k%max;  //central time bin

  // set pointers to data
  //Int_t dummy[5] ={0,0,0,0,0};
  Float_t * matrix[5]; //5x5 matrix with digits  - indexing i = 0 ..4  j = -2..2
  Float_t * resmatrix[5];
  for (Int_t di=-2;di<=2;di++){
    matrix[di+2]  =  &bins[k+di*max];
    resmatrix[di+2]  =  &fResBins[k+di*max];
  }
  //build matrix with virtual charge
  Float_t sigmay2= GetSigmaY2(j0);
  Float_t sigmaz2= GetSigmaZ2(j0);

  Float_t vmatrix[5][5];
  vmatrix[2][2] = matrix[2][0];
  c.SetType(0);
  c.SetMax((UShort_t)(vmatrix[2][2])); // write maximal amplitude
  for (Int_t di =-1;di <=1;di++)
    for (Int_t dj =-1;dj <=1;dj++){
      Float_t amp = matrix[di+2][dj];
      if ( (amp<2) && (fLoop<2)){
	// if under threshold  - calculate virtual charge
	Float_t ratio = TMath::Exp(-1.2*TMath::Abs(di)/sigmay2)*TMath::Exp(-1.2*TMath::Abs(dj)/sigmaz2);
	amp = ((matrix[2][0]-2)*(matrix[2][0]-2)/(matrix[-di+2][-dj]+2))*ratio;
	if (amp>2) amp = 2;
	vmatrix[2+di][2+dj]=amp;
	vmatrix[2+2*di][2+2*dj]=0;
	if ( (di*dj)!=0){       
	  //DIAGONAL ELEMENTS
	  vmatrix[2+2*di][2+dj] =0;
	  vmatrix[2+di][2+2*dj] =0;
	}
	continue;
      }
      if (amp<4){
	//if small  amplitude - below  2 x threshold  - don't consider other one	
	vmatrix[2+di][2+dj]=amp;
	vmatrix[2+2*di][2+2*dj]=0;  // don't take to the account next bin
	if ( (di*dj)!=0){       
	  //DIAGONAL ELEMENTS
	  vmatrix[2+2*di][2+dj] =0;
	  vmatrix[2+di][2+2*dj] =0;
	}
	continue;
      }
      //if bigger then take everything
      vmatrix[2+di][2+dj]=amp;
      vmatrix[2+2*di][2+2*dj]= matrix[2*di+2][2*dj] ;      
      if ( (di*dj)!=0){       
	  //DIAGONAL ELEMENTS
	  vmatrix[2+2*di][2+dj] = matrix[2*di+2][dj];
	  vmatrix[2+di][2+2*dj] = matrix[2+di][dj*2];
	}      
    }


  
  Float_t sumw=0;
  Float_t sumiw=0;
  Float_t sumi2w=0;
  Float_t sumjw=0;
  Float_t sumj2w=0;
  //
  for (Int_t i=-2;i<=2;i++)
    for (Int_t j=-2;j<=2;j++){
      Float_t amp = vmatrix[i+2][j+2];

      sumw    += amp;
      sumiw   += i*amp;
      sumi2w  += i*i*amp;
      sumjw   += j*amp;
      sumj2w  += j*j*amp;
    }    
  //   
  Float_t meani = sumiw/sumw;
  Float_t mi2   = sumi2w/sumw-meani*meani;
  Float_t meanj = sumjw/sumw;
  Float_t mj2   = sumj2w/sumw-meanj*meanj;
  //
  Float_t ry = mi2/sigmay2;
  Float_t rz = mj2/sigmaz2;
  
  //
  if ( ( (ry<0.6) || (rz<0.6) ) && fLoop==2) return;
  if ( (ry <1.2) && (rz<1.2) ) {
    //if cluster looks like expected 
    //+1.2 deviation from expected sigma accepted
    //    c.fMax = FitMax(vmatrix,meani,meanj,TMath::Sqrt(sigmay2),TMath::Sqrt(sigmaz2));

    meani +=i0;
    meanj +=j0;
    //set cluster parameters
    c.SetQ(sumw);
    c.SetY(meani*fPadWidth); 
    c.SetZ(meanj*fZWidth); 
    c.SetSigmaY2(mi2);
    c.SetSigmaZ2(mj2);
    AddCluster(c);
    //remove cluster data from data
    for (Int_t di=-2;di<=2;di++)
      for (Int_t dj=-2;dj<=2;dj++){
	resmatrix[di+2][dj] -= vmatrix[di+2][dj+2];
	if (resmatrix[di+2][dj]<0) resmatrix[di+2][dj]=0;
      }
    resmatrix[2][0] =0;
    return;     
  }
  //
  //unfolding when neccessary  
  //
  
  Float_t * matrix2[7]; //7x7 matrix with digits  - indexing i = 0 ..6  j = -3..3
  Float_t dummy[7]={0,0,0,0,0,0};
  for (Int_t di=-3;di<=3;di++){
    matrix2[di+3] =  &bins[k+di*max];
    if ((k+di*max)<3)  matrix2[di+3] = &dummy[3];
    if ((k+di*max)>fMaxBin-3)  matrix2[di+3] = &dummy[3];
  }
  Float_t vmatrix2[5][5];
  Float_t sumu;
  Float_t overlap;
  UnfoldCluster(matrix2,vmatrix2,meani,meanj,sumu,overlap);
  //
  //  c.fMax = FitMax(vmatrix2,meani,meanj,TMath::Sqrt(sigmay2),TMath::Sqrt(sigmaz2));
  meani +=i0;
  meanj +=j0;
  //set cluster parameters
  c.SetQ(sumu);
  c.SetY(meani*fPadWidth); 
  c.SetZ(meanj*fZWidth); 
  c.SetSigmaY2(mi2);
  c.SetSigmaZ2(mj2);
  c.SetType(Char_t(overlap)+1);
  AddCluster(c);

  //unfolding 2
  meani-=i0;
  meanj-=j0;
  if (gDebug>4)
    printf("%f\t%f\n", vmatrix2[2][2], vmatrix[2][2]);
}



void AliTPCclustererMI::UnfoldCluster(Float_t * matrix2[7], Float_t recmatrix[5][5], Float_t & meani, Float_t & meanj, 
				      Float_t & sumu, Float_t & overlap )
{
  //
  //unfold cluster from input matrix
  //data corresponding to cluster writen in recmatrix
  //output meani and meanj

  //take separatelly y and z

  Float_t sum3i[7] = {0,0,0,0,0,0,0};
  Float_t sum3j[7] = {0,0,0,0,0,0,0};

  for (Int_t k =0;k<7;k++)
    for (Int_t l = -1; l<=1;l++){
      sum3i[k]+=matrix2[k][l];
      sum3j[k]+=matrix2[l+3][k-3];
    }
  Float_t mratio[3][3]={{1,1,1},{1,1,1},{1,1,1}};
  //
  //unfold  y 
  Float_t sum3wi    = 0;  //charge minus overlap
  Float_t sum3wio   = 0;  //full charge
  Float_t sum3iw    = 0;  //sum for mean value
  for (Int_t dk=-1;dk<=1;dk++){
    sum3wio+=sum3i[dk+3];
    if (dk==0){
      sum3wi+=sum3i[dk+3];     
    }
    else{
      Float_t ratio =1;
      if (  ( ((sum3i[dk+3]+3)/(sum3i[3]-3))+1 < (sum3i[2*dk+3]-3)/(sum3i[dk+3]+3))||
	   sum3i[dk+3]<=sum3i[2*dk+3] && sum3i[dk+3]>2 ){
	Float_t xm2 = sum3i[-dk+3];
	Float_t xm1 = sum3i[+3];
	Float_t x1  = sum3i[2*dk+3];
	Float_t x2  = sum3i[3*dk+3]; 	
	Float_t w11   = TMath::Max((Float_t)(4.*xm1-xm2),(Float_t)0.000001);	  
	Float_t w12   = TMath::Max((Float_t)(4 *x1 -x2),(Float_t)0.);
	ratio = w11/(w11+w12);	 
	for (Int_t dl=-1;dl<=1;dl++)
	  mratio[dk+1][dl+1] *= ratio;
      }
      Float_t amp = sum3i[dk+3]*ratio;
      sum3wi+=amp;
      sum3iw+= dk*amp;      
    }
  }
  meani = sum3iw/sum3wi;
  Float_t overlapi = (sum3wio-sum3wi)/sum3wio;



  //unfold  z 
  Float_t sum3wj    = 0;  //charge minus overlap
  Float_t sum3wjo   = 0;  //full charge
  Float_t sum3jw    = 0;  //sum for mean value
  for (Int_t dk=-1;dk<=1;dk++){
    sum3wjo+=sum3j[dk+3];
    if (dk==0){
      sum3wj+=sum3j[dk+3];     
    }
    else{
      Float_t ratio =1;
      if ( ( ((sum3j[dk+3]+3)/(sum3j[3]-3))+1 < (sum3j[2*dk+3]-3)/(sum3j[dk+3]+3)) ||
	   (sum3j[dk+3]<=sum3j[2*dk+3] && sum3j[dk+3]>2)){
	Float_t xm2 = sum3j[-dk+3];
	Float_t xm1 = sum3j[+3];
	Float_t x1  = sum3j[2*dk+3];
	Float_t x2  = sum3j[3*dk+3]; 	
	Float_t w11   = TMath::Max((Float_t)(4.*xm1-xm2),(Float_t)0.000001);	  
	Float_t w12   = TMath::Max((Float_t)(4 *x1 -x2),(Float_t)0.);
	ratio = w11/(w11+w12);	 
	for (Int_t dl=-1;dl<=1;dl++)
	  mratio[dl+1][dk+1] *= ratio;
      }
      Float_t amp = sum3j[dk+3]*ratio;
      sum3wj+=amp;
      sum3jw+= dk*amp;      
    }
  }
  meanj = sum3jw/sum3wj;
  Float_t overlapj = (sum3wjo-sum3wj)/sum3wjo;  
  overlap = Int_t(100*TMath::Max(overlapi,overlapj)+3);  
  sumu = (sum3wj+sum3wi)/2.;
  
  if (overlap ==3) {
    //if not overlap detected remove everything
    for (Int_t di =-2; di<=2;di++)
      for (Int_t dj =-2; dj<=2;dj++){
	recmatrix[di+2][dj+2] = matrix2[3+di][dj];
      }
  }
  else{
    for (Int_t di =-1; di<=1;di++)
      for (Int_t dj =-1; dj<=1;dj++){
	Float_t ratio =1;
	if (mratio[di+1][dj+1]==1){
	  recmatrix[di+2][dj+2]     = matrix2[3+di][dj];
	  if (TMath::Abs(di)+TMath::Abs(dj)>1){
	    recmatrix[2*di+2][dj+2] = matrix2[3+2*di][dj];
	    recmatrix[di+2][2*dj+2] = matrix2[3+di][2*dj];
	  }	  
	  recmatrix[2*di+2][2*dj+2] = matrix2[3+2*di][2*dj];
	}
	else
	  {
	    //if we have overlap in direction
	    recmatrix[di+2][dj+2] = mratio[di+1][dj+1]* matrix2[3+di][dj];    
	    if (TMath::Abs(di)+TMath::Abs(dj)>1){
	      ratio =  TMath::Min((Float_t)(recmatrix[di+2][dj+2]/(matrix2[3+0*di][1*dj]+1)),(Float_t)1.);
	      recmatrix[2*di+2][dj+2] = ratio*recmatrix[di+2][dj+2];
	      //
	      ratio =  TMath::Min((Float_t)(recmatrix[di+2][dj+2]/(matrix2[3+1*di][0*dj]+1)),(Float_t)1.);
	      recmatrix[di+2][2*dj+2] = ratio*recmatrix[di+2][dj+2];
	    }
	    else{
	      ratio =  recmatrix[di+2][dj+2]/matrix2[3][0];
	      recmatrix[2*di+2][2*dj+2] = ratio*recmatrix[di+2][dj+2];
	    }
	  }
      }
  }
  if (gDebug>4) 
    printf("%f\n", recmatrix[2][2]);
  
}

Float_t AliTPCclustererMI::FitMax(Float_t vmatrix[5][5], Float_t y, Float_t z, Float_t sigmay, Float_t sigmaz)
{
  //
  // estimate max
  Float_t sumteor= 0;
  Float_t sumamp = 0;

  for (Int_t di = -1;di<=1;di++)
    for (Int_t dj = -1;dj<=1;dj++){
      if (vmatrix[2+di][2+dj]>2){
	Float_t teor = TMath::Gaus(di,y,sigmay*1.2)*TMath::Gaus(dj,z,sigmaz*1.2);
	sumteor += teor*vmatrix[2+di][2+dj];
	sumamp  += vmatrix[2+di][2+dj]*vmatrix[2+di][2+dj];
      }
    }    
  Float_t max = sumamp/sumteor;
  return max;
}

void AliTPCclustererMI::AddCluster(AliTPCclusterMI &c){
  //
  // transform cluster to the global coordinata
  // add the cluster to the array
  //
  Float_t meani = c.GetY()/fPadWidth;
  Float_t meanj = c.GetZ()/fZWidth;

  Int_t ki = TMath::Nint(meani-3);
  if (ki<0) ki=0;
  if (ki>=fMaxPad) ki = fMaxPad-1;
  Int_t kj = TMath::Nint(meanj-3);
  if (kj<0) kj=0;
  if (kj>=fMaxTime-3) kj=fMaxTime-4;
  // ki and kj shifted to "real" coordinata
  if (fRowDig) {
    c.SetLabel(fRowDig->GetTrackIDFast(kj,ki,0)-2,0);
    c.SetLabel(fRowDig->GetTrackIDFast(kj,ki,1)-2,1);
    c.SetLabel(fRowDig->GetTrackIDFast(kj,ki,2)-2,2);
  }
  
  
  Float_t s2 = c.GetSigmaY2();
  Float_t w=fParam->GetPadPitchWidth(fSector);
  
  c.SetSigmaY2(s2*w*w);
  s2 = c.GetSigmaZ2(); 
  w=fZWidth;
  c.SetSigmaZ2(s2*w*w);
  c.SetY((meani - 2.5 - 0.5*fMaxPad)*fParam->GetPadPitchWidth(fSector));
  c.SetZ(fZWidth*(meanj-3)); 
  c.SetZ(c.GetZ() - 3.*fParam->GetZSigma() + fParam->GetNTBinsL1()*fParam->GetZWidth()); // PASA delay + L1 delay
  c.SetZ(fSign*(fParam->GetZLength() - c.GetZ()));
  c.SetX(fRx);
  c.SetDetector(fSector);
  c.SetRow(fRow);

  if (ki<=1 || ki>=fMaxPad-1 || kj==1 || kj==fMaxTime-2) {
    //c.SetSigmaY2(c.GetSigmaY2()*25.);
    //c.SetSigmaZ2(c.GetSigmaZ2()*4.);
    c.SetType(-(c.GetType()+3));  //edge clusters
  }
  if (fLoop==2) c.SetType(100);

  TClonesArray * arr = fRowCl->GetArray();
  // AliTPCclusterMI * cl = 
  new ((*arr)[fNcluster]) AliTPCclusterMI(c);

  fNcluster++;
}


//_____________________________________________________________________________
void AliTPCclustererMI::Digits2Clusters()
{
  //-----------------------------------------------------------------
  // This is a simple cluster finder.
  //-----------------------------------------------------------------

  if (!fInput) { 
    Error("Digits2Clusters", "input tree not initialised");
    return;
  }
 
  if (!fOutput) {
    Error("Digits2Clusters", "output tree not initialised");
    return;
  }

  AliTPCCalPad * gainTPC = AliTPCcalibDB::Instance()->GetPadGainFactor();

  AliSimDigits digarr, *dummy=&digarr;
  fRowDig = dummy;
  fInput->GetBranch("Segment")->SetAddress(&dummy);
  Stat_t nentries = fInput->GetEntries();
  
  fMaxTime=fParam->GetMaxTBin()+6; // add 3 virtual time bins before and 3   after
    
  Int_t nclusters  = 0;

  for (Int_t n=0; n<nentries; n++) {
    fInput->GetEvent(n);
    if (!fParam->AdjustSectorRow(digarr.GetID(),fSector,fRow)) {
      cerr<<"AliTPC warning: invalid segment ID ! "<<digarr.GetID()<<endl;
      continue;
    }
    Int_t row = fRow;
    AliTPCCalROC * gainROC = gainTPC->GetCalROC(fSector);  // pad gains per given sector
    
    //
    AliTPCClustersRow *clrow= new AliTPCClustersRow();
    fRowCl = clrow;
    clrow->SetClass("AliTPCclusterMI");
    clrow->SetArray(1);

    clrow->SetID(digarr.GetID());
    fOutput->GetBranch("Segment")->SetAddress(&clrow);
    fRx=fParam->GetPadRowRadii(fSector,row);
    
    
    const Int_t kNIS=fParam->GetNInnerSector(), kNOS=fParam->GetNOuterSector();
    fZWidth = fParam->GetZWidth();
    if (fSector < kNIS) {
      fMaxPad = fParam->GetNPadsLow(row);
      fSign = (fSector < kNIS/2) ? 1 : -1;
      fPadLength = fParam->GetPadPitchLength(fSector,row);
      fPadWidth = fParam->GetPadPitchWidth();
    } else {
      fMaxPad = fParam->GetNPadsUp(row);
      fSign = ((fSector-kNIS) < kNOS/2) ? 1 : -1;
      fPadLength = fParam->GetPadPitchLength(fSector,row);
      fPadWidth  = fParam->GetPadPitchWidth();
    }
    
    
    fMaxBin=fMaxTime*(fMaxPad+6);  // add 3 virtual pads  before and 3 after
    fBins    =new Float_t[fMaxBin];
    fResBins =new Float_t[fMaxBin];  //fBins with residuals after 1 finder loop 
    memset(fBins,0,sizeof(Float_t)*fMaxBin);
    memset(fResBins,0,sizeof(Float_t)*fMaxBin);
    
    if (digarr.First()) //MI change
      do {
	Float_t dig=digarr.CurrentDigit();
	if (dig<=fParam->GetZeroSup()) continue;
	Int_t j=digarr.CurrentRow()+3, i=digarr.CurrentColumn()+3;
	Float_t gain = gainROC->GetValue(row,digarr.CurrentRow()/fParam->GetMaxTBin());
	fBins[i*fMaxTime+j]=dig/gain;
      } while (digarr.Next());
    digarr.ExpandTrackBuffer();

    FindClusters();

    fOutput->Fill();
    delete clrow;    
    nclusters+=fNcluster;    
    delete[] fBins;      
    delete[] fResBins;  
  }  

  Info("Digits2Clusters", "Number of found clusters : %d", nclusters);
}

void AliTPCclustererMI::Digits2Clusters(AliRawReader* rawReader)
{
//-----------------------------------------------------------------
// This is a cluster finder for raw data.
//-----------------------------------------------------------------

  if (!fOutput) {
    Error("Digits2Clusters", "output tree not initialised");
    return;
  }

  rawReader->Reset();
  AliTPCRawStream input(rawReader);

  fRowDig = NULL;

  Int_t nclusters  = 0;
  
  fMaxTime = fParam->GetMaxTBin() + 6; // add 3 virtual time bins before and 3 after
  const Int_t kNIS = fParam->GetNInnerSector();
  const Int_t kNOS = fParam->GetNOuterSector();
  const Int_t kNS = kNIS + kNOS;
  fZWidth = fParam->GetZWidth();
  Int_t zeroSup = fParam->GetZeroSup();

  fBins = NULL;
  Float_t** splitRows = new Float_t* [kNS*2];
  Float_t** splitRowsRes = new Float_t* [kNS*2];
  for (Int_t iSector = 0; iSector < kNS*2; iSector++)
    splitRows[iSector] = NULL;
  Int_t iSplitRow = -1;

  Bool_t next = kTRUE;
  while (next) {
    next = input.Next();

    // when the sector or row number has changed ...
    if (input.IsNewRow() || !next) {

      // ... find clusters in the previous pad row, and ...
      if (fBins) {
	if ((iSplitRow < 0) || splitRows[fSector + kNS*iSplitRow]) {
	  fRowCl = new AliTPCClustersRow;
	  fRowCl->SetClass("AliTPCclusterMI");
	  fRowCl->SetArray(1);
	  fRowCl->SetID(fParam->GetIndex(fSector, input.GetPrevRow()));
	  fOutput->GetBranch("Segment")->SetAddress(&fRowCl);

	  FindClusters();

	  fOutput->Fill();
	  delete fRowCl;    
	  nclusters += fNcluster;    
	  delete[] fBins;
	  delete[] fResBins;
	  if (iSplitRow >= 0) splitRows[fSector + kNS*iSplitRow] = NULL;

	} else if (iSplitRow >= 0) {
	  splitRows[fSector + kNS*iSplitRow] = fBins;
	  splitRowsRes[fSector + kNS*iSplitRow] = fResBins;
	}
      }

      if (!next) break;

      // ... prepare for the next pad row
      fSector = input.GetSector();
      fRow    = input.GetRow();
      Int_t iRow = input.GetRow();
      fRx = fParam->GetPadRowRadii(fSector, iRow);
    
      iSplitRow = -1;
      if (fSector < kNIS) {
	fMaxPad = fParam->GetNPadsLow(iRow);
	fSign = (fSector < kNIS/2) ? 1 : -1;
	if (iRow == 30) iSplitRow = 0;
      } else {
	fMaxPad = fParam->GetNPadsUp(iRow);
	fSign = ((fSector-kNIS) < kNOS/2) ? 1 : -1;
	if (iRow == 27) iSplitRow = 0;
	else if (iRow == 76) iSplitRow = 1;
      }
      fPadLength = fParam->GetPadPitchLength(fSector, iRow);
      fPadWidth  = fParam->GetPadPitchWidth();
    
      fMaxBin = fMaxTime*(fMaxPad+6);  // add 3 virtual pads  before and 3 after
      if ((iSplitRow < 0) || !splitRows[fSector + kNS*iSplitRow]) {
	fBins    = new Float_t[fMaxBin];
	fResBins = new Float_t[fMaxBin];  //fBins with residuals after 1 finder loop 
	memset(fBins, 0, sizeof(Float_t)*fMaxBin);
	memset(fResBins, 0, sizeof(Float_t)*fMaxBin);
      } else {
	fBins    = splitRows[fSector + kNS*iSplitRow];
	fResBins = splitRowsRes[fSector + kNS*iSplitRow];
      }
    }

    // fill fBins with digits data
    if (input.GetSignal() <= zeroSup) continue;
    Int_t i = input.GetPad() + 3;
    Int_t j = input.GetTime() + 3;
    fBins[i*fMaxTime+j] = input.GetSignal();
  }

  // find clusters in split rows that were skipped until now.
  // this can happen if the rows were not splitted 
  for (fSector = 0; fSector < kNS; fSector++)
    for (Int_t iSplit = 0; iSplit < 2; iSplit++)
      if (splitRows[fSector + kNS*iSplit]) {

	Int_t iRow = -1;
	if (fSector < kNIS) {
	  iRow = 30;
	  fMaxPad = fParam->GetNPadsLow(iRow);
	  fSign = (fSector < kNIS/2) ? 1 : -1;
	} else {
	  if (iSplit == 0) iRow = 27; else iRow = 76;
	  fMaxPad = fParam->GetNPadsUp(iRow);
	  fSign = ((fSector-kNIS) < kNOS/2) ? 1 : -1;
	}
	fRx = fParam->GetPadRowRadii(fSector, iRow);
	fPadLength = fParam->GetPadPitchLength(fSector, iRow);
	fPadWidth  = fParam->GetPadPitchWidth();

    	fMaxBin = fMaxTime*(fMaxPad+6);  // add 3 virtual pads  before and 3 after
	fBins    = splitRows[fSector + kNS*iSplit];
	fResBins = splitRowsRes[fSector + kNS*iSplit];

	fRowCl = new AliTPCClustersRow;
	fRowCl->SetClass("AliTPCclusterMI");
	fRowCl->SetArray(1);
	fRowCl->SetID(fParam->GetIndex(fSector, iRow));
	fOutput->GetBranch("Segment")->SetAddress(&fRowCl);

	FindClusters();

	fOutput->Fill();
	delete fRowCl;    
	nclusters += fNcluster;    
	delete[] fBins;
	delete[] fResBins;
      }

  delete[] splitRows;
  delete[] splitRowsRes;
  Info("Digits2Clusters", "Number of found clusters : %d\n", nclusters);
}

void AliTPCclustererMI::FindClusters()
{
  //add virtual charge at the edge   
  for (Int_t i=0; i<fMaxTime; i++){
    Float_t amp1 = fBins[i+3*fMaxTime]; 
    Float_t amp0 =0;
    if (amp1>0){
      Float_t amp2 = fBins[i+4*fMaxTime];
      if (amp2==0) amp2=0.5;
      Float_t sigma2 = GetSigmaY2(i);		
      amp0 = (amp1*amp1/amp2)*TMath::Exp(-1./sigma2);	
      if (gDebug>4) printf("\n%f\n",amp0);
    }  
    fBins[i+2*fMaxTime] = amp0;    
    amp0 = 0;
    amp1 = fBins[(fMaxPad+2)*fMaxTime+i];
    if (amp1>0){
      Float_t amp2 = fBins[i+(fMaxPad+1)*fMaxTime];
      if (amp2==0) amp2=0.5;
      Float_t sigma2 = GetSigmaY2(i);		
      amp0 = (amp1*amp1/amp2)*TMath::Exp(-1./sigma2);		
      if (gDebug>4) printf("\n%f\n",amp0);
    }        
    fBins[(fMaxPad+3)*fMaxTime+i] = amp0;           
  }

//  memcpy(fResBins,fBins, fMaxBin*2);
  memcpy(fResBins,fBins, fMaxBin);
  //
  fNcluster=0;
  //first loop - for "gold cluster" 
  fLoop=1;
  Float_t *b=&fBins[-1]+2*fMaxTime;
  Int_t crtime = Int_t((fParam->GetZLength()-AliTPCReconstructor::GetCtgRange()*fRx)/fZWidth-fParam->GetNTBinsL1()-5);

  for (Int_t i=2*fMaxTime; i<fMaxBin-2*fMaxTime; i++) {
    b++;
    if (*b<8) continue;   //threshold form maxima
    if (i%fMaxTime<crtime) {
      Int_t delta = -(i%fMaxTime)+crtime;
      b+=delta;
      i+=delta;
      continue; 
    }
     
    if (!IsMaximum(*b,fMaxTime,b)) continue;
    AliTPCclusterMI c;
    Int_t dummy=0;
    MakeCluster(i, fMaxTime, fBins, dummy,c);
    //}
  }
  //memcpy(fBins,fResBins, fMaxBin*2);
  //second  loop - for rest cluster 
  /*        
  fLoop=2;
  b=&fResBins[-1]+2*fMaxTime;
  for (Int_t i=2*fMaxTime; i<fMaxBin-2*fMaxTime; i++) {
    b++;
    if (*b<25) continue;   // bigger threshold for maxima
    if (!IsMaximum(*b,fMaxTime,b)) continue;
    AliTPCclusterMI c;
    Int_t dummy;
    MakeCluster(i, fMaxTime, fResBins, dummy,c);
    //}
  }
  */
}
