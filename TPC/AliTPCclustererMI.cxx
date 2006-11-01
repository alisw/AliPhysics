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
#include "TGraph.h"
#include "TF1.h"
#include "TRandom.h"
#include "AliMathBase.h"

#include "AliTPCClustersArray.h"
#include "AliTPCClustersRow.h"
#include "AliDigits.h"
#include "AliSimDigits.h"
#include "AliTPCParam.h"
#include "AliTPCRecoParam.h"
#include "AliRawReader.h"
#include "AliTPCRawStream.h"
#include "AliRawEventHeaderBase.h"
#include "AliRunLoader.h"
#include "AliLoader.h"
#include "Riostream.h"
#include <TTree.h>
#include "AliTPCcalibDB.h"
#include "AliTPCCalPad.h"
#include "AliTPCCalROC.h"
#include "TTreeStream.h"
#include "AliLog.h"


ClassImp(AliTPCclustererMI)



AliTPCclustererMI::AliTPCclustererMI(const AliTPCParam* par, const AliTPCRecoParam * recoParam):
  fBins(0),
  fResBins(0),
  fLoop(0),
  fMaxBin(0),
  fMaxTime(0),
  fMaxPad(0),
  fSector(-1),
  fRow(-1),
  fSign(0),
  fRx(0),
  fPadWidth(0),
  fPadLength(0),
  fZWidth(0),
  fPedSubtraction(kFALSE),
  fIsOldRCUFormat(kFALSE),
  fEventHeader(0),
  fTimeStamp(0),
  fEventType(0),
  fInput(0),
  fOutput(0),
  fRowCl(0),
  fRowDig(0),
  fParam(0),
  fNcluster(0),
  fAmplitudeHisto(0),
  fDebugStreamer(0),
  fRecoParam(0)
{
  //
  // COSNTRUCTOR
  // param     - tpc parameters for given file
  // recoparam - reconstruction parameters 
  //
  fIsOldRCUFormat = kFALSE;
  fInput =0;
  fOutput=0;
  fParam = par;
  if (recoParam) {
    fRecoParam = recoParam;
  }else{
    //set default parameters if not specified
    fRecoParam = AliTPCReconstructor::GetRecoParam();
    if (!fRecoParam)  fRecoParam = AliTPCRecoParam::GetLowFluxParam();
  }
  fDebugStreamer = new TTreeSRedirector("TPCsignal.root");
  fAmplitudeHisto = 0;
}

AliTPCclustererMI::~AliTPCclustererMI(){
  DumpHistos();
  if (fAmplitudeHisto) delete fAmplitudeHisto;
  if (fDebugStreamer) delete fDebugStreamer;
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
  //
  //  k    - Make cluster at position k  
  //  bins - 2 D array of signals mapped to 1 dimensional array - 
  //  max  - the number of time bins er one dimension
  //  c    - refernce to cluster to be filled
  //
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
  if ( (ry <1.2) && (rz<1.2) || (!fRecoParam->GetDoUnfold())) {
    //
    //if cluster looks like expected or Unfolding not switched on
    //standard COG is used
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
  if (!fRecoParam->GetBYMirror()){
    if (fSector%36>17){
      c.SetY(-(meani - 2.5 - 0.5*fMaxPad)*fParam->GetPadPitchWidth(fSector));
    }
  }
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
	Float_t gain = gainROC->GetValue(row,digarr.CurrentColumn());
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
// This is a cluster finder for the TPC raw data.
// The method assumes NO ordering of the altro channels.
// The pedestal subtraction can be switched on and off
// using an option of the TPC reconstructor
//-----------------------------------------------------------------

  if (!fOutput) {
    Error("Digits2Clusters", "output tree not initialised");
    return;
  }

  fRowDig = NULL;
  AliTPCROC * roc = AliTPCROC::Instance();
  AliTPCCalPad * gainTPC = AliTPCcalibDB::Instance()->GetPadGainFactor();
  AliTPCRawStream input(rawReader);
  fEventHeader = (AliRawEventHeaderBase*)rawReader->GetEventHeader();
  if (fEventHeader){
    fTimeStamp = fEventHeader->Get("Timestamp");  
    fEventType = fEventHeader->Get("Type");  
  }


  Int_t nclusters  = 0;
  
  fMaxTime = fParam->GetMaxTBin() + 6; // add 3 virtual time bins before and 3 after
  const Int_t kNIS = fParam->GetNInnerSector();
  const Int_t kNOS = fParam->GetNOuterSector();
  const Int_t kNS = kNIS + kNOS;
  fZWidth = fParam->GetZWidth();
  Int_t zeroSup = fParam->GetZeroSup();
  //
  //alocate memory for sector - maximal case
  //
  Float_t** allBins = NULL;
  Float_t** allBinsRes = NULL;
  Int_t nRowsMax = roc->GetNRows(roc->GetNSector()-1);
  Int_t nPadsMax = roc->GetNPads(roc->GetNSector()-1,nRowsMax-1);
  allBins = new Float_t*[nRowsMax];
  allBinsRes = new Float_t*[nRowsMax];
  for (Int_t iRow = 0; iRow < nRowsMax; iRow++) {
    //
    Int_t maxBin = fMaxTime*(nPadsMax+6);  // add 3 virtual pads  before and 3 after
    allBins[iRow] = new Float_t[maxBin];
    allBinsRes[iRow] = new Float_t[maxBin];
    memset(allBins[iRow],0,sizeof(Float_t)*maxBin);
    memset(allBinsRes[iRow],0,sizeof(Float_t)*maxBin);
  }
  //
  // Loop over sectors
  //
  for(fSector = 0; fSector < kNS; fSector++) {

    AliTPCCalROC * gainROC = gainTPC->GetCalROC(fSector);  // pad gains per given sector
 
    Int_t nRows = 0;
    Int_t nDDLs = 0, indexDDL = 0;
    if (fSector < kNIS) {
      nRows = fParam->GetNRowLow();
      fSign = (fSector < kNIS/2) ? 1 : -1;
      nDDLs = 2;
      indexDDL = fSector * 2;
    }
    else {
      nRows = fParam->GetNRowUp();
      fSign = ((fSector-kNIS) < kNOS/2) ? 1 : -1;
      nDDLs = 4;
      indexDDL = (fSector-kNIS) * 4 + kNIS * 2;
    }

    for (Int_t iRow = 0; iRow < nRows; iRow++) {
      Int_t maxPad;
      if (fSector < kNIS)
	maxPad = fParam->GetNPadsLow(iRow);
      else
	maxPad = fParam->GetNPadsUp(iRow);
      
      Int_t maxBin = fMaxTime*(maxPad+6);  // add 3 virtual pads  before and 3 after
      memset(allBins[iRow],0,sizeof(Float_t)*maxBin);
      memset(allBinsRes[iRow],0,sizeof(Float_t)*maxBin);
    }
    
    // Loas the raw data for corresponding DDLs
    rawReader->Reset();
    input.SetOldRCUFormat(fIsOldRCUFormat);
    rawReader->Select("TPC",indexDDL,indexDDL+nDDLs-1);
    Int_t digCounter=0;
    // Begin loop over altro data
    Bool_t calcPedestal = fRecoParam->GetCalcPedestal();
    Float_t gain =1;
    Int_t lastPad=-1;
    while (input.Next()) {
      digCounter++;
      if (input.GetSector() != fSector)
	AliFatal(Form("Sector index mismatch ! Expected (%d), but got (%d) !",fSector,input.GetSector()));

      
      Int_t iRow = input.GetRow();
      if (iRow < 0 || iRow >= nRows)
	AliFatal(Form("Pad-row index (%d) outside the range (%d -> %d) !",
		      iRow, 0, nRows -1));
      //pad
      Int_t iPad = input.GetPad();
      if (iPad < 0 || iPad >= nPadsMax)
	AliFatal(Form("Pad index (%d) outside the range (%d -> %d) !",
		      iPad, 0, nPadsMax-1));
      if (iPad!=lastPad){
	gain    = gainROC->GetValue(iRow,iPad);
	lastPad = iPad;
      }
      iPad+=3;
      //time
      Int_t iTimeBin = input.GetTime();
      if ( iTimeBin < 0 || iTimeBin >= fParam->GetMaxTBin())
	AliFatal(Form("Timebin index (%d) outside the range (%d -> %d) !",
		      iTimeBin, 0, iTimeBin -1));
      iTimeBin+=3;
      //signal
      Float_t signal = input.GetSignal();
      if (!calcPedestal && signal <= zeroSup) continue;      
      allBins[iRow][iPad*fMaxTime+iTimeBin] = signal/gain;
      allBins[iRow][iPad*fMaxTime+0] = 1;  // pad with signal
    } // End of the loop over altro data
    //
    //
    // Now loop over rows and perform pedestal subtraction
    if (digCounter==0) continue;
    //    if (fPedSubtraction) {
    if (calcPedestal) {
      for (Int_t iRow = 0; iRow < nRows; iRow++) {
	Int_t maxPad;
	if (fSector < kNIS)
	  maxPad = fParam->GetNPadsLow(iRow);
	else
	  maxPad = fParam->GetNPadsUp(iRow);

	for (Int_t iPad = 0; iPad < maxPad + 6; iPad++) {
	  if (allBins[iRow][iPad*fMaxTime+0] !=1) continue;  // no data
	  Float_t *p = &allBins[iRow][iPad*fMaxTime+3];
	  //Float_t pedestal = TMath::Median(fMaxTime, p);	
	  Int_t id[3] = {fSector, iRow, iPad-3};
	  Double_t rms=0;
	  Float_t pedestal = ProcesSignal(p, fMaxTime, id, rms);
	  for (Int_t iTimeBin = 0; iTimeBin < fMaxTime; iTimeBin++) {
	    allBins[iRow][iPad*fMaxTime+iTimeBin] -= pedestal;
	    if (iTimeBin < AliTPCReconstructor::GetRecoParam()->GetFirstBin())  
	      allBins[iRow][iPad*fMaxTime+iTimeBin] = 0;
	    if (iTimeBin > AliTPCReconstructor::GetRecoParam()->GetLastBin())  
	      allBins[iRow][iPad*fMaxTime+iTimeBin] = 0;
	    if (allBins[iRow][iPad*fMaxTime+iTimeBin] < zeroSup)
	      allBins[iRow][iPad*fMaxTime+iTimeBin] = 0;
	    if (allBins[iRow][iPad*fMaxTime+iTimeBin] < 3.0*rms)   // 3 sigma cut on RMS
	      allBins[iRow][iPad*fMaxTime+iTimeBin] = 0;	    
	  }
	}
      }
    }
    // Now loop over rows and find clusters
    for (fRow = 0; fRow < nRows; fRow++) {
      fRowCl = new AliTPCClustersRow;
      fRowCl->SetClass("AliTPCclusterMI");
      fRowCl->SetArray(1);
      fRowCl->SetID(fParam->GetIndex(fSector, fRow));
      fOutput->GetBranch("Segment")->SetAddress(&fRowCl);

      fRx = fParam->GetPadRowRadii(fSector, fRow);
      fPadLength = fParam->GetPadPitchLength(fSector, fRow);
      fPadWidth  = fParam->GetPadPitchWidth();
      if (fSector < kNIS)
	fMaxPad = fParam->GetNPadsLow(fRow);
      else
	fMaxPad = fParam->GetNPadsUp(fRow);
      fMaxBin = fMaxTime*(fMaxPad+6);  // add 3 virtual pads  before and 3 after

      fBins = allBins[fRow];
      fResBins = allBinsRes[fRow];

      FindClusters();

      fOutput->Fill();
      delete fRowCl;    
      nclusters += fNcluster;    
    } // End of loop to find clusters
  } // End of loop over sectors
  
  for (Int_t iRow = 0; iRow < nRowsMax; iRow++) {
    delete [] allBins[iRow];
    delete [] allBinsRes[iRow];
  }  
  delete [] allBins;
  delete [] allBinsRes;
  
  Info("Digits2Clusters", "Event\t%d\tNumber of found clusters : %d\n", rawReader->GetEventId(), nclusters);

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
  Int_t crtime = Int_t((fParam->GetZLength()-fRecoParam->GetCtgRange()*fRx)/fZWidth-fParam->GetNTBinsL1()-5);

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


Double_t AliTPCclustererMI::ProcesSignal(Float_t *signal, Int_t nchannels, Int_t id[3], Double_t &rmsOut){
  //
  // process signal on given pad - + streaming of additional information in special mode
  //
  // id[0] - sector
  // id[1] - row
  // id[2] - pad 

  //
  // ESTIMATE pedestal and the noise
  // 
  Int_t offset =100;
  const Int_t kPedMax = 100;
  Float_t  max    =  0;
  Float_t  maxPos =  0;
  Int_t    median =  -1;
  Int_t    count0 =  0;
  Int_t    count1 =  0;
  //
  UShort_t histo[kPedMax];
  memset(histo,0,kPedMax*sizeof(UShort_t));
  for (Int_t i=0; i<fMaxTime; i++){
    if (signal[i]<=0) continue;
    if (signal[i]>max) {
      max = signal[i];
      maxPos = i;
    }
    if (signal[i]>kPedMax-1) continue;
    histo[int(signal[i]+0.5)]++;
    count0++;
  }
  //
  for (Int_t i=1; i<kPedMax; i++){
    if (count1<count0*0.5) median=i;
    count1+=histo[i];
  }
  // truncated mean  
  //
  Double_t count10=histo[median] ,mean=histo[median]*median,  rms=histo[median]*median*median ;
  Double_t count06=histo[median] ,mean06=histo[median]*median,  rms06=histo[median]*median*median ;
  Double_t count09=histo[median] ,mean09=histo[median]*median,  rms09=histo[median]*median*median ;
  //
  for (Int_t idelta=1; idelta<10; idelta++){
    if (median-idelta<=0) continue;
    if (median+idelta>kPedMax) continue;
    if (count06<0.6*count1){
      count06+=histo[median-idelta];
      mean06 +=histo[median-idelta]*(median-idelta);
      rms06  +=histo[median-idelta]*(median-idelta)*(median-idelta);
      count06+=histo[median+idelta];
      mean06 +=histo[median+idelta]*(median+idelta);
      rms06  +=histo[median+idelta]*(median+idelta)*(median+idelta);
    }
    if (count09<0.9*count1){
      count09+=histo[median-idelta];
      mean09 +=histo[median-idelta]*(median-idelta);
      rms09  +=histo[median-idelta]*(median-idelta)*(median-idelta);
      count09+=histo[median+idelta];
      mean09 +=histo[median+idelta]*(median+idelta);
      rms09  +=histo[median+idelta]*(median+idelta)*(median+idelta);
    }
    if (count10<0.95*count1){
      count10+=histo[median-idelta];
      mean +=histo[median-idelta]*(median-idelta);
      rms  +=histo[median-idelta]*(median-idelta)*(median-idelta);
      count10+=histo[median+idelta];
      mean +=histo[median+idelta]*(median+idelta);
      rms  +=histo[median+idelta]*(median+idelta)*(median+idelta);
    }
  }
  mean  /=count10;
  mean06/=count06;
  mean09/=count09;
  rms    = TMath::Sqrt(TMath::Abs(rms/count10-mean*mean));
  rms06    = TMath::Sqrt(TMath::Abs(rms06/count06-mean06*mean06));
  rms09    = TMath::Sqrt(TMath::Abs(rms09/count09-mean09*mean09));
  rmsOut = rms09;
  //
  if (AliLog::GetDebugLevel("","AliTPCclustererMI")==0) return median;
  //
  UInt_t uid[3] = {UInt_t(id[0]),UInt_t(id[1]),UInt_t(id[2])};
  //
  // Dump mean signal info
  //
  (*fDebugStreamer)<<"Signal"<<
    "TimeStamp="<<fTimeStamp<<
    "EventType="<<fEventType<<
    "Sector="<<uid[0]<<
    "Row="<<uid[1]<<
    "Pad="<<uid[2]<<
    "Max="<<max<<
    "MaxPos="<<maxPos<<
    //
    "Median="<<median<<
    "Mean="<<mean<<
    "RMS="<<rms<<      
    "Mean06="<<mean06<<
    "RMS06="<<rms06<<
    "Mean09="<<mean09<<
    "RMS09="<<rms09<<
    "\n";
  //
  // fill pedestal histogram
  //
  AliTPCROC * roc = AliTPCROC::Instance();
  if (!fAmplitudeHisto){
    fAmplitudeHisto = new TObjArray(72);
  }  
  //
  if (uid[0]<roc->GetNSectors() 
      && uid[1]< roc->GetNRows(uid[0])  && 
      uid[2] <roc->GetNPads(uid[0], uid[1])){
    TObjArray  * sectorArray = (TObjArray*)fAmplitudeHisto->UncheckedAt(uid[0]);
    if (!sectorArray){
      Int_t npads =roc->GetNChannels(uid[0]);
      sectorArray = new TObjArray(npads);
      fAmplitudeHisto->AddAt(sectorArray, uid[0]);
    }
    Int_t position =  uid[2]+roc->GetRowIndexes(uid[0])[uid[1]];
    TH1F * histo = (TH1F*)sectorArray->UncheckedAt(position);
    if (!histo){
      char hname[100];
      sprintf(hname,"Amp_%d_%d_%d",uid[0],uid[1],uid[2]);
      TFile * backup = gFile;
      fDebugStreamer->GetFile()->cd();
      histo = new TH1F(hname, hname, 100, 5,100);
      //histo->SetDirectory(0);     // histogram not connected to directory -(File)
      sectorArray->AddAt(histo, position);
      if (backup) backup->cd();
    }
    for (Int_t i=0; i<nchannels; i++){
      histo->Fill(signal[i]);
    }
  }
  //
  //
  //
  Float_t kMin =fRecoParam->GetDumpAmplitudeMin();   // minimal signal to be dumped
  Float_t *dsignal = new Float_t[nchannels];
  Float_t *dtime   = new Float_t[nchannels];
  for (Int_t i=0; i<nchannels; i++){
    dtime[i] = i;
    dsignal[i] = signal[i];
  }
  //
  //
  TGraph * graph;
  Bool_t random = (gRandom->Rndm()<0.0001);
  if (max-median>kMin || rms06>2.*fParam->GetZeroSup() || random){
    graph =new TGraph(nchannels, dtime, dsignal);
    if (rms06>2.*fParam->GetZeroSup() || random)
      (*fDebugStreamer)<<"SignalN"<<    //noise pads - or random sample of pads
	"TimeStamp="<<fTimeStamp<<
	"EventType="<<fEventType<<
	"Sector="<<uid[0]<<
	"Row="<<uid[1]<<
	"Pad="<<uid[2]<<
	"Graph="<<graph<<
	"Max="<<max<<
	"MaxPos="<<maxPos<<
	//
	"Median="<<median<<
	"Mean="<<mean<<
	"RMS="<<rms<<      
	"Mean06="<<mean06<<
	"RMS06="<<rms06<<
	"Mean09="<<mean09<<
	"RMS09="<<rms09<<
	"\n";
    if (max-median>kMin &&maxPos>AliTPCReconstructor::GetRecoParam()->GetFirstBin()) 
      (*fDebugStreamer)<<"SignalB"<<     // pads with signal
	"TimeStamp="<<fTimeStamp<<
	"EventType="<<fEventType<<
	"Sector="<<uid[0]<<
	"Row="<<uid[1]<<
	"Pad="<<uid[2]<<
	"Graph="<<graph<<
	"Max="<<max<<
	"MaxPos="<<maxPos<<	
	//
	"Median="<<median<<
	"Mean="<<mean<<
	"RMS="<<rms<<      
	"Mean06="<<mean06<<
	"RMS06="<<rms06<<
	"Mean09="<<mean09<<
	"RMS09="<<rms09<<
	"\n";
    delete graph;
  }
  
  //
  //
  //  Central Electrode signal analysis  
  //
  Double_t ceQmax  =0, ceQsum=0, ceTime=0;
  Double_t cemean  = mean06, cerms=rms06 ;
  Int_t    cemaxpos= 0;
  Double_t ceThreshold=5.*cerms;
  Double_t ceSumThreshold=8.*cerms;
  const Int_t    kCemin=5;  // range for the analysis of the ce signal +- channels from the peak
  const Int_t    kCemax=5;
  for (Int_t i=nchannels-2; i>nchannels/2; i--){
    if ( (dsignal[i]-mean06)>ceThreshold && dsignal[i]>=dsignal[i+1] && dsignal[i]>=dsignal[i-1] ){
      cemaxpos=i;
      break;
    }
  }
  if (cemaxpos!=0){
    ceQmax = 0;
    Int_t cemaxpos2=0;
    for (Int_t i=cemaxpos-20; i<cemaxpos+5; i++){
      if (i<0 || i>nchannels-1) continue;
      Double_t val=dsignal[i]- cemean;
      if (val>ceQmax){
	cemaxpos2=i;
	ceQmax = val;
      }
    }
    cemaxpos = cemaxpos2;
 
    for (Int_t i=cemaxpos-kCemin; i<cemaxpos+kCemax; i++){
      if (i>0 && i<nchannels&&dsignal[i]- cemean>0){
	Double_t val=dsignal[i]- cemean;
	ceTime+=val*dtime[i];
	ceQsum+=val;
	if (val>ceQmax) ceQmax=val;
      }
    }
    if (ceQmax&&ceQsum>ceSumThreshold) {
      ceTime/=ceQsum;
      (*fDebugStreamer)<<"Signalce"<<
	"TimeStamp="<<fTimeStamp<<
	"EventType="<<fEventType<<
	"Sector="<<uid[0]<<
	"Row="<<uid[1]<<
	"Pad="<<uid[2]<<
	"Max="<<ceQmax<<
	"Qsum="<<ceQsum<<
	"Time="<<ceTime<<
	"RMS06="<<rms06<<
	//
	"\n";
    }
  }
  // end of ce signal analysis
  //

  //
  //  Gating grid signal analysis  
  //
  Double_t ggQmax  =0, ggQsum=0, ggTime=0;
  Double_t ggmean  = mean06, ggrms=rms06 ;
  Int_t    ggmaxpos= 0;
  Double_t ggThreshold=5.*ggrms;
  Double_t ggSumThreshold=8.*ggrms;

  for (Int_t i=1; i<nchannels/4; i++){
    if ( (dsignal[i]-mean06)>ggThreshold && dsignal[i]>=dsignal[i+1] && dsignal[i]>=dsignal[i-1] &&
	 (dsignal[i]+dsignal[i+1]+dsignal[i-1]-3*mean06)>ggSumThreshold){
      ggmaxpos=i;
      if (dsignal[i-1]>dsignal[i+1]) ggmaxpos=i-1;
      break;
    }
  }
  if (ggmaxpos!=0){
      for (Int_t i=ggmaxpos-1; i<ggmaxpos+3; i++){       
	  if (i>0 && i<nchannels && dsignal[i]-ggmean>0){
	      Double_t val=dsignal[i]- ggmean;
	      ggTime+=val*dtime[i];
	      ggQsum+=val;
	      if (val>ggQmax) ggQmax=val;
	  }
      }
      if (ggQmax&&ggQsum>ggSumThreshold) {
	  ggTime/=ggQsum;
	  (*fDebugStreamer)<<"Signalgg"<<
	    "TimeStamp="<<fTimeStamp<<
	    "EventType="<<fEventType<<
	      "Sector="<<uid[0]<<
	      "Row="<<uid[1]<<
	      "Pad="<<uid[2]<<
	      "Max="<<ggQmax<<
	      "Qsum="<<ggQsum<<
	      "Time="<<ggTime<<
	      "RMS06="<<rms06<<
	      //
	      "\n";
      }
  }
  // end of gg signal analysis
      

  delete [] dsignal;
  delete [] dtime;
  if (rms06>fRecoParam->GetMaxNoise()) return 1024+median; // sign noisy channel in debug mode
  return median;
}



void AliTPCclustererMI::DumpHistos(){
  //
  // Dump histogram information
  //
  if (!fAmplitudeHisto) return;
  AliTPCROC * roc = AliTPCROC::Instance();
  for (UInt_t isector=0; isector<AliTPCROC::Instance()->GetNSectors(); isector++){
    TObjArray * array = (TObjArray*)fAmplitudeHisto->UncheckedAt(isector);
    if (!array) continue;
    for (UInt_t ipad = 0; ipad <(UInt_t)array->GetEntriesFast(); ipad++){
      TH1F * histo = (TH1F*) array->UncheckedAt(ipad);
      if (!histo) continue;
      if (histo->GetEntries()<100) continue;
      histo->Fit("gaus","q");
      Float_t mean =  histo->GetMean();
      Float_t rms  =  histo->GetRMS();
      Float_t gmean = histo->GetFunction("gaus")->GetParameter(1);
      Float_t gsigma = histo->GetFunction("gaus")->GetParameter(2);
      Float_t max = histo->GetFunction("gaus")->GetParameter(0);

      // get pad number
      UInt_t row=0, pad =0;
      const UInt_t *indexes =roc->GetRowIndexes(isector);
      for (UInt_t irow=0; irow<roc->GetNRows(isector); irow++){
	if (indexes[irow]<ipad){
	  row = irow;
	  pad = ipad-indexes[irow];
	}
      }      
      Int_t rpad = pad - (AliTPCROC::Instance()->GetNPads(isector,row))/2;
      //
      (*fDebugStreamer)<<"Fit"<<
	"TimeStamp="<<fTimeStamp<<
	"EventType="<<fEventType<<
	"Sector="<<isector<<
	"Row="<<row<<
	"Pad="<<pad<<
	"RPad="<<rpad<<
	"Max="<<max<<
	"Mean="<<mean<<
	"RMS="<<rms<<      
	"GMean="<<gmean<<
	"GSigma="<<gsigma<<
	"\n";
      if (array->UncheckedAt(ipad)) fDebugStreamer->StoreObject(array->UncheckedAt(ipad));
    }
  }
}
