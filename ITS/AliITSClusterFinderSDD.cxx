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

#include <TFile.h>

#include "AliITSClusterFinderSDD.h"
#include "AliITSMapA1.h"
#include "AliITS.h"
#include "AliITSdigit.h"
#include "AliITSRawCluster.h"
#include "AliITSRecPoint.h"
#include "AliITSsegmentation.h"
#include "AliITSresponse.h"
#include "AliRun.h"



ClassImp(AliITSClusterFinderSDD)

//----------------------------------------------------------
AliITSClusterFinderSDD::AliITSClusterFinderSDD
(AliITSsegmentation *seg, AliITSresponse *response, TClonesArray *digits, TClonesArray *recp)   
{
  // constructor

    fSegmentation=seg;
    fResponse=response;
    fDigits=digits;
    fClusters=recp;
    fNclusters= fClusters->GetEntriesFast();
    SetCutAmplitude();
    SetDAnode();
    SetDTime();
    SetMinPeak();
    SetMinNCells();
    SetMaxNCells();
    SetTimeCorr();
    fMap=new AliITSMapA1(fSegmentation,fDigits,fCutAmplitude);

}

//_____________________________________________________________________________
AliITSClusterFinderSDD::AliITSClusterFinderSDD()
{
  // constructor
    fSegmentation=0;
    fResponse=0;
    fDigits=0;
    fClusters=0;
    fNclusters=0;
    fMap=0;
    SetCutAmplitude();
    SetDAnode();
    SetDTime();
    SetMinPeak();
    SetMinNCells();
    SetMaxNCells();
    SetTimeCorr();

}

//_____________________________________________________________________________
AliITSClusterFinderSDD::~AliITSClusterFinderSDD()
{
    // destructor

    if(fMap) delete fMap;

}
//__________________________________________________________________________
AliITSClusterFinderSDD::AliITSClusterFinderSDD(const AliITSClusterFinderSDD &source){
  //     Copy Constructor 
  if(&source == this) return;
  this->fClusters = source.fClusters ;
  this->fNclusters = source.fNclusters ;
  this->fMap = source.fMap ;
  this->fCutAmplitude = source.fCutAmplitude ;
  this->fDAnode = source.fDAnode ;
  this->fDTime = source.fDTime ;
  this->fTimeCorr = source.fTimeCorr ;
  this->fMinPeak = source.fMinPeak ;
  this->fMinNCells = source.fMinNCells ;
  this->fMaxNCells = source.fMaxNCells ;
  return;
}

//_________________________________________________________________________
AliITSClusterFinderSDD& 
  AliITSClusterFinderSDD::operator=(const AliITSClusterFinderSDD &source) {
  //    Assignment operator
  if(&source == this) return *this;
  this->fClusters = source.fClusters ;
  this->fNclusters = source.fNclusters ;
  this->fMap = source.fMap ;
  this->fCutAmplitude = source.fCutAmplitude ;
  this->fDAnode = source.fDAnode ;
  this->fDTime = source.fDTime ;
  this->fTimeCorr = source.fTimeCorr ;
  this->fMinPeak = source.fMinPeak ;
  this->fMinNCells = source.fMinNCells ;
  this->fMaxNCells = source.fMaxNCells ;
  return *this;
}


//_____________________________________________________________________________

void AliITSClusterFinderSDD::Find1DClusters()
{
  // find 1D clusters

    AliITS *iTS=(AliITS*)gAlice->GetModule("ITS");

   // retrieve the parameters 
   Int_t fNofMaps = fSegmentation->Npz();
   Int_t fMaxNofSamples = fSegmentation->Npx();
   Int_t fNofAnodes = fNofMaps/2;
   Int_t dummy=0;
   Float_t fTimeStep = fSegmentation->Dpx(dummy);
   Float_t fSddLength = fSegmentation->Dx();
   Float_t fDriftSpeed = fResponse->DriftSpeed();
   
   Float_t anodePitch = fSegmentation->Dpz(dummy);
   // map the signal
   fMap->SetThreshold(fCutAmplitude);
   fMap->FillMap();

   Float_t noise;
   Float_t baseline;
   fResponse->GetNoiseParam(noise,baseline);

   Float_t maxadc = fResponse->MaxAdc();    
   Float_t topValue = fResponse->MagicValue();
   Float_t norm = maxadc/topValue;

   Int_t nofFoundClusters = 0;
   Int_t i;
   Float_t **dfadc = new Float_t*[fNofAnodes];
   for(i=0;i<fNofAnodes;i++) dfadc[i] = new Float_t[fMaxNofSamples];
   Float_t fadc, fadc1, fadc2;
   Int_t j,k,idx,l,m;
   for(j=0;j<2;j++) {
     for(k=0;k<fNofAnodes;k++) {
       idx = j*fNofAnodes+k;
       // signal (fadc) & derivative (dfadc)
       dfadc[k][255]=0.;
       for(l=0; l<fMaxNofSamples; l++) {
	 fadc2=(Float_t)fMap->GetSignal(idx,l);
	 if(l>0) fadc1=(Float_t)fMap->GetSignal(idx,l-1);
	 if(l>0) dfadc[k][l-1] = fadc2-fadc1;
       } // samples
     } // anodes
    
     for(k=0;k<fNofAnodes;k++) {
       //cout << "Anode: " << k+1 << ", Wing: " << j+1 << endl;
       idx = j*fNofAnodes+k;
 
       Int_t imax = 0;
       Int_t imaxd = 0;
       Int_t it=0;
       while(it <= fMaxNofSamples-3) {
	
	 imax = it;
	 imaxd = it;
	 // maximum of signal	  
	
	 Float_t fadcmax = 0.;
	 Float_t dfadcmax = 0.;
	 Int_t lthrmina = 1;
	 //	if(it >= 60) lthrmina = 2;
	 //      if(it >= 100) lthrmina = 3;
	 Int_t lthrmint = 2;
	 //if(it >= 60) lthrmint = 3;
	 //if(it >= 100) lthrmint = 4;

	 Int_t lthra = 1;
	 Int_t lthrt = 0;
	
	 for(m=0;m<20;m++) {
	   Int_t id = it+m;
	   if(id>=fMaxNofSamples) break;
	   fadc=(float)fMap->GetSignal(idx,id);
	   if(fadc > fadcmax) { fadcmax = fadc; imax = id;}
	   if(fadc > (float)fCutAmplitude) { 
	     lthrt++; 
	   }

	   if(dfadc[k][id] > dfadcmax) {
	     dfadcmax = dfadc[k][id];
	     imaxd = id;
	   }
	 }
	 it = imaxd;
	
	 if(fMap->TestHit(idx,imax) == kEmpty) {it++; continue;}

	 // cluster charge
	 Int_t tstart = it-1;
	
	 Bool_t ilcl = 0;
	 if(lthrt >= lthrmint && lthra >= lthrmina) ilcl = 1;
	 //printf("ilcl %d\n",ilcl);

	 Float_t b,n;
	 fResponse->GetNoiseParam(n,b);

	 if(ilcl) {
	   nofFoundClusters++;
	   Int_t tstop = tstart;
	   Float_t dfadcmin = 10000.;
	   Int_t ij;
	   for(ij=0; ij<20; ij++) {
	     if(dfadc[k][it+ij] < dfadcmin) {
	       tstop = it+ij+1;
	       dfadcmin = dfadc[k][it+ij];
	     }
	   }

	   Float_t clusterCharge = 0.;
	   Float_t clusterAnode = k+0.5;
	   Float_t clusterTime = 0.;
	   Float_t clusterMult = 0.;
	   Float_t clusterPeakAmplitude = 0.;
	   Int_t its,peakpos=-1;
	   Float_t n, baseline;
	   fResponse->GetNoiseParam(n,baseline);
	   n *= norm;
	   baseline *= norm;
	   for(its=tstart; its<=tstop; its++) {
	     fadc=(float)fMap->GetSignal(idx,its);
	     if(fadc>baseline)
	       fadc-=baseline;
	     else
	       fadc=0.;
	     clusterCharge += fadc;
	     // as a matter of fact we should take the peak pos before FFT
	     // to get the list of tracks !!!
	     if(fadc > clusterPeakAmplitude) {
	       clusterPeakAmplitude = fadc;
	       //peakpos=fMap->GetHitIndex(idx,its);
	       Int_t shift=(int)(fTimeCorr/fTimeStep);
	       if(its>shift && its<(fMaxNofSamples-shift)) peakpos=fMap->GetHitIndex(idx,its+shift);
	       else peakpos=fMap->GetHitIndex(idx,its);
	       if(peakpos<0) peakpos=fMap->GetHitIndex(idx,its);
	     }
	     clusterTime += fadc*its;
	     clusterMult++;
	     if(its == tstop) {
	       // charge from ADC back to nA 
	       clusterCharge /= norm;
	       if(clusterCharge <= 0.) printf("clusterCharge %f norm %f\n",clusterCharge,norm);
	       clusterTime /= (clusterCharge/fTimeStep);   // ns
	       clusterCharge *= (fTimeStep/160.);          // keV
	       if(clusterTime > fTimeCorr) clusterTime -= fTimeCorr;   // ns
	     }
	   }
	   // cout << "Anode: " << k << ", tstart: " << tstart << ", tstop: " << tstop << ", Charge: " << clusterCharge << endl;

	   Float_t clusteranodePath = (clusterAnode - fNofAnodes/2)*anodePitch;
	   Float_t clusterDriftPath = clusterTime*fDriftSpeed;
	   clusterDriftPath = fSddLength-clusterDriftPath;

	   if(clusterCharge <= 0.) break;

	   AliITSRawClusterSDD clust(j+1,clusterAnode,clusterTime,clusterCharge,clusterPeakAmplitude,peakpos,0.,0.,clusterDriftPath,clusteranodePath,clusterMult);
	   iTS->AddCluster(1,&clust);
	   it = tstop;
	} // ilcl
	
	it++;
	
      } // while (samples)
    } // anodes
  } // detectors (2)


   //fMap->ClearMap(); 
  
  for(i=0;i<fNofAnodes;i++) delete[] dfadc[i];
  delete [] dfadc;

  return;

}

//_____________________________________________________________________________
void  AliITSClusterFinderSDD::GroupClusters()
{
  // group clusters
  Int_t dummy=0;
  Float_t fTimeStep = fSegmentation->Dpx(dummy);


  // get number of clusters for this module
  Int_t nofClusters = fClusters->GetEntriesFast();
  nofClusters -= fNclusters;

  AliITSRawClusterSDD *clusterI;
  AliITSRawClusterSDD *clusterJ;

  Int_t *label = new Int_t [nofClusters];
  Int_t i,j;
  for(i=0; i<nofClusters; i++) label[i] = 0;
  for(i=0; i<nofClusters; i++) { 
    if(label[i] != 0) continue;
    for(j=i+1; j<nofClusters; j++) { 
      if(label[j] != 0) continue;
      clusterI = (AliITSRawClusterSDD*) fClusters->At(i);
      clusterJ = (AliITSRawClusterSDD*) fClusters->At(j);
      // 1.3 good
      if(clusterI->T() < fTimeStep*60) fDAnode = 3.2;
      if(clusterI->T() < fTimeStep*10) fDAnode = 1.2;
      Bool_t pair = clusterI->Brother(clusterJ,fDAnode,fDTime);
      if(!pair) continue;
      //      clusterI->PrintInfo();
      //      clusterJ->PrintInfo();
      clusterI->Add(clusterJ);
      label[j] = 1;
      fClusters->RemoveAt(j);
    } // J clusters  
    label[i] = 1;
  } // I clusters
  fClusters->Compress();
  
  delete [] label;
  return;

}

//_____________________________________________________________________________

void AliITSClusterFinderSDD::SelectClusters()
{
  // get number of clusters for this module
  Int_t nofClusters = fClusters->GetEntriesFast();
  nofClusters -= fNclusters;

  Int_t i;
  for(i=0; i<nofClusters; i++) { 
    AliITSRawClusterSDD *clusterI = (AliITSRawClusterSDD*) fClusters->At(i);
    Int_t rmflg = 0;
    Float_t wy = 0.;
    if(clusterI->Anodes() != 0.) {
      wy = ((Float_t) clusterI->Samples())/clusterI->Anodes();
    }
    Float_t amp = clusterI->PeakAmpl();
    if(amp < fMinPeak) rmflg = 1;  
    if(wy < fMinNCells) rmflg = 1;
    //if(wy > fMaxNCells) rmflg = 1;
    if(rmflg) fClusters->RemoveAt(i);
  } // I clusters
  fClusters->Compress();
  return;

}

//_____________________________________________________________________________

void AliITSClusterFinderSDD::GetRecPoints()
{
  // get rec points

  AliITS *iTS=(AliITS*)gAlice->GetModule("ITS");

  // get number of clusters for this module
  Int_t nofClusters = fClusters->GetEntriesFast();
  nofClusters -= fNclusters;

  const Float_t kconvGeV = 1.e-6; // GeV -> KeV
  const Float_t kconv = 1.0e-4; 
  const Float_t kRMSx = 38.0*kconv; // microns->cm ITS TDR Table 1.3
  const Float_t kRMSz = 28.0*kconv; // microns->cm ITS TDR Table 1.3


  Int_t i;
  Int_t ix, iz, idx=-1;
  AliITSdigitSDD *dig=0;
  Int_t maxt=fSegmentation->Npx();
  Int_t ndigits=fDigits->GetEntriesFast();
  for(i=0; i<nofClusters; i++) { 
    AliITSRawClusterSDD *clusterI = (AliITSRawClusterSDD*)fClusters->At(i);
    if(!clusterI) Error("SDD: GetRecPoints","i clusterI ",i,clusterI);
    if(clusterI) idx=clusterI->PeakPos();
    if(idx>ndigits) Error("SDD: GetRecPoints","idx ndigits",idx,ndigits);
    // try peak neighbours - to be done 
    if(idx && idx <= ndigits) dig = (AliITSdigitSDD*)fDigits->UncheckedAt(idx);
    if(!dig) {
	// try cog
	fSegmentation->GetPadIxz(clusterI->X(),clusterI->Z(),ix,iz);
	dig = (AliITSdigitSDD*)fMap->GetHit(iz-1,ix-1);
	// if null try neighbours
	if (!dig) dig = (AliITSdigitSDD*)fMap->GetHit(iz-1,ix); 
	if (!dig) dig = (AliITSdigitSDD*)fMap->GetHit(iz-1,ix+1); 
        if (!dig) printf("SDD: cannot assign the track number!\n");
    }

    AliITSRecPoint rnew;
    rnew.SetX(clusterI->X());
    rnew.SetZ(clusterI->Z());
    rnew.SetQ(clusterI->Q());   // in KeV - should be ADC
    rnew.SetdEdX(kconvGeV*clusterI->Q());
    rnew.SetSigmaX2(kRMSx*kRMSx);
    rnew.SetSigmaZ2(kRMSz*kRMSz);
    if(dig) rnew.fTracks[0]=dig->fTracks[0];
    if(dig) rnew.fTracks[1]=dig->fTracks[1];
    if(dig) rnew.fTracks[2]=dig->fTracks[2];
    //printf("SDD: i %d track1 track2 track3 %d %d %d x y %f %f\n",i,rnew.fTracks[0],rnew.fTracks[1],rnew.fTracks[2],clusterI->X(),clusterI->Z());
    iTS->AddRecPoint(rnew);
  } // I clusters

  fMap->ClearMap();
}

//_____________________________________________________________________________

void AliITSClusterFinderSDD::FindRawClusters()
{
  // find raw clusters
    Find1DClusters();
    GroupClusters();
    SelectClusters();
    GetRecPoints();
}
