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


#include "AliITSClusterFinderSDD.h"
#include "AliITSMapA2.h"
#include "AliITSMapA1.h"
#include "AliITS.h"
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
    SetNCells();
    SetMap();
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
    fMapA2=0;
    SetCutAmplitude();
    SetDAnode();
    SetDTime();
    SetMinPeak();
    SetNCells();

}

//_____________________________________________________________________________
AliITSClusterFinderSDD::~AliITSClusterFinderSDD()
{
    // destructor

    if(fMap) delete fMap;
    if(fMapA2) delete fMapA2;

}
//__________________________________________________________________________
AliITSClusterFinderSDD::AliITSClusterFinderSDD(const AliITSClusterFinderSDD &source){
  //     Copy Constructor 
  if(&source == this) return;
  this->fClusters = source.fClusters ;
  this->fNclusters = source.fNclusters ;
  this->fMap = source.fMap ;
  this->fMapA2 = source.fMapA2 ;
  this->fCutAmplitude = source.fCutAmplitude ;
  this->fDAnode = source.fDAnode ;
  this->fDTime = source.fDTime ;
  this->fMinPeak = source.fMinPeak ;
  this->fMinNCells = source.fMinNCells ;
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
  this->fMapA2 = source.fMapA2 ;
  this->fCutAmplitude = source.fCutAmplitude ;
  this->fDAnode = source.fDAnode ;
  this->fDTime = source.fDTime ;
  this->fMinPeak = source.fMinPeak ;
  this->fMinNCells = source.fMinNCells ;
  return *this;
}

//_____________________________________________________________________________
void AliITSClusterFinderSDD::SetMap()
{
  // set map
    if(!fMapA2) fMapA2=new AliITSMapA2(fSegmentation,fDigits,(double)fCutAmplitude);
    if(!fMap) fMap=new AliITSMapA1(fSegmentation,fDigits);

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
   fMapA2->FillMap();

   // Piergiorgio's code - do not subtract baseline since we start
   // from digits and do not duplicate arrays, i.e. use fMap instead
   // of Float_t fadc[2*fNofAnodes][fMaxNofSamples];
 

  Int_t nofFoundClusters = 0;
  Int_t i;
  Float_t **dfadc = new Float_t*[fNofMaps];
  for(i=0;i<fNofMaps;i++) dfadc[i] = new Float_t[fMaxNofSamples];
  Float_t fadc, fadc1, fadc2;
  Int_t j,k,idx,l,m;
  for(j=0;j<2;j++) {
    for(k=0;k<fNofAnodes;k++) {
      idx = j*fNofAnodes+k;
      // signal (fadc) & derivative (dfadc)
      for(l=0; l<fMaxNofSamples; l++) {
	fadc2=(Float_t)fMapA2->GetSignal(idx,l);
	fadc1=(Float_t)fMapA2->GetSignal(idx,l-1);
	if(l>0) dfadc[k][l-1] = fadc2-fadc1;
      } // samples
    } // anodes
    
    for(k=0;k<fNofAnodes;k++) {
      idx = j*fNofAnodes+k;
 
      Int_t imax = 0;
      Int_t imaxd = 0;
      Int_t it=0;
      while(it <= fMaxNofSamples-3) {
	//	cout << "sample: " << it << endl;
	
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

	Int_t lthra = 0;
	Int_t lthrt = 0;
	for(m=0;m<10;m++) {
	  Int_t id = it+m;
	  if(id>=fMaxNofSamples) break;
          fadc=fMapA2->GetSignal(idx,id);
	  if(fadc > fadcmax) {
	    fadcmax = fadc;
	    if(fadc > 0) { lthra++; lthrt++; }
	    imax = id;
	  }
	  if(dfadc[k][id] > dfadcmax) {
	    dfadcmax = dfadc[k][id];
	    imaxd = id;
	  }
	}
	it = imaxd;
	// skip if no signal over threshold
	if(fMapA2->TestHit(idx,imax) == kEmpty) {it++; continue;}

	if(k>0) {
	  if(fMapA2->TestHit(idx-1,imax) != kEmpty) lthra++;
	}
	if(k<fNofAnodes-1)
	  if(fMapA2->TestHit(idx+1,imax) != kEmpty) lthra++;

	if(imax>0) {
	  if(fMapA2->TestHit(idx,imax-1) != kEmpty) lthrt++;
	}
	if(imax<fMaxNofSamples)
	  if(fMapA2->TestHit(idx,imax+1) != kEmpty) lthrt++;
	
	// cluster charge
	Int_t tstart = it-1;
	
	Bool_t ilcl = 0;
	if(lthrt >= lthrmint && lthra >= lthrmina) ilcl = 1;
	
	if(ilcl) {
	  nofFoundClusters++;
	  Int_t tstop = tstart;
	  Float_t dfadcmin = 10000.;
          Int_t ij;
	  for(ij=0; ij<10; ij++) {
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
          Int_t its;
	  for(its=tstart; its<=tstop; its++) {
            fadc=fMapA2->GetSignal(idx,its);
	    clusterCharge += fadc;
	    if(fadc > clusterPeakAmplitude) clusterPeakAmplitude = fadc;
	    clusterTime += fadc*its;
	    clusterMult++;
	    if(its == tstop) {
	       clusterTime /= (clusterCharge/fTimeStep);   // ns
	       clusterCharge *= (fTimeStep/160.);          // keV
	       if(clusterTime > 58.2) clusterTime -= 58.2;   // ns
	      /*
	      else {
		cout << "Warning: cluster with negative time " << clusterTime << ", peak ampl.: " << clusterPeakAmplitude << ", mult: " << clusterMult << ", charge: " << clusterCharge << endl;
		cout << "Anode: " << k << ", tstart: " << tstart << ", tstop: " << tstop << ", Charge: " << clusterCharge << endl;
	      }
	      */
	      }
	  }
	  // cout << "Anode: " << k << ", tstart: " << tstart << ", tstop: " << tstop << ", Charge: " << clusterCharge << endl;

	  Float_t clusteranodePath = (0.06 + clusterAnode - fNofAnodes/2)*anodePitch;
	  Float_t clusterDriftPath = clusterTime*fDriftSpeed;
	  clusterDriftPath = fSddLength-clusterDriftPath;

          if(clusterCharge < 0.) break;

	  //printf("wing clusterMult clusterAnode clusterTime %d %f %f %f \n",j+1,clusterMult, clusterAnode, clusterTime);

	  AliITSRawClusterSDD *clust = new AliITSRawClusterSDD(j+1,clusterAnode,clusterTime,clusterCharge,clusterPeakAmplitude,0.,0.,clusterDriftPath,clusteranodePath,clusterMult);
	  //fClusters->Add(point);
	  iTS->AddCluster(1,clust);
	  it = tstop;
	} // ilcl
	
	it++;
	
      } // while (samples)
    } // anodes
  } // detectors (2)


  fMapA2->ClearMap();
  
  for(i=0;i<fNofMaps;i++) delete[] dfadc[i];
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

  fMap->FillMap();

  Int_t i, ix, iz;
  for(i=0; i<nofClusters; i++) { 
    AliITSRawClusterSDD *clusterI = (AliITSRawClusterSDD*) fClusters->At(i);
    fSegmentation->GetPadIxz(clusterI->X(),clusterI->Z(),ix,iz);
    AliITSdigitSDD *dig = (AliITSdigitSDD*)fMap->GetHit(iz-1,ix-1);
    AliITSRecPoint rnew;
    rnew.SetX(clusterI->X());
    rnew.SetZ(clusterI->Z());
    rnew.SetQ(clusterI->Q());   // in KeV - should be ADC
    rnew.SetdEdX(kconvGeV*clusterI->Q());
    rnew.SetSigmaX2(kRMSx*kRMSx);
    rnew.SetSigmaZ2(kRMSz*kRMSz);
    rnew.fTracks[0]=dig->fTracks[0];
    rnew.fTracks[1]=dig->fTracks[1];
    rnew.fTracks[2]=dig->fTracks[2];
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




