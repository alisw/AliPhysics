/**************************************************************************
 * Copyright(c) 1998-2003, ALICE Experiment at CERN, All rights reserved. *
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
////////////////////////////////////////////////////////////////////////////
//            Implementation of the ITS clusterer V2 class                //
//                                                                        //
//          Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch            //
//                                                                        //
///////////////////////////////////////////////////////////////////////////


#include "AliITSClusterFinderV2SPD.h"
#include "AliITSRecPoint.h"
#include "AliITSDetTypeRec.h"
#include "AliRawReader.h"
#include "AliITSRawStreamSPD.h"
#include <TClonesArray.h>
#include "AliITSdigitSPD.h"

ClassImp(AliITSClusterFinderV2SPD)

extern AliRun *gAlice;

AliITSClusterFinderV2SPD::AliITSClusterFinderV2SPD(AliITSDetTypeRec* dettyp):AliITSClusterFinderV2(dettyp),
fLastSPD1(0),
fNySPD(256),
fNzSPD(160),
fYpitchSPD(0.0050),
fZ1pitchSPD(0.0425),
fZ2pitchSPD(0.0625),
fHwSPD(0.64),
fHlSPD(3.48){

  //Default constructor

  fLastSPD1=fDetTypeRec->GetITSgeom()->GetModuleIndex(2,1,1)-1;
  fYSPD[0]=0.5*fYpitchSPD;
  for (Int_t m=1; m<fNySPD; m++) fYSPD[m]=fYSPD[m-1]+fYpitchSPD; 
  fZSPD[0]=fZ1pitchSPD;
  for (Int_t m=1; m<fNzSPD; m++) {
    Double_t dz=fZ1pitchSPD;
    if (m==31 || m==32 || m==63  || m==64  || m==95 || m==96 || 
        m==127 || m==128) dz=fZ2pitchSPD; 
    fZSPD[m]=fZSPD[m-1]+dz;
  }
  for (Int_t m=0; m<fNzSPD; m++) {
    Double_t dz=0.5*fZ1pitchSPD;
    if (m==31 || m==32 || m==63  || m==64  || m==95 || m==96 || 
        m==127 || m==128) dz=0.5*fZ2pitchSPD; 
    fZSPD[m]-=dz;
  }

}
 

void AliITSClusterFinderV2SPD::FindRawClusters(Int_t mod){

  //Find clusters V2
  SetModule(mod);
  FindClustersSPD(fDigits);

}

void AliITSClusterFinderV2SPD::RawdataToClusters(AliRawReader* rawReader, TClonesArray** clusters){
    //------------------------------------------------------------
  // This function creates ITS clusters from raw data
  //------------------------------------------------------------
  rawReader->Reset();
  AliITSRawStreamSPD inputSPD(rawReader);
  FindClustersSPD(&inputSPD, clusters);

}

Int_t AliITSClusterFinderV2SPD::ClustersSPD(AliBin* bins, TClonesArray* digits,TClonesArray* clusters,Int_t maxBins,Int_t nzbins,Int_t iModule,Bool_t rawdata){
  
  //Cluster finder for SPD (from digits and from rawdata)

  Int_t nclu=0;
  for(Int_t iBin =0; iBin < maxBins;iBin++){
    if(bins[iBin].IsUsed()) continue;
    Int_t nBins = 0;
    Int_t idxBins[200];
    FindCluster(iBin,nzbins,bins,nBins,idxBins);
    if (nBins == 200){
      Error("ClustersSPD","SPD Too big cluster !\n"); 
      continue;
    }
    Int_t milab[10];
    for(Int_t ilab=0;ilab<10;ilab++){
      milab[ilab]=-2;
    }
    if(rawdata){
      milab[3]=fNdet[iModule];
    }
    Int_t ymin,ymax,zmin,zmax;
    if(rawdata){
      ymin = (idxBins[0] / nzbins) - 1;
      ymax = ymin;
      zmin = (idxBins[0] % nzbins) - 1;
      zmax = zmin;
    }
    else{
      AliITSdigitSPD* dig = (AliITSdigitSPD*)digits->UncheckedAt(idxBins[0]);
      ymin=dig->GetCoord2();
      ymax=ymin;
      zmin=dig->GetCoord1();
      zmax=zmin;
    }
    for (Int_t idx = 0; idx < nBins; idx++) {
      Int_t iy;
      Int_t iz; 
      if(rawdata){
	iy  = (idxBins[idx] / nzbins) - 1;
	iz  = (idxBins[idx] % nzbins) - 1;
      }
      else{
	AliITSdigitSPD* dig = (AliITSdigitSPD*)digits->UncheckedAt(idxBins[idx]);
	iy = dig->GetCoord2();
	iz = dig->GetCoord1();
      }
      if (ymin > iy) ymin = iy;
      if (ymax < iy) ymax = iy;
      if (zmin > iz) zmin = iz;
      if (zmax < iz) zmax = iz;

    }
    if(!rawdata){
      for(Int_t l=0;l<nBins;l++){
	AliITSdigitSPD* dig = (AliITSdigitSPD*)digits->UncheckedAt(idxBins[l]);
	for(Int_t dlab=0;dlab<10;dlab++){
	  Int_t digitlab = (dig->GetTracks())[dlab];
	  if(digitlab<0) continue;
	  AddLabel(milab,digitlab);
	}
	if (milab[9]>0) CheckLabels2(milab);
      }
      CheckLabels2(milab);
    }
    
    Int_t idy =0; //max 2 clusters
    if((iModule <= fLastSPD1) &&idy<3) idy=3;
    if((iModule > fLastSPD1) &&idy<4) idy=4;
    Int_t idz=3;
    for(Int_t iiz=zmin; iiz<=zmax;iiz+=idz){
      for(Int_t iiy=ymin;iiy<=ymax;iiy+=idy){

	Int_t ndigits=0;
	Float_t y=0.,z=0.,q=0.;
	for(Int_t idx=0;idx<nBins;idx++){
	  Int_t iy;
	  Int_t iz; 
	  if(rawdata){
	    iy  = (idxBins[idx] / nzbins)-1;
	    iz  = (idxBins[idx] % nzbins)-1;
	  }
	  else{
	    AliITSdigitSPD* dig = (AliITSdigitSPD*)digits->UncheckedAt(idxBins[idx]);
	    iy = dig->GetCoord2();
	    iz = dig->GetCoord1();
	  }
	  if(zmax-zmin>=idz || ymax-ymin>=idy){
	    if(TMath::Abs(iy-iiy)>0.75*idy) continue;
	    if(TMath::Abs(iz-iiz)>0.75*idz) continue;
	  }
	  ndigits++;
	  Float_t qBin;
	  if(rawdata) qBin = bins[idxBins[idx]].GetQ();
	  if(!rawdata){
	    AliITSdigitSPD* dig = (AliITSdigitSPD*)digits->UncheckedAt(idxBins[idx]);
	    qBin = (Float_t)dig->GetSignal();
	  }
	  y+= qBin * fYSPD[iy];
	  z+= qBin * fZSPD[iz];
	  q+= qBin;	
	}// for idx
	if(ndigits==0) continue;
	y /= q;
	z /= q;
	y -= fHwSPD;
	z -= fHlSPD;
	Float_t hit[5]; //y,z,sigma(y)^2, sigma(z)^2, charge
	hit[0] = -(-y+fYshift[iModule]);
	if(iModule <= fLastSPD1) hit[0] = -hit[0];
	hit[1] = -z+fZshift[iModule];
	hit[2] = fYpitchSPD*fYpitchSPD/12.;
	hit[3] = fZ1pitchSPD*fZ1pitchSPD/12.;
	hit[4] = 1.;
	if(!rawdata) milab[3]=fNdet[iModule];
	Int_t info[3] = {ymax-ymin+1,zmax-zmin+1,fNlayer[iModule]};
	if(!rawdata){
	 AliITSRecPoint cl(iModule,fDetTypeRec->GetITSgeom(),milab,hit,info);
	 fDetTypeRec->AddRecPoint(cl);
	}
        else{
	  Int_t label[4]={milab[0],milab[1],milab[2],milab[3]};
	  new (clusters->AddrAt(nclu)) 
		AliITSRecPoint(iModule,fDetTypeRec->GetITSgeom(),label, hit,info);
	} 
	nclu++;
      }// for iiy
    }// for iiz
  }//end for iBin
  return nclu;
  
}




void AliITSClusterFinderV2SPD::FindClustersSPD(AliITSRawStream* input, 
					TClonesArray** clusters) 
{
  //------------------------------------------------------------
  // Actual SPD cluster finder for raw data
  //------------------------------------------------------------
  Int_t nClustersSPD = 0;
  Int_t kNzBins = fNzSPD + 2;
  Int_t kNyBins = fNySPD + 2;
  Int_t kMaxBin = kNzBins * kNyBins;
  AliBin *binsSPD = new AliBin[kMaxBin];
  AliBin *binsSPDInit = new AliBin[kMaxBin];  
  AliBin* bins = NULL;

  // read raw data input stream
  while (kTRUE) {
    Bool_t next = input->Next();
    if (!next || input->IsNewModule()) {
      Int_t iModule = input->GetPrevModuleID();

      // when all data from a module was read, search for clusters
      if (bins) { 
	clusters[iModule] = new TClonesArray("AliITSRecPoint");
	Int_t nClusters = ClustersSPD(bins,0,clusters[iModule],kMaxBin,kNzBins,iModule,kTRUE);
	nClustersSPD += nClusters;
	bins = NULL;
      }

      if (!next) break;
      bins = binsSPD;
      memcpy(binsSPD,binsSPDInit,sizeof(AliBin)*kMaxBin);
    }

    // fill the current digit into the bins array
    Int_t index = (input->GetCoord2()+1) * kNzBins + (input->GetCoord1()+1);
    bins[index].SetIndex(index);
    bins[index].SetMask(1);
    bins[index].SetQ(1);
  }

  delete [] binsSPDInit;
  delete [] binsSPD;
  
  Info("FindClustersSPD", "found clusters in ITS SPD: %d", nClustersSPD);
}



void AliITSClusterFinderV2SPD::FindClustersSPD(TClonesArray *digits) {
  //------------------------------------------------------------
  // Actual SPD cluster finder
  //------------------------------------------------------------


  Int_t kNzBins = fNzSPD + 2;
  const Int_t kMAXBIN=kNzBins*(fNySPD+2);

  Int_t ndigits=digits->GetEntriesFast();
  AliBin *bins=new AliBin[kMAXBIN];

  Int_t k;
  AliITSdigitSPD *d=0;
  for (k=0; k<ndigits; k++) {
     d=(AliITSdigitSPD*)digits->UncheckedAt(k);
     Int_t i=d->GetCoord2()+1;   //y
     Int_t j=d->GetCoord1()+1;
     Int_t index=i*kNzBins+j;
     bins[index].SetIndex(k);
     bins[index].SetMask(1);
  }
   
  ClustersSPD(bins,digits,0,kMAXBIN,kNzBins,fModule,kFALSE); 
  delete [] bins;
}


