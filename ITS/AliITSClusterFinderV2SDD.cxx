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



#include "AliITSClusterFinderV2SDD.h"
#include "AliITSRecPoint.h"
#include "AliITSDetTypeRec.h"
#include "AliRawReader.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSDetTypeRec.h"
#include "AliITSsegmentationSDD.h"
#include <TClonesArray.h>
#include "AliITSdigitSDD.h"

ClassImp(AliITSClusterFinderV2SDD)

extern AliRun *gAlice;

AliITSClusterFinderV2SDD::AliITSClusterFinderV2SDD(AliITSDetTypeRec* dettyp):AliITSClusterFinderV2(dettyp),
fNySDD(256),
fNzSDD(256),
fYpitchSDD(0.01825),
fZpitchSDD(0.02940),
fHwSDD(3.5085),
fHlSDD(3.7632),
fYoffSDD(0.0425){

  //Default constructor


}
 

void AliITSClusterFinderV2SDD::FindRawClusters(Int_t mod){

  //Find clusters V2
  SetModule(mod);
  FindClustersSDD(fDigits);

}

void AliITSClusterFinderV2SDD::FindClustersSDD(TClonesArray *digits) {
  //------------------------------------------------------------
  // Actual SDD cluster finder
  //------------------------------------------------------------
  Int_t kNzBins = fNzSDD + 2;
  const Int_t kMAXBIN=kNzBins*(fNySDD+2);

  AliBin *bins[2];
  bins[0]=new AliBin[kMAXBIN];
  bins[1]=new AliBin[kMAXBIN];

  AliITSdigitSDD *d=0;
  Int_t i, ndigits=digits->GetEntriesFast();
  for (i=0; i<ndigits; i++) {
     d=(AliITSdigitSDD*)digits->UncheckedAt(i);
     Int_t y=d->GetCoord2()+1;   //y
     Int_t z=d->GetCoord1()+1;   //z
     Int_t q=d->GetSignal();
     if (q<3) continue;

     if (z <= fNzSDD) {
       bins[0][y*kNzBins+z].SetQ(q);
       bins[0][y*kNzBins+z].SetMask(1);
       bins[0][y*kNzBins+z].SetIndex(i);
     } else {
       z-=fNzSDD; 
       bins[1][y*kNzBins+z].SetQ(q);
       bins[1][y*kNzBins+z].SetMask(1);
       bins[1][y*kNzBins+z].SetIndex(i);
     }
  }
  
  FindClustersSDD(bins, kMAXBIN, kNzBins, digits);

  delete[] bins[0];
  delete[] bins[1];

}

void AliITSClusterFinderV2SDD::
FindClustersSDD(AliBin* bins[2], Int_t nMaxBin, Int_t nzBins, 
		TClonesArray *digits, TClonesArray *clusters) {
  //------------------------------------------------------------
  // Actual SDD cluster finder
  //------------------------------------------------------------
  Int_t ncl=0; 
  TClonesArray &cl=*clusters;
  for (Int_t s=0; s<2; s++)
    for (Int_t i=0; i<nMaxBin; i++) {
      if (bins[s][i].IsUsed()) continue;
      Int_t idx[32]; UInt_t msk[32]; Int_t npeaks=0;
      FindPeaks(i, nzBins, bins[s], idx, msk, npeaks);

      if (npeaks>30) continue;
      if (npeaks==0) continue;

      Int_t k,l;
      for (k=0; k<npeaks-1; k++){//mark adjacent peaks
        if (idx[k] < 0) continue; //this peak is already removed
        for (l=k+1; l<npeaks; l++) {
           if (idx[l] < 0) continue; //this peak is already removed
           Int_t ki=idx[k]/nzBins, kj=idx[k] - ki*nzBins;
           Int_t li=idx[l]/nzBins, lj=idx[l] - li*nzBins;
           Int_t di=TMath::Abs(ki - li);
           Int_t dj=TMath::Abs(kj - lj);
           if (di>1 || dj>1) continue;
           if (bins[s][idx[k]].GetQ() > bins[s][idx[l]].GetQ()) {
              msk[l]=msk[k];
              idx[l]*=-1;
           } else {
              msk[k]=msk[l];
              idx[k]*=-1;
              break;
           } 
        }
      }

      for (k=0; k<npeaks; k++) {
        MarkPeak(TMath::Abs(idx[k]), nzBins, bins[s], msk[k]);
      }
        
      for (k=0; k<npeaks; k++) {
         if (idx[k] < 0) continue; //removed peak
         AliITSRecPoint c(fDetTypeRec->GetITSgeom());
         MakeCluster(idx[k], nzBins, bins[s], msk[k], c);
	 //mi change
	 Int_t milab[10];
	 for (Int_t ilab=0;ilab<10;ilab++){
	   milab[ilab]=-2;
	 }
	 Int_t maxi=0,mini=0,maxj=0,minj=0;
	 //AliBin *bmax=&bins[s][idx[k]];
	 //Float_t max = TMath::Max(TMath::Abs(bmax->GetQ())/5.,3.);
	 Float_t max=3;
	 for (Int_t di=-2; di<=2;di++)
	   for (Int_t dj=-3;dj<=3;dj++){
	     Int_t index = idx[k]+di+dj*nzBins;
	     if (index<0) continue;
	     if (index>=nMaxBin) continue;
	     AliBin *b=&bins[s][index];
	     if (TMath::Abs(b->GetQ())>max){
	       if (di>maxi) maxi=di;
	       if (di<mini) mini=di;
	       if (dj>maxj) maxj=dj;
	       if (dj<minj) minj=dj;
	       //
	       if(digits) {
		 if (TMath::Abs(di)<2&&TMath::Abs(dj)<2){
		   AliITSdigitSDD* d=(AliITSdigitSDD*)digits->UncheckedAt(b->GetIndex());
		   for (Int_t itrack=0;itrack<10;itrack++){
		     Int_t track = (d->GetTracks())[itrack];
		     if (track>=0) {
		       AddLabel(milab, track); 
		     }
		   }
		 }
	       }
	     }
	   }
	 
	 /* 
	    Float_t s2 = c.GetSigmaY2()/c.GetQ() - c.GetY()*c.GetY();
	    Float_t w=par->GetPadPitchWidth(sec);
	    c.SetSigmaY2(s2);
	    if (s2 != 0.) {
	    c.SetSigmaY2(c.GetSigmaY2()*0.108);
	    if (sec<par->GetNInnerSector()) c.SetSigmaY2(c.GetSigmaY2()*2.07);
	    }	 
	    s2 = c.GetSigmaZ2()/c.GetQ() - c.GetZ()*c.GetZ();
	    w=par->GetZWidth();
	    c.SetSigmaZ2(s2);
	    
	    if (s2 != 0.) {
	    c.SetSigmaZ2(c.GetSigmaZ2()*0.169);
	    if (sec<par->GetNInnerSector()) c.SetSigmaZ2(c.GetSigmaZ2()*1.77);
	    }
	 */

         c.SetSigmaY2(0.0030*0.0030);
         c.SetSigmaZ2(0.0020*0.0020);
         c.SetDetectorIndex(fNdet[fModule]);

         Float_t y=c.GetY(),z=c.GetZ(), q=c.GetQ();
         y/=q; z/=q;
	 //
	 //Float_t s2 = c.GetSigmaY2()/c.GetQ() - y*y;
	 // c.SetSigmaY2(s2);
	 //s2 = c.GetSigmaZ2()/c.GetQ() - z*z;
         //c.SetSigmaZ2(s2);
	 //
         y=(y-0.5)*fYpitchSDD;
         y-=fHwSDD;
         y-=fYoffSDD;  //delay ?
         if (s) y=-y;

         z=(z-0.5)*fZpitchSDD;
         z-=fHlSDD;

         y=-(-y+fYshift[fModule]);
         z=  -z+fZshift[fModule];
	 //      c.SetY(y);
	 //  c.SetZ(z);
	 CorrectPosition(z,y);
	 c.SetYZ(fModule,y,z);
	 c.SetNy(maxj-minj+1);
	 c.SetNz(maxi-mini+1);
	 c.SetType(npeaks);
         c.SetQ(q/12.7);  //to be consistent with the SSD charges

         if (c.GetQ() < 20.) continue; //noise cluster
	 
	 if (digits) {	  
	   //	   AliBin *b=&bins[s][idx[k]];
	   //	   AliITSdigitSDD* d=(AliITSdigitSDD*)digits->UncheckedAt(b->GetIndex());
	   {
	     //Int_t lab[3];
	     //lab[0]=(d->GetTracks())[0];
	     //lab[1]=(d->GetTracks())[1];
	     //lab[2]=(d->GetTracks())[2];
	     //CheckLabels(lab);
	     CheckLabels2(milab); 
	     c.SetLabel(milab[0],0);
	     c.SetLabel(milab[1],1);
	     c.SetLabel(milab[2],2);
	     c.SetLayer(fNlayer[fModule]);
	   }
	 }
	 if(clusters) new (cl[ncl]) AliITSRecPoint(c); 
	 else {
	   fDetTypeRec->AddRecPoint(c);
	 }
	 ncl++;
      }
    }
}



void AliITSClusterFinderV2SDD::RawdataToClusters(AliRawReader* rawReader,TClonesArray** clusters){
    //------------------------------------------------------------
  // This function creates ITS clusters from raw data
  //------------------------------------------------------------
  rawReader->Reset();
  AliITSRawStreamSDD inputSDD(rawReader);
  FindClustersSDD(&inputSDD,clusters);

}

void AliITSClusterFinderV2SDD::FindClustersSDD(AliITSRawStream* input, 
					TClonesArray** clusters) 
{
  //------------------------------------------------------------
  // Actual SDD cluster finder for raw data
  //------------------------------------------------------------
  Int_t nClustersSDD = 0;
  Int_t kNzBins = fNzSDD + 2;
  Int_t kMaxBin = kNzBins * (fNySDD+2);
  AliBin *binsSDDInit = new AliBin[kMaxBin];
  AliBin *binsSDD1 = new AliBin[kMaxBin];
  AliBin *binsSDD2 = new AliBin[kMaxBin];
  AliBin* bins[2] = {NULL, NULL};

  // read raw data input stream
  while (kTRUE) {
    Bool_t next = input->Next();
    if (!next || input->IsNewModule()) {
      Int_t iModule = input->GetPrevModuleID();

      // when all data from a module was read, search for clusters
      if (bins[0]) { 
	clusters[iModule] = new TClonesArray("AliITSRecPoint");
	fModule = iModule;
	FindClustersSDD(bins, kMaxBin, kNzBins, NULL, clusters[iModule]);
	Int_t nClusters = clusters[iModule]->GetEntriesFast();
	nClustersSDD += nClusters;
	bins[0] = bins[1] = NULL;
	
      }

      if (!next) break;
      bins[0]=binsSDD1;
      bins[1]=binsSDD2;
      memcpy(binsSDD1,binsSDDInit,sizeof(AliBin)*kMaxBin);
      memcpy(binsSDD2,binsSDDInit,sizeof(AliBin)*kMaxBin);

    }

    // fill the current digit into the bins array
    if(input->GetSignal()>=3) {
      Int_t iz = input->GetCoord1()+1;
      Int_t side = ((iz <= fNzSDD) ? 0 : 1);
      iz -= side*fNzSDD;
      Int_t index = (input->GetCoord2()+1) * kNzBins + iz;
      bins[side][index].SetQ(input->GetSignal());
      bins[side][index].SetMask(1);
      bins[side][index].SetIndex(index);
    }
  }
  delete[] binsSDD1;
  delete[] binsSDD2;
  delete[] binsSDDInit;

  Info("FindClustersSDD", "found clusters in ITS SDD: %d", nClustersSDD);

}


//_________________________________________________________________________
void AliITSClusterFinderV2SDD::CorrectPosition(Float_t &z, Float_t&y){

  //correction of coordinates using the maps stored in the DB

  AliITSCalibrationSDD* cal = (AliITSCalibrationSDD*)GetResp(fModule);
  static const Int_t nbint = cal->GetMapTimeNBin();
  static const Int_t nbina = cal->Chips()*cal->Channels();
  Float_t stepa = (GetSeg()->Dpz(0))/10000.; //anode pitch in cm
  Float_t stept = (GetSeg()->Dx()/cal->GetMapTimeNBin()/2.)/10.;
  Float_t xdet,zdet;
  fDetTypeRec->GetITSgeom()->TrackingV2ToDetL(fModule,y,z,xdet,zdet);

  Int_t bint = TMath::Abs((Int_t)(xdet/stept));
  if(xdet>=0) bint+=(Int_t)(nbint/2.);
  if(bint>nbint) AliError("Wrong bin number!");

  Int_t bina = TMath::Abs((Int_t)(zdet/stepa));
  if(zdet>=0) bina+=(Int_t)(nbina/2.);
  if(bina>nbina) AliError("Wrong bin number!");

  Float_t devz = cal->GetMapACell(bina,bint)/10000.;
  Float_t devx = cal->GetMapTCell(bina,bint)/10000.;
  zdet+=devz;
  xdet+=devx;
  fDetTypeRec->GetITSgeom()->DetLToTrackingV2(fModule,xdet,zdet,y,z);


}
