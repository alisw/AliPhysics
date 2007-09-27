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

/* $Id$*/

////////////////////////////////////////////////////////////////////////////
//            Implementation of the ITS clusterer V2 class                //
//                                                                        //
//          Origin: Iouri Belikov, CERN, Jouri.Belikov@cern.ch            //
//                                                                        //
///////////////////////////////////////////////////////////////////////////



#include <TClonesArray.h>
#include "AliITSClusterFinderV2SDD.h"
#include "AliITSRecPoint.h"
#include "AliITSDetTypeRec.h"
#include "AliRawReader.h"
#include "AliITSRawStreamSDD.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSDetTypeRec.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSdigitSDD.h"
#include "AliITSgeomTGeo.h"
ClassImp(AliITSClusterFinderV2SDD)

AliITSClusterFinderV2SDD::AliITSClusterFinderV2SDD(AliITSDetTypeRec* dettyp):AliITSClusterFinderV2(dettyp),
fNySDD(256),
fNzSDD(256),
fZpitchSDD(0.02940),
fHwSDD(3.5085),
fHlSDD(3.7632),
fTimeOffsetSDD(0.){

  //Default constructor

  SetTimeOffset();
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
  AliITSCalibrationSDD* cal = (AliITSCalibrationSDD*)GetResp(fModule);
  AliITSresponseSDD* res  = (AliITSresponseSDD*)cal->GetResponse();
  const char *option=res->ZeroSuppOption();
  AliITSdigitSDD *d=0;
  Int_t i, ndigits=digits->GetEntriesFast();
  for (i=0; i<ndigits; i++) {
     d=(AliITSdigitSDD*)digits->UncheckedAt(i);
     Int_t y=d->GetCoord2()+1;   //y
     Int_t z=d->GetCoord1()+1;   //z
     Int_t q=d->GetSignal();
     if(!((strstr(option,"1D")) || (strstr(option,"2D")))){
       Float_t baseline = cal->GetBaseline(d->GetCoord1());
       if(q>baseline) q-=(Int_t)baseline;
       else q=0;
     }
     if(q<cal->GetThresholdAnode(d->GetCoord1())) continue;
     //if (q<3) continue;

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

  const TGeoHMatrix *mT2L=AliITSgeomTGeo::GetTracking2LocalMatrix(fModule);
  AliITSCalibrationSDD* cal = (AliITSCalibrationSDD*)GetResp(fModule);
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
         AliITSRecPoint c;
         MakeCluster(idx[k], nzBins, bins[s], msk[k], c);
	 //mi change
	 Int_t milab[10];
	 for (Int_t ilab=0;ilab<10;ilab++){
	   milab[ilab]=-2;
	 }
	 Int_t maxi=0,mini=0,maxj=0,minj=0;
	 //AliBin *bmax=&bins[s][idx[k]];
	 //Float_t max = TMath::Max(TMath::Abs(bmax->GetQ())/5.,3.);
    
	 for (Int_t di=-2; di<=2;di++){
	   for (Int_t dj=-3;dj<=3;dj++){
	     Int_t index = idx[k]+di+dj*nzBins;
	     if (index<0) continue;
	     if (index>=nMaxBin) continue;
	     AliBin *b=&bins[s][index];
	     Int_t nAnode=index%nzBins-1;
	     Int_t adcSignal=b->GetQ();
	     if(adcSignal>cal->GetThresholdAnode(nAnode)){
	       if (di>maxi) maxi=di;
	       if (di<mini) mini=di;
	       if (dj>maxj) maxj=dj;
	       if (dj<minj) minj=dj;
	     }
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


         Float_t y=c.GetY(),z=c.GetZ(), q=c.GetQ();
         y/=q; z/=q;


	 Float_t yyyy = y;

         z=(z-0.5)*fZpitchSDD;
         z-=fHlSDD;
	 Float_t zdet = z;
	 Float_t timebin = GetSeg()->Dpx(0);
	 Float_t driftTime = (yyyy-0.5)*timebin - fTimeOffsetSDD;
	 Float_t xdet = cal->GetDriftPath(driftTime,0);
	 xdet=xdet/10000.-fHwSDD;
	 if (s) xdet=-xdet;
	 
	 CorrectPosition(zdet,xdet);

         {
         Double_t loc[3]={xdet,0.,zdet},trk[3]={0.,0.,0.};
         mT2L->MasterToLocal(loc,trk);
         y=trk[1];
         z=trk[2]; 
         }
         q/=5.243;  //to have MPV 1 MIP = 86.4 KeV --> this must go to calibr.
         Float_t hit[5] = {y, z, 0.0030*0.0030, 0.0020*0.0020, q};
         Int_t  info[3] = {maxj-minj+1, maxi-mini+1, fNlayer[fModule]};
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
	   }
	 }
         milab[3]=fNdet[fModule];

         AliITSRecPoint cc(milab,hit,info);
	 cc.SetType(npeaks);

	 if(clusters) new (cl[ncl]) AliITSRecPoint(cc); 
	 else {
	   fDetTypeRec->AddRecPoint(cc);
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
  AliBin *bins[2];
  bins[0]=new AliBin[kMaxBin];
  bins[1]=new AliBin[kMaxBin];

  // read raw data input stream
  while (kTRUE) {
    Bool_t next = input->Next();
    if (!next || input->IsCompletedModule()) {
      // when all data from a module was read, search for clusters
      Int_t iModule = input->GetModuleID();
      clusters[iModule] = new TClonesArray("AliITSRecPoint");
      fModule = iModule;
      FindClustersSDD(bins, kMaxBin, kNzBins, NULL, clusters[iModule]);
      Int_t nClusters = clusters[iModule]->GetEntriesFast();
      nClustersSDD += nClusters;
      for(Int_t iBin=0;iBin<kMaxBin; iBin++){
	bins[0][iBin].Reset();
	bins[1][iBin].Reset();
      }
      if (!next) break;
    }else{
    // fill the current digit into the bins array
      AliITSCalibrationSDD* cal = (AliITSCalibrationSDD*)GetResp(input->GetModuleID());    
      AliITSresponseSDD* res  = (AliITSresponseSDD*)cal->GetResponse();
      const char *option=res->ZeroSuppOption();
      Int_t q=input->GetSignal();
      if(!((strstr(option,"1D")) || (strstr(option,"2D")))){
	Float_t baseline = cal->GetBaseline(input->GetCoord1());
	if(q>baseline) q-=(Int_t)baseline;
	else q=0;
      }
      if(q>=cal->GetThresholdAnode(input->GetCoord1())) {
	Int_t iz = input->GetCoord1()+1;
	//Int_t side = ((iz <= fNzSDD) ? 0 : 1);
	Int_t side = ((AliITSRawStreamSDD*)input)->GetChannel();
	//  iz -= side*fNzSDD;
	Int_t index = (input->GetCoord2()+1) * kNzBins + iz;
	bins[side][index].SetQ(q);
	bins[side][index].SetMask(1);
	bins[side][index].SetIndex(index);
      }
    }
  }
  delete [] bins[0];
  delete [] bins[1];

  Info("FindClustersSDD", "found clusters in ITS SDD: %d", nClustersSDD);
}


//_________________________________________________________________________
void AliITSClusterFinderV2SDD::CorrectPosition(Float_t &z, Float_t&y){

  //correction of coordinates using the maps stored in the DB

  AliITSCalibrationSDD* cal = (AliITSCalibrationSDD*)GetResp(fModule);
  static const Int_t knbint = cal->GetMapTimeNBin();
  static const Int_t knbina = cal->Chips()*cal->Channels();
  Float_t stepa = (GetSeg()->Dpz(0))/10000.; //anode pitch in cm
  Float_t stept = (GetSeg()->Dx()/cal->GetMapTimeNBin()/2.)/10.;
  
  Int_t bint = TMath::Abs((Int_t)(y/stept));
  if(y>=0) bint+=(Int_t)(knbint/2.);
  if(bint>knbint) AliError("Wrong bin number!");

  Int_t bina = TMath::Abs((Int_t)(z/stepa));
  if(z>=0) bina+=(Int_t)(knbina/2.);
  if(bina>knbina) AliError("Wrong bin number!");

  Float_t devz = cal->GetMapACell(bina,bint)/10000.;
  Float_t devx = cal->GetMapTCell(bina,bint)/10000.;
  z+=devz;
  y+=devx;


}
