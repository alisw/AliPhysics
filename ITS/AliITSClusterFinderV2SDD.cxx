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
#include "AliITSRawStreamSDDCompressed.h"
#include "AliITSCalibrationSDD.h"
#include "AliITSresponseSDD.h"
#include "AliITSDetTypeRec.h"
#include "AliITSReconstructor.h"
#include "AliITSsegmentationSDD.h"
#include "AliITSdigitSDD.h"
#include "AliITSgeomTGeo.h"

ClassImp(AliITSClusterFinderV2SDD)

AliITSClusterFinderV2SDD::AliITSClusterFinderV2SDD(AliITSDetTypeRec* dettyp):AliITSClusterFinder(dettyp)
{

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
  Int_t nAnodes = GetSeg()->NpzHalf();
  Int_t nzBins = nAnodes+2;
  Int_t nTimeBins = GetSeg()->Npx();
  Int_t nxBins = nTimeBins+2;
  const Int_t kMaxBin=nzBins*nxBins;

  AliBin *bins[2];
  bins[0]=new AliBin[kMaxBin];
  bins[1]=new AliBin[kMaxBin];
  AliITSCalibrationSDD* cal = (AliITSCalibrationSDD*)GetResp(fModule);
  if(cal==0){
    AliError(Form("Calibration object not present for SDD module %d\n",fModule));
    return;
  }

  AliITSdigitSDD *d=0;
  Int_t i, ndigits=digits->GetEntriesFast();
  for (i=0; i<ndigits; i++) {
     d=(AliITSdigitSDD*)digits->UncheckedAt(i);
     Int_t ian=d->GetCoord1();
     Int_t itb=d->GetCoord2();
     Int_t iSide=0;
     if (ian >= nAnodes) iSide=1;    
     Float_t gain=cal->GetChannelGain(ian)/fDetTypeRec->GetAverageGainSDD();
     Float_t charge=d->GetSignal(); // returns expanded signal 
                                    // (10 bit, low threshold already added)
     Float_t baseline = cal->GetBaseline(ian);
     if(charge>baseline) charge-=baseline;
     else charge=0;

     if(gain>0.){ // Bad channels have gain=0.
       charge/=gain;
       if(charge<cal->GetThresholdAnode(ian)) continue;
       Int_t q=(Int_t)(charge+0.5);
       Int_t y=itb+1;   
       Int_t z=ian+1;   
       if (z <= nAnodes){
	 bins[0][y*nzBins+z].SetQ(q);
	 bins[0][y*nzBins+z].SetMask(1);
	 bins[0][y*nzBins+z].SetIndex(i);
       } else {
	 z-=nAnodes;
	 bins[1][y*nzBins+z].SetQ(q);
	 bins[1][y*nzBins+z].SetMask(1);
	 bins[1][y*nzBins+z].SetIndex(i);
       }
     }
  }
  
  FindClustersSDD(bins, kMaxBin, nzBins, digits);

  delete[] bins[0];
  delete[] bins[1];

}

void AliITSClusterFinderV2SDD::
FindClustersSDD(AliBin* bins[2], Int_t nMaxBin, Int_t nzBins, 
		TClonesArray *digits, TClonesArray *clusters, Int_t jitter) {
  //------------------------------------------------------------
  // Actual SDD cluster finder
  //------------------------------------------------------------

  static AliITSRecoParam *repa = NULL;
  if(!repa){
    repa = (AliITSRecoParam*) AliITSReconstructor::GetRecoParam();
    if(!repa){
      repa = AliITSRecoParam::GetHighFluxParam();
      AliWarning("Using default AliITSRecoParam class");
    }
  }
  const TGeoHMatrix *mT2L=AliITSgeomTGeo::GetTracking2LocalMatrix(fModule);
  AliITSCalibrationSDD* cal = (AliITSCalibrationSDD*)GetResp(fModule);
  if(cal==0){
    AliError(Form("Calibration object not present for SDD module %d\n",fModule));
    return;
  }
  Int_t ncl=0; 
  TClonesArray &cl=*clusters;
  for (Int_t s=0; s<2; s++)
    for (Int_t i=0; i<nMaxBin; i++) {
      if(NoiseSuppress(i,s,nzBins,bins[s],cal)) continue;
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
	if(repa->GetUseUnfoldingInClusterFinderSDD()==kFALSE) msk[k]=msk[0];
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

	for (Int_t di=-5; di<=5;di++){
	  for (Int_t dj=-10;dj<=10;dj++){
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
	Int_t clSizAnode=maxi-mini+1;
	Int_t clSizTb=maxj-minj+1;
	if(repa->GetUseSDDClusterSizeSelection()){
	  if(clSizTb==1) continue; // cut common mode noise spikes
	  if(clSizAnode>3)  continue; // cut common mode noise spikes
	  if(clSizTb>10)  continue; // cut clusters on noisy anodes
	}

	AliITSresponseSDD* rsdd = fDetTypeRec->GetResponseSDD();
	Float_t y=c.GetY(),z=c.GetZ(), q=c.GetQ();
	y/=q; z/=q;
	Float_t zAnode=z-0.5;  // to have anode in range 0.-255. and centered on the mid of the pitch
	Float_t timebin=y-0.5;  // to have time bin in range 0.-255. amd centered on the mid of the bin
	if(s==1) zAnode += GetSeg()->NpzHalf();  // right side has anodes from 256. to 511.
	Float_t zdet = GetSeg()->GetLocalZFromAnode(zAnode);
	Float_t driftTimeUncorr = GetSeg()->GetDriftTimeFromTb(timebin)+jitter*rsdd->GetCarlosRXClockPeriod();
	Float_t driftTime=driftTimeUncorr-rsdd->GetTimeZero(fModule);
	Float_t driftPathMicron = cal->GetDriftPath(driftTime,zAnode);
	const Double_t kMicronTocm = 1.0e-4; 
	Float_t xdet=(driftPathMicron-GetSeg()->Dx())*kMicronTocm; // xdet is negative
	if (s==0) xdet=-xdet; // left side has positive local x
	
	Float_t corrx=0, corrz=0;
	cal->GetCorrections(zdet,xdet,corrz,corrx,GetSeg());
	zdet+=corrz;
	xdet+=corrx;

	Double_t loc[3]={xdet,0.,zdet},trk[3]={0.,0.,0.};
	mT2L->MasterToLocal(loc,trk);
	y=trk[1];
	z=trk[2]; 

	q/=rsdd->GetADC2keV();  //to have MPV 1 MIP = 86.4 KeV
	if(cal-> IsAMAt20MHz()) q*=2.; // account for 1/2 sampling freq.
	if(q<repa->GetMinClusterChargeSDD()) continue; // remove noise clusters

	Float_t hit[5] = {y, z, 0.0030*0.0030, 0.0020*0.0020, q};
	Int_t  info[3] = {clSizTb, clSizAnode, fNlayer[fModule]};
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
	cc.SetDriftTime(driftTimeUncorr);
	if(clusters) new (cl[ncl]) AliITSRecPoint(cc); 
	else {
	  fDetTypeRec->AddRecPoint(cc);
	}
	ncl++;
      }
    }
  
}
//______________________________________________________________________
void AliITSClusterFinderV2SDD::RawdataToClusters(AliRawReader* rawReader,TClonesArray** clusters){
    //------------------------------------------------------------
  // This function creates ITS clusters from raw data
  //------------------------------------------------------------
  rawReader->Reset();
  AliITSRawStream* inputSDD;
  if(fDetTypeRec->IsHLTmodeC()==kTRUE){
    inputSDD=new AliITSRawStreamSDDCompressed(rawReader);
  }else{
    inputSDD=new AliITSRawStreamSDD(rawReader);
  }

  AliITSDDLModuleMapSDD *ddlmap=(AliITSDDLModuleMapSDD*)fDetTypeRec->GetDDLModuleMapSDD();
  inputSDD->SetDDLModuleMap(ddlmap);
  for(Int_t iddl=0; iddl<AliITSDDLModuleMapSDD::GetNDDLs(); iddl++){
    for(Int_t icar=0; icar<AliITSDDLModuleMapSDD::GetNModPerDDL();icar++){
      Int_t iMod=ddlmap->GetModuleNumber(iddl,icar);
      if(iMod==-1) continue;
      AliITSCalibrationSDD* cal = (AliITSCalibrationSDD*)GetResp(iMod);
      if(cal==0){
	AliError(Form("Calibration object not present for SDD module %d\n",iMod));
	continue;
      }
      Bool_t isZeroSupp=cal->GetZeroSupp();
      if(isZeroSupp){ 
	for(Int_t iSid=0; iSid<2; iSid++) inputSDD->SetZeroSuppLowThreshold(iMod-240,iSid,cal->GetZSLowThreshold(iSid));
      }else{
	for(Int_t iSid=0; iSid<2; iSid++) inputSDD->SetZeroSuppLowThreshold(iMod-240,iSid,0);
      }
    }
  }
  FindClustersSDD(inputSDD,clusters);
  delete inputSDD;
}

void AliITSClusterFinderV2SDD::FindClustersSDD(AliITSRawStream* input, 
					TClonesArray** clusters) 
{
  //------------------------------------------------------------
  // Actual SDD cluster finder for raw data
  //------------------------------------------------------------
  Int_t nClustersSDD = 0;
  Int_t nAnodes = GetSeg()->NpzHalf();
  Int_t nzBins = nAnodes+2;
  Int_t nTimeBins = GetSeg()->Npx();
  Int_t nxBins = nTimeBins+2;
  const Int_t kMaxBin=nzBins*nxBins;
  AliBin *bins[2];
  AliBin *ddlbins[kHybridsPerDDL]; // 12 modules (=24 hybrids) of 1 DDL read "in parallel"
  for(Int_t iHyb=0;iHyb<kHybridsPerDDL;iHyb++) ddlbins[iHyb]=new AliBin[kMaxBin];
  Int_t vectModId[kModulesPerDDL];
  for(Int_t iMod=0; iMod<kModulesPerDDL; iMod++) vectModId[iMod]=-1;

  // read raw data input stream
  while (input->Next()) {
    Int_t iModule = input->GetModuleID();
    if(iModule<0){
      AliWarning(Form("Invalid SDD module number %d\n", iModule));
      continue;
    }
    Int_t iCarlos =input->GetCarlosId();
    Int_t iSide = input->GetChannel();
    Int_t iHybrid=iCarlos*2+iSide;

    if (input->IsCompletedModule()) {
      // store the module number
      vectModId[iCarlos]=iModule;
    }
    else if (input->IsCompletedDDL()) {
      // when all data from a DDL was read, search for clusters
      Int_t jitter=input->GetJitter();
      for(Int_t iMod=0; iMod<kModulesPerDDL; iMod++){
	if(vectModId[iMod]>=0){
	  fModule = vectModId[iMod];
	  clusters[fModule] = new TClonesArray("AliITSRecPoint");
	  bins[0]=ddlbins[iMod*2];   // first hybrid of the module
	  bins[1]=ddlbins[iMod*2+1]; // second hybrid of the module
	  FindClustersSDD(bins, kMaxBin, nzBins, NULL, clusters[fModule],jitter);
	  Int_t nClusters = clusters[fModule]->GetEntriesFast();
	  nClustersSDD += nClusters;
	  vectModId[iMod]=-1;
	}
	for(Int_t iBin=0;iBin<kMaxBin; iBin++){
	  ddlbins[iMod*2][iBin].Reset();
	  ddlbins[iMod*2+1][iBin].Reset();
	}
      }
    }else{
    // fill the current digit into the bins array
      if(iHybrid<0 || iHybrid>=kHybridsPerDDL){ 
	AliWarning(Form("Invalid SDD hybrid number %d on module %d\n", iHybrid,iModule));
	continue;
      }
      AliITSCalibrationSDD* cal = (AliITSCalibrationSDD*)GetResp(iModule);    
      if(cal==0){
	AliError(Form("Calibration object not present for SDD module %d\n",iModule));
	continue;
      }
      Float_t charge=input->GetSignal();
      Int_t chan=input->GetCoord1()+nAnodes*iSide;
      Float_t gain=cal->GetChannelGain(chan)/fDetTypeRec->GetAverageGainSDD();;
      Float_t baseline = cal->GetBaseline(chan);
      if(charge>baseline) charge-=baseline;
      else charge=0;
      if(gain>0.){ // Bad channels have gain=0
	charge/=gain;
	if(charge>=cal->GetThresholdAnode(chan)) {
	  Int_t q=(Int_t)(charge+0.5);
	  Int_t iz = input->GetCoord1();
	  Int_t itb = input->GetCoord2();
	  Int_t index = (itb+1) * nzBins + (iz+1);
	  if(index<kMaxBin){
	    ddlbins[iHybrid][index].SetQ(q);
	    ddlbins[iHybrid][index].SetMask(1);
	    ddlbins[iHybrid][index].SetIndex(index);
	  }else{
	    AliWarning(Form("Invalid SDD cell: Anode=%d   TimeBin=%d",iz,itb));	  
	  }
	}
      }
    }
  }
  for(Int_t iHyb=0;iHyb<kHybridsPerDDL;iHyb++) delete [] ddlbins[iHyb];
  Info("FindClustersSDD", "found clusters in ITS SDD: %d", nClustersSDD);
}

//______________________________________________________________________
Bool_t AliITSClusterFinderV2SDD::NoiseSuppress(Int_t k, Int_t sid,Int_t nzBins, AliBin* bins, AliITSCalibrationSDD* cal) const {
  // applies zero suppression using the measured noise of each anode
  // threshold values from ALICE-INT-1999-28 V10
  // returns kTRUE if the digit should eb noise suppressed, kFALSE if it should be kept
  Float_t xfactL=2.2; 
  Float_t xfactH=4.0;
  //

  Int_t iAn=(k%nzBins)-1;
  if(iAn<0 || iAn>255) return kTRUE;
  if(sid==1) iAn+=256;
  Int_t nLow=0, nHigh=0;
  Float_t noise=cal->GetNoiseAfterElectronics(iAn);
  Float_t noisem1=noise;
  if(iAn>1) noisem1=cal->GetNoiseAfterElectronics(iAn-1);
  Float_t noisep1=noise;
  if(iAn<511) noisep1=cal->GetNoiseAfterElectronics(iAn+1);
  Float_t tL=noise*xfactL;
  Float_t tH=noise*xfactH;
  Float_t tLp1=noisep1*xfactL;
  Float_t tHp1=noisep1*xfactH;
  Float_t tLm1=noisem1*xfactL;
  Float_t tHm1=noisem1*xfactH;
  Int_t cC=bins[k].GetQ();
  if(cC<=tL){
    bins[k].SetQ(0);
    bins[k].SetMask(0xFFFFFFFE);
    return kTRUE;;
  }
  nLow++; // cC is greater than tL
  if(cC>tH) nHigh++;
  Int_t sS=bins[k-1].GetQ();
  if(sS>tLm1) nLow++;
  if(sS>tHm1) nHigh++;
  Int_t nN=bins[k+1].GetQ();
  if(nN>tLp1) nLow++;
  if(nN>tHp1) nHigh++;
  Int_t eE=bins[k-nzBins].GetQ();
  if(eE>tL) nLow++;
  if(eE>tH) nHigh++;
  Int_t wW=bins[k+nzBins].GetQ();
  if(wW>tL) nLow++;
  if(wW>tH) nHigh++;
  if(nLow<3 || nHigh<1) return kTRUE;
  else return kFALSE;
}



