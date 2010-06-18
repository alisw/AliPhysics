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
#include <TBits.h>
#include "AliITSClusterFinderV2SDD.h"
#include "AliITSRecPoint.h"
#include "AliITSRecPointContainer.h"
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

AliITSClusterFinderV2SDD::AliITSClusterFinderV2SDD(AliITSDetTypeRec* dettyp):AliITSClusterFinder(dettyp),
  fNAnodes(0),
  fNTimeBins(0),
  fNZbins(0),
  fNXbins(0),
  fCutOnPeakLoose(0.),
  fCutOnPeakTight(0.),
  fMaxDrTimeForTightCut(0.)
{
  //Default constructor

  fNAnodes = GetSeg()->NpzHalf();
  fNZbins = fNAnodes+2;
  fNTimeBins = GetSeg()->Npx();
  fNXbins = fNTimeBins+2;
  AliDebug(2,Form("Cells in SDD cluster finder: Andoes=%d  TimeBins=%d",fNAnodes,fNTimeBins));
  const Int_t kMaxBin=fNZbins*fNXbins;
  for(Int_t iHyb=0;iHyb<kHybridsPerDDL;iHyb++){
   fDDLBins[iHyb]=new AliBin[kMaxBin];
  }
  SetPeakSelection(15.,30.,2000.);
}
 

//______________________________________________________________________
AliITSClusterFinderV2SDD::~AliITSClusterFinderV2SDD()
{
  //Destructor
  for(Int_t iHyb=0;iHyb<kHybridsPerDDL;iHyb++){ 
    delete [] fDDLBins[iHyb];
  }
}

//______________________________________________________________________
void AliITSClusterFinderV2SDD::FindRawClusters(Int_t mod){

  //Find clusters V2
  SetModule(mod);
  FindClustersSDD(fDigits);

}

//______________________________________________________________________
void AliITSClusterFinderV2SDD::FindClustersSDD(TClonesArray *digits) {
  //------------------------------------------------------------
  // Actual SDD cluster finder
  //------------------------------------------------------------

  const Int_t kMaxBin=fNZbins*fNXbins;
  AliBin *bins[2];
  bins[0]=new AliBin[kMaxBin];
  bins[1]=new AliBin[kMaxBin];
  TBits *anodeFired[2];
  anodeFired[0]=new TBits(fNAnodes);
  anodeFired[1]=new TBits(fNAnodes);
  anodeFired[0]->ResetAllBits();
  anodeFired[1]->ResetAllBits();
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
     if (ian >= fNAnodes) iSide=1;    
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
       if (z <= fNAnodes){
	 bins[0][y*fNZbins+z].SetQ(q);
	 bins[0][y*fNZbins+z].SetMask(1);
	 bins[0][y*fNZbins+z].SetIndex(i);
	 anodeFired[0]->SetBitNumber(ian);
       } else {
	 z-=fNAnodes;
	 bins[1][y*fNZbins+z].SetQ(q);
	 bins[1][y*fNZbins+z].SetMask(1);
	 bins[1][y*fNZbins+z].SetIndex(i);
	 anodeFired[1]->SetBitNumber(ian-fNAnodes);
       }
     }
  }
  
  FindClustersSDD(bins, anodeFired, digits);

  delete[] bins[0];
  delete[] bins[1];
  delete anodeFired[0];
  delete anodeFired[1];

}

//______________________________________________________________________
void AliITSClusterFinderV2SDD::
FindClustersSDD(AliBin* bins[2], TBits* anodeFired[2],  
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
  const Int_t kMaxBin=fNZbins*fNXbins;
  Int_t ncl=0; 
  TClonesArray &cl=*clusters;
  for (Int_t s=0; s<2; s++){
    for(Int_t iAnode=0; iAnode<GetSeg()->NpzHalf(); iAnode++){
      if(anodeFired[s]->TestBitNumber(iAnode)==kFALSE) continue;
      for(Int_t iTimeBin=0; iTimeBin<GetSeg()->Npx(); iTimeBin++){
	Int_t index=(iTimeBin+1)*fNZbins+(iAnode+1);
	if (bins[s][index].IsUsed()) continue;
	if(NoiseSuppress(index,s,bins[s],cal)) continue;
	Int_t idx[32]; UInt_t msk[32]; Int_t npeaks=0;
	FindPeaks(index, fNZbins, bins[s], idx, msk, npeaks);

	if (npeaks>30) continue;
	if (npeaks==0) continue;

	Int_t k,l;
	Int_t nClust;
	if(repa->GetUseUnfoldingInClusterFinderSDD()){
	  for (k=0; k<npeaks-1; k++){//mark adjacent peaks	    
	    if (idx[k] < 0) continue; //this peak is already removed
	    for (l=k+1; l<npeaks; l++) {
	      if (idx[l] < 0) continue; //this peak is already removed
	      Int_t ki=idx[k]/fNZbins, kj=idx[k] - ki*fNZbins;
	      Int_t li=idx[l]/fNZbins, lj=idx[l] - li*fNZbins;
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
	  nClust=npeaks;
	}else{
	  for (k=1; k<npeaks; k++) msk[k]=msk[0];
	  nClust=1;
	}
	Float_t maxADC=0;
	for (k=0; k<npeaks; k++) {
	  if(idx[k]>0. && bins[s][idx[k]].GetQ() > maxADC) maxADC=bins[s][idx[k]].GetQ();
	  MarkPeak(TMath::Abs(idx[k]), fNZbins, bins[s], msk[k]);
	}
	if(maxADC<fCutOnPeakLoose) continue;

	for (k=0; k<nClust; k++) {
	  if (idx[k] < 0) continue; //removed peak
	  AliITSRecPoint c;
	  MakeCluster(idx[k], fNZbins, bins[s], msk[k], c);
	  //mi change
	  Int_t milab[10];
	  for (Int_t ilab=0;ilab<10;ilab++){
	    milab[ilab]=-2;
	  }
	  
	  if(digits) {
	    for (Int_t di=-2; di<=2;di++){
	      for (Int_t dj=-2;dj<=2;dj++){
		index = idx[k]+di+dj*fNZbins;
		if (index<0) continue;
		if (index>=kMaxBin) continue;
		AliBin *b=&bins[s][index];
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
	  

	  Int_t clSizAnode=fZmax-fZmin+1;
	  Int_t clSizTb=fXmax-fXmin+1;	  
	  if(repa->GetUseSDDClusterSizeSelection()){
	    if(clSizTb==1) continue; // cut common mode noise spikes
	    if(clSizAnode>5)  continue; // cut common mode noise spikes
	    if(clSizTb>10)  continue; // cut clusters on noisy anodes
	    if(cal-> IsAMAt20MHz() && clSizTb>8)  continue; // cut clusters on noisy anodes
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
	  if(driftTime<fMaxDrTimeForTightCut && maxADC<fCutOnPeakTight) continue;

	  Float_t driftSpeed = cal->GetDriftSpeedAtAnode(zAnode) + rsdd->GetDeltaVDrift(fModule,zAnode>255);
	  Float_t driftPathMicron = driftTime*driftSpeed;
	  const Double_t kMicronTocm = 1.0e-4; 
	  Float_t xdet=(driftPathMicron-GetSeg()->Dx())*kMicronTocm; // xdet is negative
	  if (s==0) xdet=-xdet; // left side has positive local x
	  
	  if(repa->GetUseSDDCorrectionMaps()){
	    Float_t corrx=0, corrz=0;
	    cal->GetCorrections(zdet,xdet,corrz,corrx,GetSeg());
	    zdet+=corrz;
	    xdet+=corrx;
	  }
	  
	  Double_t loc[3]={xdet,0.,zdet},trk[3]={0.,0.,0.};
	  mT2L->MasterToLocal(loc,trk);
	  y=trk[1];
	  z=trk[2]; 
	  
	  q/=rsdd->GetADCtokeV(fModule);
	  q+=(driftTime*rsdd->GetChargevsTime()); // correction for zero supp.
	  if(cal-> IsAMAt20MHz()) q*=2.; // account for 1/2 sampling freq.
	  if(q<repa->GetMinClusterChargeSDD()) continue; // remove noise clusters
	  
	  Float_t hit[5] = {y, z, 0.0030*0.0030, 0.0020*0.0020, q};
	  Int_t  info[3] = {clSizTb, clSizAnode, fNlayer[fModule]};
	  if (digits) CheckLabels2(milab);
	  milab[3]=fNdet[fModule];
	  AliITSRecPoint cc(milab,hit,info);
	  cc.SetType(nClust*100+npeaks);
	  cc.SetDriftTime(driftTimeUncorr);
	  cc.SetDriftSide(s);
	  cc.SetChargeRatio(maxADC);
	  if(clusters) new (cl[ncl]) AliITSRecPoint(cc); 
	  else {
	    fDetTypeRec->AddRecPoint(cc);
	  }
	  ncl++;
	}
      }
    }
  }
  AliDebug(2,Form("Clusters found on SDD module %d (unfolding %d) = %d\n",fModule,repa->GetUseUnfoldingInClusterFinderSDD(),ncl));

} 
//______________________________________________________________________
void AliITSClusterFinderV2SDD::RawdataToClusters(AliRawReader* rawReader){
    //------------------------------------------------------------
  // This function creates ITS clusters from raw data
  //------------------------------------------------------------

  AliITSRawStream* inputSDD=AliITSRawStreamSDD::CreateRawStreamSDD(rawReader);
  AliDebug(1,Form("%s is used",inputSDD->ClassName()));

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
  FindClustersSDD(inputSDD);
  delete inputSDD;
}

void AliITSClusterFinderV2SDD::FindClustersSDD(AliITSRawStream* input) 
{
  //------------------------------------------------------------
  // Actual SDD cluster finder for raw data
  //------------------------------------------------------------
  AliITSRecPointContainer* rpc = AliITSRecPointContainer::Instance();
  Int_t nClustersSDD = 0;
  AliBin *bins[2];
  TBits* anodeFired[2];
  TBits* ddlAnodeFired[kHybridsPerDDL];
  for(Int_t iHyb=0;iHyb<kHybridsPerDDL;iHyb++){
    ddlAnodeFired[iHyb]=new TBits(fNAnodes);
    ddlAnodeFired[iHyb]->ResetAllBits();
  }
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
	  TClonesArray* clusters = rpc->UncheckedGetClusters(fModule);
	  bins[0]=fDDLBins[iMod*2];   // first hybrid of the module
	  bins[1]=fDDLBins[iMod*2+1]; // second hybrid of the module
	  anodeFired[0]=ddlAnodeFired[iMod*2];
	  anodeFired[1]=ddlAnodeFired[iMod*2+1];
	  FindClustersSDD(bins, anodeFired, NULL, clusters,jitter);
	  Int_t nClusters = clusters->GetEntriesFast();
	  nClustersSDD += nClusters;
	  vectModId[iMod]=-1;
	}
	for (Int_t s=0; s<2; s++){
	  Int_t indexHyb=iMod*2+s;
	  for(Int_t iAnode=0; iAnode<GetSeg()->NpzHalf(); iAnode++){
	    if(ddlAnodeFired[indexHyb]->TestBitNumber(iAnode)==kFALSE) continue;
	    for(Int_t iTimeBin=0; iTimeBin<GetSeg()->Npx(); iTimeBin++){
	      Int_t index=(iTimeBin+1)*fNZbins+(iAnode+1);
	      fDDLBins[indexHyb][index].Reset();
	    }
	  }
	}
	ddlAnodeFired[iMod*2]->ResetAllBits();
	ddlAnodeFired[iMod*2+1]->ResetAllBits();
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
      Int_t chan=input->GetCoord1()+fNAnodes*iSide;
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
	  Int_t index = (itb+1) * fNZbins + (iz+1);
	  if((itb < fNTimeBins) && (iz < fNAnodes)) {
	    fDDLBins[iHybrid][index].SetQ(q);
	    fDDLBins[iHybrid][index].SetMask(1);
	    fDDLBins[iHybrid][index].SetIndex(index);
	    ddlAnodeFired[iHybrid]->SetBitNumber(iz);
	  }else{
	    AliWarning(Form("Invalid SDD cell: Anode=%d   TimeBin=%d",iz,itb));	  
	  }
	}
      }
    }
  }
  for(Int_t iHyb=0;iHyb<kHybridsPerDDL;iHyb++){ 
   delete ddlAnodeFired[iHyb];
  }
  AliDebug(1,Form("found clusters in ITS SDD: %d", nClustersSDD));
}

//______________________________________________________________________
Bool_t AliITSClusterFinderV2SDD::NoiseSuppress(Int_t k, Int_t sid, AliBin* bins, AliITSCalibrationSDD* cal) const {
  // applies zero suppression using the measured noise of each anode
  // threshold values from ALICE-INT-1999-28 V10
  // returns kTRUE if the digit should eb noise suppressed, kFALSE if it should be kept
  Float_t xfactL=2.2; 
  Float_t xfactH=4.0;
  //

  Int_t iAn=(k%fNZbins)-1;
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
  Int_t eE=bins[k-fNZbins].GetQ();
  if(eE>tL) nLow++;
  if(eE>tH) nHigh++;
  Int_t wW=bins[k+fNZbins].GetQ();
  if(wW>tL) nLow++;
  if(wW>tH) nHigh++;
  if(nLow<2 || nHigh<1) return kTRUE;
  else return kFALSE;
}



