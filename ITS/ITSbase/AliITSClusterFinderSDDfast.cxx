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
//         Implementation of the ITS SDD fast clusterer  class            //
//                                                                        //
//   Origin: Simone Capodicasa, Universita e INFN, capodica@to.infn.it    //
//                                                                        //
////////////////////////////////////////////////////////////////////////////

#include <vector>
#include <TClonesArray.h>
#include <TBits.h>
#include <TMath.h>
#include <TH2F.h>
#include <TFile.h>
#include "AliITSClusterFinderSDDfast.h"
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

ClassImp(AliITSClusterFinderSDDfast)

AliITSClusterFinderSDDfast::AliITSClusterFinderSDDfast(AliITSDetTypeRec* dettyp):AliITSClusterFinder(dettyp),
  fNAnodes(0),
  fNTimeBins(0),
  fNZbins(0),
  fNXbins(0),
  fDDLBins(),
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
  SetPeakSelection(15.,30.,2000.);
  fDDLBins.resize(kHybridsPerDDL);
}


//______________________________________________________________________
AliITSClusterFinderSDDfast::~AliITSClusterFinderSDDfast()
{
  //Default destructor
}

//______________________________________________________________________
void AliITSClusterFinderSDDfast::FindRawClusters(Int_t mod){

  //Find clusters 
  SetModule(mod);
  FindClustersSDD(fDigits);

}

//______________________________________________________________________
void AliITSClusterFinderSDDfast::FindClustersSDD(TClonesArray *digits){

  std::vector<int> bins0;
  std::vector<int> bins1;
  const Int_t kMapDim=fNZbins*fNXbins/32;
  Int_t map0[kMapDim];
  Int_t map1[kMapDim];
  for(Int_t j=0;j<kMapDim;++j){
    map0[j]=map1[j]=0;
  }
  AliITSCalibrationSDD* cal = (AliITSCalibrationSDD*)GetResp(fModule);
  if(cal==0){
    AliError(Form("Calibration object not present for SDD module %d\n",fModule));
    return;
  }

  AliITSdigitSDD *d=0;
  Int_t i, ndigits=digits->GetEntriesFast();
  for (i=0; i<ndigits; i++){
    d=(AliITSdigitSDD*)digits->UncheckedAt(i);
    Int_t ian=d->GetCoord1();
    Int_t itb=d->GetCoord2();
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
      Int_t iindex=y*fNZbins+z;
      Float_t noise=cal->GetNoiseAfterElectronics(ian)*2.2; // applies zero suppression using the measured noise of each anode. Threshold values from ALICE-INT-1999-28 V10
      if (z<=fNAnodes){
        if(q>noise){
          bins0.push_back(iindex);
          bins0.push_back(q);
          bins0.push_back(0);
          bins0.push_back(i);
          map0[iindex/32]|=(1<<(iindex%32));
        }
      }
      else{
        z-=fNAnodes;
        if(q>noise){
          iindex=y*fNZbins+z;
          bins1.push_back(iindex);
          bins1.push_back(q);
          bins1.push_back(0);
          bins1.push_back(i);
          map1[iindex/32]|=(1<<(iindex%32));
        }
      }
    }
  }
  FindClustersSDD(bins0, bins1, map0, map1, digits);
}

//______________________________________________________________________
void AliITSClusterFinderSDDfast::FindClustersSDD(std::vector<int>& bins0, std::vector<int>& bins1, const Int_t map0[], const Int_t map1[], TClonesArray *digits, TClonesArray *clusters, Int_t jitter){

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

  TClonesArray &cl=*clusters;
  Int_t nrp=0;
  for (Int_t s=0; s<2; s++){
    Int_t *bins;
    unsigned int binssize;
    const Int_t* map;
    if(s==0){
      binssize=bins0.size();
      bins = &bins0[0];
      map=map0;
    }
    if(s==1){
      binssize=bins1.size();
      bins=&bins1[0];
      map=map1;
    }

    const Int_t rresto=fNZbins-1;
    Int_t cid=0;
    for(std::vector<int>::size_type i=0;i<binssize;i+=4){
      if(!bins[i+2])
	bins[i+2]=++cid;
      Int_t me=bins[i];
      Int_t resto=me%fNZbins;
      if(resto==rresto){
        Int_t idxs[1]={me+fNZbins};
        for(Int_t k=0;k<1;++k)
	  if(map[idxs[k]/32]&(1<<(idxs[k]%32)))
	    for(std::vector<int>::size_type j=i+4;j<binssize;j+=4)
	      if(bins[j]==idxs[k]){
		bins[j+2]=bins[i+2];
		break;
	      }
      }
      else{
        Int_t idxs[2]={me+1,me+fNZbins};
        for(int k=0;k<2;++k)
	  if(map[idxs[k]/32]&(1<<(idxs[k]%32)))
	    for(std::vector<int>::size_type j=i+4;j<binssize;j+=4)
	      if(bins[j]==idxs[k]){
		bins[j+2]=bins[i+2];
		break;
	      }
      }
    }
    for(std::vector<int>::size_type i=0;i<binssize;i+=4){
      Int_t me=bins[i];
      Int_t resto=me%fNZbins;
      if(resto==fNZbins-1){
        Int_t idxs[1]={me+fNZbins};
        Int_t myid=bins[i+2];
        for(Int_t l=0;l<1;++l){
          if(map[idxs[l]/32]&(1<<(idxs[l]%32)))
	    for(std::vector<int>::size_type j=i+4;j<binssize;j+=4){
	      if(bins[j]==idxs[l]){
		Int_t hisid=bins[j+2];
		if(myid!=hisid){
		  for(std::vector<int>::size_type k=2;k<binssize;k+=4)
		    if(bins[k]==hisid)
		      bins[k]=myid;
		}
		break;
	      }
	    }
        }
      }
      else{
        Int_t idxs[2]={me+1,me+fNZbins};
        Int_t myid=bins[i+2];
        for(Int_t l=0;l<2;++l){
          if(map[idxs[l]/32]&(1<<(idxs[l]%32)))
	    for(std::vector<int>::size_type j=i+4;j<binssize;j+=4){
	      if(bins[j]==idxs[l]){
		Int_t hisid=bins[j+2];
		if(myid!=hisid){
		  for(std::vector<int>::size_type k=2;k<binssize;k+=4)
		    if(bins[k]==hisid)
		      bins[k]=myid;
		}
		break;
	      }
	    }
        }
      }
    }

    Int_t recp[cid][12];
    for(Int_t i=0;i<cid;++i)
      for(Int_t j=0;j<12;++j)
	recp[i][j]=0;

    Int_t kplab[cid][10];
    for(Int_t i=0;i<cid;++i)
      for(Int_t j=0;j<10;++j)
	kplab[i][j]=-2;

    for(std::vector<int>::size_type i=0;i<binssize;i+=4){
      Int_t q=bins[i+1];
      Int_t me=bins[i+2]-1;
      Int_t z=bins[i]%fNZbins;
      Int_t x=bins[i]/fNZbins;
      recp[me][0]+=q;                   //sumq
      recp[me][1]+=z*q;                 //sumz
      recp[me][2]+=x*q;                 //sumx

#ifdef CSBASEDERROR
      recp[me][3]+=z*z*q;               //sigmaZ2
      recp[me][4]+=x*x*q;               //sigmaX2
#endif

      if(recp[me][5]==0){
        recp[me][6]=z;
        recp[me][7]=z;
        recp[me][8]=x;
        recp[me][9]=x;
        recp[me][10]=q;
        recp[me][11]=bins[i];
      }
      else{
        if(recp[me][6]<z) recp[me][6]=z;
        if(recp[me][7]>z) recp[me][7]=z;
        if(recp[me][8]<x) recp[me][8]=x;
        if(recp[me][9]>x) recp[me][9]=x;
        if(recp[me][10]<q){
          recp[me][10]=q;
          recp[me][11]=bins[i];
        }
      }

      if(digits){
        Int_t kplab2[10];
        for(Int_t ilab=0;ilab<10;++ilab)
	  kplab2[ilab]=kplab[me][ilab];
        AliITSdigitSDD* d=(AliITSdigitSDD*)digits->UncheckedAt(bins[i+3]);
        for (Int_t itrack=0;itrack<10;itrack++){
          Int_t track = (d->GetTracks())[itrack];
          if (track>=0) {
            AddLabel(kplab2, track);
          }
        }
        for(Int_t ilab=0;ilab<10;++ilab)
	  kplab[me][ilab]=kplab2[ilab];
      }
      ++recp[me][5];                    //nPiInClu
    }

    for(Int_t i=0;i<cid;++i){
      if(recp[i][5]==0) continue;
      if(recp[i][5]==1) continue;

      Float_t q=recp[i][0];

      Int_t clSizAnode=recp[i][6]-recp[i][7]+1;
      Int_t clSizTb=recp[i][8]-recp[i][9]+1;
      if(repa->GetUseSDDClusterSizeSelection()){
        if(clSizTb==1) continue; // cut common mode noise spikes
        if(clSizAnode>5) continue; // cut common mode noise spikes
        if(clSizTb>10) continue; // cut clusters on noisy anodes
        if(cal->IsAMAt20MHz() && clSizTb>8) continue; // cut clusters on noisy anodes
      }

      Float_t zz=(Float_t)recp[i][1]/q;
      Float_t xx=(Float_t)recp[i][2]/q;

      AliITSresponseSDD* rsdd = fDetTypeRec->GetResponseSDD();

      Float_t zAnode=zz-0.5;  // to have anode in range 0.-255. and centered on the mid of the pitch
      Float_t timebin=xx-0.5;  // to have time bin in range 0.-255. amd centered on the mid of the bin

      if(s==1) zAnode+=fNAnodes;  // right side has anodes from 256. to 511.

      Float_t zdet=GetSeg()->GetLocalZFromAnode(zAnode);
      Float_t driftTimeUncorr=GetSeg()->GetDriftTimeFromTb(timebin)+jitter*rsdd->GetCarlosRXClockPeriod();
      Float_t driftTime=driftTimeUncorr-rsdd->GetTimeZero(fModule);

      if(driftTime<fMaxDrTimeForTightCut && recp[i][10]<fCutOnPeakTight) continue;
      if(recp[i][10]<fCutOnPeakLoose) continue;

      Float_t driftSpeed=cal->GetDriftSpeedAtAnode(zAnode) + rsdd->GetDeltaVDrift(fModule,zAnode>255);
      Float_t driftPathMicron=driftTime*driftSpeed;
      const Double_t kMicronTocm=1.0e-4;
      Float_t xdet=(driftPathMicron-GetSeg()->Dx())*kMicronTocm; // xdet is negative
      if(s==0) xdet=-xdet; // left side has positive local x

      if(repa->GetUseSDDCorrectionMaps()){
        Float_t corrx=0, corrz=0;
        cal->GetCorrections(zdet,xdet,corrz,corrx,GetSeg());
        zdet+=corrz;
        xdet+=corrx;
      }

      Double_t loc[3]={xdet,0.,zdet},trk[3]={0.,0.,0.};
      mT2L->MasterToLocal(loc,trk);
      xx=trk[1];
      zz=trk[2];

      Double_t dEdxslope;
      if(digits) dEdxslope=rsdd->GetADCvsDriftTime(fModule,kTRUE);
      else dEdxslope=rsdd->GetADCvsDriftTime(fModule);
      q+=(driftTime*dEdxslope); // correction for zero supp.
      q/=rsdd->GetADCtokeV(fModule);
      if(cal-> IsAMAt20MHz()) q*=2.; // account for 1/2 sampling freq.
      if(q<repa->GetMinClusterChargeSDD()) continue; // remove noise clusters

#ifdef  CSBASEDERROR
      Float_t hit[6]={xx,zz,recp[i][3],recp[i][4],q,0.};
#else
      Float_t hit[6]={xx, zz, 0.0030*0.0030, 0.0020*0.0020, q, 0.};
#endif

      Int_t  info[3]={clSizTb, clSizAnode, fNlayer[fModule]};

      Int_t kplab2[10];
      if(digits){
        for(Int_t ilab=0;ilab<10;++ilab)
	  if(kplab[i][ilab]!=-2)
	    kplab2[ilab]=kplab[i][ilab];
	  else
	    kplab2[ilab]=-2;
      }
      else{
        for(Int_t ilab=10;ilab--;) kplab2[ilab]=-2;
        if(fRawID2ClusID) kplab2[0]=fNClusters+1; // RS: store clID+1 as a reference to the cluster
      }
      if(digits) CheckLabels2(kplab2);
      kplab2[3]=fNdet[fModule];
      AliITSRecPoint cc(kplab2,hit,info);
      cc.SetType(101);
      cc.SetDriftTime(driftTimeUncorr);
      cc.SetDriftSide(s);
      cc.SetChargeRatio(recp[i][10]);
      if(clusters) new (cl[nrp]) AliITSRecPoint(cc);
      else {
        fDetTypeRec->AddRecPoint(cc);
      }
      ++nrp;
      ++fNClusters;
    }
  }
  AliDebug(2,Form("Clusters found on SDD module %d (unfolding %d) = %d\n",fModule,repa->GetUseUnfoldingInClusterFinderSDD(),nrp));
}

//______________________________________________________________________
void AliITSClusterFinderSDDfast::RawdataToClusters(AliRawReader* rawReader){
  //------------------------------------------------------------
  // This function creates ITS clusters from raw data
  //------------------------------------------------------------
  fNClusters = 0; //RS
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

void AliITSClusterFinderSDDfast::FindClustersSDD(AliITSRawStream* input){

  AliITSRecPointContainer* rpc = AliITSRecPointContainer::Instance();
  Int_t nClustersSDD=0;

  const Int_t kMapDim=fNZbins*fNXbins/32;
  Int_t mapsDDL[kHybridsPerDDL][kMapDim];
  for(Int_t i=0;i<kHybridsPerDDL;++i)
    for(Int_t j=0;j<kMapDim;++j)
      mapsDDL[i][j]=0;

  Int_t vectModId[kModulesPerDDL];
  for(Int_t iMod=0; iMod<kModulesPerDDL; iMod++) vectModId[iMod]=-1;
  // read raw data input stream
  int countRW = 0; //RS
  if (fRawID2ClusID) fRawID2ClusID->Reset(); //RS if array was provided, we shall store the rawID -> ClusterID

  while (input->Next()) {
    Int_t iModule = input->GetModuleID();
    if(iModule<0){
      AliWarning(Form("Invalid SDD module number %d\n", iModule));
      continue;
    }
    Int_t iCarlos=input->GetCarlosId();
    Int_t iSide=input->GetChannel();
    Int_t iHybrid=iCarlos*2+iSide;

    if(input->IsCompletedModule()){
      vectModId[iCarlos]=iModule; // store the module number
    }
    else if(input->IsCompletedDDL()){
      // when all data from a DDL was read, search for clusters
      Int_t jitter=input->GetJitter();
      for(Int_t iMod=0;iMod<kModulesPerDDL;iMod++){
        if(vectModId[iMod]>=0){
          fModule = vectModId[iMod];
          TClonesArray* clusters=rpc->UncheckedGetClusters(fModule);
	  std::vector<int> bins0;
	  std::vector<int> bins1;
          bins0=fDDLBins[iMod*2];
          bins1=fDDLBins[iMod*2+1];
	  Int_t map0[kMapDim];
	  Int_t map1[kMapDim];
	  for(Int_t j=0;j<kMapDim;++j){
	    map0[j]=map1[j]=0;
	  }
          for(Int_t i=0;i<kMapDim;++i){
            map0[i]=mapsDDL[iMod*2][i];
            map1[i]=mapsDDL[iMod*2+1][i];
          }

          FindClustersSDD(bins0, bins1, map0, map1, NULL, clusters,jitter);

          Int_t nClusters = clusters->GetEntriesFast();
          nClustersSDD += nClusters;
          vectModId[iMod]=-1;
        }
        for (Int_t s=0; s<2; s++){
          Int_t indexHyb=iMod*2+s;
	  for(std::vector<int>::size_type i=0;i<fDDLBins[indexHyb].size();++i)
	    fDDLBins[indexHyb].erase(fDDLBins[indexHyb].begin(),fDDLBins[indexHyb].end());
	  for(Int_t j=0;j<kMapDim;++j)
	    mapsDDL[indexHyb][j]=0;
        }
      }
    }
    else{ // fill the current digit into the bins array
      if(iHybrid<0 || iHybrid>=kHybridsPerDDL){
        AliWarning(Form("Invalid SDD hybrid number %d on module %d\n", iHybrid,iModule));
        continue;
      }
      AliITSCalibrationSDD* cal=(AliITSCalibrationSDD*)GetResp(iModule);
      if(cal==0){
        AliError(Form("Calibration object not present for SDD module %d\n",iModule));
        continue;
      }
      Float_t charge=input->GetSignal();
      Int_t chan=input->GetCoord1()+fNAnodes*iSide;
      Float_t gain=cal->GetChannelGain(chan)/fDetTypeRec->GetAverageGainSDD();;
      Float_t baseline=cal->GetBaseline(chan);
      if(charge>baseline) charge-=baseline;
      else charge=0;
      if(gain>0.){ // Bad channels have gain=0
        charge/=gain;
        if(charge>=cal->GetThresholdAnode(chan)){
          Int_t q=(Int_t)(charge+0.5);
          Int_t iz = input->GetCoord1();
          Int_t itb = input->GetCoord2();
          Int_t ian=iz;
          if(iSide==1)
            ian+=256;
          Float_t noise=cal->GetNoiseAfterElectronics(ian)*2.2;  // applies zero suppression using the measured noise of each anode. Threshold values from ALICE-INT-1999-28 V10
          Int_t index=(itb+1)*fNZbins+(iz+1);
          if((itb<fNTimeBins) && (iz<fNAnodes)){
            if(q<noise) continue;
            fDDLBins[iHybrid].push_back(index);
            fDDLBins[iHybrid].push_back(q);
            fDDLBins[iHybrid].push_back(0);
            fDDLBins[iHybrid].push_back(countRW);
            mapsDDL[iHybrid][index/32]|=(1<<(index%32));
          }
          else{
            static bool show_info = !(getenv("HLT_ONLINE_MODE") && strcmp(getenv("HLT_ONLINE_MODE"), "on") == 0);
	    static int nErrors = 0;
            if (show_info || nErrors++ < 10)
	    {
        	AliWarning(Form("Invalid SDD cell: Carlos=%d,Module=%d,Anode=%d,TimeBin=%d",iCarlos,iModule,iz,itb));
    	    }
          }
        }
      }
    }
    countRW++; //RS
  }
  AliDebug(1,Form("found clusters in ITS SDD: %d", nClustersSDD));
}
