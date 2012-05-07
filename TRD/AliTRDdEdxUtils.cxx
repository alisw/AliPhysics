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

//
// class to calculate TRD dEdx
// xx
// xx
// xx
// xx
//
//  Xianguo Lu 
//  lu@physi.uni-heidelberg.de
//  Xianguo.Lu@cern.ch
//  
//


#include "TF1.h"
#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "THn.h"
#include "THnSparse.h"
#include "TMath.h"
#include "TMatrixD.h"
#include "TMinuit.h"
#include "TObjArray.h"
#include "TRandom3.h"
#include "TStopwatch.h"
#include "TVectorD.h"

#include "TTreeStream.h"

#include "AliCDBId.h"
#include "AliCDBMetaData.h"
#include "AliCDBStorage.h"
#include "AliESDEvent.h"
#include "AliESDfriendTrack.h"
#include "AliESDtrack.h"
#include "AliTRDcalibDB.h"
#include "AliTRDCalROC.h"
#include "AliTRDtrackV1.h"


#include "AliTRDdEdxUtils.h"

#define EPSILON 1e-12

THnF * AliTRDdEdxUtils::fgHistGain=0x0;
THnF * AliTRDdEdxUtils::fgHistT0=0x0;
THnF * AliTRDdEdxUtils::fgHistVd=0x0;
TObjArray * AliTRDdEdxUtils::fgHistPHQ=new TObjArray(8);

TString AliTRDdEdxUtils::fgCalibFile;
TObjArray * AliTRDdEdxUtils::fgObjGain = 0x0;
TObjArray * AliTRDdEdxUtils::fgObjT0 = 0x0;
TObjArray * AliTRDdEdxUtils::fgObjVd = 0x0;
TObjArray * AliTRDdEdxUtils::fgObjPHQ = 0x0;

Int_t AliTRDdEdxUtils::fgNchamber = -999;
Double_t AliTRDdEdxUtils::fgChamberQ[6];
Double_t AliTRDdEdxUtils::fgChamberTmean[6];

Double_t AliTRDdEdxUtils::fgTrackTmean = -999;

Bool_t   AliTRDdEdxUtils::fgPadGainOn = kTRUE;
Bool_t   AliTRDdEdxUtils::fgExBOn = kTRUE; 
Double_t AliTRDdEdxUtils::fgQScale = 45;
Double_t AliTRDdEdxUtils::fgQ0Frac = 0.3;
Double_t AliTRDdEdxUtils::fgQ1Frac = 0.5;
Double_t AliTRDdEdxUtils::fgTimeBinCountCut = 0.0; 
Int_t    AliTRDdEdxUtils::fgCalibTPCnclsCut = 70;

//===================================================================================
//                                   Math and Histogram
//===================================================================================
void AliTRDdEdxUtils::GetCDFCuts(const TH1D *hh, const Int_t ncut, Double_t cuts[], const Double_t cdfs[], const Double_t thres)
{
  //
  //counts of hh is sorted
  //

  for(Int_t ii=0; ii<ncut; ii++){
    cuts[ii] = -999;
  }

  Int_t nsel = 0; 
  const Int_t nbin = hh->GetNbinsX();
  Double_t datas[nbin];
  for(Int_t ii=1; ii<=nbin; ii++){
    const Double_t res = hh->GetBinContent(ii);
    if(res<thres){
      continue;
    }

    datas[nsel] = res;
    nsel++;
  }
  if(!nsel)
    return;

  Int_t id[nsel];
  TMath::Sort(nsel, datas, id, kFALSE);

  for(Int_t ii=0; ii<ncut; ii++){
    const Double_t icdf = cdfs[ii];
    if(icdf<0 || icdf>1){
      printf("AliTRDdEdxUtils::GetCDFCuts error cdfs[%d] %15f out of range!\n", ii, icdf); exit(1);
    }
    cuts[ii] = datas[id[Int_t(icdf*nsel)]];
  }
}

Double_t AliTRDdEdxUtils::GetMeanRMS(const Double_t nn, const Double_t sum, const Double_t w2s, Double_t * grms, Double_t * gerr)
{
  //
  //calculate mean (with error) and rms from sum, w2s, nn
  //if nn=0, mean, error, and rms all = 0
  //

  Double_t tmean = 0, trms = 0, terr = 0;

  if(nn>EPSILON){
    tmean = sum/nn;

    const Double_t arg = w2s/nn-tmean*tmean;
    if(TMath::Abs(arg)<EPSILON){
      trms = 0;
    }
    else{
      if( arg <0 ){
        printf("AliTRDdEdxUtils::GetMeanRMS error negative sqrt argument!! %e -- %e %e %f\n", arg, w2s, sum, nn); exit(1);
      }
      
      trms = TMath::Sqrt(arg);
    }

    terr = trms/TMath::Sqrt(nn);
  }

  if(grms){
    (*grms) = trms;
  }

  if(gerr){
    (*gerr) = terr;
  }

  return tmean;
}

Double_t AliTRDdEdxUtils::TruncatedMean(const Int_t nx, const Double_t xdata[], const Double_t lowfrac, const Double_t highfrac, Double_t * grms, Double_t * gerr, Double_t *wws)
{
  //
  //calculate truncated mean
  //return <x*w>_{low-high according to x}
  //

  /*
  //test->
  for(Int_t ii=0; ii<nx; ii++){
    printf("test %d/%d %f\n", ii, nx, xdata[ii]);
  }
  //<--test
  */

  Int_t index[nx];
  TMath::Sort(nx, xdata, index, kFALSE);

  Int_t nused = 0;
  Double_t sum = 0;
  Double_t w2s = 0;
  const Int_t istart = Int_t (nx*lowfrac);
  const Int_t istop  = Int_t (nx*highfrac);

  //=,< correct, because when low=0, high=1 it is correct
  for(Int_t ii=istart; ii<istop; ii++){
    Double_t weight = 1;
    if(wws){
      weight = wws[index[ii]];
    }
    const Double_t sx = xdata[index[ii]]*weight;

    sum += sx;
    w2s += sx*sx;

    nused++;
    //printf("test in loop %d/%d %f %f %f\n", ii, nused, sx, sum, w2s);
    
  }

  return GetMeanRMS(nused, sum, w2s, grms, gerr);
}

Double_t AliTRDdEdxUtils::TruncatedMean(const TH1 *hh, const Double_t lowfrac, const Double_t highfrac, Double_t * grms, Double_t * gerr)
{
  //
  //do truncation on histogram
  //
  //if hh is scaled, be sure Sumw2 is called before scaling!! then mean, rms and err will all be correct
  
  //with under-/over-flow
  Double_t npreTrunc = 0;
  for(Int_t itmp=0; itmp<=hh->GetNbinsX()+1; itmp++){
    const Double_t be = hh->GetBinError(itmp);
    const Double_t bc = hh->GetBinContent(itmp);
    if(be<EPSILON){
      if(bc>EPSILON){
        printf("AliTRDdEdxUtils::TruncatedMean (hist) error %e %e %d\n", bc, be, itmp); exit(1);
      }
      continue;
    }
    npreTrunc += bc*bc/be/be;
  }

  const Double_t nstart = npreTrunc*lowfrac;
  const Double_t nstop = npreTrunc*highfrac;

  //with Double_t this should also handle normalized hist
  Double_t ntot = 0;
  Double_t nused = 0;
  Double_t sum = 0;
  Double_t w2s = 0;
  for(Int_t itmp=0; itmp<=hh->GetNbinsX()+1; itmp++){
    const Double_t be = hh->GetBinError(itmp);
    const Double_t bc = hh->GetBinContent(itmp);
    if(be<EPSILON){
      if(bc>EPSILON){
        printf("AliTRDdEdxUtils::TruncatedMean (hist) error %e %e %d\n", bc, be, itmp); exit(1);
      }
      continue;
    }
    const Double_t weight = bc*bc/be/be;
    ntot+=weight;
    //<= correct, because when high=1, nstop has to be included
    if(ntot>nstart && ntot<=nstop){

      const Double_t val = hh->GetBinCenter(itmp);
      sum += weight*val;
      w2s += weight*val*val;
    
      nused += weight;

      //printf("test %d %f %f --- %f %f -- %f %f\n", itmp, weight, val, sum, w2s, nused, nsample);
    }
    else if(ntot>nstop){
      if(itmp>=hh->GetNbinsX()){
        printf("AliTRDdEdxUtils::TruncatedMean warning hist range too small %s %f %f %d %d, %15f %15f %15f; nused w2s sum set to 0\n", hh->GetName(), hh->GetBinLowEdge(1), hh->GetBinLowEdge(itmp), itmp, hh->GetNbinsX(), hh->GetBinContent(hh->GetNbinsX())/hh->Integral(0,hh->GetNbinsX()+1), hh->GetBinContent(hh->GetNbinsX()), hh->Integral(0,hh->GetNbinsX()+1)); //exit(1);
        nused = 0;
        w2s = sum = 0;
      }
      break;
    }
  }

  return GetMeanRMS(nused, sum, w2s, grms, gerr);
}

void AliTRDdEdxUtils::FitSlicesY(const TH2D *hh, TH1D *&hnor, TH1D *&hmpv, TH1D *&hwid, TH1D *&hres, const Double_t thres, const Double_t lowfrac, const Double_t highfrac)
{
  //
  //fit slices of hh using truncation
  //

  const Int_t x0 = hh->GetXaxis()->GetFirst();
  const Int_t x1 = hh->GetXaxis()->GetLast();
  const Int_t y0 = hh->GetYaxis()->GetFirst();
  const Int_t y1 = hh->GetYaxis()->GetLast();

  const Int_t nx = hh->GetNbinsX();
  const Int_t ny = hh->GetNbinsY();
  const Double_t xmin = hh->GetXaxis()->GetXmin();
  const Double_t xmax = hh->GetXaxis()->GetXmax();
  const Double_t ymin = hh->GetYaxis()->GetXmin();
  const Double_t ymax = hh->GetYaxis()->GetXmax();

  hnor = new TH1D(Form("%s_amp",hh->GetName()), "", nx, xmin, xmax);
  hmpv = new TH1D(Form("%s_mpv",hh->GetName()), "", nx, xmin, xmax);
  hwid = new TH1D(Form("%s_wid",hh->GetName()), "", nx, xmin, xmax);
  hres = new TH1D(Form("%s_res",hh->GetName()), "", nx, xmin, xmax);

  for(Int_t ix=x0; ix<=x1; ix++){
    //to speed up
    const Double_t rawcount = hh->Integral(ix,ix,0, ny+1);
    if(rawcount<EPSILON){
      continue;
    }

    TH1D *htmp = new TH1D(Form("FitSlicesY_%s_%d", hh->GetName(), ix),"",ny, ymin, ymax);
    Double_t ntot = 0;
    for(Int_t iy=y0; iy<=y1; iy++){
      const Double_t be = hh->GetBinError(ix,iy);
      const Double_t bc = hh->GetBinContent(ix, iy);

      if(be<EPSILON){
        if(bc>EPSILON){
          printf("AliTRDdEdxUtils::FitSlicesY error %d %d %e %e\n", ix, iy, be, bc); exit(1);
        }
        continue;
      }

      htmp->SetBinContent(iy, bc);
      htmp->SetBinError(iy, be);

      ntot += (bc/be)*(bc/be);

      //if(be) printf("test %d %d : %f %f %f\n", ix, iy, bc, be, pow(bc/be,2));
    }

    hnor->SetBinContent(ix, ntot);
    hnor->SetBinError(  ix, 0);
    
    if(ntot<thres || htmp->GetRMS()<EPSILON){
      delete htmp;
      continue;
    }

    //test htmp->Draw();
    Double_t trms = -999, terr = -999;
    const Double_t tmean = TruncatedMean(htmp, lowfrac, highfrac, &trms, &terr);

    hmpv->SetBinContent(ix, tmean);
    hmpv->SetBinError(  ix, terr);

    hwid->SetBinContent(ix, trms);
    hwid->SetBinError(  ix, 0);

    hres->SetBinContent(ix, tmean>EPSILON ? trms/tmean:0);
    hres->SetBinError(  ix, 0);

    delete htmp;
  }

  TH1 *hhs[]={hnor, hmpv, hwid, hres};
  const TString yt[]={"N", "MPV", "#sigma", "#sigma/MPV"};
  const Int_t nh = sizeof(hhs)/sizeof(TH1*);
  for(Int_t ii=0; ii<nh; ii++){
    hhs[ii]->SetYTitle(Form("%s of %s", yt[ii].Data(), hh->GetYaxis()->GetTitle()));
    hhs[ii]->SetXTitle(hh->GetXaxis()->GetTitle());
    hhs[ii]->GetYaxis()->SetTitleOffset(hh->GetYaxis()->GetTitleOffset());
    hhs[ii]->SetTitle(hh->GetTitle());
  }
}

//===================================================================================
//                                TRD Analysis Fast Tool
//===================================================================================

Int_t AliTRDdEdxUtils::GetNtracklet(const AliESDEvent *esd)
{
  //
  //number of trd tracklet in one esd event
  //
  const Int_t ntrk0 = esd->GetNumberOfTracks();
  Int_t ntrdv1=0;
  for(Int_t ii=0; ii<ntrk0; ii++){
    ntrdv1 += esd->GetTrack(ii)->GetTRDntracklets();
  }
  return ntrdv1;
}

AliTRDtrackV1 * AliTRDdEdxUtils::GetTRDtrackV1(const AliESDtrack * esdtrack)
{
  //
  //Get TRD friend track
  //

  AliESDfriendTrack *  friendtrk = (AliESDfriendTrack *)esdtrack->GetFriendTrack();
  if(!friendtrk){
    //printf("xlulog AliAnalysisTaskCosmicTRD::GetTRDtrack no friend!!\n"); exit(1);
    return 0x0;
  }

  TObject *calibObject=0x0;
  AliTRDtrackV1 * trdtrack=0x0;
  for(Int_t l=0; (calibObject=friendtrk->GetCalibObject(l)); l++) {
    if( (trdtrack=dynamic_cast<AliTRDtrackV1*>(calibObject)) )
      break;
  }

  return trdtrack;
}

Bool_t AliTRDdEdxUtils::IsInSameStack(const AliTRDtrackV1 *trdtrack)
{
  //
  // to check if all tracklets are in the same stack, useful in cosmic
  //

  TVectorD secs(18), stks(5);

  for(Int_t ilayer = 0; ilayer < 6; ilayer++){
    AliTRDseedV1 *tracklet=trdtrack->GetTracklet(ilayer);
    if(!tracklet)
      continue;
    
    const Int_t det = tracklet->GetDetector();
    const Int_t isector = AliTRDgeometry::GetSector(det);
    const Int_t istack  = AliTRDgeometry::GetStack(det);

    secs[isector] = 1;
    stks[istack]  = 1;
 }

  if(secs.Sum()!=1 || stks.Sum()!=1){
    return kFALSE;
  }
  else 
    return kTRUE;
}

Bool_t AliTRDdEdxUtils::GetFirstSectorStackMomentum(const AliTRDtrackV1 *trdtrack, Int_t & isec, Int_t & istk, Double_t & mom)
{
  //
  //as function name
  //
  isec = istk = -999;
  mom = -999;

  for(Int_t ilayer = 0; ilayer < 6; ilayer++){
    AliTRDseedV1 *tracklet=trdtrack->GetTracklet(ilayer);
    if(!tracklet)
      continue;
    
    const Int_t det = tracklet->GetDetector();
    isec = AliTRDgeometry::GetSector(det);
    istk = AliTRDgeometry::GetStack(det);

    mom = tracklet->GetMomentum();

    break;
  }

  if(isec<0)
    return kFALSE;
  else 
    return kTRUE;
}

//===================================================================================
//                                Calibration
//===================================================================================
Double_t AliTRDdEdxUtils::GetCalibTPCscale(const Int_t tpcncls, const Double_t tpcsig)
{
  //
  //the scale used in calibration
  //

  if(tpcncls < CalibTPCnclsCut())
    return -999;

  if(tpcsig<EPSILON)
    return -999;

  return tpcsig/120;

}

Int_t AliTRDdEdxUtils::GetPHQIterator(const Bool_t kinvq, const Double_t mag, const Int_t charge)
{
  //
  //iterator for calib obj and hist
  //
  return kinvq*4 + (mag>0)*2 + (charge>0); 
}

TObjArray * AliTRDdEdxUtils::GetObjPHQ()
{
  //
  //return fgObjPHQ, initialized if null
  //

  if(!fgObjPHQ){
    fgObjPHQ = new TObjArray(8);
  }

  return fgObjPHQ;
}

TObjArray * AliTRDdEdxUtils::GetObjPHQ(const Bool_t kinvq, const Double_t mag, const Int_t charge)
{
  //
  //return calib obj
  //
  if(!fgObjPHQ){
    printf("AliTRDdEdxUtils::GetObjPHQ(kinvq, mag, charge) error fgObjPHQ null!!\n"); exit(1);
  }

  return (TObjArray*) fgObjPHQ->At(GetPHQIterator(kinvq, mag, charge));
}

THnF * AliTRDdEdxUtils::GetHistPHQ(const Bool_t kinvq, const Double_t mag, const Int_t charge)
{
  //
  //return calib hist
  //
  return (THnF*) fgHistPHQ->At(GetPHQIterator(kinvq, mag, charge));
}

TString AliTRDdEdxUtils::GetPHQName(const Bool_t kobj, const Int_t iter)
{
  //
  //get name of calib obj/hist of PHQ
  //
  return Form("TRDCalib%sPHQ%d", kobj?"Obj":"Hist", iter);
}

void AliTRDdEdxUtils::DeleteCalibObj()
{
  //
  //delete calib obj
  //
  delete fgObjGain;
  delete fgObjT0;
  delete fgObjVd;
  
  fgObjGain = 0x0;
  fgObjT0 = 0x0;
  fgObjVd = 0x0;

  if(fgObjPHQ){
    fgObjPHQ->SetOwner();
    delete fgObjPHQ;
    fgObjPHQ = 0x0;
  }
}

Bool_t AliTRDdEdxUtils::GenerateDefaultPHQOCDB(const TString path)
{
  //
  //generate default OCDB object PHQ, do like
  //AliTRDdEdxUtils::GenerateDefaultPHQOCDB("local://./")
  //

  TObjArray * arr8 = new TObjArray(8);
  arr8->SetOwner();

  for(Int_t ii=0; ii<8; ii++){
    TObjArray * arr1 = new TObjArray(1);
    arr1->SetOwner();
    arr1->SetName(GetPHQName(1, ii));

    const Int_t nbins = NTRDtimebin();
    TVectorD * vec = new TVectorD(nbins);
    for(Int_t jj=0; jj<nbins; jj++){
      (*vec)[jj] = 1;
    }
    arr1->Add(vec);
    arr8->Add(arr1);
  }

  AliCDBMetaData *metaData= new AliCDBMetaData();
  metaData->SetObjectClassName("TObjArray");
  metaData->SetResponsible("Raphaelle Bailhache and Xianguo Lu");

  AliCDBId id1("TRD/Calib/PHQ", 0, 999999999);
  AliCDBStorage * gStorage = AliCDBManager::Instance()->GetStorage(path);
  gStorage->Put(arr8, id1, metaData);

  delete metaData;
  delete arr8;

  return kTRUE;
}

void AliTRDdEdxUtils::IniCalibObj()
{
  //
  //set CalibObj from file, clone to static calib obj
  //

  DeleteCalibObj();
  
  TFile *cfile=new TFile(fgCalibFile);
  if(!cfile->IsOpen()){
    printf("AliTRDdEdxUtils::IniCalibObj error fgCalibFile not open! %s\n", fgCalibFile.Data());exit(1);
  }

  printf("\nAliTRDdEdxUtils::IniCalibObj file: %s\n", fgCalibFile.Data());

  //---
  const TString objnames[] ={"TRDCalibObjGain", "TRDCalibObjT0", "TRDCalibObjVd"};
  TObjArray ** gobjs[]={           &fgObjGain,        &fgObjT0,        &fgObjVd};

  const Int_t nobj = sizeof(objnames)/sizeof(TString);
  for(Int_t iobj=0; iobj<nobj; iobj++){
    TObjArray *tmpo=0x0;
    cfile->GetObject(objnames[iobj], tmpo);
    if(!tmpo){
      printf("AliTRDdEdxUtils::IniCalibObj error obj %s not found!\n", objnames[iobj].Data()); exit(1);
    }

    (*gobjs[iobj])=(TObjArray*)tmpo->Clone();
    (*gobjs[iobj])->SetOwner();
  }

  fgObjPHQ = new TObjArray(8);
  for(Int_t iter=0; iter<8; iter++){
    const TString objn = GetPHQName(1, iter);
    TObjArray *tmpo=0x0;
    cfile->GetObject(objn, tmpo);
    if(!tmpo){
      printf("AliTRDdEdxUtils::IniCalibObj error obj %s not found!\n", objn.Data()); exit(1);
    }

    TObjArray *obji=(TObjArray*) tmpo->Clone();
    obji->SetOwner();
    fgObjPHQ->AddAt(obji, iter);
  }

  //---
  
  cfile->Close();
  delete cfile;
}

void AliTRDdEdxUtils::DeleteCalibHist()
{
  //
  //delete calib hist
  //
  delete fgHistGain;
  delete fgHistT0;
  delete fgHistVd;

  fgHistGain = 0x0;
  fgHistT0 = 0x0;
  fgHistVd = 0x0;

  //fgHistPHQ owns the hists
  fgHistPHQ->SetOwner();
  fgHistPHQ->Clear();
}

void AliTRDdEdxUtils::IniCalibHist(TList *list, const Bool_t kPHQonly)
{
  //
  //initialize calib hist, list should not own the hist, or list->Clear/delete hist should not be called
  //

  DeleteCalibHist();

  Int_t nbin[2];
  const Double_t xmin[2]={0, 0};
  Double_t xmax[2];

  nbin[0]=NTRDtimebin(); nbin[1]= 11250; xmax[0]=nbin[0]; xmax[1]=20; 
  for(Int_t iter=0; iter<8; iter++){
    const TString hn = GetPHQName(0, iter);
    THnF *hi = new THnF(hn, "", 2, nbin, xmin, xmax);
    //fgHistPHQ owns the hists
    fgHistPHQ->AddAt(hi, iter);
    list->Add(hi);
  }

  if(kPHQonly)
    return;

  nbin[0]=NTRDchamber(); nbin[1]= 11250; xmax[0]=nbin[0]; xmax[1]=20;                 fgHistGain = new THnF("TRDCalibHistGain", "", 2, nbin, xmin, xmax);
  nbin[0]=NTRDchamber(); nbin[1]= 11250; xmax[0]=nbin[0]; xmax[1]=AliTRDseedV1::kNtb; fgHistT0   = new THnF("TRDCalibHistT0",   "", 2, nbin, xmin, xmax);
  nbin[0]=NTRDchamber(); nbin[1]= 11250; xmax[0]=nbin[0]; xmax[1]=AliTRDseedV1::kNtb; fgHistVd   = new THnF("TRDCalibHistVd",   "", 2, nbin, xmin, xmax);

  list->Add(fgHistGain);
  list->Add(fgHistT0);
  list->Add(fgHistVd);
}

Bool_t AliTRDdEdxUtils::ReadCalibHist(const TString filename, const TString listname)
{
  //
  //used in AliTRDPreprocessorOffline
  //read in calib hist from file, only for PHQ
  //
  DeleteCalibHist();

  //maybe already open by others... don't close
  TFile fcalib(filename);
  
  TObjArray * array = (TObjArray*)fcalib.Get(listname);

  for(Int_t iter=0; iter<8; iter++){
    const TString hn = GetPHQName(0, iter);
    THnF * tmph=0x0;
    if(array){
      tmph = (THnF *) array->FindObject(hn);
    }
    else{
      tmph = (THnF *) fcalib.Get(hn);
    }
    if(!tmph){
      printf("AliTRDdEdxUtils::ReadCalibHist warning calib hist not found! %s %s\n", filename.Data(), listname.Data());
      fcalib.ls();
      if(array){
        array->ls();
      }
      return kFALSE;
    }
    THnF *hi = (THnF*)tmph->Clone();
    fgHistPHQ->AddAt(hi, iter);
  }

  return kTRUE;
}

void AliTRDdEdxUtils::FillCalibHist(const Int_t ncls, const TVectorD *arrayQ, const TVectorD *arrayX, THnF * hcalib, const Double_t scale)
{
  //
  //fill calibration hist
  //
  if(!hcalib){printf("AliTRDdEdxUtils::FillCalibHist errro hcalib null!!\n"); exit(1);}

  for(Int_t ii=0; ii<ncls; ii++){
    const Double_t dq = (*arrayQ)[ii];
    const Double_t xx = (*arrayX)[ii];

    const Double_t qmax = hcalib->GetAxis(1)->GetXmax() -0.5 * hcalib->GetAxis(1)->GetBinWidth(1);
    const Double_t xmin = hcalib->GetAxis(0)->GetXmin();
    const Double_t xmax = hcalib->GetAxis(0)->GetXmax();

    if(xx>=xmax || xx<xmin){
      printf("AliTRDdEdxUtils::FillCalibHist error x overflow or underflow! %s %15f %15f %15f\n", hcalib->GetName(),  xx, xmin, xmax); exit(1);
    }

    const Double_t var[]={xx, TMath::Min(dq, qmax)/scale};
    hcalib->Fill(var);
  }
}

void AliTRDdEdxUtils::FillCalibHist(const AliTRDtrackV1 *trdv1, const Bool_t kinvq, const Double_t mag, const Int_t charge, const Double_t scale) 
{
  //
  //get cluster Q and fill calib hist, if kinvq = kTRUE, 1/Q is filled
  //

  THnF * hcalib = AliTRDdEdxUtils::GetHistPHQ(kinvq, mag, charge);

  TVectorD arrayQ(200), arrayX(200);
  const Int_t ncls = AliTRDdEdxUtils::GetArrayClusterQ(kinvq, &arrayQ, &arrayX, trdv1);
  FillCalibHist(ncls, &arrayQ, &arrayX, hcalib, kinvq ? 1/scale : scale);

  static Int_t kprint = 100;
  if(kprint<0){
    printf("\nAliTRDdEdxUtils::FillCalibHist summary: \n");
    printf("\nkinvq= %d;\n", kinvq);
    for(Int_t iq=0; iq<ncls; iq++){
      printf("arrayX[%3d] = %15.0f; arrayQ[%3d] =  %15f;\n", iq, arrayX[iq], iq, arrayQ[iq]);
    }
    printf("\n");
  }
  kprint++;
}

Int_t AliTRDdEdxUtils::ApplyCalib(const Int_t nc0, TVectorD *arrayQ, TVectorD *arrayX, const TObjArray *cobj)
{
  //
  //apply calibration on arrayQ
  //
  if(!cobj){ printf("AliTRDdEdxUtils::ApplyCalib error gain array null!!\n"); exit(1);}

  TVectorD tmpq(arrayQ->GetNrows());
  TVectorD tmpx(arrayX->GetNrows());
  Int_t ncls = 0;

  const TVectorD * gain = (TVectorD*) cobj->At(0); 
  for(Int_t ii=0; ii<nc0; ii++){
    const Double_t dq = (*arrayQ)[ii];
    const Int_t xx = (Int_t)(*arrayX)[ii];
    const Double_t gg = (*gain)[xx];

    if(gg<EPSILON){
      continue;
    }

    tmpq[ncls] = dq*gg;
    tmpx[ncls] = xx;
    ncls++;
  }

  (*arrayQ)=tmpq;
  (*arrayX)=tmpx;

  return ncls;
}

void AliTRDdEdxUtils::GetPHCountMeanRMS(const TH1D *hnor, TH1D *&hmean)
{
  //
  //calculate from the ph calib hist the (mean-3sigma) ph-count in the chamber, save in the TH1D output
  //
  const Int_t ndet = 540;
  TObjArray *obj=new TObjArray(ndet);
  obj->SetOwner();
  for(Int_t ii=0; ii<ndet; ii++){
    obj->Add(new TVectorD(AliTRDseedV1::kNtb));
  }

  //ibin = binlowedge of bin(ibin+1) = the number fills this bin
  for(Int_t ibin=0; ibin<hnor->GetNbinsX(); ibin++){
    const Double_t stat = hnor->GetBinContent(ibin+1);
    if(stat<EPSILON){
      continue;
    }

    const Int_t idet = ToDetector(ibin);
    const Int_t itb  = ToTimeBin(ibin);
    TVectorD *vec=(TVectorD *)obj->At(idet);
    (*vec)[itb] = stat;
  }

  hmean = new TH1D(Form("%sdetmean", hnor->GetName()), "", hnor->GetNbinsX(), hnor->GetXaxis()->GetXmin(), hnor->GetXaxis()->GetXmax());
  for(Int_t ibin=0; ibin<hnor->GetNbinsX(); ibin++){
    const Int_t idet = ToDetector(ibin);
    const TVectorD *vec=(TVectorD *)obj->At(idet);

    Int_t nonzero = 0;
    for(Int_t ii=0; ii<vec->GetNrows(); ii++){
      if((*vec)[ii]>EPSILON){
        nonzero++;
      }
    }

    Double_t mean = 0;
    const Double_t lowfrac = 0.6;
    //if there are too many 0's, reject this chamber by settig mean=rms=0
    if(nonzero> (AliTRDseedV1::kNtb*(1-lowfrac)) ){
      //only highest (1-lowfrac)*31 timebins are used to estimate the mean and rms! important! otherwise the 0' will make rms very large!
      mean = TruncatedMean(AliTRDseedV1::kNtb, vec->GetMatrixArray(), lowfrac, 1);
    }

    hmean->SetBinContent(ibin+1, mean);
  }

  delete obj;
}

void AliTRDdEdxUtils::CalibOutput(const TList *lin, Int_t run)
{
  //
  //produce calibration objects
  //

  TString objnames("TRDCalibHistGain TRDCalibHistT0 TRDCalibHistVd ");
  for(Int_t iter=0; iter<8; iter++){
    objnames+= GetPHQName(0, iter)+" ";
  }

  TList * lout = new TList;
  lout->SetOwner();

  TTreeSRedirector *calibStream = new TTreeSRedirector(Form("TRDCalibStream_%010d.root", run));
    
  const Int_t nh=lin->GetEntries();
  for(Int_t ii=0; ii<nh; ii++){
    const THnF *hh=(THnF*)lin->At(ii);
    const TString hname = hh->GetName();
    if(!objnames.Contains(hname))
      continue;

    TObjArray * cobj0 = GetCalibObj(hh, run, lout, calibStream);
    lout->Add(cobj0);
  }

  //lout->ls();

  //=============================================================
  //=============================================================
  
  TFile *fout=new TFile(Form("TRDCalibObj_%010d.root", run),"recreate");
  fout->cd();
  const Int_t nout=lout->GetEntries();
  for(Int_t ii=0; ii<nout; ii++){
    const TString oname = lout->At(ii)->GetName();
    if(oname.Contains("Obj")){
      TObjArray * cobj = (TObjArray*) lout->At(ii);
      cobj->Write(oname, TObjArray::kSingleKey);
    }
  }
  fout->Save();
  fout->Close();
  delete fout;

  fout=new TFile(Form("TRDCalibList_%010d.root", run),"recreate");
  fout->cd();
  lin->Write();
  lout->Write();
  fout->Save();
  fout->Close();
  delete fout;
  
  delete calibStream;

  /*
    http://root.cern.ch/root/html/TH1.html
    When an histogram is created, a reference to it is automatically added to the list of in-memory objects for the current file or directory. This default behaviour can be changed by:
    
    h->SetDirectory(0);          for the current histogram h
    TH1::AddDirectory(kFALSE);   sets a global switch disabling the reference
    
    When the histogram is deleted, the reference to it is removed from the list of objects in memory. When a file is closed, all histograms in memory associated with this file are automatically deleted. 
  */
  delete lout;
}

TObjArray* AliTRDdEdxUtils::GetCalibObj(const THnF *hh, Int_t run, TList *lout, TTreeSRedirector *calibStream)
{
  //
  //produce calibration objects
  //

  const TString hname = hh->GetName();
  const Bool_t kinvq = TString(hname(hname.First('Q')+1,1)).Atoi()&4;

  //----------------------------------------
  //               Define nbin, tag, cobj0
  //----------------------------------------
  Int_t nbin =-999;
  if(hname.Contains("Gain") || hname.Contains("T0") || hname.Contains("Vd")){
    nbin = NTRDchamber();
  }
  else if(hname.Contains("PHQ")){
    nbin = NTRDtimebin();
  }
  else{
    printf("AliTRDdEdxUtils::GetCalibObj error wrong hname!! %s\n", hname.Data()); exit(1);
  }
    
  TString tag(hname);
  tag.ReplaceAll("Hist","Obj");
  
  TObjArray * cobj0 = new TObjArray(1);
  cobj0->SetOwner();
  cobj0->SetName(tag);
  cobj0->Add(new TVectorD(nbin));
  
  //----------------------------------------
  //               Define lowFrac, highFrac
  //----------------------------------------
  Double_t lowFrac = -999, highFrac = -999;
  if(hname.Contains("Gain") || (hname.Contains("PHQ") && !kinvq) ){
    lowFrac = 0.01; highFrac = Q0Frac();
  }
  else if(hname.Contains("PHQ") && kinvq){
    lowFrac = Q1Frac(); highFrac = 0.99;
  }
  else{
    lowFrac = 0.01;
    highFrac = 0.99;
  }
  
  //----------------------------------------
  //               Get analysis result
  //----------------------------------------
  TH1::AddDirectory(kFALSE);//important!
  TH1D *hnor=0x0, *hmpv=0x0, *hres=0x0, *hwid=0x0, *htrdphmean = 0x0;//if(!lout), these have to be deleted
  TH2D *hpj = hh->Projection(1,0);
  FitSlicesY(hpj, hnor, hmpv, hwid, hres, 0, lowFrac, highFrac);
  if(hname.Contains("PHQ")){
    GetPHCountMeanRMS(hnor, htrdphmean);
    if(lout){
      lout->Add(htrdphmean);
    }
  }
  delete hpj;
  
  if(lout){
    lout->Add(hnor);
    lout->Add(hmpv);
    lout->Add(hwid);
    lout->Add(hres);
  }
  
  //----------------------------------------
  //               Define Counter
  //----------------------------------------
  TVectorD *countDet=0x0;
  TObjArray *countSSL=0x0;

  if(hname.Contains("PHQ") && !kinvq){
    countDet=new TVectorD(540);
    countSSL=new TObjArray(90);//SectorStackLayer
    countSSL->SetOwner();
    for(Int_t ii=0; ii<90; ii++){
      countSSL->Add(new TVectorD(6));
    }
  }

  //----------------------------------------
  //               Fill cobj0
  //----------------------------------------

  //ibin = binlowedge of bin(ibin+1) = the number fills this bin
  for(Int_t ibin=0; ibin<nbin; ibin++){
    Double_t gnor = hnor->GetBinContent(ibin+1);
    Double_t gmpv = hmpv->GetBinContent(ibin+1);
    Double_t gwid = hwid->GetBinContent(ibin+1);
    Double_t gres = hres->GetBinContent(ibin+1);

    //--- set additional cut by kpass
    Bool_t kpass = kTRUE;
    Double_t gtrdphmean = -999;
    if(htrdphmean){
      gtrdphmean = htrdphmean->GetBinContent(ibin+1);
      //chamber no statistics (e.g. too many 0's), not usual, not seen in run 143237
      if(gtrdphmean<EPSILON){
        kpass = kFALSE;
      }
      if(gnor<TimeBinCountCut()*gtrdphmean){
        kpass = kFALSE;
      }
    }
    
    //--- set calibration constant p0
    Double_t p0= 0;
    
    //reason for gmpv=0:
    //1)gnor<=3; truncation in hist: (0, 0.6*ntot=1.8 with ntot=3]={1}, in hist entries can pile up so that ntot=2, or 3, and (ntot>nstart && ntot<=nstop) is skipped;
    //2)TruncatedMean(hist) out of range (only for Q0, not Q1).
    
    if(gmpv>EPSILON && kpass){ 
      if(tag.Contains("T0")){
        p0 = gmpv;
      }
      else{
        p0 = 1/gmpv;
      }
      //printf("outcalibobj%s %d %15.6e\n", tag.Data(), ibin, p0);
    }

    (*( (TVectorD*)cobj0->At(0) ))[ibin] = p0;

    //--- save optional record
    if(p0>EPSILON && countDet && countSSL){
      const Int_t idet = ToDetector(ibin);
      (*countDet)[idet]=1;
      
      const Int_t isector = ToSector(ibin);
      const Int_t istack = ToStack(ibin);
      const Int_t ilayer = ToLayer(ibin);
      TVectorD * vecsectorstack = (TVectorD*)countSSL->At(istack*18+isector);
      (*vecsectorstack)[ilayer]=1;
    }
    
    if(calibStream){
      (*calibStream)<<tag<<
        "run="<<run<<
        "p0="<<p0<<
        
        "nor="<<gnor<<
        "mpv="<<gmpv<<
        "wid="<<gwid<<
        "res="<<gres<<
        "gtrdphmean="<<gtrdphmean<<
        
        "ibin="<<ibin<<
        "\n";
    }
  }
  
  //----------------------------------------
  //               Status Report
  //----------------------------------------
  if(countDet && countSSL){
    TVectorD count2Dstack(90);
    for(Int_t ii=0; ii<90; ii++){
      TVectorD * vecsectorstack = (TVectorD*)countSSL->At(ii);
      const Int_t nlayer = (Int_t)vecsectorstack->Sum();
      if(nlayer==6){
        count2Dstack[ii]=1;
      }
    }

    printf("\nAliTRDdEdxUtils::GetCalibObj Summary run: %d name: %s entries: %.0f ndetector: %03.0f n2dstack %02.0f\n\n", run, hname.Data(), hh->GetEntries(), countDet->Sum(), count2Dstack.Sum());
  }

  //----------------------------------------
  //               Clean Up
  //----------------------------------------
  
  TH1D **hhs[]={&hnor, &hmpv, &hwid, &hres, &htrdphmean};
  const Int_t nhh=sizeof(hhs)/sizeof(TH1D**);
  for(Int_t ihh=0; ihh<nhh; ihh++){
    if(!lout){
      delete (*hhs[ihh]);
    }
  }
  
  delete countDet;
  delete countSSL;

  //----------------------------------------

  return cobj0;
}

//===================================================================================
//                                   dEdx calculation
//===================================================================================
Double_t AliTRDdEdxUtils::ToyCook(const Bool_t kinvq, Int_t &ncluster, TVectorD *arrayQ, TVectorD *arrayX, const TObjArray *cobj)
{
  //
  //template for cookdedx
  //
  if(cobj){
    if(arrayQ && arrayX){
      ncluster = ApplyCalib(ncluster, arrayQ, arrayX, cobj);
    }
    else{
      printf("AliTRDdEdxUtils::ToyCook arrayQ arrayX null, applycalib can not be applied!\n"); exit(1);
    }
  }

  Double_t lowFrac =-999, highFrac = -999;
  if(kinvq){
    lowFrac = Q1Frac(); highFrac = 0.99;
  }
  else{
    lowFrac = 0.01; highFrac = Q0Frac();
  }

  Double_t meanQ = TruncatedMean(ncluster, arrayQ->GetMatrixArray(), lowFrac, highFrac);
  if(kinvq){
    if(meanQ>EPSILON){
      meanQ = 1/meanQ;
    }
  }

  return meanQ;
}

Double_t AliTRDdEdxUtils::CombineddEdx(const Bool_t kinvq, Int_t &concls, TVectorD *coarrayQ, TVectorD *coarrayX, const Int_t tpcncls, const TVectorD *tpcarrayQ, const TVectorD *tpcarrayX, const Int_t trdncls, const TVectorD *trdarrayQ, const TVectorD *trdarrayX)
{
  //
  //combine tpc and trd dedx
  //

  for(Int_t iq=0; iq<tpcncls; iq++){
    (*coarrayQ)[iq]=(*tpcarrayQ)[iq];
    if(tpcarrayX && trdarrayX && coarrayX){
      (*coarrayX)[iq]=(*tpcarrayX)[iq];
    }
  }
  for(Int_t iq=0; iq<trdncls; iq++){
    (*coarrayQ)[tpcncls+iq]=(*trdarrayQ)[iq];
    if(tpcarrayX && trdarrayX && coarrayX){
      (*coarrayX)[tpcncls+iq]=159+(*trdarrayX)[iq];
    }
  }

  concls=trdncls+tpcncls;

  const Double_t coQ = ToyCook(kinvq, concls, coarrayQ, coarrayX);

  return coQ;
}


//===================================================================================
//                                   dEdx Getter and Setter
//===================================================================================
Double_t AliTRDdEdxUtils::GetAngularCorrection(const AliTRDseedV1 *seed)
{
  //
  //return angular normalization factor
  //

  return TMath::Sqrt(1+seed->GetYref(1)*seed->GetYref(1)+seed->GetZref(1)*seed->GetZref(1));
}

Double_t AliTRDdEdxUtils::GetPadGain(const Int_t det, const Int_t icol, const Int_t irow)
{
  //
  //get pad calibration
  //
  AliTRDcalibDB* calibration = AliTRDcalibDB::Instance();
  if (!calibration) {
    printf("AliTRDdEdxUtils::GetPadCalib No AliTRDcalibDB instance available\n"); exit(1);
  }
  AliTRDCalROC * calGainFactorROC = calibration->GetGainFactorROC(det);
  if(!calGainFactorROC){
    printf("AliTRDdEdxUtils::GetPadCalib no calGainFactorROC!\n"); exit(1);
  }

  Double_t padgain = -999;
  if( icol >= 0 && 
      icol < calGainFactorROC->GetNcols() && 
      irow >=0 && 
      irow < calGainFactorROC->GetNrows()){
    padgain = calGainFactorROC->GetValue(icol, irow);
    if(padgain<EPSILON){
      printf("AliTRDdEdxUtils::GetPadGain padgain 0! %f %f -- %d %d %d -- %d %d\n", padgain, EPSILON, det, icol, irow, calGainFactorROC->GetNcols(), calGainFactorROC->GetNrows()); exit(1);
    }
  }
  else{
    //printf("\nAliTRDdEdxUtils::GetPadGain warning!! indices out of range %d %d %d -- %d %d\n\n", det, icol, irow, calGainFactorROC->GetNcols(), calGainFactorROC->GetNrows() );  
  }

  return padgain;
}

Double_t AliTRDdEdxUtils::GetRNDClusterQ(AliTRDcluster *cl)
{
  //
  //get cluter q from GetRawQ, apply baseline and Kr pad-calibration
  //

  const Int_t det     = cl->GetDetector();
  const Int_t pad3col = cl->GetPadCol();
  const Int_t padrow  = cl->GetPadRow();

  const Double_t baseline = 10;

  Double_t rndqsum = 0;
  for(Int_t ii=0; ii<7; ii++){
    if(cl->GetSignals()[ii] < EPSILON){//bad pad marked by electronics
      continue;
    }

    const Int_t icol = pad3col+(ii-3);
    const Double_t padgain = GetPadGain(det, icol, padrow);
    if(padgain<0){//indices out of range, pad3col near boundary case
      continue;
    }

    const Double_t rndsignal = (cl->GetSignals()[ii] - baseline)/(IsPadGainOn()? padgain : 1);

    //sum it anyway even if signal below baseline, as long as the total is positive
    rndqsum += rndsignal;
  }

  return rndqsum;
}

Double_t AliTRDdEdxUtils::GetClusterQ(const Bool_t kinvq, const AliTRDseedV1 * seed, const Int_t itb)
{
  //
  //get cluster charge
  //
  Double_t dq = 0;
  AliTRDcluster *cl = 0x0;
      
  //GetRNDClusterQ(cl)>0 ensures that the total sum of q is above baseline*NsignalPhysical. 
  cl = seed->GetClusters(itb);                    if(cl && GetRNDClusterQ(cl)>0 ) dq+= GetRNDClusterQ(cl);//cl->GetRawQ();
  cl = seed->GetClusters(itb+AliTRDseedV1::kNtb); if(cl && GetRNDClusterQ(cl)>0 ) dq+= GetRNDClusterQ(cl);//cl->GetRawQ();

  dq /= GetAngularCorrection(seed);
  
  dq /= QScale();
      
  if(kinvq){
    if(dq>EPSILON){
      dq = 1/dq;
    }
  }

  return dq;
}

Int_t AliTRDdEdxUtils::GetArrayClusterQ(const Bool_t kinvq, TVectorD *arrayQ, TVectorD *arrayX, const AliTRDtrackV1 *trdtrack, Int_t timeBin0, Int_t timeBin1, Int_t tstep)
{
  //
  //return nclustter
  //(if kinvq, return 1/q array), size of array must be larger than 31*6
  //
  if(!arrayQ || arrayQ->GetNrows()< (AliTRDseedV1::kNtb*AliTRDtrackV1::kNplane)){
    printf("AliTRDdEdxUtils::GetArrayClusterQ error arrayQ null or size too small! %d\n", arrayQ? arrayQ->GetNrows() : -999); exit(1);
  }
  if(!arrayX || arrayX->GetNrows()< (AliTRDseedV1::kNtb*AliTRDtrackV1::kNplane)){
    printf("AliTRDdEdxUtils::GetArrayClusterQ error arrayX null or size too small! %d\n", arrayX? arrayX->GetNrows() : -999); exit(1);
  }

  const Int_t mintb = 0;
  const Int_t maxtb = AliTRDseedV1::kNtb-1;
  if(timeBin0<mintb) timeBin0=mintb;
  if(timeBin1>maxtb) timeBin1=maxtb;
  if(tstep<=0) tstep=1;

  //============
  Int_t tbN=0;
  Double_t tbQ[200];
  Int_t tbBin[200];
    
  for(Int_t ichamber=0; ichamber < AliTRDtrackV1::kNplane; ichamber++){
    const AliTRDseedV1 * seed = trdtrack->GetTracklet(ichamber);
    if(!seed)
      continue;
    
    const Int_t det = seed->GetDetector();

    for(Int_t itb=timeBin0; itb<=timeBin1; itb+=tstep){
      const Double_t dq = GetClusterQ(kinvq, seed, itb);
      if(dq<EPSILON)
        continue;

      const Int_t gtb = det * AliTRDseedV1::kNtb + itb;

      tbQ[tbN]=dq;
      tbBin[tbN]=gtb;
      tbN++;
    }
  }

  Int_t ncls = 0;
  for(Int_t iq=0; iq<tbN; iq++){
    if(tbQ[iq]<EPSILON)
      continue;

    (*arrayQ)[ncls] = tbQ[iq];
    (*arrayX)[ncls] = tbBin[iq];

    ncls++;
  }

  static Int_t kprint = 100;
  if(kprint<0){
    printf("\nAliTRDdEdxUtils::GetArrayClusterQ raw cluster-Q\n");
    for(Int_t iq=0; iq<ncls; iq++){
      const Int_t ichamber =  ToLayer((*arrayX)[iq]);
      const AliTRDseedV1 * seed = trdtrack->GetTracklet(ichamber);
      if(!seed){
        printf("error seed null!!\n"); exit(1);
      }
      const Double_t rawq =  (*arrayQ)[iq] * 45. * GetAngularCorrection(seed);
      printf("esdid=%d; chamber=%d; timebin=%d; rawq= %.3f; myq[%d]= %e;\n", trdtrack->GetESDid(), ichamber, ToTimeBin((*arrayX)[iq]), rawq, iq, (*arrayQ)[iq]);
    }
    printf("\n");
  }
  kprint++;

  return ncls;
}

Int_t AliTRDdEdxUtils::UpdateArrayX(const Int_t ncls, TVectorD* arrayX)
{
  //
  //arrayX det*Ntb+itb -> itb
  //

  TVectorD countChamber(6);
  for(Int_t ii = 0; ii<ncls; ii++){
    const Int_t xx = (Int_t)(*arrayX)[ii];
    const Int_t idet = ToDetector(xx);
    
    const Double_t ich = AliTRDgeometry::GetLayer(idet);
    const Double_t itb = ToTimeBin(xx);
    (*arrayX)[ii] = ich*AliTRDseedV1::kNtb+itb;

    countChamber[ich] = 1;
  }

  const Double_t nch = countChamber.Sum();
  return (Int_t) nch;
}

void AliTRDdEdxUtils::SetChamberQT(const AliTRDtrackV1 *trdtrack, const Int_t kcalib, THnF * hgain, THnF * ht0, THnF * hvd)
{
  //
  //CookdEdx at TRD track level, use chamber info, related calibrations: chamber-gain; T0, Vd based on raw PH distribution
  //

  static Int_t kprint = 100;

  fgNchamber = 0;
  for(Int_t ichamber=0; ichamber < AliTRDtrackV1::kNplane; ichamber++){
    //initialize output, default values: 0, so that summation and weighting will automatically discard default quantities
    fgChamberQ[ichamber] = fgChamberTmean[ichamber] = 0;

    const AliTRDseedV1 *seed = trdtrack->GetTracklet(ichamber);
    if (!seed) 
      continue;

    const Int_t idet = seed->GetDetector();

    //-------------------------------------------------------------------------

    Double_t qsum = 0, qtsum = 0, w2sum = 0;
    for(Int_t itb=0; itb<AliTRDseedV1::kNtb; itb++){
      const Double_t dq = GetClusterQ(0, seed, itb);
      if(dq<EPSILON)
        continue;

      qsum += dq;
      qtsum += dq*itb; 
      w2sum += dq*itb*itb;
    }
    if(qsum<EPSILON)
      continue;

    //-------------------------------------------------------------------------

    Double_t tbm, tbr = 0;
    tbm = GetMeanRMS(qsum, qtsum, w2sum, &tbr);

    qsum /= 1.25e3/45.;

    if(hgain){ 
      const Double_t var[]={idet, qsum};      
      hgain->Fill(var);
    }
    if(ht0){
      const Double_t var[]={idet, tbm};
      ht0->Fill(var);
    }
    if(hvd){
      const Double_t var[]={idet, tbr};
      hvd->Fill(var);
    }

    Double_t gain = 1, t0 = 0, vd = 1;
    if(kcalib){
      if(!fgObjGain) {printf("AliTRDdEdxUtils::SetChamberQT error Gain array null!!\n"); exit(1);}
      if(!  fgObjT0) {printf("AliTRDdEdxUtils::SetChamberQT error T0   array null!!\n"); exit(1);}
      if(!  fgObjVd) {printf("AliTRDdEdxUtils::SetChamberQT error Vd   array null!!\n"); exit(1);}

      const TVectorD * gainvec = (TVectorD*) fgObjGain->At(0); gain = (*gainvec)[idet];
      const TVectorD *   t0vec = (TVectorD*)   fgObjT0->At(0);   t0 = (*  t0vec)[idet];
      const TVectorD *   vdvec = (TVectorD*)   fgObjVd->At(0);   vd = (*  vdvec)[idet];
    }
    if(kprint<0){
      printf("\nAliTRDdEdxUtils::CookdEdxV2\n");
      printf("idet = %d;\n", idet);
      printf("gain = %15f; t0 = %15f; vd = %15f;\n", gain, t0, vd);
      printf("\n");
    }

    qsum *= gain;
    tbm = (tbm-t0)*vd;

    if(qsum<EPSILON)
      continue;

    //-------------------------------------------------------------------------

    //should have non-zero value, initialized with default 0 (except for calibrated tbm, may be very close to 0)
    fgChamberQ[ichamber]  = qsum;
    fgChamberTmean[ichamber] = tbm;  
    fgNchamber++;
  }
  
  if(kprint<0){
    printf("\nAliTRDdEdxUtils::CookdEdxV2 summary:\n");

    printf("\nfgNchamber = %d\n", fgNchamber);
    for(Int_t ich=0; ich<AliTRDtrackV1::kNplane; ich++){
      printf("fgChamberTmean[%d] = %15f; fgChamberQ[%d] = %15f;\n", ich, fgChamberTmean[ich], ich, fgChamberQ[ich]);
    }
  }

  fgTrackTmean = -999;
  if(fgNchamber){
    fgTrackTmean = 0;
    for(Int_t ich=0; ich<AliTRDtrackV1::kNplane; ich++){
      fgTrackTmean += fgChamberTmean[ich];
    }
    fgTrackTmean /= fgNchamber;
  }

  if(kprint<0){
    printf("\nAliTRDdEdxUtils::CookdEdxV2\n");
    printf("GetTrackTmean() %15f\n", GetTrackTmean());
    printf("\n");
  }
  kprint++;

  return;
}


//===================================================================================
//                                 dEdx Parameterization
//===================================================================================

Double_t AliTRDdEdxUtils::Q0MeanTRDpp(const Double_t bg)
{
  //
  //truncated Mean Q_{xx} in TRD
  //
 
  Double_t par[8];
  //03132012161150
  //opt: ppQ0
par[0]=   2.397001e-01;
par[1]=   1.334697e+00;
par[2]=   6.967470e+00;
par[3]=   9.055289e-02;
par[4]=   9.388760e+00;
par[5]=   9.452742e-04;
par[6]=  -1.866091e+00;
par[7]=   1.403545e+00;

  ///u/xlu/.task/CommondEdx/myAnaData/Optimum/check11/Merged/LHC10e_plot/Separation/see2.log:hhtype2Q0b2c2 scale        0.428543 at ltbg        0.650000
  return   0.428543*MeandEdxTR(&bg, par);
}

Double_t AliTRDdEdxUtils::Q0MeanTRDPbPb(const Double_t bg)
{
  //
  //truncated Mean Q_{xx} in TRD
  //
 
  Double_t par[8];

  //03132012161259
  //opt: PbPbQ0
par[0]=   1.844912e-01;
par[1]=   2.509702e+00;
par[2]=   6.744031e+00;
par[3]=   7.355123e-02;
par[4]=   1.166023e+01;
par[5]=   1.736186e-04;
par[6]=  -1.716063e+00;
par[7]=   1.611366e+00;

  ///u/xlu/.task/CommondEdx/myAnaData/Optimum/check11/Merged/LHC10e_plot/Separation/see4.log:hhtype4Q0b2c2 scale        0.460994 at ltbg        0.650000  
  return   0.460994*MeandEdxTR(&bg, par);
}

Double_t AliTRDdEdxUtils::Q1MeanTRDpp(const Double_t bg)
{
  //
  //truncated Mean 1/(1/Q)_{xx} in TRD
  //
 
  Double_t par[8];

  //So 4. Mär 13:30:51 CET 2012
  //opt: trdppQ1
  par[0]=   2.434646e-01;
  par[1]=   1.400211e+00;
  par[2]=   6.937471e+00;
  par[3]=   7.758118e-02;
  par[4]=   1.097372e+01;
  par[5]=   4.297518e-04;
  par[6]=  -1.806266e+00;
  par[7]=   1.543811e+00;

  //hhtype2Q1b2c2 scale        0.418629 at ltbg        0.650000

  return  0.418629*MeandEdxTR(&bg, par);
}

Double_t AliTRDdEdxUtils::Q1MeanTRDPbPb(const Double_t bg)
{
  //
  //truncated Mean 1/(1/Q)_{xx} in TRD
  //
 
  Double_t par[8];

  //So 4. Mär 13:30:52 CET 2012
  //opt: trdPbPbQ1
  par[0]=   2.193660e-01;
  par[1]=   2.051864e+00;
  par[2]=   6.825112e+00;
  par[3]=   6.151693e-02;
  par[4]=   1.390343e+01;
  par[5]=   6.010032e-05;
  par[6]=  -1.676324e+00;
  par[7]=   1.838873e+00;

  //hhtype4Q1b2c2 scale        0.457988 at ltbg        0.650000

  return  0.457988*MeandEdxTR(&bg, par);
}

Double_t AliTRDdEdxUtils::QMeanTPC(const Double_t bg)
{
  //
  //bethe bloch in TPC
  //

  Double_t par[5];
  //Mi 15. Feb 14:48:05 CET 2012
  //train_2012-02-13_1214.12001, tpcsig
  par[0]=       4.401269;
  par[1]=       9.725370;
  par[2]=       0.000178;
  par[3]=       1.904962;
  par[4]=       1.426576;

  return MeandEdx(&bg, par);
}

Double_t AliTRDdEdxUtils::MeandEdxTR(const Double_t * xx, const Double_t * pin)
{
  //
  //ALEPH+LOGISTIC parametrization for dEdx+TR, in unit of MIP
  //npar = 8 = 3+5
  //
  Double_t ptr[4]={0};
  for(int ii=0; ii<3; ii++){
    ptr[ii+1]=pin[ii];
  }
  return MeanTR(xx,ptr) + MeandEdx(xx,&(pin[3]));
}

Double_t AliTRDdEdxUtils::MeanTR(const Double_t * xx, const Double_t * par)
{
  //
  //ALEPH+LOGISTIC parametrization for dEdx+TR, in unit of MIP
  //npar = 4
  //

  const Double_t bg = xx[0];
  const Double_t gamma = sqrt(1+bg*bg);

  const Double_t p0 = TMath::Abs(par[1]);
  const Double_t p1 = TMath::Abs(par[2]);
  const Double_t p2 = TMath::Abs(par[3]);

  const Double_t zz = TMath::Log(gamma);
  const Double_t tryield = p0/( 1 + TMath::Exp(-p1*(zz-p2)) );

  return par[0]+tryield;
}

Double_t AliTRDdEdxUtils::MeandEdx(const Double_t * xx, const Double_t * par)
{
  //
  //ALEPH parametrization for dEdx
  //npar = 5
  //

  const Double_t bg = xx[0];
  const Double_t beta = bg/TMath::Sqrt(1.+ bg*bg);

  const Double_t p0 = TMath::Abs(par[0]);
  const Double_t p1 = TMath::Abs(par[1]);
  const Double_t p2 = TMath::Abs(par[2]);
  const Double_t p3 = TMath::Abs(par[3]);
  const Double_t p4 = TMath::Abs(par[4]);

  const Double_t aa = TMath::Power(beta, p3);
  const Double_t bb = TMath::Log( p2 + TMath::Power(1./bg, p4) );

  //printf("test----- %f %f -- %f %f %f %f %f --- %f %f %f\n", bg, beta, p0, p1, p2, p3, p4, p0/aa, aa, bb);

  return (p1-aa-bb)*p0/aa;
}

Double_t AliTRDdEdxUtils::ToLogx(FFunc func, const Double_t * xx, const Double_t * par)
{
  //
  //f(x)-> f(y) with y=log10(x)
  //
  const Double_t x2[]={TMath::Power(10, xx[0])};
  return func(x2, par);
}

//===================================================================================
//                                Detector, Data and Control Constant 
//===================================================================================
Int_t  AliTRDdEdxUtils::ToDetector(const Int_t gtb)
{
  //
  //gtb = det*Ntb+itb
  //
  return gtb/AliTRDseedV1::kNtb;
}

Int_t  AliTRDdEdxUtils::ToTimeBin(const Int_t gtb)
{ 
  //
  //gtb = det*Ntb+itb
  //
  return gtb%AliTRDseedV1::kNtb;
}

Int_t  AliTRDdEdxUtils::ToSector(const Int_t gtb)
{
  //
  //return sector
  //
  return AliTRDgeometry::GetSector(ToDetector(gtb));
}

Int_t  AliTRDdEdxUtils::ToStack(const Int_t gtb)
{
  //
  //return stack
  //
  return AliTRDgeometry::GetStack(ToDetector(gtb));
}

Int_t  AliTRDdEdxUtils::ToLayer(const Int_t gtb)
{
  //
  //return layer
  //
  return AliTRDgeometry::GetLayer(ToDetector(gtb));
}

TString AliTRDdEdxUtils::GetRunType(const Int_t run)
{
  //
  //return run type
  //

  TString type;
  if(run>=121527 && run<= 126460)//LHC10d
    type="pp2010LHC10d";
  else if(run>=126461 && run<= 130930)//LHC10e
    type="pp2010LHC10e";
  else if(run>=136782 && run <= 139846)//LHC10h
    type="PbPb2010LHC10h";
  else if(run>= 143224 && run<= 143237)//2011Feb
    type="cosmic2011Feb";
  else if(run>= 150587 && run<= 154930){
    type="cosmic2011MayJun";

    TString runstr(Form("%d",run));
    const TString listrun1kg("154601 154602 154629 154634 154636 154639 154643");
    if(listrun1kg.Contains(runstr)){
      type+="1kG";
    }
    else{
      type+="5kG";
    }      
  }
  else{
    type="unknown";
  }

  type.ToUpper();
  return type;
}

void AliTRDdEdxUtils::PrintControl()
{
  //
  //print out control variable
  //
  printf("\nAliTRDdEdxUtils::PrintControl Q0Frac %.1f, Q1Frac %.1f, TimeBinCountCut %.2f, CalibTPCnclsCut %d, IsExBOn %d, IsPadGainOn %d, QScale %.2f\n\n", Q0Frac(), Q1Frac(), TimeBinCountCut(), CalibTPCnclsCut(), IsExBOn(), IsPadGainOn(), QScale());
}
//===================================================================================
//===================================================================================
