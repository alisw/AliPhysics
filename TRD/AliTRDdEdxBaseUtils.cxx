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
// TRD dEdx base utils
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

#include "AliTRDdEdxBaseUtils.h"

#define EPSILON 1e-12

Double_t AliTRDdEdxBaseUtils::fgQ0Frac = 0.3;
Double_t AliTRDdEdxBaseUtils::fgQ1Frac = 0.5;
Double_t AliTRDdEdxBaseUtils::fgTimeBinCountCut = 0.0; 
Int_t    AliTRDdEdxBaseUtils::fgCalibTPCnclsCut = 70;
Bool_t   AliTRDdEdxBaseUtils::fgExBOn = kTRUE; 
Bool_t   AliTRDdEdxBaseUtils::fgPadGainOn = kTRUE;
Double_t AliTRDdEdxBaseUtils::fgQScale = 45;

//===================================================================================
//                                   Math and Histogram
//===================================================================================
void AliTRDdEdxBaseUtils::BinLogX(TAxis *axis) 
{
  //
  // Method for the correct logarithmic binning of histograms
  // copied and modified from AliTPCcalibBase

  const Int_t bins = axis->GetNbins();

  const Double_t from = axis->GetXmin();
  const Double_t to = axis->GetXmax();
  if (from<EPSILON) return;
  Double_t *new_bins = new Double_t[bins + 1];

  new_bins[0] = from;
  const Double_t factor = pow(to/from, 1./bins);

  for (int i = 1; i <= bins; i++) {
   new_bins[i] = factor * new_bins[i-1];
  }
  axis->Set(bins, new_bins);
  delete [] new_bins;
}

void AliTRDdEdxBaseUtils::GetCDFCuts(const TH1D *hh, const Int_t ncut, Double_t cuts[], const Double_t cdfs[], const Double_t thres)
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
      printf("AliTRDdEdxBaseUtils::GetCDFCuts error cdfs[%d] %15f out of range!\n", ii, icdf); exit(1);
    }
    cuts[ii] = datas[id[Int_t(icdf*nsel)]];
  }
}

Double_t AliTRDdEdxBaseUtils::GetMeanRMS(const Double_t nn, const Double_t sum, const Double_t w2s, Double_t * grms, Double_t * gerr)
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
        printf("AliTRDdEdxBaseUtils::GetMeanRMS error negative sqrt argument!! %e -- %e %e %f\n", arg, w2s, sum, nn); exit(1);
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

Double_t AliTRDdEdxBaseUtils::TruncatedMean(const Int_t nx, const Double_t xdata[], const Double_t lowfrac, const Double_t highfrac, Double_t * grms, Double_t * gerr, Double_t *wws)
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

Double_t AliTRDdEdxBaseUtils::TruncatedMean(const TH1 *hh, const Double_t lowfrac, const Double_t highfrac, Double_t * grms, Double_t * gerr)
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
        printf("AliTRDdEdxBaseUtils::TruncatedMean (hist) error %e %e %d\n", bc, be, itmp); exit(1);
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
        printf("AliTRDdEdxBaseUtils::TruncatedMean (hist) error %e %e %d\n", bc, be, itmp); exit(1);
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
        printf("AliTRDdEdxBaseUtils::TruncatedMean warning hist range too small %s %f %f %d %d, %15f %15f %15f; nused w2s sum set to 0\n", hh->GetName(), hh->GetBinLowEdge(1), hh->GetBinLowEdge(itmp), itmp, hh->GetNbinsX(), hh->GetBinContent(hh->GetNbinsX())/hh->Integral(0,hh->GetNbinsX()+1), hh->GetBinContent(hh->GetNbinsX()), hh->Integral(0,hh->GetNbinsX()+1)); //exit(1);
        nused = 0;
        w2s = sum = 0;
      }
      break;
    }
  }

  return GetMeanRMS(nused, sum, w2s, grms, gerr);
}

void AliTRDdEdxBaseUtils::FitSlicesY(const TH2D *hh, TH1D *&hnor, TH1D *&hmpv, TH1D *&hwid, TH1D *&hres, const Double_t thres, const Double_t lowfrac, const Double_t highfrac)
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
          printf("AliTRDdEdxBaseUtils::FitSlicesY error %d %d %e %e\n", ix, iy, be, bc); exit(1);
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

Int_t AliTRDdEdxBaseUtils::GetNtracklet(const AliESDEvent *esd)
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

AliTRDtrackV1 * AliTRDdEdxBaseUtils::GetTRDtrackV1(const AliESDtrack * esdtrack)
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

Bool_t AliTRDdEdxBaseUtils::IsInSameStack(const AliTRDtrackV1 *trdtrack)
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

Bool_t AliTRDdEdxBaseUtils::GetFirstSectorStackMomentum(const AliTRDtrackV1 *trdtrack, Int_t & isec, Int_t & istk, Double_t & mom)
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
//                                Detector and Data Constant 
//===================================================================================
Int_t  AliTRDdEdxBaseUtils::ToDetector(const Int_t gtb)
{
  //
  //gtb = det*Ntb+itb
  //
  return gtb/AliTRDseedV1::kNtb;
}

Int_t  AliTRDdEdxBaseUtils::ToTimeBin(const Int_t gtb)
{ 
  //
  //gtb = det*Ntb+itb
  //
  return gtb%AliTRDseedV1::kNtb;
}

Int_t  AliTRDdEdxBaseUtils::ToSector(const Int_t gtb)
{
  //
  //return sector
  //
  return AliTRDgeometry::GetSector(ToDetector(gtb));
}

Int_t  AliTRDdEdxBaseUtils::ToStack(const Int_t gtb)
{
  //
  //return stack
  //
  return AliTRDgeometry::GetStack(ToDetector(gtb));
}

Int_t  AliTRDdEdxBaseUtils::ToLayer(const Int_t gtb)
{
  //
  //return layer
  //
  return AliTRDgeometry::GetLayer(ToDetector(gtb));
}

TString AliTRDdEdxBaseUtils::GetRunType(const Int_t run)
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

void AliTRDdEdxBaseUtils::PrintControl()
{
  //
  //print out control variable
  //
  printf("\nAliTRDdEdxBaseUtils::PrintControl Q0Frac %.1f, Q1Frac %.1f, TimeBinCountCut %.2f, CalibTPCnclsCut %d, IsExBOn %d, IsPadGainOn %d, QScale %.2f\n\n", Q0Frac(), Q1Frac(), TimeBinCountCut(), CalibTPCnclsCut(), IsExBOn(), IsPadGainOn(), QScale());
}

//===================================================================================
//                                 dEdx Parameterization
//===================================================================================

Double_t AliTRDdEdxBaseUtils::Q0MeanTRDpp(const Double_t bg)
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

Double_t AliTRDdEdxBaseUtils::Q0MeanTRDPbPb(const Double_t bg)
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

Double_t AliTRDdEdxBaseUtils::Q1MeanTRDpp(const Double_t bg)
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

Double_t AliTRDdEdxBaseUtils::Q1MeanTRDPbPb(const Double_t bg)
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

Double_t AliTRDdEdxBaseUtils::QMeanTPC(const Double_t bg)
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

Double_t AliTRDdEdxBaseUtils::MeandEdxTR(const Double_t * xx, const Double_t * pin)
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

Double_t AliTRDdEdxBaseUtils::MeanTR(const Double_t * xx, const Double_t * par)
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

Double_t AliTRDdEdxBaseUtils::MeandEdx(const Double_t * xx, const Double_t * par)
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

Double_t AliTRDdEdxBaseUtils::ToLogx(FFunc func, const Double_t * xx, const Double_t * par)
{
  //
  //f(x)-> f(y) with y=log10(x)
  //
  const Double_t x2[]={TMath::Power(10, xx[0])};
  return func(x2, par);
}

