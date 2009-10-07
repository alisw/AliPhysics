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

/* $Id: AliTRDpidRefMakerLQ.cxx 34163 2009-08-07 11:28:51Z cblume $ */

///////////////////////////////////////////////////////////////////////////////
//
//
//  TRD calibration class for building reference data for PID
//  - 2D reference histograms (responsible A.Bercuci) 
//  - 3D reference histograms (not yet implemented) (responsible A.Bercuci)
//  - Neural Network (responsible A.Wilk)
//
//   Origin
//   Alex Bercuci  (A.Bercuci@gsi.de)
//
///////////////////////////////////////////////////////////////////////////////

#include <TFile.h>
#include <TROOT.h>
#include <TTree.h>
#include <TMath.h>
#include <TEventList.h>
#include <TH2D.h>
#include <TH2I.h>
#include <TPrincipal.h>
#include <TVector3.h>
#include <TLinearFitter.h>
#include <TCanvas.h>
#include <TEllipse.h>
#include <TMarker.h>

#include "AliLog.h"
#include "../../STAT/TKDPDF.h"
#include "AliTRDpidRefMakerLQ.h"
#include "../Cal/AliTRDCalPID.h"
#include "AliTRDseedV1.h"
#include "AliTRDcalibDB.h"
#include "AliTRDgeometry.h"

ClassImp(AliTRDpidRefMakerLQ)

//__________________________________________________________________
AliTRDpidRefMakerLQ::AliTRDpidRefMakerLQ()
  :AliTRDpidRefMaker("pidRefMakerLQ", "PID(LQ) Reference Maker")
  ,fPbin(-1)
  ,fSbin(-1)
  ,fResults(0x0)
{
  //
  // AliTRDpidRefMakerLQ default constructor
  //
}

//__________________________________________________________________
AliTRDpidRefMakerLQ::~AliTRDpidRefMakerLQ()
{
  //
  // AliTRDCalPIDQRef destructor
  //
  if(fResults) {
    fResults->Delete();
    delete fResults;
  }
}

//________________________________________________________________________
void AliTRDpidRefMakerLQ::CreateOutputObjects()
{
  // Create histograms
  // Called once

  AliTRDpidRefMaker::CreateOutputObjects();

  // save dE/dx references
  TH2 *h2 = 0x0;
  for(Int_t ip=AliTRDCalPID::kNMom; ip--;){ 
    TObjArray *arr = new TObjArray(AliPID::kSPECIES);
    arr->SetName(Form("Pbin%02d", ip));
    for(Int_t is=AliPID::kSPECIES; is--;) {
      h2 = new TH2I(Form("h%s%d", AliPID::ParticleShortName(is), ip), Form("%s ref. dEdx @ Pbin[%d]", AliPID::ParticleName(is), ip), 50, 5., 10., 50, 5., 10.);
      h2->GetXaxis()->SetTitle("log(dE/dx_{am}) [au]");
      h2->GetYaxis()->SetTitle("log(dE/dx_{dr}) [au]");
      h2->GetZaxis()->SetTitle("#");
      arr->AddAt(h2, is);
    }
    fContainer->AddAt(arr, 1+ip);
  }

  fData = new TTree(GetName(), Form("Reference data for %s", GetName()));
  fData->Branch("s", &fSbin, "l/b");
  fData->Branch("p", &fPbin, "p/b");
  fData->Branch("dEdx" , fdEdx , Form("dEdx[%d]/F", GetNslices()));
}


//________________________________________________________________________
Float_t* AliTRDpidRefMakerLQ::CookdEdx(AliTRDseedV1 *trklt)
{
  trklt->CookdEdx(AliTRDpidUtil::kLQslices);
  const Float_t *dedx = trklt->GetdEdx();
  if(dedx[0]+dedx[1] <= 0.) return 0x0;
  if(dedx[2] <= 0.) return 0x0;

  fdEdx[0] = TMath::Log(dedx[0]+dedx[1]);
  fdEdx[1] = TMath::Log(dedx[2]);
  return fdEdx;
}

#include "../Cal/AliTRDCalPIDLQ.h"

//__________________________________________________________________
TObject* AliTRDpidRefMakerLQ::GetOCDBEntry(Option_t *opt)
{
  TDirectoryFile *d = 0x0;
  if(!TFile::Open(Form("TRD.Calib%s.root", GetName()))) return 0x0;
  if(!(d=(TDirectoryFile*)gFile->Get(Form("PDF_%s", opt)))) return 0x0;
  AliTRDCalPIDLQ *cal = new AliTRDCalPIDLQ("pidLQ", "LQ TRD PID object");
  cal->LoadPDF(d);
  return cal;
}

//__________________________________________________________________
Bool_t AliTRDpidRefMakerLQ::GetRefFigure(Int_t ifig)
{
  if(ifig<0 || ifig>AliTRDCalPID::kNMom-1){ 
    AliError("Ref fig requested outside definition.");
    return kFALSE;
  }
  if(!fResults){ 
    AliError("No results processed.");
    return kFALSE;
  }
  if(!gPad){
    AliWarning("Please provide a canvas to draw results.");
    return kFALSE;
  }

  TObjArray *arr = (TObjArray*)fResults->At(ifig);
  gPad->Divide(3, 2, 1.e-5, 1.e-5); 
  TList *l=gPad->GetListOfPrimitives(); 
  for(Int_t is=0; is<AliPID::kSPECIES; is++){
    ((TVirtualPad*)l->At(is))->cd();
    ((TH2*)arr->At(is))->Draw("cont4z");
  }

  return kTRUE;
}


//________________________________________________________________________
void AliTRDpidRefMakerLQ::Fill()
{
  // get particle type
  fSbin = TMath::LocMax(AliPID::kSPECIES, fPID);
  // get momentum bin
  fPbin = AliTRDpidUtil::GetMomentumBin(fP);
  // fill data
  fData->Fill();
  // fill monitor
  ((TH2*)fContainer->At(0))->Fill(fSbin, fPbin);
  TH2* h2 = (TH2*)((TObjArray*)fContainer->At(1+fPbin))->At(fSbin);
  h2->Fill(fdEdx[0], fdEdx[1]);
}

//________________________________________________________________________
Bool_t AliTRDpidRefMakerLQ::PostProcess()
{
// Analyse merged dedx = f(p) distributions.
//   - select momentum - species bins
//   - rotate to principal components
//   - locally interpolate with TKDPDF
//   - save interpolation to monitoring histograms
//   - write pdf to file for loading to OCDB
// 


  TFile *fCalib = TFile::Open(Form("TRD.Calib%s.root", GetName()), "update");
  fData = dynamic_cast<TTree*>(gFile->Get(GetName()));
  if (!fData) {
    AliError("Tree not available");
    return kFALSE;
  }
  TDatime d;
  TDirectoryFile *pdfs = new TDirectoryFile(Form("PDF_%d", d.GetDate()), "PDFs for LQ TRD-PID", "", gFile);
  pdfs->Write();
  AliDebug(2, Form("Data[%d]", fData->GetEntries()));

  // save a volatile PDF representation
  gROOT->cd();
  TH2 *h2 = 0x0;
  fResults = new TObjArray(AliTRDCalPID::kNMom);
  for(Int_t ip=AliTRDCalPID::kNMom; ip--;){ 
    TObjArray *arr = new TObjArray(AliPID::kSPECIES);
    arr->SetName(Form("Pbin%02d", ip));
    for(Int_t is=AliPID::kSPECIES; is--;) {
      h2 = new TH2D(Form("i%s%d", AliPID::ParticleShortName(is), ip), Form("Interpolated %s dEdx @ Pbin[%d]", AliPID::ParticleName(is), ip), 50, -3., 3., 50, -3., 3.);
      h2->GetXaxis()->SetTitle("log(dE/dx^{*}_{am}) [au]");
      h2->GetYaxis()->SetTitle("log(dE/dx^{*}_{dr}) [au]");
      h2->GetZaxis()->SetTitle("#");
      arr->AddAt(h2, is);
    }
    fResults->AddAt(arr, ip);
  }
  pdfs->cd();

  TCanvas *cc = new TCanvas("cc", "", 500, 500);

  Float_t *data[] = {0x0, 0x0};
  TPrincipal principal(2, "ND");
  for(Int_t ip=AliTRDCalPID::kNMom; ip--; ){ 
    for(Int_t is=AliPID::kSPECIES; is--;) {
      principal.Clear();
      Int_t n = fData->Draw("dEdx[0]:dEdx[1]", Form("p==%d&&s==%d", ip, is), "goff");

      // estimate bucket statistics
      Int_t nb(kMinBuckets), // number of buckets
            ns(Int_t(Float_t(n)/nb));    //statistics/bucket
            
// if(Float_t(n)/nb < 220.) ns = 200; // 7% stat error
//       else if(Float_t(n)/nb < 420.) ns = 400; // 5% stat error

      AliDebug(2, Form("pBin[%d] sBin[%d] N[%d] n[%d] nb[%d]", ip, is, n, ns, nb));
      if(ns<Int_t(kMinStat)){
        AliWarning(Form("Not enough entries [%d] for %s[%d].", ns, AliPID::ParticleShortName(is), ip));
        continue;
      }
      // allocate storage
      data[0] = new Float_t[n];data[1] = new Float_t[n];

      // fill and make principal
      Double_t *v1 = fData->GetV1(), *v2 = fData->GetV2();
      while(n--){
        Double_t dedx[] = {v1[n], v2[n]};
        principal.AddRow(dedx);
      }
      principal.MakePrincipals();

//       // calculate covariance ellipse
//       eValues  = principal.GetEigenValues();
//       x0  = 0.;
//       y0  = 0.;
//       rx  = 3.5*sqrt((*eValues)[0]);
//       ry  = 3.5*sqrt((*eValues)[1]);

      // rotate to principal axis
      const Double_t *xx; Double_t rxy[2];
      while((xx = principal.GetRow(++n))){
        principal.X2P(xx, rxy);
        data[0][n]=rxy[0]; data[1][n]=rxy[1];
      }

      // build PDF
      TKDPDF pdf(n, 2, ns, data);
      pdf.SetStore();
      pdf.SetAlpha(5.);
      //pdf.GetStatus();
      Float_t *c, v, ve; Double_t r, e;
      for(Int_t in=pdf.GetNTNodes(); in--;){
        pdf.GetCOGPoint(in, c, v, ve);
        rxy[0] = (Double_t)c[0];rxy[1] = (Double_t)c[1];
        pdf.Eval(rxy, r, e, kTRUE);
      }
      // visual on-line monitoring
      pdf.DrawBins(0,1,-5.,5.,-8.,8.);cc->Modified(); cc->Update(); cc->SaveAs(Form("pdf_%s%02d.gif", AliPID::ParticleShortName(is), ip));
      cc->SaveAs(Form("%s_%s%02d.gif", GetName(), AliPID::ParticleShortName(is), ip));

      // save a discretization of the PDF for result monitoring
      TH2 *h2s = (TH2D*)((TObjArray*)fResults->At(ip))->At(is);
      TAxis *ax = h2s->GetXaxis(), *ay = h2s->GetYaxis();
      for(int ix=1; ix<=ax->GetNbins(); ix++){
        rxy[0] = ax->GetBinCenter(ix);
        for(int iy=1; iy<=ay->GetNbins(); iy++){
          rxy[1] = ay->GetBinCenter(iy);
      
          Double_t r,e;
          pdf.Eval(rxy, r, e);
          if(r<0. || e/r>.15) continue; // 15% relative error
          //printf("x[%2d] x[%2d] r[%f] e[%f]\n", ix, iy, r, e);
          h2s->SetBinContent(ix, iy, r);
        }
      }

      // write results to output array
      pdf.Write(Form("%s[%d]", AliPID::ParticleShortName(is), ip));

      delete [] data[0]; delete [] data[1];
    }
  }
  pdfs->Write();
  fCalib->Close(); delete fCalib;

  return kTRUE; // testing protection
}


