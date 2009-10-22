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
//#include <TVector3.h>
//#include <TLinearFitter.h>
#include <TCanvas.h>
//#include <TEllipse.h>
//#include <TMarker.h>

#include "AliLog.h"
#include "../../STAT/TKDPDF.h"
#include "AliTRDpidRefMakerLQ.h"
#include "../Cal/AliTRDCalPID.h"
#include "../Cal/AliTRDCalPIDLQ.h"
#include "AliTRDseedV1.h"
#include "AliTRDcalibDB.h"
#include "AliTRDgeometry.h"

ClassImp(AliTRDpidRefMakerLQ)

//__________________________________________________________________
AliTRDpidRefMakerLQ::AliTRDpidRefMakerLQ()
  :AliTRDpidRefMaker("PIDrefMakerLQ", "PID(LQ) Reference Maker")
  ,fPIDbin(0xff)
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
      h2 = new TH2D(Form("h%s%d", AliPID::ParticleShortName(is), ip), Form("%s ref. dEdx @ Pbin[%d]", AliPID::ParticleName(is), ip), 50, 5., 10., 50, 5., 10.);
      h2->GetXaxis()->SetTitle("log(dE/dx_{am}) [au]");
      h2->GetYaxis()->SetTitle("log(dE/dx_{dr}) [au]");
      h2->GetZaxis()->SetTitle("#");
      arr->AddAt(h2, is);
    }
    fContainer->AddAt(arr, 1+ip);
  }

  fData = new TTree(GetName(), Form("Reference data for %s", GetName()));
  fData->Branch("ps", &fPIDbin, "p/b");
  fData->Branch("dEdx" , fdEdx , Form("dEdx[%d]/F", GetNslices()));
}


//________________________________________________________________________
Float_t* AliTRDpidRefMakerLQ::CookdEdx(AliTRDseedV1 *trklt)
{
// Fill dEdx array for multidim LQ PID

  trklt->CookdEdx(AliTRDpidUtil::kLQslices);
  const Float_t *dedx = trklt->GetdEdx();
  if(dedx[0]+dedx[1] <= 0.) return 0x0;
  if(dedx[2] <= 0.) return 0x0;

  fdEdx[0] = TMath::Log(dedx[0]+dedx[1]);
  fdEdx[1] = TMath::Log(dedx[2]);
  return fdEdx;
}

//__________________________________________________________________
TObject* AliTRDpidRefMakerLQ::GetOCDBEntry(Option_t *opt)
{
// Steer loading of OCDB LQ PID

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
// Steering reference picture

  if(ifig<0 || ifig>AliTRDCalPID::kNMom-1){ 
    AliError("Ref fig requested outside definition.");
    return kFALSE;
  }
  if(!gPad){
    AliWarning("Please provide a canvas to draw results.");
    return kFALSE;
  }

  TObjArray *arr = (TObjArray*)fContainer->At(ifig);
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
  Int_t fSbin = TMath::LocMax(AliPID::kSPECIES, fPID);
  // get momentum bin
  Int_t fPbin = AliTRDpidUtil::GetMomentumBin(fP);
  fPIDbin = ((fSbin&0xf)<<4)|(fPbin&0xf);
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
    AliError("Calibration data not available");
    return kFALSE;
  }
  TObjArray *o = 0x0;
  if(!(o = (TObjArray*)gFile->Get(Form("Moni%s", GetName())))){
    AliWarning("Missing monitoring container.");
    return kFALSE;
  }
  fContainer = (TObjArray*)o->Clone("monitor");

  TDatime d;
  TDirectoryFile *pdfs = new TDirectoryFile(Form("PDF_%d", d.GetDate()), "PDFs for LQ TRD-PID", "", gFile);
  pdfs->Write();
  AliDebug(2, Form("Data[%d]", fData->GetEntries()));
  pdfs->cd();

  //TCanvas *cc = new TCanvas("cc", "", 500, 500);

  Float_t *data[] = {0x0, 0x0};
//  TPrincipal principal(2, "ND");
//  for(Int_t ip=AliTRDCalPID::kNMom; ip--; ){ 
  for(Int_t ip=5; ip>0; ip-=10){ 
    for(Int_t is=AliPID::kSPECIES; is--;) {
//      principal.Clear();
      Int_t bin = (is<<4) | ip;
      Int_t n = fData->Draw("dEdx[0]:dEdx[1]", Form("ps==%d", bin), "goff");

      // estimate bucket statistics
      Int_t nb(kMinBuckets), // number of buckets
            ns(Int_t(Float_t(n)/nb));    //statistics/bucket
            
// if(Float_t(n)/nb < 220.) ns = 200; // 7% stat error
//       else if(Float_t(n)/nb < 420.) ns = 400; // 5% stat error

      AliDebug(2, Form("pBin[%d] sBin[%d] n[%d] ns[%d] nb[%d]", ip, is, n, ns, nb));
      if(ns<Int_t(kMinStat)){
        AliWarning(Form("Not enough entries [%d] for %s[%d].", n, AliPID::ParticleShortName(is), ip));
        continue;
      }
      // allocate storage
      data[0] = new Float_t[n];data[1] = new Float_t[n];
      // fill
      Double_t *v0 = fData->GetV1(), *v1 = fData->GetV2();
      for(Int_t in=n; in--;){
        data[0][in]=v0[in]; data[1][in]=v1[in];
      }

      // build PDF
      TKDPDF pdf(n, 2, ns, data);
      pdf.SetCOG(kFALSE);
      pdf.SetWeights();
      pdf.SetStore();
      pdf.SetAlpha(5.);
      pdf.GetStatus();
      Float_t *c, v, ve; Double_t r, e, rxy[2];
      for(Int_t in=pdf.GetNTNodes(); in--;){
        pdf.GetCOGPoint(in, c, v, ve);
        rxy[0] = (Double_t)c[0];rxy[1] = (Double_t)c[1];
        pdf.Eval(rxy, r, e, kTRUE);
      }
//       // visual on-line monitoring
//       pdf.DrawProjection();cc->Modified(); cc->Update(); cc->SaveAs(Form("pdf_%s%02d.gif", AliPID::ParticleShortName(is), ip));
//       cc->SaveAs(Form("%s_%s%02d.gif", GetName(), AliPID::ParticleShortName(is), ip));

      // save a discretization of the PDF for result monitoring
      TH2 *h2s = (TH2D*)((TObjArray*)fContainer->At(ip))->At(is);
      TAxis *ax = h2s->GetXaxis(), *ay = h2s->GetYaxis();
      h2s->Clear();
      for(int ix=1; ix<=ax->GetNbins(); ix++){
        rxy[0] = ax->GetBinCenter(ix);
        for(int iy=1; iy<=ay->GetNbins(); iy++){
          rxy[1] = ay->GetBinCenter(iy);
      
          Double_t rr,ee;
          pdf.Eval(rxy, rr, ee, kFALSE);
          if(rr<0. || ee/rr>.15) continue; // 15% relative error
          //printf("x[%2d] x[%2d] r[%f] e[%f]\n", ix, iy, rr, ee);
          h2s->SetBinContent(ix, iy, rr);
        }
      }

      // write results to output array
      //pdf.GetStatus();
      pdf.Write(Form("%s[%d]", AliPID::ParticleShortName(is), ip));

      delete [] data[0]; delete [] data[1];
    }
  }
  pdfs->Write();
  fCalib->Close(); delete fCalib;

  return kTRUE; // testing protection
}


