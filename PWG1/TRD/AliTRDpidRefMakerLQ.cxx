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
//
//   Origin
//   Alex Bercuci  (A.Bercuci@gsi.de)
//
///////////////////////////////////////////////////////////////////////////////

#include <TSystem.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TMath.h>
#include <TEventList.h>
#include <TH2D.h>
#include <TH2I.h>
#include <TCanvas.h>

#include "AliLog.h"
#include "../STAT/TKDPDF.h"
#include "../STAT/TKDInterpolator.h"
#include "AliTRDpidRefMakerLQ.h"
#include "Cal/AliTRDCalPID.h"
#include "Cal/AliTRDCalPIDLQ.h"
#include "AliTRDseedV1.h"
#include "AliTRDcalibDB.h"
#include "AliTRDgeometry.h"
#include "info/AliTRDpidInfo.h"

ClassImp(AliTRDpidRefMakerLQ)

//__________________________________________________________________
AliTRDpidRefMakerLQ::AliTRDpidRefMakerLQ() 
  : AliTRDpidRefMaker("PIDrefMakerLQ", "PID(LQ) Reference Maker")
  , fPDF(NULL)
{
  //
  // AliTRDpidRefMakerLQ default constructor
  //
  DefineOutput(2, TObjArray::Class());
}

//__________________________________________________________________
AliTRDpidRefMakerLQ::~AliTRDpidRefMakerLQ()
{
  //
  // AliTRDCalPIDQRef destructor
  //
  if(fPDF){
    fPDF->Write("PDF_2DLQ", TObject::kSingleKey);
    fPDF->Delete();
    delete fPDF;
  }
}

// //________________________________________________________________________
void AliTRDpidRefMakerLQ::CreateOutputObjects()
{
  // Create histograms
  // Called once

  //OpenFile(0, "RECREATE");
  fContainer = Histos();

  OpenFile(2, "RECREATE");
  fPDF = new TObjArray(AliTRDCalPID::kNMom*AliPID::kSPECIES);
  fPDF->SetOwner();fPDF->SetName("PDF_2DLQ");
}


//__________________________________________________________________
TObjArray* AliTRDpidRefMakerLQ::Histos()
{
  // Create histograms

  if(fContainer) return fContainer;

  fContainer = new TObjArray(AliTRDCalPID::kNMom);
  //fContainer->SetOwner(kTRUE);
  //fContainer->SetName(Form("Moni%s", GetName()));

  // save dE/dx references
  TH2 *h2 = 0x0;
  for(Int_t ip=AliTRDCalPID::kNMom; ip--;){ 
    TObjArray *arr = new TObjArray(AliPID::kSPECIES);
    arr->SetName(Form("Pbin%02d", ip)); arr->SetOwner();
    for(Int_t is=AliPID::kSPECIES; is--;) {
      h2 = new TH2D(Form("h%s%d", AliPID::ParticleShortName(is), ip), Form("%s @ Pbin[%d]", AliPID::ParticleName(is), ip), 50, 7., 12., 50, 6.5, 11.);
      h2->GetXaxis()->SetTitle("log(dE/dx_{am}) [au]");
      h2->GetYaxis()->SetTitle("log(dE/dx_{dr}) [au]");
      h2->GetZaxis()->SetTitle("#");
      arr->AddAt(h2, is);
    }
    fContainer->AddAt(arr, ip);
  }
  fNRefFigures=AliTRDCalPID::kNMom;
  return fContainer;
}


//__________________________________________________________________
TObject* AliTRDpidRefMakerLQ::GetOCDBEntry(Option_t */*opt*/)
{
// Steer loading of OCDB LQ PID

  if(gSystem->AccessPathName(Form("TRD.Calib%s.root", GetName()), kReadPermission)){
    AliError(Form("File TRD.Calib%s.root not readable", GetName()));
    return NULL;
  }
  AliTRDCalPIDLQ *cal = new AliTRDCalPIDLQ("pidLQ", "LQ TRD PID object");
  cal->LoadReferences(Form("TRD.Calib%s.root", GetName()));
  return cal;
}

//__________________________________________________________________
Bool_t AliTRDpidRefMakerLQ::GetRefFigure(Int_t ifig)
{
// Steering reference picture

  if(ifig<0 || ifig>=fNRefFigures){ 
    AliError("Ref fig requested outside definition.");
    return kFALSE;
  }
  if(!gPad){
    AliWarning("Please provide a canvas to draw results.");
    return kFALSE;
  }

  TObjArray *arr(NULL);TList *l(NULL);TH2 *h2(NULL);
  switch(ifig){
  case 0: // PDG plot
    if(!(h2=(TH2*)fContainer->At(0))){
      AliError("Abundance Plot missing.");
      return kFALSE;
    }
    h2->Draw("cont4z");
    break;
  default:
    if(!(arr = (TObjArray*)fContainer->At(ifig))){
      AliError(Form("2D container @ pBin[%d] missing.", ifig-1));
      return kFALSE;
    }
    gPad->Divide(3, 2, 1.e-5, 1.e-5);l=gPad->GetListOfPrimitives(); 
    for(Int_t is=0; is<AliPID::kSPECIES; is++){
      ((TVirtualPad*)l->At(is))->cd();
      if(!(h2=(TH2*)arr->At(is))){
        AliError(Form("2D for %s @ pBin[%d] missing.", AliPID::ParticleShortName(is), ifig-1));
        return kFALSE;
      }
      h2->Draw("cont4z");
    }
  }
  return kTRUE;
}


//________________________________________________________________________
Bool_t AliTRDpidRefMakerLQ::Load(const Char_t */*fname*/)
{
  const Char_t *name("PIDrefMaker");
  if(gSystem->AccessPathName(Form("TRD.Calib%s.root", name), kReadPermission)){
    AliError(Form("File TRD.Calib%s.root not readable", name));
    return kFALSE;
  }
  if(!TFile::Open(Form("TRD.Calib%s.root", name))){
    AliError(Form("File TRD.Calib%s.root corrupted", name));
    return kFALSE;
  }
  if (!(fData = dynamic_cast<TTree*>(gFile->Get(name)))) {
    AliError(Form("Tree %s not available", name));
    return kFALSE;
  }
  LinkPIDdata();

  TObjArray *o(NULL);
  if(!(o = (TObjArray*)gFile->Get(Form("Moni%s", GetName())))){
    AliWarning(Form("Monitor container Moni%s not available.", name));
    return kFALSE;
  }
  fContainer = (TObjArray*)o->Clone("monitor");
  fNRefFigures=AliTRDCalPID::kNMom;


  if(!TFile::Open(Form("TRD.Calib%s.root", GetName()), "UPDATE")){
    AliError(Form("File TRD.Calib%s.root corrupted", GetName()));
    return kFALSE;
  }
  if(!(o = (TObjArray*)gFile->Get(Form("PDF_2DLQ")))) {
    AliInfo("PDF container PDF_2DLQ not available. Create.");
    fPDF = new TObjArray(AliTRDCalPID::kNMom*AliPID::kSPECIES);
    fPDF->SetOwner();fPDF->SetName("PDF_2DLQ");
  } else fPDF = (TObjArray*)o->Clone("PDF_2DLQ");

  return kTRUE;
}


//________________________________________________________________________
void AliTRDpidRefMakerLQ::Exec(Option_t */*opt*/)
{
// Load PID data into local data storage

  AliTRDpidInfo *pid(NULL);
  const AliTRDpidInfo::AliTRDpidData *data(NULL);
  Char_t s(-1);
  for(Int_t itrk=fInfo->GetEntriesFast(); itrk--;){
    if(!(pid=(AliTRDpidInfo*)fInfo->At(itrk))) continue;
    if((s=pid->GetPID())<0) continue;
    for(Int_t itrklt=pid->GetNtracklets();itrklt--;){
      data=pid->GetData(itrklt);
      Int_t ip(data->Momentum());
      

      Double_t dedx[] = {0., 0.};
      if(!AliTRDCalPIDLQ::CookdEdx(data->fdEdx, dedx)) continue;
      ((TH2*)((TObjArray*)fContainer->At(ip))->At(s))->Fill(dedx[0], dedx[1]);
    }
  }

  PostData(0, fContainer);
  PostData(2, fPDF);
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

  TCanvas *fMonitor(NULL);
  // allocate working storage
  const Int_t kWS(AliPID::kSPECIES*AliTRDCalPID::kNMom);
  Float_t *data[2*kWS];
  for(Int_t i(0); i<2*kWS; i++) data[i]=new Float_t[kMaxStat];
  Int_t ndata[kWS]; memset(ndata, 0, kWS*sizeof(Int_t));

  AliDebug(1, Form("Loading data[%d]", fData->GetEntries()));
  for(Int_t itrk=0; itrk < fData->GetEntries(); itrk++){
    if(!(fData->GetEntry(itrk))) continue;
    Int_t sbin(fPIDbin);
    for(Int_t ily=fPIDdataArray->fNtracklets; ily--;){
      Int_t pbin(fPIDdataArray->fData[ily].fPLbin & 0xf);
      
      Double_t dedx[] = {0., 0.};
      if(!AliTRDCalPIDLQ::CookdEdx(fPIDdataArray->fData[ily].fdEdx, dedx)) continue;
      Int_t idx=AliTRDCalPIDLQ::GetModelID(pbin,sbin);
      if(ndata[idx]==kMaxStat) continue;

      // store data
      data[idx][ndata[idx]] = dedx[0];
      data[idx+kWS][ndata[idx]] = dedx[1];
      ndata[idx]++;
    }
  }
  for(Int_t ip=0;ip<AliTRDCalPID::kNMom;ip++){ 
    for(Int_t is=AliPID::kSPECIES; is--;) {
      // estimate bucket statistics
      Int_t idx(AliTRDCalPIDLQ::GetModelID(ip,is)),
            nb(kMinBuckets), // number of buckets
	ns((Int_t)((Float_t)(ndata[idx])/nb));    //statistics/bucket
            
      AliDebug(2, Form("pBin[%d] sBin[%d] n[%d] ns[%d] nb[%d]", ip, is, ndata[idx], ns, nb));
      if(ns<Int_t(kMinStat)){
        AliWarning(Form("Not enough entries [%d] for %s[%d].", ndata[idx], AliPID::ParticleShortName(is), ip));
        continue;
      }

      // build helper PDF
      Float_t *ldata[2]={data[idx], data[kWS+idx]};
      TKDPDF *pdf = new TKDPDF(ndata[idx], 2, ns, ldata);
      pdf->SetCOG(kFALSE);
      pdf->SetWeights();
      //pdf->SetStore();
      pdf->SetAlpha(5.);
      //pdf.GetStatus();
      fPDF->AddAt(pdf, idx);
      Int_t in=pdf->GetNTNodes(); Float_t par[6], *pp=&par[0];
      while(in--){
        const TKDNodeInfo *nn = pdf->GetNodeInfo(in);
        nn->GetCOG(pp);
        Double_t p[] = {par[0], par[1]}, r,e;
        pdf->Eval(p,r,e,1);
      }
/*
      Int_t nnodes = pdf.GetNTNodes(),
            nside = Int_t(0.05*nnodes),
            nzeros = 4*(nside+1);
      printf("nnodes[%d] nside[%d] nzeros[%d]\n", nnodes, nside, nzeros);
    

      // Build interpolator on the pdf skeleton
      TKDInterpolator interpolator(2, nnodes+nzeros); 
      for(Int_t in=nnodes; in--;)
        interpolator.SetNode(in, *pdf.GetNodeInfo(in));
      TKDNodeInfo *nodes = new TKDNodeInfo[nzeros], *node = &nodes[0];
      Float_t ax0min, ax0max, ax1min, ax1max;
      pdf.GetRange(0, ax0min, ax0max); Float_t dx = (ax0max-ax0min)/nside;
      pdf.GetRange(1, ax1min, ax1max); Float_t dy = (ax1max-ax1min)/nside;
      printf("x=[%f %f] y[%f %f]\n", ax0min, ax0max, ax1min, ax1max);

      Int_t jn = nnodes; 
      SetZeroes(&interpolator, node, nside, jn, ax0min, dx, ax1min, -dy, 'x');
      SetZeroes(&interpolator, node, nside, jn, ax1min, dy, ax0max, dx, 'y');
      SetZeroes(&interpolator, node, nside, jn, ax0max,-dx, ax1max, dy, 'x');
      SetZeroes(&interpolator, node, nside, jn ,ax1max, -dy, ax0min, -dx, 'y');
      delete [] nodes;
      Int_t in=nnodes; Float_t par[6], *pp=&par[0];
      while(in--){
        const TKDNodeInfo *nn = interpolator.GetNodeInfo(in);
        nn->GetCOG(pp);
        //printf("evaluate for node[%d] @ [%f %f]\n",in, par[0], par[1]);
        Double_t p[] = {par[0], par[1]}, r,e;
        interpolator.Eval(p,r,e,1);
      }
*/
      // visual on-line monitoring
      if(HasOnlineMonitor()){
        if(!fMonitor) fMonitor = new TCanvas("cc", "PDF 2D LQ", 500, 500);
        pdf->DrawProjection();
        fMonitor->Modified(); fMonitor->Update(); 
        fMonitor->SaveAs(Form("pdf_%s%02d.gif", AliPID::ParticleShortName(is), ip));
      }

      //fContainer->ls();
      // save a discretization of the PDF for result monitoring
      Double_t rxy[]={0.,0.};
      TH2 *h2s = (TH2D*)((TObjArray*)fContainer->At(ip))->At(is);
      TAxis *ax = h2s->GetXaxis(), *ay = h2s->GetYaxis();
      h2s->Clear();
      for(int ix=1; ix<=ax->GetNbins(); ix++){
        rxy[0] = ax->GetBinCenter(ix);
        for(int iy=1; iy<=ay->GetNbins(); iy++){
          rxy[1] = ay->GetBinCenter(iy);
      
          Double_t rr,ee;
          pdf->Eval(rxy, rr, ee, kFALSE);
          if(rr<0. || ee/rr>.15) continue; // 15% relative error
          //printf("x[%2d] x[%2d] r[%f] e[%f]\n", ix, iy, rr, ee);
          h2s->SetBinContent(ix, iy, rr);
        }
      }
    }
  }
  for(Int_t i(0); i<2*kWS; i++) delete [] data[i];
  return kTRUE;
}


//__________________________________________________________________
void AliTRDpidRefMakerLQ::SetZeroes(TKDInterpolator *interpolator, TKDNodeInfo *node, Int_t n, Int_t& idx, Float_t x, Float_t dx, Float_t y, Float_t dy, const Char_t opt)
{
// Set extra nodes to ensure boundary conditions
  
  printf("SetZeroes(%c)\n", opt);
  Float_t par[6], val[] = {0., 1.};
  Int_t a[6];
  if(opt=='x'){
    a[0]=0; a[1]=1; a[2]=2; a[3]=3; a[4]=4; a[5]=5;
  } else if(opt=='y'){
    a[0]=1; a[1]=0; a[2]=4; a[3]=5; a[4]=2; a[5]=3;
  } else return;
  Float_t tmp;
  par[a[1]] = y;
  par[a[4]] = y; par[a[5]] = y+dy;
  if(dy<0.){tmp=par[a[4]]; par[a[4]]=par[a[5]]; par[a[5]]=tmp;}
  for(Int_t in=n; in--; node++, idx++, x+=dx){
    par[a[0]] = x+.5*dx;
    par[a[2]] = x;  par[a[3]] = x+dx;
    if(dx<0.){tmp=par[a[2]]; par[a[2]]=par[a[3]]; par[a[3]]=tmp;}
    node->SetNode(2, par, val);
    printf("\n\tnode[%d]\n", idx); node->Print();
    interpolator->SetNode(idx, *node);
  }
  par[a[0]] = x;
  par[a[2]] = x;  par[a[3]] = x+dx;
  if(dx<0.){tmp=par[a[2]]; par[a[2]]=par[a[3]]; par[a[3]]=tmp;}
  node->SetNode(2, par, val);
  printf("\n\tnode[%d]\n", idx); node->Print();
  interpolator->SetNode(idx, *node);node++;idx++;
}

