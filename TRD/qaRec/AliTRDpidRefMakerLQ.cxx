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
  :AliTRDpidRefMaker("PidRefMakerLQ", "PID(LQ) Reference Maker")
  ,fPbin(-1)
  ,fSbin(-1)
{
  //
  // AliTRDpidRefMakerLQ default constructor
  //

  memset(fH2dEdx, 0x0, AliPID::kSPECIES*sizeof(TH2*));
}

//__________________________________________________________________
AliTRDpidRefMakerLQ::~AliTRDpidRefMakerLQ()
{
  //
  // AliTRDCalPIDQRef destructor
  //

  for(int ispec=0; ispec<AliPID::kSPECIES; ispec++){
    if(fH2dEdx[ispec]) delete fH2dEdx[ispec];
  }	
}

//________________________________________________________________________
void AliTRDpidRefMakerLQ::CreateOutputObjects()
{
  // Create histograms
  // Called once

  AliTRDpidRefMaker::CreateOutputObjects();

  fData = new TTree(GetName(), Form("Reference data for %s", GetName()));
  fData->Branch("s", &fSbin, "l/b");
  fData->Branch("p", &fPbin, "p/b");
  fData->Branch("dEdx" , fdEdx , Form("dEdx[%d]/F", GetNslices()));
}


//________________________________________________________________________
Float_t* AliTRDpidRefMakerLQ::GetdEdx(AliTRDseedV1 *trklt)
{
  trklt->CookdEdx(AliTRDpidUtil::kLQslices);
  const Float_t *dedx = trklt->GetdEdx();
  if(dedx[0]+dedx[1] <= 0.) return 0x0;
  if(dedx[2] <= 0.) return 0x0;

  fdEdx[0] = TMath::Log(dedx[0]+dedx[1]);
  fdEdx[1] = TMath::Log(dedx[2]);
  return fdEdx;
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
}

//________________________________________________________________________
Bool_t AliTRDpidRefMakerLQ::PostProcess()
{
  TFile::Open(Form("TRD.Calib%s.root", GetName()));
  fData = dynamic_cast<TTree*>(gFile->Get(GetName()));
  if (!fData) {
    AliError("Tree not available");
    return kFALSE;
  }
  AliDebug(2, Form("Data[%d]", fData->GetEntries()));

  Float_t *data[] = {0x0, 0x0};
  TPrincipal principal(2, "ND");

  TCanvas *cc = new TCanvas("cc", "", 500, 500);

  for(Int_t ip=AliTRDCalPID::kNMom; ip--;){ 
    for(Int_t is=AliPID::kSPECIES; is--;) {
      principal.Clear();
      Int_t n = fData->Draw("dEdx[0]:dEdx[1]", Form("p==%d&&s==%d", ip, is), "goff");
      AliDebug(2, Form("pBin[%d] sBin[%d] n[%d]", ip, is, n));
      if(n<10){
        AliWarning(Form("Not enough data for %s[%d].", AliPID::ParticleShortName(is), ip));
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
        //hProj->Fill(rxy[0], rxy[1]);
      }

      // estimate acceptable statistical error per bucket

      // estimate number of buckets

      // build PDF
      TKDPDF pdf(n, 2, 10, data);
      printf("PDF nodes[%d]\n", pdf.GetNTNodes());
      pdf.SetStore();
      pdf.SetAlpha(5.);
      Float_t *c, v, ve; Double_t r, e;
      for(Int_t in=pdf.GetNTNodes(); in--;){
        pdf.GetCOGPoint(in, c, v, ve);
        rxy[0] = (Double_t)c[0];rxy[1] = (Double_t)c[1];
        printf("%2d x[%f] y[%f]\n", in, rxy[0], rxy[1]);
        pdf.Eval(rxy, r, e, kTRUE);
      }
      pdf.DrawBins(0,1,-6,6,-6,6);cc->Modified(); cc->Update();
      AliDebug(2, Form("n[%d]", n));

      //pdf.Write(Form("%s[%d]", AliPID::ParticleShortName(is), ip));
      

      delete [] data[0]; delete [] data[1];
    }
  }

  return kTRUE; // testing protection
}


//__________________________________________________________________
void	AliTRDpidRefMakerLQ::Reset()
{
  //
  // Reset reference histograms
  //

  for(int ispec=0; ispec<AliPID::kSPECIES; ispec++){
    if(fH2dEdx[ispec]) fH2dEdx[ispec]->Reset();
  }	
}

//__________________________________________________________________
void  AliTRDpidRefMakerLQ::SaveReferences(const Int_t mom, const char *fn)
{
  //
  // Save the reference histograms
  //

  TFile *fSave = 0x0;
  TListIter it((TList*)gROOT->GetListOfFiles());
  Bool_t kFOUND = kFALSE;
  TDirectory *pwd = gDirectory;
  while((fSave=(TFile*)it.Next()))
    if(strcmp(fn, fSave->GetName())==0){
      kFOUND = kTRUE;
      break;
    }
  if(!kFOUND) fSave = TFile::Open(fn, "RECREATE");
  fSave->cd();

  // save dE/dx references
  TH2 *h2 = 0x0;
  for(int ispec=0; ispec<AliPID::kSPECIES; ispec++){
    h2 = (TH2D*)fH2dEdx[ispec]->Clone(Form("h2dEdx%s%d", AliTRDCalPID::GetPartSymb(ispec), mom));
    h2->SetTitle(Form("2D dEdx for particle %s @ %d", AliTRDCalPID::GetPartName(ispec), mom));
    h2->GetXaxis()->SetTitle("dE/dx_{TRD}^{amplif} [au]");
    h2->GetYaxis()->SetTitle("dE/dx_{TRD}^{drift} [au]");
    h2->GetZaxis()->SetTitle("Entries");
    h2->Write();
  }

  pwd->cd();
}

