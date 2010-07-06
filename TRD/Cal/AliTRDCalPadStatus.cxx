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

/* $Id$ */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  TRD calibration class for the single pad status                          //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TStyle.h>
#include <TCanvas.h>
#include <TH1F.h>
#include <TH2F.h>

#include "AliTRDCalPadStatus.h"
#include "AliTRDgeometry.h"
#include "AliTRDpadPlane.h"
#include "AliTRDCalSingleChamberStatus.h"  // test

ClassImp(AliTRDCalPadStatus)

//_____________________________________________________________________________
AliTRDCalPadStatus::AliTRDCalPadStatus()
  :TNamed()
{
  //
  // AliTRDCalPadStatus default constructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    fROC[idet] = 0;
  }

}

//_____________________________________________________________________________
AliTRDCalPadStatus::AliTRDCalPadStatus(const Text_t *name, const Text_t *title)
  :TNamed(name,title)
{
  //
  // AliTRDCalPadStatus constructor
  //

  for (Int_t isec = 0; isec < kNsect; isec++) {
    for (Int_t ipla = 0; ipla < kNplan; ipla++) {
      for (Int_t icha = 0; icha < kNcham; icha++) {
        Int_t idet = AliTRDgeometry::GetDetector(ipla,icha,isec);
        fROC[idet] = new AliTRDCalSingleChamberStatus(ipla,icha,144);
      }
    }
  }

}

//_____________________________________________________________________________
AliTRDCalPadStatus::AliTRDCalPadStatus(const AliTRDCalPadStatus &c)
  :TNamed(c)
{
  //
  // AliTRDCalPadStatus copy constructor
  //

  ((AliTRDCalPadStatus &) c).Copy(*this);

}

//_____________________________________________________________________________
AliTRDCalPadStatus::~AliTRDCalPadStatus()
{
  //
  // AliTRDCalPadStatus destructor
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    if (fROC[idet]) {
      delete fROC[idet];
      fROC[idet] = 0;
    }
  }

}

//_____________________________________________________________________________
AliTRDCalPadStatus &AliTRDCalPadStatus::operator=(const AliTRDCalPadStatus &c)
{
  //
  // Assignment operator
  //

  if (this != &c) ((AliTRDCalPadStatus &) c).Copy(*this);
  return *this;

}

//_____________________________________________________________________________
void AliTRDCalPadStatus::Copy(TObject &c) const
{
  //
  // Copy function
  //

  for (Int_t idet = 0; idet < kNdet; idet++) {
    if (fROC[idet]) {
      fROC[idet]->Copy(*((AliTRDCalPadStatus &) c).fROC[idet]);
    }
  }

  TObject::Copy(c);

}

//_____________________________________________________________________________
Bool_t AliTRDCalPadStatus::CheckStatus(Int_t d, Int_t col, Int_t row, Int_t bitMask) const
{
  //
  // Checks the pad status
  //

  AliTRDCalSingleChamberStatus *roc = GetCalROC(d);
  if (!roc) {
    return kFALSE;
  }
  else {
    return (roc->GetStatus(col, row) & bitMask) ? kTRUE : kFALSE;
  }

}

//_____________________________________________________________________________
AliTRDCalSingleChamberStatus* AliTRDCalPadStatus::GetCalROC(Int_t p, Int_t c, Int_t s) const
{ 
  //
  // Returns the readout chamber of this pad
  //

  return fROC[AliTRDgeometry::GetDetector(p,c,s)];   

}

//_____________________________________________________________________________
TH1F *AliTRDCalPadStatus::MakeHisto1D()
{
  //
  // Make 1D histo
  //

  char  name[1000];
  sprintf(name,"%s Pad 1D",GetTitle());
  TH1F * his = new TH1F(name,name,6, -0.5,5.5);
  his->GetXaxis()->SetBinLabel(1,"Good");
  his->GetXaxis()->SetBinLabel(2,"Masked");
  his->GetXaxis()->SetBinLabel(3,"PadBridgedLeft");
  his->GetXaxis()->SetBinLabel(4,"PadBridgedRight");
  his->GetXaxis()->SetBinLabel(5,"ReadSecond");
  his->GetXaxis()->SetBinLabel(6,"NotConnected");

  for (Int_t idet = 0; idet < kNdet; idet++) 
    {
      if (fROC[idet])
	{
	  for (Int_t ichannel=0; ichannel<fROC[idet]->GetNchannels(); ichannel++)
	    {
	      Int_t status = (Int_t) fROC[idet]->GetStatus(ichannel);
	      if(status==2)  status= 1;
	      if(status==4)  status= 2;
	      if(status==8)  status= 3;
	      if(status==16) status= 4;
	      if(status==32) status= 5;
	      his->Fill(status);
	    }
	}
    }

  return his;

}

//_____________________________________________________________________________
TH2F *AliTRDCalPadStatus::MakeHisto2DSmPl(Int_t sm, Int_t pl)
{
  //
  // Make 2D graph
  //

  gStyle->SetPalette(1);
  AliTRDgeometry *trdGeo = new AliTRDgeometry();
  AliTRDpadPlane *padPlane0 = trdGeo->GetPadPlane(pl,0);
  Double_t row0    = padPlane0->GetRow0();
  Double_t col0    = padPlane0->GetCol0();

  char  name[1000];
  sprintf(name,"%s Pad 2D sm %d pl %d",GetTitle(),sm,pl);
  TH2F * his = new TH2F( name, name, 88,-TMath::Abs(row0),TMath::Abs(row0)
                                   ,148,-TMath::Abs(col0),TMath::Abs(col0));

  // Where we begin
  Int_t offsetsmpl = 30*sm+pl;

  for (Int_t k = 0; k < kNcham; k++){
    Int_t det = offsetsmpl+k*6;
    if (fROC[det]){
      AliTRDCalSingleChamberStatus * calRoc = fROC[det];
      for (Int_t icol=0; icol<calRoc->GetNcols(); icol++){
	for (Int_t irow=0; irow<calRoc->GetNrows(); irow++){
	  Int_t binz     = 0;
	  Int_t kb       = kNcham-1-k;
	  Int_t krow     = calRoc->GetNrows()-1-irow;
	  Int_t kcol     = calRoc->GetNcols()-1-icol;
	  if(kb > 2) binz = 16*(kb-1)+12+krow+1+2*(kb+1);
	  else binz = 16*kb+krow+1+2*(kb+1); 
	  Int_t biny = kcol+1+2;
	  Float_t value = calRoc->GetStatus(icol,irow);
	  his->SetBinContent(binz,biny,value);
	}
      }
      for(Int_t icol = 1; icol < 147; icol++){
	for(Int_t l = 0; l < 2; l++){
	  Int_t binz     = 0;
	  Int_t kb       = kNcham-1-k;
	  if(kb > 2) binz = 16*(kb-1)+12+1+2*(kb+1)-(l+1);
	  else binz = 16*kb+1+2*(kb+1)-(l+1); 
	  his->SetBinContent(binz,icol,50.0);
	}
      }
    }
  }
  for(Int_t icol = 1; icol < 147; icol++){
    his->SetBinContent(88,icol,50.0);
    his->SetBinContent(87,icol,50.0);
  }
  for(Int_t irow = 1; irow < 89; irow++){
    his->SetBinContent(irow,1,50.0);
    his->SetBinContent(irow,2,50.0);
    his->SetBinContent(irow,147,50.0);
    his->SetBinContent(irow,148,50.0);
  }

  his->SetXTitle("z (cm)");
  his->SetYTitle("y (cm)");
  his->SetMaximum(50);
  his->SetMinimum(0.0);
  his->SetStats(0);

  return his;

}

//_____________________________________________________________________________
void AliTRDCalPadStatus::PlotHistos2DSm(Int_t sm, const Char_t *name)
{
  //
  // Make 2D graph
  //

  gStyle->SetPalette(1);
  TCanvas *c1 = new TCanvas(name,name,50,50,600,800);
  c1->Divide(3,2);
  c1->cd(1);
  MakeHisto2DSmPl(sm,0)->Draw("colz");
  c1->cd(2);
  MakeHisto2DSmPl(sm,1)->Draw("colz");
  c1->cd(3);
  MakeHisto2DSmPl(sm,2)->Draw("colz");
  c1->cd(4);
  MakeHisto2DSmPl(sm,3)->Draw("colz");
  c1->cd(5);
  MakeHisto2DSmPl(sm,4)->Draw("colz");
  c1->cd(6);
  MakeHisto2DSmPl(sm,5)->Draw("colz");

}
