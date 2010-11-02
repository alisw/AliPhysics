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

/*   Origin: Alberica Toia, CERN, Alberica.Toia@cern.ch                   */

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//  class to determine centrality percentiles from 1D distributions          // 
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#include <TNamed.h>
#include <TH1D.h>
#include <TString.h>
#include <TFile.h>
#include <TMath.h>
#include <TROOT.h>
#include <TH2F.h>
#include <TF1.h>
#include <TStyle.h>
#include <TGraphErrors.h>
#include <vector>
#include "AliCentralityBy1D.h"


ClassImp(AliCentralityBy1D)  
 
//______________________________________________________________________________


AliCentralityBy1D::AliCentralityBy1D() {
  // standard constructor
}

AliCentralityBy1D::~AliCentralityBy1D() {
  // destructor
}

void AliCentralityBy1D::AddHisto(TString name) {
  histnames.push_back(name);
}

void AliCentralityBy1D::SetPercentileFile(TString filename) {
  outrootfilename = filename;
}

void AliCentralityBy1D::SetPercentileCrossSection(Float_t xsec) {
  percentXsec = xsec;
}

void AliCentralityBy1D::MakePercentiles(TString infilename) {
  TH1D *thist;
  
  TFile *outrootfile;
  
  // open inrootfile, outrootfile
  inrootfile  = new TFile(infilename);
  outrootfile = new TFile(outrootfilename,"RECREATE");
  
  // loop over all distribution names
  
  vector<TString>::const_iterator hni;
  for(hni=histnames.begin(); hni!=histnames.end(); hni++) {
    thist = MakePercentHisto(*hni);
    SaveHisto(thist,outrootfile);
    delete thist; //??
  }
  // close inrootfile, outrootfile
  inrootfile->Close();
  outrootfile->Close();
  
}

TH1D * AliCentralityBy1D::MakePercentHisto(TString hdistributionName) {
  TH1D *htemp  = (TH1D*) (inrootfile->Get(hdistributionName)); 
  TH1D *hpercent  = (TH1D*) htemp->Clone("hpercent");
  hpercent->SetName(hdistributionName.Append("_percentile"));
  hpercent->Reset();

  for (int ibin=1; ibin<=htemp->GetNbinsX(); ibin++) {
    
    hpercent->SetBinContent(ibin, percentXsec *
			    htemp->Integral(ibin,htemp->GetNbinsX())  / 
			    htemp->Integral(1,htemp->GetNbinsX()));
    
  }
  
  delete htemp;
  return hpercent;
  
}

void AliCentralityBy1D::SaveHisto(TH1D *hist, TFile *outrootfile) {
  outrootfile->cd();
  hist->Write();
}



