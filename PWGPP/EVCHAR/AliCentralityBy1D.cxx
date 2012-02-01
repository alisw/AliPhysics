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
#include <TH1F.h>
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


AliCentralityBy1D::AliCentralityBy1D():
  fInrootfilename(0),
  fOutrootfilename(0),
  fHistnames()
{
  // standard constructor
}

AliCentralityBy1D::~AliCentralityBy1D() {
  // destructor
}

void AliCentralityBy1D::MakePercentiles() {
  TH1F *htemp;  
  TFile *inrootfile;
  TFile *outrootfile;

  // open inrootfile, outrootfile
  std::cout << "input file "  << fInrootfilename  << std::endl;
  std::cout << "output file " << fOutrootfilename << std::endl;
  inrootfile  = new TFile(fInrootfilename,"OPEN");
  //outrootfile = new TFile(fOutrootfilename,"RECREATE");
  outrootfile = new TFile(fOutrootfilename,"UPDATE");

  // loop over all distribution names  
   std::vector<TString>::const_iterator hni;
   for(hni=fHistnames.begin(); hni!=fHistnames.end(); hni++) {
     htemp  = (TH1F*) (inrootfile->Get(*hni)); 
     if (!htemp) {
       TList *list  = (TList*) (inrootfile->Get("CentralityStat")); 
       htemp  = (TH1F*) (list->FindObject(*hni));
     } 

     TH1F *hpercent  = (TH1F*) htemp->Clone("hpercent");
     TString name=htemp->GetName();   
     name.Append("_percentile");
     hpercent->SetNameTitle(name.Data(),name.Data());
     hpercent->Reset();

     int start_bin=htemp->FindBin(fMultLowBound);
     for (int ibin=1; ibin<=htemp->GetNbinsX(); ibin++) {
      
       if (ibin>=start_bin)
   	hpercent->SetBinContent(ibin, fPercentXsec *
   				htemp->Integral(ibin,htemp->GetNbinsX())  / 
   				htemp->Integral(start_bin,htemp->GetNbinsX()));
       else
   	hpercent->SetBinContent(ibin, 100);
     }    

     SaveHisto(htemp,hpercent,outrootfile);

  }
  
  // close inrootfile, outrootfile
  inrootfile->Close();
  outrootfile->Close();
  
}

void AliCentralityBy1D::SaveHisto(TH1F *hist1, TH1F *hist2, TFile *outrootfile) {
  outrootfile->cd();
  hist1->Write();
  hist2->Write();

  // int n=12;
  // Float_t x1[n];
  // Float_t centrality[]= {0.,5.,10.,20.,30.,40.,50.,60.,70.,80.,90.,100.};
  // for (int i=n-1; i>=0; i--) {
  //   x1[i] = hist2->GetBinCenter(hist2->FindLastBinAbove(centrality[i]));
  //   cout << x1[i] << ",";
  // }
  // cout << endl;
  
}



