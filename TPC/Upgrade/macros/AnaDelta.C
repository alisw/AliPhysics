#include <TStyle.h>
#include <TROOT.h>
#include <TAxis.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <THn.h>
#include <TObjArray.h>
#include <TObject.h>
#include <TString.h>
#include <TTreeStream.h>

/*

.L $ALICE_ROOT/TPC/Upgrade/macros/AnaDelta.C+g


*/

void AnaDelta(TString file, TString outDir=".")
{
  //
  //
  //

  gStyle->SetOptFit();
  
  TTreeSRedirector stream(Form("%s/deltas.root",outDir.Data()));
  gROOT->cd();
  
  TFile f(file);
  THn *hn=(THn*)f.Get("hn");

  TAxis *ar   = hn->GetAxis(0);
  TAxis *aphi = hn->GetAxis(1);
  TAxis *az   = hn->GetAxis(2);

  
  for (Int_t iz=0; iz<az->GetNbins(); ++iz) {
    az->SetRange(iz+1,iz+1);
    TObjArray arrFits;
    arrFits.SetName(Form("z_%02d",iz));
    arrFits.SetOwner();
    
    for (Int_t ir=0; ir<ar->GetNbins(); ++ir) {
      ar->SetRange(ir+1,ir+1);
      for (Int_t iphi=0; iphi<aphi->GetNbins(); ++iphi) {
        aphi->SetRange(iphi+1,iphi+1);
      
        TH1 *hProj = hn->Projection(3);
        if (hProj->GetEntries()<1) {
          delete hProj;
          continue;
        }

        TF1 fg("fg","gaus",-2,2);
        Double_t cr   = ar->GetBinCenter(ir+1);
        Double_t cphi = aphi->GetBinCenter(iphi+1);
        Double_t cz   = az->GetBinCenter(iz+1);
        hProj->SetNameTitle(Form("h_%02d_%02d_%d02",iz+1, iphi+1, ir+1),
                         Form("z,phi,r: %02d,%02d,%d02 (%.2f, %.2f, %.2f)",iz+1,iphi+1,ir+1, cr, cphi, cz ));
        hProj->Fit(&fg,"QR");
        arrFits.Add(hProj);

        Double_t mean     = fg.GetParameter(1);
        Double_t meanErr  = fg.GetParError(1);
        Double_t sigma    = fg.GetParameter(2);
        Double_t sigmaErr = fg.GetParError(2);
        Int_t    entries  = hProj->GetEntries();
        Double_t chi2ndf  = fg.GetChisquare()/fg.GetNDF();
        
        stream << "d" <<
        "ir="        << ir       <<
        "iphi="      << iphi     <<
        "iz="        << iz       <<
        "cr="        << cr       <<
        "cphi="      << cphi     <<
        "cz="        << cz       <<
        "mean="      << mean     <<
        "meanErr="   << meanErr  <<
        "sigma="     << sigma    <<
        "sigmaErr="  << sigmaErr <<
        "entries="   << entries  <<
        "chi2ndf="   << chi2ndf  <<
        "\n";
      }
    }
    stream.GetFile()->cd();
    arrFits.Write(0x0,TObject::kSingleKey);
    gROOT->cd();
  }

  delete hn;
}

