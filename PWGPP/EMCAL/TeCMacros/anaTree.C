#if !defined(__CINT__) || defined(__MAKECINT__)
#include <AliEMCALGeometry.h>
#include <TCanvas.h>
#include <TClonesArray.h>
#include <TDatime.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TGrid.h>
#include <TH2F.h>
#include <TLegend.h>
#include <TMap.h>
#include <TNtuple.h>
#include <TObject.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TSystem.h>
#include "createTree.C"
#endif

void anaTree(
              const char *ifile     ="treefile.root",
              const char *ofile     ="outhist.root",
              Bool_t applyDurLimit  = kFALSE,         // bolean to switch on duration limit
              Float_t durMin        = 0,              // minimum length of a run in minutes
              Float_t durMax        = 10000           // maximum length of a run in minutes
            )
{
  AliEMCALGeometry *g=AliEMCALGeometry::GetInstance("EMCAL_COMPLETE12SMV1_DCAL_8SM");
  const Int_t kSM=g->GetNumberOfSuperModules();
  const Int_t kNcells=g->GetNCells();
  const Int_t gain = 1;

  TCalInfo *info = 0;
  TFile *in = TFile::Open(ifile,"read");
  TTree *tt = (TTree*)in->Get("tcal");
  tt->SetBranchAddress("event",&info);
  tt->Branch("event", &info, 32000, 99);

  TProfile *gLedVsT[20];
  TProfile *gLedMonVsT[20];
//   TProfile *gRatVsT[20];
  for (Int_t i=0;i<20;++i) {
    gLedVsT[i] = new TProfile("","Led info;T;",10000,10,50);
    gLedVsT[i]->SetMarkerStyle(20);
    gLedVsT[i]->SetMarkerSize(1.2);
    gLedVsT[i]->SetMarkerColor(1);
    gLedVsT[i]->SetLineColor(1);
    gLedVsT[i]->SetLineWidth(3);
    gLedVsT[i]->Sumw2();
    gLedVsT[i]->SetName(Form("ledsm%d",i));
    gLedMonVsT[i] = new TProfile("","Led info;T;",10000,10,50);
    gLedMonVsT[i]->SetMarkerStyle(20);
    gLedMonVsT[i]->SetMarkerSize(1.2);
    gLedMonVsT[i]->SetMarkerColor(1);
    gLedMonVsT[i]->SetLineColor(1);
    gLedMonVsT[i]->SetLineWidth(3);
    gLedMonVsT[i]->Sumw2();
    gLedMonVsT[i]->SetName(Form("ledmonsm%d",i));
//     gRatVsT[i] = new TProfile("","Led/LedMon;T",10000,10,50);
//     gRatVsT[i]->SetMarkerStyle(20);
//     gRatVsT[i]->SetMarkerSize(1.2);
//     gRatVsT[i]->SetMarkerColor(1);
//     gRatVsT[i]->SetLineColor(1);
//     gRatVsT[i]->SetLineWidth(3);
//     gRatVsT[i]->Sumw2();
//     gRatVsT[i]->SetName(Form("ledovermonsm%d",i));
  }

  TProfile *gLedCellVsT[kNcells+1];
  TProfile *gLedMonCellVsT[kNcells+1];
//   TProfile *gRatCellVsT[kNcells+1];

  for (Int_t j=0;j<kNcells+1;++j) {
    gLedCellVsT[j] = new TProfile("",Form("Led info cell ID%i ;T;",j),10000,10,50);
    gLedCellVsT[j]->SetMarkerStyle(20);
    gLedCellVsT[j]->SetMarkerSize(1.2);
    gLedCellVsT[j]->SetMarkerColor(1);
    gLedCellVsT[j]->SetLineColor(1);
    gLedCellVsT[j]->SetLineWidth(3);
    gLedCellVsT[j]->Sumw2();
    gLedCellVsT[j]->SetName(Form("ledCell%d",j));
    gLedMonCellVsT[j] = new TProfile("",Form("LedMon info cell ID%i ;T;",j),10000,10,50);
    gLedMonCellVsT[j]->SetMarkerStyle(20);
    gLedMonCellVsT[j]->SetMarkerSize(1.2);
    gLedMonCellVsT[j]->SetMarkerColor(1);
    gLedMonCellVsT[j]->SetLineColor(1);
    gLedMonCellVsT[j]->SetLineWidth(3);
    gLedMonCellVsT[j]->Sumw2();
    gLedMonCellVsT[j]->SetName(Form("ledMonCell%d",j));
//     gRatCellVsT[j] = new TProfile("",Form("Led/LedMon cell ID%i ;T;",j),10000,10,50);
//     gRatCellVsT[j]->SetMarkerStyle(20);
//     gRatCellVsT[j]->SetMarkerSize(1.2);
//     gRatCellVsT[j]->SetMarkerColor(1);
//     gRatCellVsT[j]->SetLineColor(1);
//     gRatCellVsT[j]->SetLineWidth(3);
//     gRatCellVsT[j]->Sumw2();
//     gRatCellVsT[j]->SetName(Form("ledovermonCell%d",j));
  }


  Int_t Nev=tt->GetEntries();
  for (Int_t i=0;i<Nev;++i) {
    tt->GetEvent(i);
    Float_t deltaTime = ((Float_t)info->fLastTime-(Float_t)info->fFirstTime)/60.;         // run duration in minutes
    cout << info->fRunNo << "\t" << info->fFirstTime << "\t"<<info->fLastTime << "\t"<< deltaTime<<  endl;
    if (applyDurLimit){
      if (deltaTime < durMin || deltaTime > durMax){
        cout << "INFO: skipped run due to mismatch in run length" << endl;
        continue;
      }
    }
    TClonesArray &cells = info->fCells;
    for (Int_t j=0;j<cells.GetEntries();++j) {
      TCalCell *cell = static_cast<TCalCell*>(cells.At(j));
      Int_t cellID  = cell->fId;
      Int_t sm      = cell->fSM;
      Double_t ledM = cell->fLedM;
      Double_t ledR = cell->fLedR;
      Double_t monM = cell->fMonM;
      Double_t monR = cell->fMonR;
      Double_t locT = cell->fLocT;
      Double_t smT  = cell->fSMT;

      Double_t T = smT; // use local or SM T, change here
      if ((T<5)||(T>45))
        continue;
      if ((ledM<=0)||(ledR<=0))
        continue;
      if ((monM<=0)||(monR<=0))
        continue;
      Double_t w = ledR;
      TProfile *h     = gLedVsT[sm];
      h->Fill(T,ledM,w);
      gLedCellVsT[cellID]->Fill(smT,ledM,w);
      Double_t w3 = monR;
      TProfile *h3     = gLedMonVsT[sm];
      h3->Fill(T,monM,w3);
      gLedMonCellVsT[cellID]->Fill(smT,monM,w3);
//       Double_t w2=TMath::Sqrt(ledR*ledR+monM*monM);
//       TProfile *h2 = gRatVsT[sm];
//       h2->Fill(T,ledM/monM,w2);
//       gRatCellVsT[cellID]->Fill(smT,ledM/monM,w2);
    }
  }

  TFile *out = TFile::Open(ofile,"recreate");
  if (!out)
    out = TFile::Open("dummyfile.root","update");
  for (Int_t i=0;i<20;++i) {
    gLedVsT[i]->Write();
    gLedMonVsT[i]->Write();
//     gRatVsT[i]->Write();
  }
  Int_t smcurr          = -1;
  Int_t cellIDFirstInSM = 0;
  for (Int_t j=0;j<kNcells+1;++j) {
    Int_t iTower = -1, iIphi = -1,  iIeta = -1, sm=-1;
    g->GetCellIndex(j,sm,iTower,iIphi,iIeta);
    if (smcurr != sm){
      smcurr = sm;
      out->mkdir(Form("cellsInSM%d",sm));
      out->cd(Form("cellsInSM%d",sm));
      cout << "SM: \t" << sm-1 << "\t"<< cellIDFirstInSM << "\t" << j << "\t" << j -cellIDFirstInSM<< endl;
      cellIDFirstInSM = j;
    }

    if (gLedCellVsT[j]->GetEntries() > 0)
      gLedCellVsT[j]->Write();
    if (gLedMonCellVsT[j]->GetEntries() > 0)
      gLedMonCellVsT[j]->Write();
//     if (gRatCellVsT[j]->GetEntries() > 0)
//       gRatCellVsT[j]->Write();
  }
}
