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

void anaTree(const char *ifile="treefile.root", const char *ofile="outhist.root")
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
  TProfile *gRatVsT[20];
  for (Int_t i=0;i<20;++i) {
    gLedVsT[i] = new TProfile("","Led info;T;",1000,15,35);
    gLedVsT[i]->SetMarkerStyle(20);
    gLedVsT[i]->SetMarkerSize(1.2);
    gLedVsT[i]->SetMarkerColor(i+1);
    gLedVsT[i]->SetLineColor(i+1);
    gLedVsT[i]->SetLineWidth(3);
    gLedVsT[i]->SetName(Form("ledsm%d",i));
    gRatVsT[i] = new TProfile("","Led/LedMon;T",1000,15,35);
    gRatVsT[i]->SetMarkerStyle(20);
    gRatVsT[i]->SetMarkerSize(1.2);
    gRatVsT[i]->SetMarkerColor(i+1);
    gRatVsT[i]->SetLineColor(i+1);
    gRatVsT[i]->SetLineWidth(3);
    gRatVsT[i]->SetName(Form("ledovermonsm%d",i));
  }

  Int_t Nev=tt->GetEntries();
  for (Int_t i=0;i<Nev;++i) {
    tt->GetEvent(i);
    cout << info->fRunNo << endl;
    TClonesArray &cells = info->fCells;
    for (Int_t j=0;j<cells.GetEntries();++j) {
      TCalCell *cell = static_cast<TCalCell*>(cells.At(j));
      Int_t sm    = cell->fSM;
      Double_t ledM = cell->fLedM;
      Double_t ledR = cell->fLedR;
      Double_t monM = cell->fMonM;
      Double_t monR = cell->fMonR;
      Double_t locT = cell->fLocT;
      if ((locT<5)||(locT>45))
	continue;
      if ((ledM<=0)||(ledR<=0))
	continue;
      Double_t w = ledR;
      TProfile *h = gLedVsT[sm];
      h->Fill(locT,ledM,w);
      if ((monM<=0)||(monR<=0))
	continue;
      Double_t w2=TMath::Sqrt(ledR*ledR+monM*monM);
      TProfile *h2 = gRatVsT[sm];
      h2->Fill(locT,ledM/monM,w2);
    }
  }

  TFile *out = TFile::Open(ofile,"recreate");
  if (!out) 
    out = TFile::Open("dummyfile.root","update");
  for (Int_t i=0;i<20;++i) {
    gLedVsT[i]->Write();
    gRatVsT[i]->Write();
  }
}
