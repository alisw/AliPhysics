void AliTPCDigitsDisplay(int sec, int row,
                 int max_t_chan=500, float min_t=0., float max_t=500.,
                 int max_p_chan=200, float min_p=0., float max_p=200.)
{
// Dynamically link some shared libs
   if (gClassTable->GetID("AliRun") < 0) {
      gROOT->LoadMacro("loadlibs.C");
      loadlibs();
   } else {
      delete gAlice;
      gAlice=0;
   }

// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("galice.root");
   if (f) f->Close();
   f = new TFile("galice.root");

   TClonesArray *fDigits=new TClonesArray("AliTPCdigit",10000);
   TTree *t=(TTree*)f->Get("TreeD0_Param1");
   t->GetBranch("Digits")->SetAddress(&fDigits);
   Int_t sectors_by_rows=(Int_t)t->GetEntries();
   for (Int_t n=0; n<sectors_by_rows; n++) {
      if (!t->GetEvent(n)) continue;
      AliTPCdigit *dig=(AliTPCdigit*)fDigits->UncheckedAt(0);

      if (sec  < dig->fSector) break;
      if (sec != dig->fSector) continue;
      if (row != dig->fPadRow) continue;

      char s[80];
      sprintf(s,"Sector %d   Row %d\n",sec,row);
      TH2F *h = new TH2F("h",s,max_t_chan,min_t,max_t,max_p_chan,min_p,max_p);
      Int_t ndigits=fDigits->GetEntriesFast();
      for (Int_t ndig=0; ndig<ndigits; ndig++) {
         dig=(AliTPCdigit*)fDigits->UncheckedAt(ndig);
         if (dig->fSignal < 0) continue; //cluster finder threshold
         h->Fill(dig->fTime,dig->fPad,dig->fSignal);
      }
      h->SetMaximum(200);
      gStyle->SetOptStat(0);
      TCanvas *c1=new TCanvas("c1","TPC digits display",0,0,1110,680);
      TPad *p1=new TPad("p1","",0,0,1,0.5);
      p1->Draw();
      TPad *p2=new TPad("p2","",0,0.5,1,1);
      p2->Draw();
      p2->cd();
      h->Draw("lego");
      p1->cd();
      h->Draw("colz");
   }
}

