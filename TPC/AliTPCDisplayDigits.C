Int_t AliTPCDisplayDigits(int sec, int row, int lab=-1,
                 int max_t_chan=500, float min_t=0., float max_t=500.,
                 int max_p_chan=150, float min_p=0., float max_p=150.)
{
   cerr<<"Displaying digits...\n";

// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *f = TFile::Open("galice.root");
   if (!f->IsOpen()) {cerr<<"Can't open galice.root !\n"; return 1;}

   AliTPCParam *par=(AliTPCParam *)f->Get("75x40_100x60");

   char s[80];
   sprintf(s,"Sector %d   Row %d\n",sec,row);
   TH2F *h = new TH2F("h",s,max_t_chan,min_t,max_t,max_p_chan,min_p,max_p);

   TTree *t=(TTree*)f->Get("TreeD_75x40_100x60");
   AliSimDigits dummy, *digit=&dummy;
   t->GetBranch("Segment")->SetAddress(&digit);
   for (int i=0; i<t->GetEntries(); i++) {
       t->GetEvent(i);
       int ss,rr;
       par->AdjustSectorRow(digit->GetID(),ss,rr);
       if (ss==sec && rr==row) goto ok;
   }
   return 2;

ok:
   int imax=0, jmax=0, qmax=0;
   digit->First();
   do {
      Short_t dig=digit->CurrentDigit();
      int i=digit->CurrentRow(), j=digit->CurrentColumn();
      if (lab >= 0) {
         int lab0=digit->GetTrackID(i,j,0);
         int lab1=digit->GetTrackID(i,j,1);
         int lab2=digit->GetTrackID(i,j,2);
         if (lab0!=lab) if (lab1!=lab) if (lab2!=lab) continue;
         if (dig>qmax) {imax=i; jmax=j; qmax=dig;}
         cerr<<lab0<<' '<<lab1<<' '<<lab2<<endl;
      }
      h->Fill(i,j,dig);
   } while (digit->Next());
   if (qmax>0){cerr<<"Peak (time,pad,q) : "<<imax<<' '<<jmax<<' '<<qmax<<endl;}

   h->SetMaximum(100);
   gStyle->SetOptStat(0);
   TCanvas *c1=new TCanvas("ddisplay","TPC digits display",0,0,1110,680);
   TPad *p1=new TPad("p1","",0,0,1,0.5);
   p1->Draw();
   TPad *p2=new TPad("p2","",0,0.5,1,1);
   p2->Draw();
   p2->cd();
   h->DrawCopy("lego");
   p1->cd();
   h->DrawCopy("colz");

   c1->Modified(); c1->Update(); 

   f->Close();
   return 0;
}

