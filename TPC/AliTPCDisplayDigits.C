/// \file AliTPCDisplayDigits.C
///
/// \author I.Belikov, CERN, Jouri.Belikov@cern.ch

#ifndef __CINT__
#include <Riostream.h>
#include <TFile.h>
#include <TTree.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TH2.h>

#include "AliTPCParam.h"
#include "AliSimDigits.h"
#endif

Int_t AliTPCDisplayDigits(Int_t eventn, int sec, int row, int lab=-1,
                 int max_t_chan=500, float min_t=0., float max_t=500.,
                 int max_p_chan=150, float min_p=0., float max_p=150.)
{
   cerr<<"Displaying digits...\n";

// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *f = TFile::Open("rfio:galice.root");
   if (!f->IsOpen()) {cerr<<"Can't open galice.root !\n"; return 1;}

   
   
   AliTPCParam *par=(AliTPCParam *)f->Get("75x40_100x60_150x60");
   if (!par) { cerr<<"TPC parameters have not been found !\n"; return 2; }

   Char_t s[80];
   sprintf(s,"Sector %d   Row %d\n",sec,row);
   TH2F *h = new TH2F("h",s,max_t_chan,min_t,max_t,max_p_chan,min_p,max_p);
   
   Char_t  name[100];
   sprintf(name,"TreeD_75x40_100x60_150x60_%d",eventn);

   TTree *t=(TTree*)f->Get(name);
   if (!t) { cerr<<"TPC digits have not been found !\n"; return 3; }
   AliSimDigits dummy, *digit=&dummy;
   t->GetBranch("Segment")->SetAddress(&digit);
   Int_t sbr=(Int_t)t->GetEntries();
   for (Int_t i=0; i<sbr; i++) {
       if (!t->GetEvent(i)) continue;
       Int_t s,r;
       par->AdjustSectorRow(digit->GetID(),s,r);
       if (s==sec && r==row) goto ok;
   }
   return 4;
   
ok:

   Int_t imax=0, jmax=0, qmax=0;
   digit->First();
   do {
     //Short_t dig=digit->CurrentDigit();
      Int_t i=digit->CurrentRow(), j=digit->CurrentColumn();
      Short_t dig=digit->GetDigit(i,j);
      if (lab >= 0) {
         Int_t lab0=digit->GetTrackID(i,j,0);
         Int_t lab1=digit->GetTrackID(i,j,1);
         Int_t lab2=digit->GetTrackID(i,j,2);
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

