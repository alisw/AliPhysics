/****************************************************************************
 *           Very important, delicate and rather obscure macro.             *
 *                                                                          *
 *               Creates list of "findable" cascades,                       *
 *             calculates efficiency, resolutions etc.                      *
 *                                                                          *
 *  Origin: Christian Kuhn, IReS, Strasbourg, christian.kuhn@ires.in2p3.fr  *
 ****************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <Riostream.h>
  #include <fstream.h>

  #include "AliRun.h"
  #include "AliHeader.h"
  #include "AliRunLoader.h"
  #include "AliITSLoader.h"

  #include "TH1.h"
  #include "TFile.h"
  #include "TTree.h"
  #include "TObjArray.h"
  #include "TStyle.h"
  #include "TCanvas.h"
  #include "TLine.h"
  #include "TText.h"
  #include "TParticle.h"
  #include "TStopwatch.h"
  #include "TPDGCode.h"

  #include "AliRun.h"
  #include "AliPDG.h"
  #include "AliCascadeVertex.h"
#endif

struct GoodCascade {
  Int_t nlab,plab;   // V0's daughter labels
  Int_t blab;        // Bachelor label
  Int_t code;
  Float_t px,py,pz;
  Float_t x,y,z;
};
Int_t good_cascades(GoodCascade *gt, Int_t max);

extern AliRun *gAlice;

Int_t AliCascadeComparison(Int_t code=3312) {
  //code= 3312; //kXiMinus 
  //code=-3312; //kXiPlusBar
  //code= 3334; //kOmegaMinus
  //code=-3334; //kOmegaPlusBar

   cerr<<"Doing comparison...\n";

   TStopwatch timer;

   const Double_t cascadeWindow=0.05, cascadeWidth=0.015; 
   Double_t ptncut=0.12, ptpcut=0.33, kine0cut=0.003;
   Double_t ptbcut=0.11, kinecut=0.002;
   Double_t cascadeMass=1.32131;
   switch (code) {
   case kXiMinus:
        break;
   case kXiPlusBar:
        ptncut=0.33; ptpcut=0.12;
        break;
   case kOmegaMinus: 
        cascadeMass=1.67245;
        kine0cut=0.001;
        ptbcut=0.22; kinecut=0.006;
        break; 
   case kOmegaPlusBar:
        cascadeMass=1.67245; 
        kine0cut=0.001;
        ptncut=0.33; ptpcut=0.12;
        ptbcut=0.22; kinecut=0.006;
        break;
   default: cerr<<"Invalid PDG code !\n"; return 1;
   }

   /*** create some histograms ***/
   TH1F *hp=new TH1F("hp","Angular Resolution",30,-30.,30.); //phi resolution 
   hp->SetXTitle("(mrad)"); hp->SetFillColor(2);
   TH1F *hl=new TH1F("hl","Lambda Resolution",30,-30,30);
   hl->SetXTitle("(mrad)"); hl->SetFillColor(1); hl->SetFillStyle(3013); 
   TH1F *hpt=new TH1F("hpt","Relative Pt Resolution",30,-10.,10.); 
   hpt->SetXTitle("(%)"); hpt->SetFillColor(2); 

   TH1F *hx=new TH1F("hx","Position Resolution (X,Y)",30,-3.,3.); //x res. 
   hx->SetXTitle("(mm)"); hx->SetFillColor(6);
   TH1F *hy=new TH1F("hy","Position Resolution (Y)",30,-3.,3.);   //y res
   hy->SetXTitle("(mm)"); hy->SetFillColor(1); hy->SetFillStyle(3013);
   TH1F *hz=new TH1F("hz","Position Resolution (Z)",30,-3.,3.);   //z res. 
   hz->SetXTitle("(mm)"); hz->SetFillColor(6);

   Double_t pmin=0.2, pmax=4.2; Int_t nchan=20;
   TH1F *hgood=new TH1F("hgood","Good Cascades",nchan,pmin,pmax);    
   TH1F *hfound=new TH1F("hfound","Found Cascades",nchan,pmin,pmax);
   TH1F *hfake=new TH1F("hfake","Fake Cascades",nchan,pmin,pmax);
   TH1F *hg=new TH1F("hg","Efficiency for Good Cascades",nchan,pmin,pmax);
   hg->SetLineColor(4); hg->SetLineWidth(2);
   TH1F *hf=new TH1F("hf","Probability of Fake Cascades",nchan,pmin,pmax);
   hf->SetFillColor(1); hf->SetFillStyle(3013); hf->SetLineWidth(2);

   Double_t mmin=cascadeMass-cascadeWindow, mmax=cascadeMass+cascadeWindow;
   TH1F *cs =new TH1F("cs","Cascade Effective Mass",40, mmin, mmax);
   cs->SetXTitle("(GeV)");
   cs->SetLineColor(4); cs->SetLineWidth(4);
   TH1F *csf =new TH1F("csf","Fake Cascade Effective Mass",40, mmin, mmax);
   csf->SetXTitle("(GeV)"); csf->SetFillColor(6);

   if (gAlice) {
      delete gAlice->GetRunLoader();
      delete gAlice; 
      gAlice=0;
   }   
   AliRunLoader *rl = AliRunLoader::Open("galice.root");
   if (!rl) {
       cerr<<"AliV0Comparison.C :Can't start sesion !\n";
       return 1;
   }
   AliITSLoader* itsl = (AliITSLoader*)rl->GetLoader("ITSLoader");
   if (itsl == 0x0) {
       cerr<<"AliV0Comparison.C : Can not find the ITSLoader\n";
       delete rl;
       return 2;
   }

   /*** Load reconstructed cascades ***/
   TObjArray carray(1000);
   itsl->LoadCascades();
   TTree *xTree=itsl->TreeX();
   TBranch *branch=xTree->GetBranch("cascades");
   Int_t nentr=(Int_t)xTree->GetEntries();
   for (Int_t i=0; i<nentr; i++) {
       AliCascadeVertex *iovertex=new AliCascadeVertex; 
       branch->SetAddress(&iovertex);
       xTree->GetEvent(i);
       carray.AddLast(iovertex);
   }

   /*** Check if the file with the "good" cascades exists ***/
   GoodCascade gc[100];
   Int_t ngood=0;
   ifstream in("good_cascades");
   if (in) {
      cerr<<"Reading good cascades...\n";
      while (in>>gc[ngood].nlab>>gc[ngood].plab>>
	         gc[ngood].blab>>gc[ngood].code>>
                 gc[ngood].px>>gc[ngood].py>>gc[ngood].pz>>
                 gc[ngood].x >>gc[ngood].y >>gc[ngood].z) {
         ngood++;
         cerr<<ngood<<'\r';
         if (ngood==100) {
            cerr<<"Too many good cascades !\n";
            break;
         }
      }
      if (!in.eof()) cerr<<"Read error (good_cascades) !\n";
   } else {
     /*** generate a file with the "good" cascades ***/
      cerr<<"Marking good cascades (this will take a while)...\n";
      ngood=good_cascades(gc,100);
      ofstream out("good_cascades");
      if (out) {
         for (Int_t ngd=0; ngd<ngood; ngd++)            
            out<<gc[ngd].nlab<<' '<<gc[ngd].plab<<' '<<
	         gc[ngd].blab<<' '<<gc[ngd].code<<' '<<
                 gc[ngd].px<<' '<<gc[ngd].py<<' '<<gc[ngd].pz<<' '<<
                 gc[ngd].x <<' '<<gc[ngd].y <<' '<<gc[ngd].z <<endl;
      } else cerr<<"Can not open file (good_cascades) !\n";
      out.close();
   }

   Double_t pxg=0.,pyg=0.,ptg=0.;
   Int_t nlab=-1, plab=-1, blab=-1;
   Int_t i;
   for (i=0; i<nentr; i++) {
       AliCascadeVertex *cascade=(AliCascadeVertex*)carray.UncheckedAt(i);
       nlab=TMath::Abs(cascade->GetNindex()); 
       plab=TMath::Abs(cascade->GetPindex());
       blab=TMath::Abs(cascade->GetBindex());

       /** Kinematical cuts **/
       Double_t pxn,pyn,pzn; cascade->GetNPxPyPz(pxn,pyn,pzn); 
       Double_t ptn=TMath::Sqrt(pxn*pxn + pyn*pyn);
       if (ptn < ptncut) continue;
       Double_t pxp,pyp,pzp; cascade->GetPPxPyPz(pxp,pyp,pzp); 
       Double_t ptp=TMath::Sqrt(pxp*pxp + pyp*pyp);
       if (ptp < ptpcut) continue;
       Double_t pxb,pyb,pzb; cascade->GetBPxPyPz(pxb,pyb,pzb); 
       Double_t ptb=TMath::Sqrt(pxb*pxb + pyb*pyb);
       if (ptb < ptbcut) continue;
       Double_t kine0;
       Double_t kine=cascade->ChangeMassHypothesis(kine0,code);
       if (TMath::Abs(kine0)>kine0cut) continue;
       //if (TMath::Abs(kine)>kinecut) continue;

       Double_t mass=cascade->GetEffMass();
       cs->Fill(mass);
       csf->Fill(mass);

       if (TMath::Abs(mass-cascadeMass)>cascadeWidth) continue;

       Int_t j;
       for (j=0; j<ngood; j++) {
          if (gc[j].code != cascade->GetPdgCode()) continue;
          if (gc[j].nlab == nlab)
          if (gc[j].plab == plab)
          if (gc[j].blab == blab) break;
       }

       Double_t px,py,pz; cascade->GetPxPyPz(px,py,pz);
       Double_t pt=TMath::Sqrt(px*px+py*py);

       if (j==ngood) {
          hfake->Fill(pt);
          cerr<<"Fake cascade: ("<<nlab<<","<<plab<<","<<blab<<")\n";
          continue;
       }
       csf->Fill(mass,-1);

       pxg=gc[j].px; pyg=gc[j].py; ptg=TMath::Sqrt(pxg*pxg+pyg*pyg);
       Double_t phig=TMath::ATan2(pyg,pxg), phi=TMath::ATan2(py,px);
       Double_t lamg=TMath::ATan2(gc[j].pz,ptg), lam=TMath::ATan2(pz,pt);
       hp->Fill((phi - phig)*1000.);
       hl->Fill((lam - lamg)*1000.);
       hpt->Fill((1/pt - 1/ptg)/(1/ptg)*100.);

       Double_t x,y,z; cascade->GetXYZ(x,y,z);
       hx->Fill((x-gc[j].x)*10);
       hy->Fill((y-gc[j].y)*10);
       hz->Fill((z-gc[j].z)*10);

       hfound->Fill(ptg);
       gc[j].nlab=-1;

   }
   for (i=0; i<ngood; i++) {
      if (gc[i].code != code) continue;
      pxg=gc[i].px; pyg=gc[i].py; ptg=TMath::Sqrt(pxg*pxg+pyg*pyg);
      hgood->Fill(ptg);
      nlab=gc[i].nlab; plab=gc[i].plab; blab=gc[i].blab;
      if (nlab < 0) continue;
     cerr<<"Cascade ("<<nlab<<','<<plab<<","<<blab<<") has not been found !\n";
   }

   carray.Delete();

   Stat_t ng=hgood->GetEntries();
   Stat_t nf=hfound->GetEntries();

   cerr<<"Number of found cascades: "<<nentr<<" ("<<nf<<" in the peak)\n";
   cerr<<"Number of good cascades: "<<ng<<endl;

   if (ng!=0) 
      cerr<<"Integral efficiency is about "<<nf/ng*100.<<" %\n";

   gStyle->SetOptStat(111110);
   gStyle->SetOptFit(1);

   TCanvas *c1=new TCanvas("c1","",0,0,580,610);
   c1->Divide(2,2);

   c1->cd(1);
   gPad->SetFillColor(42); gPad->SetFrameFillColor(10); 
   //hp->Fit("gaus");
   hp->Draw();
   hl->Draw("same"); c1->cd();

   c1->cd(2);
   gPad->SetFillColor(42); gPad->SetFrameFillColor(10);
   //hpt->Fit("gaus"); c1->cd();
   hpt->Draw(); c1->cd();

   c1->cd(3);
   gPad->SetFillColor(42); gPad->SetFrameFillColor(10); 
   //hx->Fit("gaus"); 
   hx->Draw(); 
   hy->Draw("same"); c1->cd();

   c1->cd(4);
   gPad->SetFillColor(42); gPad->SetFrameFillColor(10);
   //hz->Fit("gaus");
   hz->Draw();


   TCanvas *c2=new TCanvas("c2","",600,0,580,610);
   c2->Divide(1,2);

   c2->cd(1);
   gPad->SetFillColor(42); gPad->SetFrameFillColor(10);
   hfound->Sumw2(); hgood->Sumw2(); hfake->Sumw2();
   hg->Divide(hfound,hgood,1,1.,"b");
   hf->Divide(hfake,hgood,1,1.,"b");
   hg->SetMaximum(1.4);
   hg->SetYTitle("Cascade reconstruction efficiency");
   hg->SetXTitle("Pt (GeV/c)");
   hg->Draw();

   TLine *line1 = new TLine(pmin,1.0,pmax,1.0); line1->SetLineStyle(4);
   line1->Draw("same");
   TLine *line2 = new TLine(pmin,0.9,pmax,0.9); line2->SetLineStyle(4);
   line2->Draw("same");

   hf->SetFillColor(1);
   hf->SetFillStyle(3013);
   hf->SetLineColor(2);
   hf->SetLineWidth(2);
   hf->Draw("histsame");
   TText *text = new TText(0.461176,0.248448,"Fake cascades");
   text->SetTextSize(0.05);
   text->Draw();
   text = new TText(0.453919,1.11408,"Good cascades");
   text->SetTextSize(0.05);
   text->Draw();


   c2->cd(2);
   gPad->SetFillColor(42); gPad->SetFrameFillColor(10);
   //cs->Fit("gaus","","",cascadeMass-cascadeWidth,cascadeMass+cascadeWidth);
   cs->Draw();
   csf->Draw("same");
   Double_t max=cs->GetMaximum();
   TLine *line3 = 
      new TLine(cascadeMass-cascadeWidth,0.,cascadeMass-cascadeWidth,max);
   line3->Draw("same");
   TLine *line4 = 
      new TLine(cascadeMass+cascadeWidth,0.,cascadeMass+cascadeWidth,max);
   line4->Draw("same");

   timer.Stop(); timer.Print();

   delete rl;

   return 0;
}


Int_t good_cascades(GoodCascade *gc, Int_t max) {
   Int_t nc=0;
   /*** Get information about the cuts ***/
   Double_t r2min=0.9*0.9;
   Double_t r2max=2.9*2.9;

   /*** Get labels of the "good" tracks ***/
   Double_t dd; Int_t id, label[15000], ngt=0;
   ifstream in("good_tracks_its");
   if (!in) {
     cerr<<"Can't open the file good_tracks_its \n";
     return nc;
   }
   while (in>>label[ngt]>>id>>dd>>dd>>dd>>dd>>dd>>dd) {
     ngt++;
     if (ngt>=15000) {
       cerr<<"Too many good ITS tracks !\n";
       return nc;
     }
   }   
   if (!in.eof()) {
      cerr<<"Read error (good_tracks_its) !\n";
      return nc;
   }

   /*** Get an access to the kinematics ***/
   AliRunLoader *rl =
        AliRunLoader::GetRunLoader(AliConfig::fgkDefaultEventFolderName);
   if (rl == 0x0) {
  ::Fatal("AliCascadeComparison.C::good_cascades","Can not find Run Loader !");
   }

   AliITSLoader* itsl = (AliITSLoader*)rl->GetLoader("ITSLoader");
   if (itsl == 0x0) {
       cerr<<"AliITSComparisonV2.C : Can not find TPCLoader\n";
       delete rl;
       return 1;
   }
   rl->LoadgAlice();
   rl->LoadHeader();
   rl->LoadKinematics();
   Int_t np = rl->GetHeader()->GetNtrack();

   while (np--) {
      cerr<<np<<'\r';
      TParticle *cp=gAlice->GetMCApp()->Particle(np);

      /*** only these cascades are "good" ***/
      Int_t code=cp->GetPdgCode();
      if (code!=kXiMinus)    if (code!=kXiPlusBar)
      if (code!=kOmegaMinus) if (code!=kOmegaPlusBar) continue; 

      /*** daughter tracks must be "good" ***/
      Int_t v0lab=cp->GetFirstDaughter(), blab=cp->GetLastDaughter();
      if (v0lab==blab) continue;
      if (v0lab<0) continue;
      if (blab<0) continue;

      TParticle *p0=gAlice->GetMCApp()->Particle(v0lab);
      TParticle *bp=gAlice->GetMCApp()->Particle(blab);
      if ((p0->GetPdgCode()!=kLambda0) && (p0->GetPdgCode()!=kLambda0Bar)) {
         TParticle *p=p0; p0=bp; bp=p;
         Int_t i=v0lab; v0lab=blab; blab=i;         
         if ((p0->GetPdgCode()!=kLambda0) && (p0->GetPdgCode()!=kLambda0Bar))
                                                                   continue;
      }

      /** is the bachelor "good" ? **/
      Int_t i;
      for (i=0; i<ngt; i++) if (label[i]==blab) break;
      if (i==ngt) continue; 

      /** is the V0 "good" ? **/
      Int_t plab=p0->GetFirstDaughter(), nlab=p0->GetLastDaughter();
      if (nlab==plab) continue;
      if (nlab<0) continue;
      if (plab<0) continue;

      for (i=0; i<ngt; i++) if (label[i]==nlab) break;
      if (i==ngt) continue;
      for (i=0; i<ngt; i++) if (label[i]==plab) break;
      if (i==ngt) continue;

      /*** fiducial volume ***/
      Double_t x=bp->Vx(), y=bp->Vy(), r2=x*x+y*y; //bachelor
      if (r2<r2min) continue;
      if (r2>r2max) continue;
      TParticle *pp=gAlice->GetMCApp()->Particle(plab);
      x=pp->Vx(); y=pp->Vy(); r2=x*x+y*y;                      //V0
      if (r2<r2min) continue;
      if (r2>r2max) continue;

      if (gAlice->GetMCApp()->Particle(plab)->GetPDG()->Charge() < 0.) {
         i=plab; plab=nlab; nlab=i;
      }

      gc[nc].code=code;
      gc[nc].plab=plab;   gc[nc].nlab=nlab; gc[nc].blab=blab;
      gc[nc].px=cp->Px(); gc[nc].py=cp->Py(); gc[nc].pz=cp->Pz();
      gc[nc].x=bp->Vx(); gc[nc].y=bp->Vy(); gc[nc].z=bp->Vz();
      nc++;

   }

   delete rl;

   return nc;
}
