/****************************************************************************
 *           Very important, delicate and rather obscure macro.             *
 *                                                                          *
 *                  Creates list of "findable" V0s,                         *
 *             calculates efficiency, resolutions etc.                      *
 *                                                                          *
 *   Origin: I.Belikov, IReS, Strasbourg, Jouri.Belikov@cern.ch             *
 ****************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
  #include "Riostream.h"
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

  #include "AliRun.h"
  #include "AliPDG.h"
  #include "AliV0vertex.h"
#endif

struct GoodVertex {
  Int_t nlab,plab;
  Int_t code;
  Float_t px,py,pz;
  Float_t x,y,z;
};
Int_t good_vertices(GoodVertex *gt, Int_t max);

extern AliRun *gAlice;

Int_t AliV0Comparison(Int_t code=310) { //Lambda=3122, LambdaBar=-3122
   cerr<<"Doing comparison...\n";

   TStopwatch timer;

   const Double_t V0window=0.05;
   Double_t ptncut=0.13, ptpcut=0.13, kinecut=0.03; 
   Double_t V0mass=0.497672, V0width=0.020;
   switch (code) {
   case kK0Short: 
        break;
   case kLambda0:    
        V0mass=1.115683; V0width=0.015; ptpcut=0.50; kinecut=0.002; 
        break;
   case kLambda0Bar: 
        V0mass=1.115683; V0width=0.015; ptncut=0.50; kinecut=0.002;
        break;
   default: cerr<<"Invalid PDG code !\n"; return 1;
   }

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
   TH1F *hgood=new TH1F("hgood","Good Vertices",nchan,pmin,pmax);    
   TH1F *hfound=new TH1F("hfound","Found Vertices",nchan,pmin,pmax);
   TH1F *hfake=new TH1F("hfake","Fake Vertices",nchan,pmin,pmax);
   TH1F *hg=new TH1F("hg","Efficiency for Good Vertices",nchan,pmin,pmax);
   hg->SetLineColor(4); hg->SetLineWidth(2);
   TH1F *hf=new TH1F("hf","Probability of Fake Vertices",nchan,pmin,pmax);
   hf->SetFillColor(1); hf->SetFillStyle(3013); hf->SetLineWidth(2);

   Double_t mmin=V0mass-V0window, mmax=V0mass+V0window;
   TH1F *v0s =new TH1F("v0s","V0s Effective Mass",40, mmin, mmax);
   v0s->SetXTitle("(GeV)");
   v0s->SetLineColor(4); v0s->SetLineWidth(4);
   TH1F *v0sf =new TH1F("v0sf","Fake V0s Effective Mass",40, mmin, mmax);
   v0sf->SetXTitle("(GeV)"); v0sf->SetFillColor(6);


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

   /*** Load reconstructed vertices ***/
   TObjArray varray(1000);
   itsl->LoadV0s();
   TTree *vTree=itsl->TreeV0();
   TBranch *branch=vTree->GetBranch("vertices");
   Int_t nentr=(Int_t)vTree->GetEntries();
   for (Int_t i=0; i<nentr; i++) {
       AliV0vertex *iovertex=new AliV0vertex; branch->SetAddress(&iovertex);
       vTree->GetEvent(i);
       varray.AddLast(iovertex);
   }

   /*** Check if the file with the "good" vertices exists ***/
   GoodVertex gv[1000];
   Int_t ngood=0;
   ifstream in("good_vertices");
   if (in) {
      cerr<<"Reading good vertices...\n";
      while (in>>gv[ngood].nlab>>gv[ngood].plab>>gv[ngood].code>>
                 gv[ngood].px>>gv[ngood].py>>gv[ngood].pz>>
                 gv[ngood].x >>gv[ngood].y >>gv[ngood].z) {
         ngood++;
         cerr<<ngood<<'\r';
         if (ngood==1000) {
            cerr<<"Too many good vertices !\n";
            break;
         }
      }
      if (!in.eof()) cerr<<"Read error (good_vertices) !\n";
   } else {
     /*** generate a file with the "good" vertices ***/
      cerr<<"Marking good vertices (this will take a while)...\n";
      ngood=good_vertices(gv,1000);
      ofstream out("good_vertices");
      if (out) {
         for (Int_t ngd=0; ngd<ngood; ngd++)            
	    out<<gv[ngd].nlab<<' '<<gv[ngd].plab<<' '<<gv[ngd].code<<' '<<
                 gv[ngd].px<<' '<<gv[ngd].py<<' '<<gv[ngd].pz<<' '<<
                 gv[ngd].x <<' '<<gv[ngd].y <<' '<<gv[ngd].z <<endl;
      } else cerr<<"Can not open file (good_vertices) !\n";
      out.close();
   }


   Double_t pxg=0.,pyg=0.,ptg=0.;
   Int_t nlab=-1, plab=-1;
   Int_t i;
   for (i=0; i<nentr; i++) {
       AliV0vertex *vertex=(AliV0vertex*)varray.UncheckedAt(i);
       nlab=TMath::Abs(vertex->GetNindex());
       plab=TMath::Abs(vertex->GetPindex());

       /** Kinematical cuts **/
       Double_t pxn,pyn,pzn; vertex->GetNPxPyPz(pxn,pyn,pzn); 
       Double_t ptn=TMath::Sqrt(pxn*pxn + pyn*pyn);
       if (ptn < ptncut) continue;
       Double_t pxp,pyp,pzp; vertex->GetPPxPyPz(pxp,pyp,pzp); 
       Double_t ptp=TMath::Sqrt(pxp*pxp + pyp*pyp);
       if (ptp < ptpcut) continue;
       Double_t kine=vertex->ChangeMassHypothesis(code);
       //if (TMath::Abs(kine)>kinecut) continue;

       Double_t mass=vertex->GetEffMass();
       v0s->Fill(mass);
       v0sf->Fill(mass);

       Int_t j;
       for (j=0; j<ngood; j++) {
          if (gv[j].code != vertex->GetPdgCode()) continue;
          if (gv[j].nlab == nlab)
          if (gv[j].plab == plab) break;
       }

       if (TMath::Abs(mass-V0mass)>V0width) continue;

       Double_t px,py,pz; vertex->GetPxPyPz(px,py,pz);
       Double_t pt=TMath::Sqrt(px*px+py*py);

       if (j==ngood) {
          hfake->Fill(pt);
          cerr<<"Fake vertex: ("<<nlab<<","<<plab<<")\n";
          continue;
       }
       v0sf->Fill(mass,-1);

       pxg=gv[j].px; pyg=gv[j].py; ptg=TMath::Sqrt(pxg*pxg+pyg*pyg);
       Double_t phig=TMath::ATan2(pyg,pxg), phi=TMath::ATan2(py,px);
       Double_t lamg=TMath::ATan2(gv[j].pz,ptg), lam=TMath::ATan2(pz,pt);
       hp->Fill((phi - phig)*1000.);
       hl->Fill((lam - lamg)*1000.);
       hpt->Fill((1/pt - 1/ptg)/(1/ptg)*100.);

       Double_t x,y,z; vertex->GetXYZ(x,y,z);
       hx->Fill((x-gv[j].x)*10);
       hy->Fill((y-gv[j].y)*10);
       hz->Fill((z-gv[j].z)*10);

       hfound->Fill(ptg);
       gv[j].nlab=-1;

   }
   for (i=0; i<ngood; i++) {
      if (gv[i].code != code) continue;
      pxg=gv[i].px; pyg=gv[i].py; ptg=TMath::Sqrt(pxg*pxg+pyg*pyg);
      hgood->Fill(ptg);
      nlab=gv[i].nlab; plab=gv[i].plab;
      if (nlab < 0) continue;
      cerr<<"Vertex ("<<nlab<<','<<plab<<") has not been found !\n";
   }

   varray.Delete();

   Stat_t ng=hgood->GetEntries();
   Stat_t nf=hfound->GetEntries();

   cerr<<"Number of found vertices: "<<nentr<<" ("<<nf<<" in the peak)\n";
   cerr<<"Number of good vertices: "<<ng<<endl;

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
   hg->SetYTitle("V0 reconstruction efficiency");
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
   TText *text = new TText(0.461176,0.248448,"Fake vertices");
   text->SetTextSize(0.05);
   text->Draw();
   text = new TText(0.453919,1.11408,"Good vertices");
   text->SetTextSize(0.05);
   text->Draw();


   c2->cd(2);
   gPad->SetFillColor(42); gPad->SetFrameFillColor(10);
   //v0s->Fit("gaus","","",V0mass-V0width,V0mass+V0width);
   v0s->Draw();
   v0sf->Draw("same");
   Double_t max=v0s->GetMaximum();
   TLine *line3 = new TLine(V0mass-V0width,0.,V0mass-V0width,max);
   line3->Draw("same");
   TLine *line4 = new TLine(V0mass+V0width,0.,V0mass+V0width,max);
   line4->Draw("same");

   timer.Stop(); timer.Print();

   delete rl;

   return 0;
}

Int_t good_vertices(GoodVertex *gv, Int_t max) {
   Int_t nv=0;
   /*** Get information about the cuts ***/
   Double_t r2min=0.9*0.9;
   Double_t r2max=2.9*2.9;

   /*** Get labels of the "good" tracks ***/
   Double_t dd; Int_t id, label[15000], ngt=0;
   ifstream in("good_tracks_its");
   if (!in) {
     cerr<<"Can't open the file good_tracks_its \n";
     return nv;
   }
   while (in>>label[ngt]>>id>>dd>>dd>>dd>>dd>>dd>>dd) {
     ngt++;
     if (ngt>=15000) {
       cerr<<"Too many good ITS tracks !\n";
       return nv;
     }
   }   
   if (!in.eof()) {
      cerr<<"Read error (good_tracks_its) !\n";
      return nv;
   }

   /*** Get an access to the kinematics ***/
   AliRunLoader *rl = 
        AliRunLoader::GetRunLoader(AliConfig::fgkDefaultEventFolderName);
   if (rl == 0x0) {
     ::Fatal("AliV0Comparison.C::good_vertices","Can not find Run Loader !");
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
      TParticle *p0=gAlice->GetMCApp()->Particle(np);

      /*** only these V0s are "good" ***/
      Int_t code=p0->GetPdgCode();
      if (code!=kK0Short)
      if (code!=kLambda0)
      if (code!=kLambda0Bar) continue; 

      /*** daughter tracks should be "good" ***/
      Int_t plab=p0->GetFirstDaughter(), nlab=p0->GetLastDaughter();
      if (nlab==plab) continue;
      if (nlab<0) continue;
      if (plab<0) continue;
      Int_t i;
      for (i=0; i<ngt; i++) if (label[i]==nlab) break;
      if (i==ngt) continue;
      for (i=0; i<ngt; i++) if (label[i]==plab) break;
      if (i==ngt) continue;

      /*** fiducial volume ***/
      TParticle *p=gAlice->GetMCApp()->Particle(nlab);
      Double_t x=p->Vx(), y=p->Vy(), z=p->Vz(), r2=x*x+y*y;
      if (r2<r2min) continue;
      if (r2>r2max) continue;
       
      if (gAlice->GetMCApp()->Particle(plab)->GetPDG()->Charge() < 0.) {
         i=plab; plab=nlab; nlab=i;
      }
      
      gv[nv].code=code;
      gv[nv].plab=plab; gv[nv].nlab=nlab;
      gv[nv].px=p0->Px(); gv[nv].py=p0->Py(); gv[nv].pz=p0->Pz();
      gv[nv].x=x; gv[nv].y=y; gv[nv].z=z;
      nv++;
   }
 
   delete rl;

   return nv;
}
