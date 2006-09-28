/****************************************************************************
 *           Very important, delicate and rather obscure macro.             *
 *                                                                          *
 *                  Creates list of "findable" V0s,                         *
 *             calculates efficiency, resolutions etc.                      *
 *                                                                          *
 *   Origin: I.Belikov, IReS, Strasbourg, Jouri.Belikov@cern.ch             *
 ****************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <TMath.h>
  #include <TError.h>
  #include <Riostream.h>
  #include <TH1F.h>
  #include <TTree.h>
  #include <TParticle.h>
  #include <TCanvas.h>
  #include <TLine.h>
  #include <TText.h>
  #include <TBenchmark.h>
  #include <TStyle.h>
  #include <TFile.h>
  #include <TROOT.h>

  #include "AliStack.h"
  #include "AliHeader.h"
  #include "AliTrackReference.h"
  #include "AliRunLoader.h"
  #include "AliRun.h"
  #include "AliESD.h"
#endif

Int_t GoodV0s(const Char_t *dir=".");

extern AliRun *gAlice;
extern TBenchmark *gBenchmark;
extern TROOT *gROOT;

static Int_t allgood=0;
static Int_t allfound=0;

Int_t AliV0Comparison(Int_t code=310, const Char_t *dir=".") { 
   //Lambda=3122, LambdaBar=-3122
   gBenchmark->Start("AliV0Comparison");

   ::Info("AliV0Comparison.C","Doing comparison...");


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

   TH1F *hp=(TH1F*)gROOT->FindObject("hp");
   if (!hp) hp=new TH1F("hp","PHI Resolution",30,-30.,30.);
   hp->SetXTitle("(mrad)"); hp->SetFillColor(2);

   TH1F *hl=(TH1F*)gROOT->FindObject("hl");
   if (!hl) hl=new TH1F("hl","LAMBDA Resolution",30,-30,30);
   hl->SetXTitle("(mrad)"); hl->SetFillColor(1); hl->SetFillStyle(3013);
 
   TH1F *hpt=(TH1F*)gROOT->FindObject("hpt");
   if (!hpt) hpt=new TH1F("hpt","Relative Pt Resolution",30,-10.,10.); 
   hpt->SetXTitle("(%)"); hpt->SetFillColor(2); 

   TH1F *hx=(TH1F*)gROOT->FindObject("hx");
   if (!hx) hx=new TH1F("hx","Position Resolution (X)",30,-3.,3.);
   hx->SetXTitle("(mm)"); hx->SetFillColor(6);

   TH1F *hy=(TH1F*)gROOT->FindObject("hy");
   if (!hy) hy=new TH1F("hy","Position Resolution (Y)",30,-3.,3.); 
   hy->SetXTitle("(mm)"); hy->SetFillColor(1); hy->SetFillStyle(3013);

   TH1F *hz=(TH1F*)gROOT->FindObject("hz");
   if (!hz) hz=new TH1F("hz","Position Resolution (Z)",30,-3.,3.);
   hz->SetXTitle("(mm)"); hz->SetFillColor(6);


   Double_t pmin=0.2, pmax=4.2; Int_t nchan=20;
   TH1F *hgood=(TH1F*)gROOT->FindObject("hgood");    
   if (!hgood) hgood=new TH1F("hgood","Good Vertices",nchan,pmin,pmax);
    
   TH1F *hfound=(TH1F*)gROOT->FindObject("hfound");
   if (!hfound) hfound=new TH1F("hfound","Found Vertices",nchan,pmin,pmax);

   TH1F *hfake=(TH1F*)gROOT->FindObject("hfake");
   if (!hfake) hfake=new TH1F("hfake","Fake Vertices",nchan,pmin,pmax);

   TH1F *hg=(TH1F*)gROOT->FindObject("hg");
   if (!hg) hg=new TH1F("hg","Efficiency for Good Vertices",nchan,pmin,pmax);
   hg->SetLineColor(4); hg->SetLineWidth(2);

   TH1F *hf=(TH1F*)gROOT->FindObject("hf");
   if (!hf) hf=new TH1F("hf","Probability of Fake Vertices",nchan,pmin,pmax);
   hf->SetFillColor(1); hf->SetFillStyle(3013); hf->SetLineWidth(2);


   Double_t mmin=V0mass-V0window, mmax=V0mass+V0window;
   TH1F *v0s=(TH1F*)gROOT->FindObject("v0s");
   if (!v0s) v0s=new TH1F("v0s","V0s Effective Mass",40, mmin, mmax);
   v0s->SetXTitle("(GeV)"); v0s->SetLineColor(4); v0s->SetLineWidth(4);

   TH1F *v0sf=(TH1F*)gROOT->FindObject("v0sf");
   if (!v0sf) v0sf=new TH1F("v0sf","Fake V0s Effective Mass",40, mmin, mmax);
   v0sf->SetXTitle("(GeV)"); v0sf->SetFillColor(6);


   Char_t fname[100];
   sprintf(fname,"%s/GoodV0s.root",dir);

   TFile *refFile=TFile::Open(fname,"old");
   if (!refFile || !refFile->IsOpen()) {
   ::Info("AliV0Comparison.C","Marking good V0s (will take a while)...");
     if (GoodV0s(dir)) {
        ::Error("AliV0Comparison.C","Can't generate the reference file !");
        return 1;
     }
   }
   refFile=TFile::Open(fname,"old");
   if (!refFile || !refFile->IsOpen()) {
     ::Error("AliV0Comparison.C","Can't open the reference file !");
     return 1;
   }   
  
   TTree *v0Tree=(TTree*)refFile->Get("v0Tree");
   if (!v0Tree) {
     ::Error("AliV0Comparison.C","Can't get the reference tree !");
     return 2;
   }
   TBranch *pbranch=v0Tree->GetBranch("positive");
   if (!pbranch) {
     ::Error("AliV0Comparison.C","Can't get the positive daughter branch !");
     return 3;
   }
   TClonesArray dummy("AliTrackReference",1000), *prefs=&dummy;
   pbranch->SetAddress(&prefs);

   TBranch *nbranch=v0Tree->GetBranch("negative");
   if (!nbranch) {
     ::Error("AliV0Comparison.C","Can't get the negative daughter branch !");
     return 4;
   }
   TClonesArray dumm("AliTrackReference",1000), *nrefs=&dumm;
   nbranch->SetAddress(&nrefs);

   
   sprintf(fname,"%s/AliESDs.root",dir);
   TFile *ef=TFile::Open(fname);
   if ((!ef)||(!ef->IsOpen())) {
      sprintf(fname,"%s/AliESDv0.root",dir);
      ef=TFile::Open(fname);
      if ((!ef)||(!ef->IsOpen())) {
         ::Error("AliV0Comparison.C","Can't open AliESDv0.root !");
         return 5;
      }
   }
   AliESD* event = new AliESD;
   TTree* esdTree = (TTree*) ef->Get("esdTree");
   if (!esdTree) {
      ::Error("AliV0Comparison.C", "no ESD tree found");
      return 6;
   }
   esdTree->SetBranchAddress("ESD", &event);


   //******* Loop over events *********
   Int_t e=0;
   while (esdTree->GetEvent(e)) {
      cout<<endl<<endl<<"********* Processing event number: "<<e<<"*******\n";
 
      Int_t nentr=event->GetNumberOfV0s();
      allfound+=nentr;

      if (v0Tree->GetEvent(e++)==0) {
         cerr<<"No reconstructable V0s !\n";
         continue;
      }

      Int_t ngood=prefs->GetEntriesFast(),ng=0; 

      Double_t pxg=0.,pyg=0.,pzg=0.,ptg=0.;
      Int_t nlab=-1, plab=-1;
      Int_t i;
      for (i=0; i<nentr; i++) {
          AliESDv0 *vertex=event->GetV0(i);

          if (vertex->GetOnFlyStatus()) continue;

          Int_t nidx=TMath::Abs(vertex->GetNindex());
          Int_t pidx=TMath::Abs(vertex->GetPindex());

          AliESDtrack *ntrack=event->GetTrack(nidx);
          AliESDtrack *ptrack=event->GetTrack(pidx);

          nlab=TMath::Abs(ntrack->GetLabel());
          plab=TMath::Abs(ptrack->GetLabel());

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

          AliTrackReference *nref=0, *pref=0;
          Int_t j;
          for (j=0; j<ngood; j++) {
              nref=(AliTrackReference*)nrefs->UncheckedAt(j);
              pref=(AliTrackReference*)prefs->UncheckedAt(j);
              if (nref->Label() == nlab)
              if (pref->Label() == plab) break;
          }

          if (TMath::Abs(mass-V0mass)>V0width) continue;

          Double_t px,py,pz; vertex->GetPxPyPz(px,py,pz);
          Double_t pt=TMath::Sqrt(px*px+py*py);

          if (j==ngood) {
             hfake->Fill(pt);
             cout<<"Fake vertex: ("<<nlab<<","<<plab<<")\n";
             continue;
          }
          v0sf->Fill(mass,-1);

          pxg=nref->Px()+pref->Px(); pyg=nref->Py()+pref->Py();
          pzg=nref->Pz()+pref->Pz(); 
          ptg=TMath::Sqrt(pxg*pxg+pyg*pyg);
          Double_t phig=TMath::ATan2(pyg,pxg), phi=TMath::ATan2(py,px);
          Double_t lamg=TMath::ATan2(pzg,ptg), lam=TMath::ATan2(pz,pt);
          hp->Fill((phi - phig)*1000.);
          hl->Fill((lam - lamg)*1000.);
          hpt->Fill((1/pt - 1/ptg)/(1/ptg)*100.);

          Double_t x,y,z; vertex->GetXYZ(x,y,z);
          hx->Fill((x - nref->X())*10);
          hy->Fill((y - nref->Y())*10);
          hz->Fill((z - nref->Z())*10);

          hfound->Fill(ptg);
          nref->SetLabel(-1);

      }
      for (i=0; i<ngood; i++) {
         AliTrackReference *nref=(AliTrackReference*)nrefs->UncheckedAt(i);
         AliTrackReference *pref=(AliTrackReference*)prefs->UncheckedAt(i);
         Int_t pdg=(Int_t)nref->GetLength();//this is the mother's PDG !
         if (code!=pdg) continue;
         ng++;
         pxg=nref->Px()+pref->Px(); pyg=nref->Py()+pref->Py(); 
         ptg=TMath::Sqrt(pxg*pxg+pyg*pyg);
         hgood->Fill(ptg);
         nlab=nref->Label(); plab=pref->Label();
         if (nlab < 0) continue;
         cout<<"Vertex ("<<nlab<<','<<plab<<") has not been found !\n";
      }
      allgood+=ng;

      cout<<"Number of found V0s : "<<nentr<<endl;
      cout<<"Number of \"good\" V0s : "<<ng<<endl;

      prefs->Clear();
      nrefs->Clear();

   } //**** End of the loop over events

   delete event;
   ef->Close();
   
   delete v0Tree;
   refFile->Close();

   Stat_t ng=hgood->GetEntries(), nf=hfound->GetEntries();
   if (ng!=0) cout<<"Integral efficiency is about "<<nf/ng*100.<<" %\n";
   cout<<"Total number of found V0s: "<<allfound<<" ("<<nf<<" in the peak)\n";
   cout<<"Total number of \"good\" V0s: "<<allgood<<endl;

   gStyle->SetOptStat(111110);
   gStyle->SetOptFit(1);

   TCanvas *c1=(TCanvas*)gROOT->FindObject("c1");
   if (!c1) {
      c1=new TCanvas("c1","",0,0,580,610);
      c1->Divide(2,2);
   }

   Int_t minc=33; 

   c1->cd(1);
   gPad->SetFillColor(42); gPad->SetFrameFillColor(10); 
   if (hp->GetEntries()<minc) hp->Draw(); else hp->Fit("gaus");
   hl->Draw("same"); c1->cd();

   c1->cd(2);
   gPad->SetFillColor(42); gPad->SetFrameFillColor(10);
   if (hpt->GetEntries()<minc) hpt->Draw(); else hpt->Fit("gaus");
   c1->cd();

   c1->cd(3);
   gPad->SetFillColor(42); gPad->SetFrameFillColor(10); 
   if (hx->GetEntries()<minc) hx->Draw(); else hx->Fit("gaus"); 
   hy->Draw("same"); c1->cd();

   c1->cd(4);
   gPad->SetFillColor(42); gPad->SetFrameFillColor(10);
   if (hz->GetEntries()<minc) hz->Draw(); else hz->Fit("gaus");

   c1->Update();   

   TCanvas *c2=(TCanvas*)gROOT->FindObject("c2");
   if (!c2) {
      c2=new TCanvas("c2","",600,0,580,610);
      c2->Divide(1,2);
   }

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
   if (v0s->GetEntries()<minc) v0s->Draw();
   else v0s->Fit("gaus","","",V0mass-V0width,V0mass+V0width);
   v0sf->Draw("same");
   Double_t max=v0s->GetMaximum();
   TLine *line3 = new TLine(V0mass-V0width,0.,V0mass-V0width,max);
   line3->Draw("same");
   TLine *line4 = new TLine(V0mass+V0width,0.,V0mass+V0width,max);
   line4->Draw("same");

   c2->Update();

   TFile fc("AliV0Comparison.root","RECREATE");
   c1->Write();
   c2->Write();
   fc.Close();

   gBenchmark->Stop("AliV0Comparison");
   gBenchmark->Show("AliV0Comparison");

   return 0;
}


Int_t GoodV0s(const Char_t *dir) {
   if (gAlice) { 
       delete gAlice->GetRunLoader();
       delete gAlice;//if everything was OK here it is already NULL
       gAlice = 0x0;
   }

   Char_t fname[100];
   sprintf(fname,"%s/galice.root",dir);

   AliRunLoader *rl = AliRunLoader::Open(fname,"COMPARISON");
   if (!rl) {
      ::Error("GoodTracksITS","Can't start session !");
      return 1;
   }

   rl->LoadgAlice();
   rl->LoadHeader();
   rl->LoadKinematics();


   Int_t nev=rl->GetNumberOfEvents();
   ::Info("GoodV0s","Number of events : %d\n",nev);  

   sprintf(fname,"%s/GoodTracksITS.root",dir);
   TFile *itsFile=TFile::Open(fname);
   if ((!itsFile)||(!itsFile->IsOpen())) {
       ::Error("GoodV0s","Can't open the GoodTracksITS.root !");
       delete rl;
       return 5; 
   }
   TClonesArray dum("AliTrackReference",1000), *itsRefs=&dum;
   TTree *itsTree=(TTree*)itsFile->Get("itsTree");
   if (!itsTree) {
       ::Error("GoodV0s","Can't get the ITS reference tree !");
       delete rl;
       return 6;
   }
   TBranch *itsBranch=itsTree->GetBranch("ITS");
   if (!itsBranch) {
      ::Error("GoodV0s","Can't get the ITS reference branch !");
      delete rl;
      return 7;
   }
   itsBranch->SetAddress(&itsRefs);


   sprintf(fname,"%s/GoodV0s.root",dir);
   TFile *v0File=TFile::Open(fname,"recreate");
   TClonesArray dummy("AliTrackReference",1000), *nrefs=&dummy;
   TClonesArray dumm("AliTrackReference",1000), *prefs=&dumm;
   TTree v0Tree("v0Tree","Tree with info about the reconstructable V0s");
   v0Tree.Branch("negative",&nrefs);
   v0Tree.Branch("positive",&prefs);


   /*** Some information about the cuts ***/
   Double_t r2min=0.9*0.9;
   Double_t r2max=2.9*2.9;



   //********  Loop over generated events
   for (Int_t e=0; e<nev; e++) {
      rl->GetEvent(e);  v0File->cd();

      Int_t np = rl->GetHeader()->GetNtrack();
      cout<<"Event "<<e<<" Number of particles: "<<np<<endl;

      itsTree->GetEvent(e);
      Int_t nk=itsRefs->GetEntriesFast();

      AliStack *stack=rl->Stack();

      AliTrackReference *nref=0, *pref=0;

      Int_t nv=0;
      while (np--) {
	//cerr<<np<<'\r';
         TParticle *p0=stack->Particle(np);

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
	 TParticle * part = stack->Particle(plab);
	 if (part) {
	   TParticlePDG * partPDG = part->GetPDG();
	     if (partPDG && partPDG->Charge() < 0.) {
	       i=plab; plab=nlab; nlab=i;
	     }
	 }

         for (i=0; i<nk; i++) {
	     nref=(AliTrackReference*)itsRefs->UncheckedAt(i);
             if (nref->Label()==nlab) break;
         }
         if (i==nk) continue;
         for (i=0; i<nk; i++) {
	     pref=(AliTrackReference*)itsRefs->UncheckedAt(i);
             if (pref->Label()==plab) break;
         }
         if (i==nk) continue;

         /*** fiducial volume ***/
         TParticle *p=stack->Particle(nlab);
         Double_t x=p->Vx(), y=p->Vy(), z=p->Vz(), r2=x*x+y*y;
         if (r2<r2min) continue;
         if (r2>r2max) continue;
       
         Int_t pdg=p0->GetPdgCode();
         nref->SetLength(pdg);  //This will the V0's PDG !

         new((*nrefs)[nv]) AliTrackReference(*nref);
         new((*prefs)[nv]) AliTrackReference(*pref);

         nv++;
      }
      itsRefs->Clear();

      v0Tree.Fill();
      nrefs->Clear(); prefs->Clear();

   }  //**** end of the loop over generated events

   v0Tree.Write();
   v0File->Close();

   delete itsTree;
   itsFile->Close();

   delete rl;
   return 0;
}
