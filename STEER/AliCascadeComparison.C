/****************************************************************************
 *           Very important, delicate and rather obscure macro.             *
 *                                                                          *
 *               Creates list of "findable" cascades,                       *
 *             calculates efficiency, resolutions etc.                      *
 *                                                                          *
 *  Origin: Christian Kuhn, IReS, Strasbourg, christian.kuhn@ires.in2p3.fr  *
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
  #include "AliESDEvent.h"
  #include "AliESDcascade.h"
#else
  const Int_t kXiMinus    = 3312;
  const Int_t kXiPlusBar  = -3312;
  const Int_t kOmegaMinus = 3334;
  const Int_t kOmegaPlusBar = -3334;
#endif

Int_t GoodCascades(const Char_t *dir=".");

extern AliRun *gAlice;
extern TBenchmark *gBenchmark;
extern TROOT *gROOT;

static Int_t allgood=0;
static Int_t allfound=0;

Int_t AliCascadeComparison(Int_t code=3312, const Char_t *dir=".") {
  //code= 3312; //kXiMinus 
  //code=-3312; //kXiPlusBar
  //code= 3334; //kOmegaMinus
  //code=-3334; //kOmegaPlusBar
   gBenchmark->Start("AliCascadeComparison");

   ::Info("AliCascadeComparison.C","Doing comparison...");

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


   TH1F *hp=(TH1F*)gROOT->FindObject("hp");
   if (!hp) hp=new TH1F("hp","Angular Resolution",30,-30.,30.);
   hp->SetXTitle("(mrad)"); hp->SetFillColor(2);

   TH1F *hl=(TH1F*)gROOT->FindObject("hl");
   if (!hl) hl=new TH1F("hl","Lambda Resolution",30,-30,30);
   hl->SetXTitle("(mrad)"); hl->SetFillColor(1); hl->SetFillStyle(3013);
 
   TH1F *hpt=(TH1F*)gROOT->FindObject("hpt");
   if (!hpt) hpt=new TH1F("hpt","Relative Pt Resolution",30,-10.,10.); 
   hpt->SetXTitle("(%)"); hpt->SetFillColor(2); 

   TH1F *hx=(TH1F*)gROOT->FindObject("hx");
   if (!hx) hx=new TH1F("hx","Position Resolution (X,Y)",30,-3.,3.);
   hx->SetXTitle("(mm)"); hx->SetFillColor(6);

   TH1F *hy=(TH1F*)gROOT->FindObject("hy");
   if (!hy) hy=new TH1F("hy","Position Resolution (Y)",30,-3.,3.);  
   hy->SetXTitle("(mm)"); hy->SetFillColor(1); hy->SetFillStyle(3013);

   TH1F *hz=(TH1F*)gROOT->FindObject("hz");
   if (!hz) hz=new TH1F("hz","Position Resolution (Z)",30,-3.,3.);
   hz->SetXTitle("(mm)"); hz->SetFillColor(6);


   Double_t pmin=0.2, pmax=4.2; Int_t nchan=20;
   TH1F *hgood=(TH1F*)gROOT->FindObject("hgood");    
   if (!hgood) hgood=new TH1F("hgood","Good Cascades",nchan,pmin,pmax);
    
   TH1F *hfound=(TH1F*)gROOT->FindObject("hfound");
   if (!hfound) hfound=new TH1F("hfound","Found Cascades",nchan,pmin,pmax);

   TH1F *hfake=(TH1F*)gROOT->FindObject("hfake");
   if (!hfake) hfake=new TH1F("hfake","Fake Cascades",nchan,pmin,pmax);

   TH1F *hg=(TH1F*)gROOT->FindObject("hg");
   if (!hg) hg=new TH1F("hg","Efficiency for Good Cascades",nchan,pmin,pmax);
   hg->SetLineColor(4); hg->SetLineWidth(2);

   TH1F *hf=(TH1F*)gROOT->FindObject("hf");
   if (!hf) hf=new TH1F("hf","Probability of Fake Cascades",nchan,pmin,pmax);
   hf->SetFillColor(1); hf->SetFillStyle(3013); hf->SetLineWidth(2);


   Double_t mmin=cascadeMass-cascadeWindow, mmax=cascadeMass+cascadeWindow;
   TH1F *cs=(TH1F*)gROOT->FindObject("cs");
   if (!cs) cs=new TH1F("cs","Cascade Effective Mass",40, mmin, mmax);
   cs->SetXTitle("(GeV)"); cs->SetLineColor(4); cs->SetLineWidth(4);

   TH1F *csf=(TH1F*)gROOT->FindObject("csf");
   if (!csf) csf=new TH1F("csf","Fake Cascade Effective Mass",40, mmin, mmax);
   csf->SetXTitle("(GeV)"); csf->SetFillColor(6);


   Char_t fname[100];
   sprintf(fname,"%s/GoodCascades.root",dir);

   TFile *refFile=TFile::Open(fname,"old");
   if (!refFile || !refFile->IsOpen()) {
   ::Info("AliCascadeComparison.C","Marking good cascades (will take a while)...");
     if (GoodCascades(dir)) {
        ::Error("AliCascadesComparison.C","Can't generate the reference file !");
        return 1;
     }
   }
   refFile=TFile::Open(fname,"old");
   if (!refFile || !refFile->IsOpen()) {
     ::Error("AliCascadeComparison.C","Can't open the reference file !");
     return 1;
   }   
  
   TTree *csTree=(TTree*)refFile->Get("csTree");
   if (!csTree) {
     ::Error("AliCascadeComparison.C","Can't get the reference tree !");
     return 2;
   }
   TBranch *pbranch=csTree->GetBranch("positive");
   if (!pbranch) {
     ::Error("AliCascadeComparison.C","Can't get the positive daughter branch !");
     return 3;
   }
   TClonesArray dummy("AliTrackReference",1000), *prefs=&dummy;
   pbranch->SetAddress(&prefs);

   TBranch *nbranch=csTree->GetBranch("negative");
   if (!nbranch) {
     ::Error("AliCascadeComparison.C","Can't get the negative daughter branch !");
     return 4;
   }
   TClonesArray dumm("AliTrackReference",1000), *nrefs=&dumm;
   nbranch->SetAddress(&nrefs);

   TBranch *bbranch=csTree->GetBranch("bachelor");
   if (!nbranch) {
     ::Error("AliCascadeComparison.C","Can't get the bachelor branch !");
     return 4;
   }
   TClonesArray dum("AliTrackReference",1000), *brefs=&dum;
   bbranch->SetAddress(&brefs);


   
   sprintf(fname,"%s/AliESDs.root",dir);
   TFile *ef=TFile::Open(fname);
   if ((!ef)||(!ef->IsOpen())) {
      sprintf(fname,"%s/AliESDcascade.root",dir);
      ef=TFile::Open(fname);
      if ((!ef)||(!ef->IsOpen())) {
         ::Error("AliCascadeComparison.C","Can't open AliESDcascade.root !");
         return 5;
      }
   }
   AliESDEvent* event = new AliESDEvent();
   TTree* esdTree = (TTree*) ef->Get("esdTree");
   if (!esdTree) {
      ::Error("AliCascadeComparison.C", "no ESD tree found");
      return 6;
   }
   event->ReadFromTree(esdTree);


   //******* Loop over events *********
   Int_t e=0;
   while (esdTree->GetEvent(e)) {
      cout<<endl<<endl<<"********* Processing event number: "<<e<<"*******\n";
 
      Int_t nentr=event->GetNumberOfCascades();
      allfound+=nentr;


      if (csTree->GetEvent(e++)==0) {
         cerr<<"No reconstructable cascades !\n";
         continue;
      }

      Int_t ngood=prefs->GetEntriesFast(),ng=0; 

      Double_t pxg=0.,pyg=0.,pzg=0.,ptg=0.;
      Int_t nlab=-1, plab=-1, blab=-1;
      Int_t i;
      for (i=0; i<nentr; i++) {
          AliESDcascade *cascade=event->GetCascade(i);

          Int_t nidx=TMath::Abs(cascade->GetNindex());
          Int_t pidx=TMath::Abs(cascade->GetPindex());
          Int_t bidx=TMath::Abs(cascade->GetBindex());

          AliESDtrack *ntrack=event->GetTrack(nidx);
          AliESDtrack *ptrack=event->GetTrack(pidx);
          AliESDtrack *btrack=event->GetTrack(bidx);

          nlab=TMath::Abs(ntrack->GetLabel()); 
          plab=TMath::Abs(ptrack->GetLabel());
          blab=TMath::Abs(btrack->GetLabel());

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

          AliTrackReference *nref=0, *pref=0, *bref=0;
          Int_t j;
          for (j=0; j<ngood; j++) {
              bref=(AliTrackReference*)brefs->UncheckedAt(j);
              nref=(AliTrackReference*)nrefs->UncheckedAt(j);
              pref=(AliTrackReference*)prefs->UncheckedAt(j);
              if (bref->Label() == blab)
              if (nref->Label() == nlab)
              if (pref->Label() == plab) break;
          }

          if (TMath::Abs(mass-cascadeMass)>cascadeWidth) continue;

          Double_t px,py,pz; cascade->GetPxPyPz(px,py,pz);
          Double_t pt=TMath::Sqrt(px*px+py*py);

          if (j==ngood) {
             hfake->Fill(pt);
             cout<<"Fake cascade: ("<<nlab<<","<<plab<<","<<blab<<")\n";
             continue;
          }
          csf->Fill(mass,-1);

          pxg=bref->Px()+nref->Px()+pref->Px(); 
          pyg=bref->Px()+nref->Py()+pref->Py(); 
          pzg=nref->Pz()+pref->Pz(); 
          ptg=TMath::Sqrt(pxg*pxg+pyg*pyg);
          Double_t phig=TMath::ATan2(pyg,pxg), phi=TMath::ATan2(py,px);
          Double_t lamg=TMath::ATan2(pzg,ptg), lam=TMath::ATan2(pz,pt);
          hp->Fill((phi - phig)*1000.);
          hl->Fill((lam - lamg)*1000.);
          hpt->Fill((1/pt - 1/ptg)/(1/ptg)*100.);

          Double_t x,y,z; cascade->GetXYZ(x,y,z);
          hx->Fill((x-nref->X())*10);
          hy->Fill((y-nref->Y())*10);
          hz->Fill((z-nref->Z())*10);

          hfound->Fill(ptg);
          nref->SetLabel(-1);

      }
      for (i=0; i<ngood; i++) {
         AliTrackReference *bref=(AliTrackReference*)brefs->UncheckedAt(i);
         AliTrackReference *nref=(AliTrackReference*)nrefs->UncheckedAt(i);
         AliTrackReference *pref=(AliTrackReference*)prefs->UncheckedAt(i);
         Int_t pdg=(Int_t)nref->GetLength();  //this is the cascade's PDG !
         if (code!=pdg) continue;
         ng++;
         pxg=bref->Px()+nref->Px()+pref->Px(); 
         pyg=bref->Px()+nref->Py()+pref->Py(); 
         ptg=TMath::Sqrt(pxg*pxg+pyg*pyg);
         hgood->Fill(ptg);
         nlab=nref->Label(); plab=pref->Label(); blab=bref->Label();
         if (nlab < 0) continue;
         cout<<"Cascade ("<<nlab<<','<<plab<<","<<blab<<") has not been found !\n";
      }
      allgood+=ng;

      cout<<"Number of found cascades : "<<nentr<<endl;
      cout<<"Number of \"good\" cascades : "<<ng<<endl;

      brefs->Clear();
      prefs->Clear();
      nrefs->Clear();

   } //**** End of the loop over events

   delete event;
   delete esdTree;
   ef->Close();

   delete csTree;
   refFile->Close();

   Stat_t ngg=hgood->GetEntries(), nf=hfound->GetEntries();
   if (ngg!=0) cout<<"Integral efficiency is about "<<nf/ngg*100.<<" %\n";
   cout<<
   "Total number of found cascades: "<<allfound<<" ("<<nf<<" in the peak)\n";
   cout<<"Total number of \"good\" cascades: "<<allgood<<endl;


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
   if (cs->GetEntries()<minc) cs->Draw();
   else cs->Fit("gaus","","",cascadeMass-cascadeWidth,cascadeMass+cascadeWidth);
   csf->Draw("same");
   Double_t max=cs->GetMaximum();
   TLine *line3 = 
      new TLine(cascadeMass-cascadeWidth,0.,cascadeMass-cascadeWidth,max);
   line3->Draw("same");
   TLine *line4 = 
      new TLine(cascadeMass+cascadeWidth,0.,cascadeMass+cascadeWidth,max);
   line4->Draw("same");

   c2->Update();

   TFile fc("AliCascadeComparison.root","RECREATE");
   c1->Write();
   c2->Write();
   fc.Close();

   gBenchmark->Stop("AliCascadeComparison");
   gBenchmark->Show("AliCascadeComparison");


   return 0;
}


Int_t GoodCascades(const Char_t *dir) {
   if (gAlice) {
      delete AliRunLoader::GetRunLoader();
      delete gAlice; 
      gAlice=0;
   }   

   Char_t fname[100];
   sprintf(fname,"%s/galice.root",dir);

   AliRunLoader *rl = AliRunLoader::Open(fname,"COMPARISON");
   if (!rl) {
       ::Error("GoodCascades","Can't start session !");
       return 1;
   }

   rl->LoadgAlice();
   rl->LoadHeader();
   rl->LoadKinematics();


   Int_t nev=rl->GetNumberOfEvents();
   ::Info("GoodCascades","Number of events : %d\n",nev);  

 
   sprintf(fname,"%s/GoodTracksITS.root",dir);
   TFile *itsFile=TFile::Open(fname);
   if ((!itsFile)||(!itsFile->IsOpen())) {
       ::Error("GoodCAscades","Can't open the GoodTracksITS.root !");
       delete rl;
       return 5; 
   }
   TClonesArray dm("AliTrackReference",1000), *itsRefs=&dm;
   TTree *itsTree=(TTree*)itsFile->Get("itsTree");
   if (!itsTree) {
       ::Error("GoodCascades","Can't get the ITS reference tree !");
       delete rl;
       return 6;
   }
   TBranch *itsBranch=itsTree->GetBranch("ITS");
   if (!itsBranch) {
      ::Error("GoodCascades","Can't get the ITS reference branch !");
      delete rl;
      return 7;
   }
   itsBranch->SetAddress(&itsRefs);


   sprintf(fname,"%s/GoodCascades.root",dir);
   TFile *csFile=TFile::Open(fname,"recreate");
   TClonesArray dummy("AliTrackReference",1000), *nrefs=&dummy;
   TClonesArray dumm("AliTrackReference",1000), *prefs=&dumm;
   TClonesArray dum("AliTrackReference",1000), *brefs=&dum;
   TTree csTree("csTree","Tree with info about the reconstructable cascades");
   csTree.Branch("negative",&nrefs);
   csTree.Branch("positive",&prefs);
   csTree.Branch("bachelor",&brefs);


   // *** Get information about the cuts ***
   Double_t r2min=0.9*0.9;
   Double_t r2max=2.9*2.9;


   //********  Loop over generated events
   for (Int_t e=0; e<nev; e++) {
      rl->GetEvent(e);  csFile->cd();

      Int_t np = rl->GetHeader()->GetNtrack();
      cout<<"Event "<<e<<" Number of particles: "<<np<<endl;

      itsTree->GetEvent(e);
      Int_t nk=itsRefs->GetEntriesFast();

      AliStack *stack=rl->Stack();

      AliTrackReference *nref=0, *pref=0, *bref=0;

      Int_t nc=0;
      while (np--) {
	 //cerr<<np<<'\r';
         TParticle *cp=stack->Particle(np);

         // *** only these cascades are "good" ***
         Int_t code=cp->GetPdgCode();
         if (code!=kXiMinus)    if (code!=kXiPlusBar)
         if (code!=kOmegaMinus) if (code!=kOmegaPlusBar) continue; 

         // *** daughter tracks must be "good" ***
         Int_t v0lab=cp->GetFirstDaughter(), blab=cp->GetLastDaughter();
         if (v0lab==blab) continue;
         if (v0lab<0) continue;
         if (blab<0) continue;

         TParticle *p0=stack->Particle(v0lab);
         TParticle *bp=stack->Particle(blab);
         Int_t i;
         if ((p0->GetPdgCode()!=kLambda0) && (p0->GetPdgCode()!=kLambda0Bar)) {
            TParticle *p=p0; p0=bp; bp=p;
            i=v0lab; v0lab=blab; blab=i;         
            if ((p0->GetPdgCode()!=kLambda0)&&(p0->GetPdgCode()!=kLambda0Bar))
                                                                   continue;
         }

         // ** is the bachelor "good" ? **
         for (i=0; i<nk; i++) {
	     bref=(AliTrackReference*)itsRefs->UncheckedAt(i);
             if (bref->Label()==blab) break;
         }
         if (i==nk) continue; 

         // ** is the V0 "good" ? **
         Int_t plab=p0->GetFirstDaughter(), nlab=p0->GetLastDaughter();
         if (nlab==plab) continue;
         if (nlab<0) continue;
         if (plab<0) continue;
         if (stack->Particle(plab)->GetPDG()->Charge() < 0.) {
            i=plab; plab=nlab; nlab=i;
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


         // *** fiducial volume ***
         Double_t x=bp->Vx(), y=bp->Vy(), r2=x*x+y*y; //bachelor
         if (r2<r2min) continue;
         if (r2>r2max) continue;
         TParticle *pp=stack->Particle(plab);
         x=pp->Vx(); y=pp->Vy(); r2=x*x+y*y;          //V0
         if (r2<r2min) continue;
         if (r2>r2max) continue;

         Int_t pdg=cp->GetPdgCode();
         nref->SetLength(pdg);  //This will the cascade's PDG !

         new((*nrefs)[nc]) AliTrackReference(*nref);
         new((*prefs)[nc]) AliTrackReference(*pref);
         new((*brefs)[nc]) AliTrackReference(*bref);

         nc++;
      }
      itsRefs->Clear();

      csTree.Fill();
      nrefs->Clear(); prefs->Clear(); brefs->Clear();

   } //**** end of the loop over generated events

   csTree.Write();
   csFile->Close();

   delete itsTree;
   itsFile->Close();

   delete rl;
   return 0;
}
