/****************************************************************************
 *  This macro calculates the efficiency of pileup reconstruction.          *
 *  Works only for events generated with the AliGenPileup generator.        *
 *                                                                          *
 *  Before running, load the ITSU libraries:                                *
 *  gSystem->Load("libITSUpgradeBase");gSystem->Load("libITSUpgradeRec");   *
 *                                                                          *
 *  Definifions:                                                            *
 *  1) Reconstructable track: physical primary, charged, pT > pTmin         * 
 *  2) Reconstructable vertex: has at least nMin reconstructable tracks     *
 *  3) Associated vertex: has at least nAssMin correctly associated tracks  *
 *  4) Good associated vertex: the fraction of correctly associated tracks  *
 *        is at least fracMin                                               *
 *  5) Fake associated vertex: not a good associated vertex                 *
 *  6) Efficiency:  the ratio of 4) over 3)                                 *
 *  7) Fake rate:   the ratio of 5) over 3)                                 *
 *                                                                          *
 *           Origin: I.Belikov, IPHC, Iouri.Belikov@iphc.cnrs.fr            *
 ****************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <Riostream.h>
  #include <TMath.h>
  #include <TTree.h>
  #include <TParticle.h>
  #include <TParticlePDG.h>
  #include <TCanvas.h>
  #include <TFile.h>
  #include <TLine.h>
  #include <TROOT.h>
  #include <TStyle.h>
  #include <TLegend.h>

  #include "AliStack.h"
  #include "AliHeader.h"
  #include "AliGenCocktailEventHeader.h"
  #include "AliRunLoader.h"
  #include "AliRun.h"
  #include "AliESDEvent.h"
  #include "AliESDtrack.h"
#endif

//**** Parameters used in the definitions
const Float_t pTmin=0.2;   // Minimal pT for a reconstructable track
const Int_t nMin=3;        // Minimal N of reconstructable tracks per vertex
const Int_t nAssMin=2;     // Minimal number of correctly associated tracks
const Float_t fracMin=0.8; // Minimal fraction of correctly associated tracks
const Float_t tWin=30e-6;  // Time-acceptance window for "good" MC vertices

extern AliRun *gAlice;
extern TROOT *gROOT;
extern TStyle *gStyle;

Int_t AliITSUComparisonPileup(const Char_t *dir=".") {
   ::Info("AliITSUComparisonPileup.C","Doing comparison...");
   Int_t GoodPileupVertices(const Char_t *dir=".");

   // **** Book histogramms   
   Int_t nb=35;
   Float_t min=0, max=70.;
   TH2F *h2spd=(TH2F*)gROOT->FindObject("h2spd");
   if (!h2spd) 
     h2spd=new TH2F("h2spd","SPD vertices;Number of good vertices;Number of reconstructed vertices",nb,min,max, nb,min,max);
   h2spd->SetLineColor(2);
   TH2F *h2trk=(TH2F*)gROOT->FindObject("h2trk");
   if (!h2trk) 
     h2trk=new TH2F("h2trk","TRK vertices;Good vertices;Reconstructed vertices",
     nb,min,max, nb,min,max);
   h2trk->SetLineColor(4);

   nb=100;
   min=-0.03; max=0.03;
   TH1F *hzspd=(TH1F*)gROOT->FindObject("hzspd");
   if (!hzspd) 
     hzspd=new TH1F("hzspd","SPD resolution in Z;#DeltaZ (cm);",nb,min,max);
   hzspd->SetLineColor(2);
   TH1F *hztrk=(TH1F*)gROOT->FindObject("hztrk");
   if (!hztrk) 
     hztrk=new TH1F("hztrk","TRK resolution in Z;#DeltaZ (cm);",nb,min,max);
   hztrk->SetLineColor(4);

   
   nb=30;
   min=-10.; max=10.; 
   TH1F *hgood=(TH1F*)gROOT->FindObject("hgood");
   if (!hgood) 
     hgood=new TH1F("hgood",";Z (cm);",nb,min,max);

   TH1F *hfoundspd=(TH1F*)gROOT->FindObject("hfoundspd");
   if (!hfoundspd) 
    hfoundspd=new TH1F("hfoundspd",";Z (cm);",nb,min,max);
   TH1F *heffspd=(TH1F*)gROOT->FindObject("heffspd");
   if (!heffspd) 
      heffspd=new TH1F("heffspd","SPD efficiency + fake rate;Z position of a prim. vertex (cm);Efficiency",nb,min,max);
   heffspd->SetLineColor(2);
   heffspd->Sumw2();

   TH1F *hfoundtrk=(TH1F*)gROOT->FindObject("hfoundtrk");
   if (!hfoundtrk) 
    hfoundtrk=new TH1F("hfoundtrk",";Z (cm);",nb,min,max);
   TH1F *hefftrk=(TH1F*)gROOT->FindObject("hefftrk");
   if (!hefftrk) 
      hefftrk=new TH1F("hefftrk","TRK efficiency;Z (cm);Efficiency",nb,min,max);
   hefftrk->SetLineColor(4);
   hefftrk->Sumw2();

   TH1F *hfaketrk=(TH1F*)gROOT->FindObject("hfaketrk");
   if (!hfaketrk) 
    hfaketrk=new TH1F("hfaketrk",";Z (cm);",nb,min,max);
   TH1F *heffaketrk=(TH1F*)gROOT->FindObject("heffaketrk");
   if (!heffaketrk) 
      heffaketrk=new TH1F("heffaketrk","TRK fake rate;Z (cm);Fake rate",nb,min,max);
   heffaketrk->SetLineColor(4);
   heffaketrk->SetFillColor(590);
   heffaketrk->Sumw2();



   nb=51;
   min=-0.5; max=50.5;
   TH1F *hngood=(TH1F*)gROOT->FindObject("hngood");
   if (!hngood) 
      hngood=new TH1F("hngood",";Z (cm);",nb,min,max);
   TH1F *hnfoundtrk=(TH1F*)gROOT->FindObject("hnfoundtrk");
   if (!hnfoundtrk) 
      hnfoundtrk=new TH1F("hnfoundtrk",";Z (cm);",nb,min,max);
   TH1F *hnefftrk=(TH1F*)gROOT->FindObject("hnefftrk");
   if (!hnefftrk) 
      hnefftrk=new TH1F("hnefftrk","TRK efficiency;Number of tracks;Efficiency",nb,min,max);
   hnefftrk->SetLineColor(4);
   hnefftrk->Sumw2();

   TH1F *hnfaketrk=(TH1F*)gROOT->FindObject("hnfaketrk");
   if (!hnfaketrk) 
      hnfaketrk=new TH1F("hnfaketrk",";Z (cm);",nb,min,max);
   TH1F *hneffaketrk=(TH1F*)gROOT->FindObject("hneffaketrk");
   if (!hneffaketrk) 
      hneffaketrk=new TH1F("hneffaketrk","TRK fake rate;Number of tracks;Efficiency",nb,min,max);
   hneffaketrk->SetLineColor(4);
   hneffaketrk->SetFillColor(590);
   hneffaketrk->Sumw2();


   // **** Generate a rerefence file with reconstructable vertices
   Char_t fname[100];
   sprintf(fname,"%s/GoodPileupVertices.root",dir);
   TFile *refFile=TFile::Open(fname,"old");
   if (!refFile || !refFile->IsOpen()) {
      ::Info("AliITSUComparisonPileup.C",
      "Marking good pileup vertices (will take a while)...");
      if (GoodPileupVertices(dir)) {
     ::Error("AliITSUComparisonPileup.C","Can't generate the reference file !");
         return 1;
      }
   }
   refFile=TFile::Open(fname,"old");
   if (!refFile || !refFile->IsOpen()) {
     ::Error("AliITSUComparisonPileup.C","Can't open the reference file !");
     return 1;
   }   
   TTree *refTree=(TTree*)refFile->Get("refTree");
   if (!refTree) {
     ::Error("AliITSUComparisonPileup.C","Can't get the reference tree !");
     return 2;
   }
   TBranch *branch=refTree->GetBranch("Vertices");
   if (!branch) {
     ::Error("AliITSUComparisonPileup.C","Can't get the vertex branch !");
     return 3;
   }
   TClonesArray dummy("AliESDVertex",100), *refs=&dummy;
   branch->SetAddress(&refs);    


   // **** Open the ESD 
   sprintf(fname,"%s/AliESDs.root",dir);
   TFile *ef=TFile::Open(fname);
   if ((!ef)||(!ef->IsOpen())) {
      ::Error("AliITSUComparisonPileup.C","Can't open AliESDs.root !");
      return 4;
   }
   AliESDEvent* event = new AliESDEvent();
   TTree* esdTree = (TTree*) ef->Get("esdTree");
   if (!esdTree) {
      ::Error("AliITSUComparisonPileup.C", "no ESD tree found");
      return 6;
   }
   event->ReadFromTree(esdTree);


   //******* Loop over reconstructed events *********
   Int_t ntrk=0, ntrkcor=0, ntrkwro=0;
   Int_t e=0;
   while (esdTree->GetEvent(e)) {
     cout<<endl<<endl<<"********* Processing event number: "<<e<<"*******\n";
     Int_t nn=event->GetNumberOfTracks();
     ntrk += nn;

     TClonesArray *verticesSPD=event->GetPileupVerticesSPD();
     Int_t nfoundSPD=verticesSPD->GetEntries(); 
     TClonesArray *verticesTRK=event->GetPileupVerticesTracks();
     Int_t nfoundTRK=verticesTRK->GetEntries(); 

     if (refTree->GetEvent(e++)==0) {
        cerr<<"No reconstructable vertices for this event !\n";
        continue;
     }
     Int_t ngood=refs->GetEntriesFast(); 
     cout<<"Found SPD vertices: "<<nfoundSPD<<
           "  Reconstructable vertices: "<<ngood<<endl;

     h2spd->Fill(ngood,nfoundSPD);
     h2trk->Fill(ngood,nfoundTRK);

     Int_t ncor=0, nwro=0;
     for (Int_t g=0; g<ngood; g++) {
         Int_t Associate(const AliESDVertex *g, const AliESDVertex *f, 
            const AliESDEvent *esd); 
         const AliESDVertex *vtxg=(AliESDVertex*)refs->UncheckedAt(g);
         Double_t zg=vtxg->GetZ();
         Double_t ng=vtxg->GetNIndices();
         hgood->Fill(zg);
         hngood->Fill(ng);

         AliESDVertex *vtxf=0;
         Double_t zf=0.;
         Int_t f=0;
         for (; f<nfoundSPD; f++) {
             vtxf=(AliESDVertex*)verticesSPD->UncheckedAt(f);
             if (!vtxf->GetStatus()) continue;
             if (Associate(vtxg,vtxf,event)==0) continue; 
             break;
         }
         if (f>=nfoundSPD) {
	     vtxf=(AliESDVertex *)event->GetPrimaryVertexSPD();
             if (!vtxf->GetStatus()) goto trk;
             if (Associate(vtxg,vtxf,event)==0) goto trk; 
	 }

         zf=vtxf->GetZ();
         hfoundspd->Fill(zg);
         hzspd->Fill(zf-zg);

     trk:
         Int_t n=0;
         for (f=0; f<nfoundTRK; f++) {
             vtxf=(AliESDVertex*)verticesTRK->UncheckedAt(f);
             if (!vtxf->GetStatus()) continue;
             n=Associate(vtxg,vtxf,event);
	     if (n < nAssMin) continue;
             break;
         }
         if (f>=nfoundTRK) {
	     vtxf=(AliESDVertex*)event->GetPrimaryVertexTracks();
             if (!vtxf->GetStatus()) continue;
             n=Associate(vtxg,vtxf,event);
             if (n < nAssMin) continue;
	 }

         ncor+=n;
         nwro+=(vtxf->GetNIndices()-n); 
         zf=vtxf->GetZ();

         if (Float_t(n)/vtxf->GetNIndices() > fracMin) {
	    hfoundtrk->Fill(zg);
	    hnfoundtrk->Fill(ng);
	 } else {
	    hfaketrk->Fill(zg);
	    hnfaketrk->Fill(ng);
	 }
         hztrk->Fill(zf-zg);

         vtxf->SetNContributors(0); // Mark this vertex as already associated

     }
     // Increase the counter of tracks (not)associated with verices
     ntrkcor += ncor;
     ntrkwro += nwro;

   } //***** End of the loop over reconstructed events

   delete event;
   delete esdTree;
   ef->Close();

   refFile->Close();
   
   cout<<"\nTotal number of found tracks: "<<ntrk<<endl;
   cout<<"Number of tracks correctly associated with vertices: "<<ntrkcor<<endl;
   cout<<"Number of tracks wrongly associated with vertices: "<<ntrkwro<<endl;
   if (ntrk != 0) {
     cout<<"Correctly associated/Found:\t"<<Float_t(ntrkcor)/ntrk<<endl;
     cout<<"Wrongly associated/Found:\t"<<Float_t(ntrkwro)/ntrk<<endl;
     cout<<"Not associated/Found:\t\t"<<1.-Float_t(ntrkwro+ntrkcor)/ntrk<<endl;
   }

   gStyle->SetOptStat(0);
   gStyle->SetOptTitle(0);
   TCanvas *c1=new TCanvas("c1","",0,0,700,1000);
   c1->Divide(1,2);
   c1->cd(1);
   gPad->SetGridx(); gPad->SetGridy();
   h2spd->Draw("box");
   h2trk->Draw("boxsame");
   gPad->BuildLegend(0.13,0.65,0.46,0.86)->SetFillColor(0);
   TLine *l=new TLine(0,0,70,70);
   l->Draw("same");

   c1->cd(2);
   gPad->SetGridx(); gPad->SetGridy();
   hztrk->Draw();
   hzspd->Draw("same");
   gPad->BuildLegend(0.13,0.65,0.46,0.86)->SetFillColor(0);


   TCanvas *c2=new TCanvas("c2","",0,0,700,1000);
   c2->Divide(1,2);
   c2->cd(1);
   gPad->SetGridx(); gPad->SetGridy();
   heffspd->Divide(hfoundspd,hgood,1,1,"b");
   heffspd->SetMinimum(0.); heffspd->SetMaximum(1.2);
   heffspd->Draw("ehist");
   hefftrk->Divide(hfoundtrk,hgood,1,1,"b");
   hefftrk->Draw("ehistsame");
   heffaketrk->Divide(hfaketrk,hgood,1,1,"b");
   heffaketrk->Draw("ehistsame");
   gPad->BuildLegend(0.13,0.65,0.46,0.86)->SetFillColor(0);

   c2->cd(2);
   hnefftrk->Divide(hnfoundtrk,hngood,1,1,"b");
   hnefftrk->SetMinimum(0.); hnefftrk->SetMaximum(1.2);
   hneffaketrk->Divide(hnfaketrk,hngood,1,1,"b");
   gPad->SetGridx(); gPad->SetGridy();
   hnefftrk->Draw("ehist");
   hneffaketrk->Draw("ehistsame");
   gPad->BuildLegend(0.13,0.65,0.46,0.86)->SetFillColor(0);


   TFile fc("AliITSUComparisonPileup.root","RECREATE");
   c1->Write();
   c2->Write();
   fc.Close();

   return 0;
}

Int_t 
Associate(const AliESDVertex *g,const AliESDVertex *f,const AliESDEvent *esd) { 
   UShort_t *idxg=g->GetIndices(); Int_t ng=g->GetNIndices(); 
   UShort_t *idxf=f->GetIndices(); Int_t nf=f->GetNIndices();

   if (nf==0) { 
   // SPD vertex
       Double_t zg=g->GetZ();
       Double_t zf=f->GetZ();
       if (TMath::Abs(zf-zg)>2e-2) return 0;
       return 1;
   }
   // TRK vertex
   Int_t nass=0;
   for (Int_t i=0; i<ng; i++) {
       UShort_t labg=idxg[i];
       for (Int_t j=0; j<nf; j++) {
           const AliESDtrack *t=esd->GetTrack(idxf[j]);
           UShort_t labf=TMath::Abs(t->GetLabel());
           if (labg != labf) continue;
           nass++;
           break; 
       }
   } 

   return nass;
}

Int_t GoodPileupVertices(const Char_t *dir) {
  Bool_t FindContributors(Float_t tz, AliStack *stack, UShort_t *idx, Int_t &n);
   if (gAlice) { 
       delete AliRunLoader::Instance();
       delete gAlice;//if everything was OK here it is already NULL
       gAlice = 0x0;
   }

   Char_t fname[100];
   sprintf(fname,"%s/galice.root",dir);

   AliRunLoader *rl = AliRunLoader::Open(fname,"COMPARISON");
   if (!rl) {
      ::Error("GoodPileupVertices","Can't start session !");
      return 1;
   }

   rl->LoadgAlice();
   rl->LoadHeader();
   rl->LoadKinematics();


   Int_t nev=rl->GetNumberOfEvents();
   ::Info("GoodPileupVertices","Number of events : %d\n",nev);  

   sprintf(fname,"%s/GoodPileupVertices.root",dir);
   TFile *refFile=TFile::Open(fname,"recreate");
   TClonesArray dummy("AliESDVertex",100), *refs=&dummy;
   TTree refTree("refTree","Tree with the reconstructable vertices");
   refTree.Branch("Vertices",&refs);

   //********  Loop over generated events 
   for (Int_t e=0; e<nev; e++) {
     rl->GetEvent(e);  refFile->cd();
     AliStack* stack = rl->Stack();

     AliHeader *ah=rl->GetHeader();
     AliGenCocktailEventHeader *cock=
            (AliGenCocktailEventHeader*)ah->GenEventHeader();
     TList *headers=cock->GetHeaders();
     const Int_t nvtx=headers->GetEntries();
     const Int_t np=ah->GetNtrack();
     cout<<"Event "<<e<<" Number of vertices, particles: "
	 <<nvtx<<' '<<np<<endl;

     Int_t nv=0;
     for (Int_t v=0; v<nvtx; v++) {
         AliGenEventHeader *h=(AliGenEventHeader *)headers->At(v);
         TArrayF vtx(3); h->PrimaryVertex(vtx);
         Float_t t=h->InteractionTime();
         if (TMath::Abs(t)>tWin) continue;
         UShort_t *idx=new UShort_t[np];
         Int_t ntrk=0;
         if (!FindContributors(t,stack,idx,ntrk)) {delete[] idx; continue;}
         AliESDVertex *vertex=new ((*refs)[nv]) AliESDVertex();
         vertex->SetXv(vtx[0]);
         vertex->SetYv(vtx[1]);
         vertex->SetZv(vtx[2]);
         vertex->SetNContributors(ntrk);
         vertex->SetIndices(ntrk,idx);
         delete[] idx;
         nv++;
     }
     refTree.Fill();
     refs->Clear();
   } //*** end of the loop over generated events

   refTree.Write();
   refFile->Close();

   delete rl;
   return 0;
}

Bool_t FindContributors(Float_t tz, AliStack *stack, UShort_t *idx, Int_t &n) {
  Int_t ntrk=0;
  Int_t np=stack->GetNtrack();
  for (Int_t i=0; i<np; i++) {
      if (!stack->IsPhysicalPrimary(i)) continue;
      TParticle *part=stack->Particle(i);
      if (!part) continue;
      TParticlePDG *partPDG = part->GetPDG();
      if (!partPDG) continue;
      if (TMath::Abs(partPDG->Charge())<1e-10) continue;
      Float_t dt=0.5*(tz-part->T())/(tz+part->T());
      if (TMath::Abs(dt)>1e-5) continue;
      idx[n++]=i;
      if (part->Pt() > pTmin) ntrk++;
  }
  return (ntrk<nMin) ? kFALSE : kTRUE;
}
