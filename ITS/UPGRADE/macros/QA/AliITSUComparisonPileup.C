/****************************************************************************
 *  This macro calculates the efficiency of pileup reconstruction.          *
 *  Works only for events generated with the AliGenPileup generator.        *
 *                                                                          *
 *  Before running, load the ITSU libraries:                                *
 *  gSystem->Load("libITSUpgradeBase");gSystem->Load("libITSUpgradeRec");   *
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

  #include "AliStack.h"
  #include "AliHeader.h"
  #include "AliGenCocktailEventHeader.h"
  #include "AliRunLoader.h"
  #include "AliRun.h"
  #include "AliESDEvent.h"
  #include "AliESDtrack.h"
#endif

Int_t GoodPileupVertices(const Char_t *dir=".");

extern AliRun *gAlice;
extern TROOT *gROOT;

Int_t AliITSUComparisonPileup(const Char_t *dir=".") {
   ::Info("AliITSUComparisonPileup.C","Doing comparison...");

   // **** Book histogramms   
   Int_t nb=20;
   Float_t min=0, max=30.;
   TH2F *h2spd=(TH2F*)gROOT->FindObject("h2spd");
   if (!h2spd) 
      h2spd=new TH2F("h2spd",";Good vertices;Reconstructed vertices",
      nb,min,max, nb,min,max);
   h2spd->SetLineColor(2);
   TH2F *h2trk=(TH2F*)gROOT->FindObject("h2trk");
   if (!h2trk) 
      h2trk=new TH2F("h2trk",";Good vertices;Reconstructed vertices",
      nb,min,max, nb,min,max);
   h2trk->SetLineColor(4);

   nb=100;
   min=-0.03; max=0.03;
   TH1F *hzspd=(TH1F*)gROOT->FindObject("hzspd");
   if (!hzspd) 
     hzspd=new TH1F("hzspd","Resolution in Z;#DeltaZ (cm);",nb,min,max);
   hzspd->SetLineColor(2);
   TH1F *hztrk=(TH1F*)gROOT->FindObject("hztrk");
   if (!hztrk) 
     hztrk=new TH1F("hztrk","Resolution in Z;#DeltaZ (cm);",nb,min,max);
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
      heffspd=new TH1F("heffspd","Efficiency;Z (cm);Efficiency",nb,min,max);
   heffspd->SetLineColor(2);
   TH1F *hfoundtrk=(TH1F*)gROOT->FindObject("hfoundtrk");
   if (!hfoundtrk) 
    hfoundtrk=new TH1F("hfoundtrk",";Z (cm);",nb,min,max);
   TH1F *hefftrk=(TH1F*)gROOT->FindObject("hefftrk");
   if (!hefftrk) 
      hefftrk=new TH1F("hefftrk",";Z (cm);Efficiency",nb,min,max);
   hefftrk->SetLineColor(4);


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
   Int_t e=0;
   while (esdTree->GetEvent(e)) {
     cout<<endl<<endl<<"********* Processing event number: "<<e<<"*******\n";
     TClonesArray *verticesSPD=event->GetPileupVerticesSPD();
     Int_t nfoundSPD=verticesSPD->GetEntries(); 
     TClonesArray *verticesTRK=event->GetPileupVerticesTracks();
     Int_t nfoundTRK=verticesTRK->GetEntries(); 

     if (refTree->GetEvent(e++)==0) {
        cerr<<"No reconstructable vertices for this event !\n";
        continue;
     }
     Int_t ngood=refs->GetEntriesFast(); 
     cout<<"Found SPD: "<<nfoundSPD<<"  good: "<<ngood<<endl;

     h2spd->Fill(ngood,nfoundSPD);
     h2trk->Fill(ngood,nfoundTRK);

     for (Int_t g=0; g<ngood; g++) {
         const AliESDVertex *vtxg=(AliESDVertex*)refs->UncheckedAt(g);
         Double_t zg=vtxg->GetZv();
         hgood->Fill(zg);

         const AliESDVertex *vtxf=0;
         Double_t zf=0.;
         Int_t f=0;
         for (; f<nfoundSPD; f++) {
             vtxf=(AliESDVertex*)verticesSPD->UncheckedAt(f);
             if (!vtxf->GetStatus()) continue;
             zf=vtxf->GetZv();
             if (TMath::Abs(zf-zg)>2e-2) continue;
             break;
         }
         if (f>=nfoundSPD) {
	     vtxf=event->GetPrimaryVertexSPD();
             if (!vtxf->GetStatus()) goto trk;
             zf=vtxf->GetZv();
             if (TMath::Abs(zf-zg)>2e-2) goto trk;
	 }
         hfoundspd->Fill(zg);
         hzspd->Fill(zf-zg);

     trk:
         for (f=0; f<nfoundTRK; f++) {
             vtxf=(AliESDVertex*)verticesTRK->UncheckedAt(f);
             if (!vtxf->GetStatus()) continue;
             zf=vtxf->GetZv();
             if (TMath::Abs(zf-zg)>2e-2) continue;
             break;
         }
         if (f>=nfoundTRK) {
	     vtxf=event->GetPrimaryVertexTracks();
             if (!vtxf->GetStatus()) continue;
             zf=vtxf->GetZv();
             if (TMath::Abs(zf-zg)>2e-2) continue;
	 }
	 hfoundtrk->Fill(zg);
         hztrk->Fill(zf-zg);

     }

   } //***** End of the loop over reconstructed events

   delete event;
   delete esdTree;
   ef->Close();

   refFile->Close();
   
   TCanvas *c1=new TCanvas("c1","",0,0,700,1000);
   c1->Divide(1,3);
   c1->cd(1);
   gPad->SetGridx(); gPad->SetGridy();
   h2spd->Draw("box");
   h2trk->Draw("boxsame");
   TLine *l=new TLine(0,0,30,30);
   l->Draw("same");

   c1->cd(2);
   gPad->SetGridx(); gPad->SetGridy();
   heffspd->Divide(hfoundspd,hgood,1,1,"b");
   heffspd->SetMinimum(0.); heffspd->SetMaximum(1.);
   heffspd->Draw("hist");
   hefftrk->Divide(hfoundtrk,hgood,1,1,"b");
   hefftrk->Draw("histsame");

   c1->cd(3);
   gPad->SetGridx(); gPad->SetGridy();
   hzspd->Draw();
   hztrk->Draw("same");

   TFile fc("AliITSUComparisonPileup.root","RECREATE");
   c1->Write();
   fc.Close();

   return 0;
}



Int_t GoodPileupVertices(const Char_t *dir) {
  Int_t FindContributors(Float_t tz, AliStack *stack, Int_t ntrk);
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
         Int_t ntrk=FindContributors(t,stack,np);
         if (ntrk<3) continue;
         AliESDVertex *vertex=new ((*refs)[nv]) AliESDVertex();
         vertex->SetXv(vtx[0]);
         vertex->SetYv(vtx[1]);
         vertex->SetZv(vtx[2]);
         vertex->SetNContributors(ntrk);
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

Int_t FindContributors(Float_t tz, AliStack *stack, Int_t np) {
  Int_t ntrk=0;
  for (Int_t i=0; i<np; i++) {
      if (!stack->IsPhysicalPrimary(i)) continue;
      TParticle *part=stack->Particle(i);
      if (!part) continue;
      TParticlePDG *partPDG = part->GetPDG();
      if (!partPDG) continue;
      if (TMath::Abs(partPDG->Charge())<1e-10) continue;
      if (part->Pt()<1) continue;
      Float_t dt=0.5*(tz-part->T())/(tz+part->T());
      if (TMath::Abs(dt)>1e-5) continue;
      ntrk++;
  }
  return ntrk;
}
