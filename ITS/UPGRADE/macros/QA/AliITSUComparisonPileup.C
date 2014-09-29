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
  #include <TCanvas.h>
  #include <TFile.h>
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
   TH2F *h2=(TH2F*)gROOT->FindObject("h2");
   if (!h2) h2=new TH2F("h2",";Number of good vertices;Number of reconstructed vertices",100,0.,100.,10,0.,10.);


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


   //******* Loop over events *********
   Int_t e=0;
   while (esdTree->GetEvent(e)) {
     cout<<endl<<endl<<"********* Processing event number: "<<e<<"*******\n";
     TClonesArray *vertices=event->GetPileupVerticesTracks();
     Int_t nrec=vertices->GetEntriesFast(); 

     if (refTree->GetEvent(e++)==0) {
        cerr<<"No reconstructable vertices for this event !\n";
        continue;
     }
     Int_t ngood=refs->GetEntriesFast(); 
     cout<<"reconstructed: "<<nrec<<"  good: "<<ngood<<endl;

     h2->Fill(ngood,nrec);

   } //***** End of the loop over events

   delete event;
   delete esdTree;
   ef->Close();

   refFile->Close();
   
   TCanvas *c1=new TCanvas("c1","",0,0,750,500);
   h2->Draw("box");

   TFile fc("AliITSUComparisonCooked.root","RECREATE");
   c1->Write();
   fc.Close();

   return 0;
}



Int_t GoodPileupVertices(const Char_t *dir) {
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

     AliHeader *ah=rl->GetHeader();
     AliGenCocktailEventHeader *cock=
            (AliGenCocktailEventHeader*)ah->GenEventHeader();
     TList *headers=cock->GetHeaders();
     const Int_t nvtx=headers->GetEntries();
     cout<<"Event "<<e<<" Number of vertices: "<<nvtx<<endl;

     Int_t nv=0;
     for (Int_t v=0; v<nvtx; v++) {
         AliGenEventHeader *h=(AliGenEventHeader *)headers->At(v);
         TArrayF vtx(3); h->PrimaryVertex(vtx);
         //if (...) continue; // Check if this vertex is reconstructable
         AliESDVertex *vertex=new ((*refs)[nv]) AliESDVertex();
         vertex->SetXv(vtx[0]);
         vertex->SetYv(vtx[1]);
         vertex->SetZv(vtx[2]);
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


