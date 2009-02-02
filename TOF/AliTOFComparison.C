/****************************************************************************
 *      This macro estimates efficiency of matching with the TOF.           *
 *      TOF "Good" tracks are those originating from the primary vertex,    *
 *      being "good" in the ITS and having at least one digit in the TOF.   * 
 *         (To get the list of "good" tracks one should first run           *
 *          AliTPCComparison.C and AliITSComparisonV2.C macros)             *
 *           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                 *
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
  #include "AliESDtrack.h"

  #include "AliTOFcluster.h"
  #include "AliLoader.h"

  #include "TClonesArray.h"
#endif

Int_t GoodTracksTOF(const Char_t *dir=".");

extern AliRun *gAlice;
extern TBenchmark *gBenchmark;
extern TROOT *gROOT;

static Int_t allgood=0;
static Int_t allmatched=0;
static Int_t allmismatched=0;

Int_t AliTOFComparison(const Char_t *dir=".") {
   gBenchmark->Start("AliTOFComparison");

   ::Info("AliTOFComparison.C","Doing comparison...");
   

   Double_t pmin=0.2;
   Double_t pmax=3.0;

   TH1F *hgood=(TH1F*)gROOT->FindObject("hgood");    
   if (!hgood) hgood=new TH1F("hgood","Good tracks",30,pmin,pmax);
    
   TH1F *hfound=(TH1F*)gROOT->FindObject("hfound");
   if (!hfound) hfound=new TH1F("hfound","Matched tracks",30,pmin,pmax);

   TH1F *hfake=(TH1F*)gROOT->FindObject("hfake");
   if (!hfake) hfake=new TH1F("hfake","Mismatched tracks",30,pmin,pmax);

   TH1F *hgp=(TH1F*)gROOT->FindObject("hgp");
   if (!hgp) hgp=new TH1F("hgp","",30,pmin,pmax);
   hgp->SetLineColor(4); hgp->SetLineWidth(2);

   TH1F *hfp=(TH1F*)gROOT->FindObject("hfp");
   if (!hfp) hfp=new TH1F("hfp","Probability of mismatching",30,pmin,pmax);
   hfp->SetFillColor(1); hfp->SetFillStyle(3013); hfp->SetLineWidth(2);

   TH1F *hgoo=(TH1F*)gROOT->FindObject("hgoo");    
   if (!hgoo) hgoo=new TH1F("hgoo","Good tracks",30,-1,1);
    
   TH1F *hfoun=(TH1F*)gROOT->FindObject("hfoun");
   if (!hfoun) hfoun=new TH1F("hfoun","Matched tracks",30,-1,1);

   TH1F *hfak=(TH1F*)gROOT->FindObject("hfak");
   if (!hfak) hfak=new TH1F("hfak","Mismatched tracks",30,-1,1);

   TH1F *hgl=(TH1F*)gROOT->FindObject("hgl");
   if (!hgl) hgl=new TH1F("hgl","",30,-1,1);
   hgl->SetLineColor(4); hgl->SetLineWidth(2);

   TH1F *hfl=(TH1F*)gROOT->FindObject("hfl");
   if (!hfl) hfl=new TH1F("hfl","Probability of mismatching",30,-1,1);
   hfl->SetFillColor(1); hfl->SetFillStyle(3013); hfl->SetLineWidth(2);



   Char_t fname[100];
   sprintf(fname,"%s/GoodTracksTOF.root",dir);

   TFile *refFile=TFile::Open(fname,"old");
   if (!refFile || !refFile->IsOpen()) {
   ::Info("AliTOFComparison.C","Marking good tracks (will take a while)...");
     if (GoodTracksTOF(dir)) {
        ::Error("AliTOFComparison.C","Can't generate the reference file !");
        return 1;
     }
   }
   refFile=TFile::Open(fname,"old");
   if (!refFile || !refFile->IsOpen()) {
     ::Error("AliTOFComparison.C","Can't open the reference file !");
     return 1;
   }   

   TTree *tofTree=(TTree*)refFile->Get("tofTree");
   if (!tofTree) {
     ::Error("AliTOFComparison.C","Can't get the reference tree !");
     return 2;
   }
   TBranch *branch=tofTree->GetBranch("TOF");
   if (!branch) {
     ::Error("AliTOFComparison.C","Can't get the TOF branch !");
     return 3;
   }
   TClonesArray dummy("AliTrackReference",1000), *refs=&dummy;
   branch->SetAddress(&refs);


   if (gAlice) { 
     delete AliRunLoader::Instance();
     delete gAlice;//if everything was OK here it is already NULL
     gAlice = 0x0;
   }
   sprintf(fname,"%s/galice.root",dir);
   AliRunLoader *rl = AliRunLoader::Open(fname,"COMPARISON");
   if (rl == 0x0) {
      cerr<<"Can not open session"<<endl;
      return 1;
   }
   AliLoader* tofl = rl->GetLoader("TOFLoader");
   if (tofl == 0x0) {
      cerr<<"Can not get the TOF loader"<<endl;
      return 2;
   }
   tofl->LoadRecPoints("read");


   sprintf(fname,"%s/AliESDs.root",dir);
   TFile *ef=TFile::Open(fname);
   if ((!ef)||(!ef->IsOpen())) {
      ::Error("AliTOFComparison.C","Can't open AliESDs.root !");
      delete rl;
      return 4;
   }
   AliESDEvent* event = new AliESDEvent();
   TTree* esdTree = (TTree*) ef->Get("esdTree");
   if (!esdTree) {
      ::Error("AliTOFComparison.C", "no ESD tree found");
      return 5;
   }
   event->ReadFromTree(esdTree);



   //******* Loop over events *********
   Int_t e=0;
   while (esdTree->GetEvent(e)) {
     cout<<endl<<endl<<"********* Processing event number: "<<e<<"*******\n";

     rl->GetEvent(e);

     TTree *clsTree=tofl->TreeR();
     if (!clsTree) {
        cerr<<"Can't get the TOF cluster tree !\n";
        return 3;
     } 
     TBranch *branch=clsTree->GetBranch("TOF");
     if (!branch) { 
        cerr<<"Can't get the branch with the TOF digits !\n";
        return 4;
     }
     TClonesArray dummy("AliTOFcluster",10000), *clusters=&dummy;   
     branch->SetAddress(&clusters);

     clsTree->GetEvent(0);
     Int_t nd=clusters->GetEntriesFast();
     cerr<<"Number of the TOF clusters: "<<nd<<endl;



     Int_t ntrk=event->GetNumberOfTracks();
     cerr<<"Number of ESD tracks : "<<ntrk<<endl; 


     if (tofTree->GetEvent(e++)==0) {
        cerr<<"No reconstructable tracks !\n";
        continue;
     }

     Int_t ngood=refs->GetEntriesFast(); 

     Int_t matched=0;
     Int_t mismatched=0;
     for (Int_t i=0; i<ngood; i++) {
	AliTrackReference *ref=(AliTrackReference*)refs->UncheckedAt(i); 
        Int_t lab=ref->Label();
        Float_t ptg=TMath::Sqrt(ref->Px()*ref->Px() + ref->Py()*ref->Py());

        Double_t tgl=ref->Pz()/ptg; //tan(lambda)

         if (ptg>pmin) { hgood->Fill(ptg); hgoo->Fill(tgl); }

         Int_t j;
	 AliESDtrack *t=0;
         for (j=0; j<ntrk; j++) {
             AliESDtrack *tt=event->GetTrack(j);
             if (lab!=TMath::Abs(tt->GetLabel())) continue;
             t=tt;
             //if ((tt->GetStatus()&AliESDtrack::kTOFpid) == 0) continue;
             if (tt->GetTOFsignal() < 0) continue;
             UInt_t idx=tt->GetTOFcluster();
             if ((Int_t)idx>=nd) {
	       cerr<<"Wrong cluster index ! "<<idx<<endl;
               return 5;
             }
             AliTOFcluster *cls=(AliTOFcluster*)clusters->At(idx);
	     if (cls) {
	       if (cls->GetLabel(0)!=lab)
		 if (cls->GetLabel(1)!=lab)
		   if (cls->GetLabel(2)!=lab) {
		     mismatched++; 
		     if (ptg>pmin) { hfake->Fill(ptg); hfak->Fill(tgl); } 
		     break;
		   }
	       if (ptg>pmin) { hfound->Fill(ptg); hfoun->Fill(tgl); }
	       matched++;
	       break;
	     }
         }
         if (j==ntrk) {
	    cerr<<"Not matched: "<<lab<<"   ";
            if (t) {
               cerr<<(t->GetStatus()&AliESDtrack::kITSout)<<' '
	           <<(t->GetStatus()&AliESDtrack::kTPCout)<<' '
	           <<(t->GetStatus()&AliESDtrack::kTRDout)<<' '
	           <<(t->GetStatus()&AliESDtrack::kTIME);
	    } else cerr<<"No ESD track !";
            cerr<<endl;
         }
     }
     cout<<"Number of good tracks: "<<ngood<<endl;
     cout<<"Number of matched tracks: "<<matched<<endl;
     cout<<"Number of mismatched tracks: "<<mismatched<<endl;

     allgood+=ngood; allmatched+=matched; allmismatched+=mismatched;

     refs->Clear();
   } //***** End of the loop over events

   delete event;
   delete esdTree;
   ef->Close();
   
   delete tofTree;
   refFile->Close();

   if (allgood!=0) cerr<<"\n\nEfficiency: "<<Float_t(allmatched)/allgood<<endl;
   cout<<"Total number of good tracks: "<<allgood<<endl;
   cout<<"Total number of matched tracks: "<<allmatched<<endl;
   cout<<"Total number of mismatched tracks: "<<allmismatched<<endl;

   TCanvas *c1=(TCanvas*)gROOT->FindObject("c1");
   if (!c1) {
      c1=new TCanvas("c1","",0,0,600,900);
      c1->Divide(1,2);
   }
   hfound->Sumw2(); hgood->Sumw2(); hfake->Sumw2();
   hgp->Divide(hfound,hgood,1,1.,"b");
   hfp->Divide(hfake,hgood,1,1.,"b");
   hgp->SetMaximum(1.4);
   hgp->SetYTitle("Matching efficiency");
   hgp->SetXTitle("Pt (GeV/c)");

   hfoun->Sumw2(); hgoo->Sumw2(); hfak->Sumw2();
   hgl->Divide(hfoun,hgoo,1,1.,"b");
   hfl->Divide(hfak,hgoo,1,1.,"b");
   hgl->SetMaximum(1.4);
   hgl->SetYTitle("Matching efficiency");
   hgl->SetXTitle("Tan(lambda)");

   c1->cd(1);

   hgp->Draw();
   hfp->Draw("histsame");
   TLine *line1 = new TLine(pmin,1.0,pmax,1.0); line1->SetLineStyle(4);
   line1->Draw("same");
   TLine *line2 = new TLine(pmin,0.9,pmax,0.9); line2->SetLineStyle(4);
   line2->Draw("same");

   c1->cd(2);

   hgl->Draw();
   hfl->Draw("histsame");
   TLine *line3 = new TLine(-1,1.0,1,1.0); line3->SetLineStyle(4);
   line3->Draw("same");
   TLine *line4 = new TLine(-1,0.9,1,0.9); line4->SetLineStyle(4);
   line4->Draw("same");

   c1->Update();
   
   TFile fc("AliTOFComparison.root","RECREATE");
   c1->Write();
   fc.Close();

   gBenchmark->Stop("AliTOFComparison");
   gBenchmark->Show("AliTOFComparison");

   delete rl;
   return 0;
}

Int_t GoodTracksTOF(const Char_t *dir) {
   if (gAlice) { 
       delete AliRunLoader::Instance();
       delete gAlice;//if everything was OK here it is already NULL
       gAlice = 0x0;
   }

   Char_t fname[100];
   sprintf(fname,"%s/galice.root",dir);

   AliRunLoader *rl = AliRunLoader::Open(fname,"COMPARISON");
   if (!rl) {
      ::Error("GoodTracksTOF","Can't start session !");
      return 1;
   }

   rl->LoadgAlice();
   rl->LoadHeader();
   rl->LoadKinematics();


   AliLoader* tofl = rl->GetLoader("TOFLoader");
   if (tofl == 0x0) {
      ::Error("GoodTracksTOF","Can not get the TOF loader !");
      delete rl;
      return 2;
   }
   tofl->LoadRecPoints("read");

   Int_t nev=rl->GetNumberOfEvents();
   ::Info("GoodTracksTOF","Number of events : %d\n",nev);  

   sprintf(fname,"%s/GoodTracksITS.root",dir);
   TFile *itsFile=TFile::Open(fname);
   if ((!itsFile)||(!itsFile->IsOpen())) {
       ::Error("GoodTracksTOF","Can't open the GoodTracksITS.root !");
       delete rl;
       return 5; 
   }
   TClonesArray dum("AliTrackReference",1000), *itsRefs=&dum;
   TTree *itsTree=(TTree*)itsFile->Get("itsTree");
   if (!itsTree) {
       ::Error("GoodTracksTOF","Can't get the ITS reference tree !");
       delete rl;
       return 6;
   }
   TBranch *itsBranch=itsTree->GetBranch("ITS");
   if (!itsBranch) {
      ::Error("GoodTracksTOF","Can't get the ITS reference branch !");
      delete rl;
      return 7;
   }
   itsBranch->SetAddress(&itsRefs);


   sprintf(fname,"%s/GoodTracksTOF.root",dir);
   TFile *tofFile=TFile::Open(fname,"recreate");
   TClonesArray dummy("AliTrackReference",1000), *tofRefs=&dummy;
   TTree tofTree("tofTree","Tree with info about the reconstructable TOF tracks");
   tofTree.Branch("TOF",&tofRefs);


   //********  Loop over generated events 
   for (Int_t e=0; e<nev; e++) {

     rl->GetEvent(e); tofFile->cd();

     Int_t np = rl->GetHeader()->GetNtrack();
     cout<<"Event "<<e<<" Number of particles: "<<np<<endl;

     //******** Fill the "good" masks
     Int_t *good=new Int_t[np]; Int_t k; for (k=0; k<np; k++) good[k]=0;

     TTree *cTree=tofl->TreeR();
     if (!cTree) {
        ::Error("GoodTracksTOF","Can't get the TOF cluster tree !");
        delete rl;
        return 8;
     } 
     TBranch *branch=cTree->GetBranch("TOF");
     if (!branch) { 
        ::Error("GoodTracksTOF","Can't get the branch with the TOF digits !");
        return 9;
     }
     TClonesArray dummy("AliTOFcluster",10000), *clusters=&dummy;
     branch->SetAddress(&clusters);
   
     cTree->GetEvent(0);
     Int_t nd=clusters->GetEntriesFast();

     for (Int_t i=0; i<nd; i++) {
       AliTOFcluster *c=(AliTOFcluster*)clusters->UncheckedAt(i);
       Int_t l0=c->GetLabel(0);
          if (l0>=np) {cerr<<"Wrong label: "<<l0<<endl; continue;}
       Int_t l1=c->GetLabel(1);
          if (l1>=np) {cerr<<"Wrong label: "<<l1<<endl; continue;}
       Int_t l2=c->GetLabel(2);
          if (l2>=np) {cerr<<"Wrong label: "<<l2<<endl; continue;}
       if (l0>=0) good[l0]++; 
       if (l1>=0) good[l1]++; 
       if (l2>=0) good[l2]++;
     }
     clusters->Clear();


     //****** select tracks which are "good" enough
     AliStack* stack = rl->Stack();

     itsTree->GetEvent(e);
     Int_t nk=itsRefs->GetEntriesFast();

     Int_t nt=0;
     for (k=0; k<nk; k++) {
        AliTrackReference *itsRef=(AliTrackReference *)itsRefs->UncheckedAt(k);
        Int_t lab=itsRef->Label();
        if (good[lab] == 0) continue;
        TParticle *p = (TParticle*)stack->Particle(lab);
        if (p == 0x0) {
           cerr<<"Can not get particle "<<lab<<endl;
           continue;
        }

        if (TMath::Abs(p->Vx())>0.1) continue;
        if (TMath::Abs(p->Vy())>0.1) continue;
        //if (TMath::Abs(p->Vz())>0.1) continue;

        new((*tofRefs)[nt]) AliTrackReference(*itsRef);
        nt++;
     }
     itsRefs->Clear();

     tofTree.Fill();
     tofRefs->Clear();

     delete[] good;

   } //*** end of the loop over generated events

   tofTree.Write();
   tofFile->Close();

   delete itsTree;
   itsFile->Close();

   delete rl;
   return 0;
}
