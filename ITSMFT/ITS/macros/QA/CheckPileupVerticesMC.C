//
// This macro checks a few basic properties of MC pileup vertices
//

#if !defined(__CINT__) || defined(__MAKECINT__)
   #include <Riostream.h>
   #include <TFile.h>
   #include <TTree.h>
   #include <TH1F.h>
   #include <TCanvas.h>
   #include <TMath.h>
   #include <TClonesArray.h>
   #include <TParticle.h>
   #include <TParticlePDG.h>

   #include "AliRunLoader.h"
   #include "AliHeader.h"
   #include "AliStack.h"
   #include "AliGenCocktailEventHeader.h"
#endif



Int_t FindContributors(Float_t t, AliStack *stack) {
  Int_t ntrk=0;
  Int_t np=stack->GetNtrack();
  for (Int_t i=0; i<np; i++) {
      if (!stack->IsPhysicalPrimary(i)) continue;
      TParticle *part=stack->Particle(i);
      if (!part) continue;
      TParticlePDG *partPDG = part->GetPDG();
      if (!partPDG) continue;
      if (TMath::Abs(partPDG->Charge())<1e-10) continue;
      if (part->Pt()<0.5) continue;
      Float_t dt=0.5*(t - part->T())/(t + part->T());
      if (TMath::Abs(dt)>1e-5) continue;
      ntrk++;
  }
  return ntrk;
}


void CheckPileupVerticesMC() {
    TH1F *hnv=new TH1F("hnv","Number of vertices",100,0,300);
    TH1F *hnt=new TH1F("hnt","Number of tracks per vertex",100,0,200);
    TH1F *ht=new TH1F("ht","Interaction time;t (sec);", 100,-0.05e-3,0.05e-3);
    TH1F *hdt=new TH1F("hdt","Difference in time between interactions; #Deltat (sec)",100,0,0.01e-3);
    TH1F *hz=new TH1F("hz","Z position of vertices; Z (cm)",100,-20,20);
    TH1F *hdz=new TH1F("hdz","Difference in Z position between vertices; #DeltaZ (cm)",100,0,2.5);

    AliRunLoader *rl = AliRunLoader::Open("galice.root","something");
    rl->LoadHeader();
    rl->LoadKinematics();

    Int_t nEvents=rl->GetNumberOfEvents();
    for (Int_t i=0; i<nEvents; i++) {
        cout<<"Event "<<i<<" out of "<<nEvents<<endl;
        rl->GetEvent(i);
        AliStack *stack=rl->Stack();    

        AliHeader *ah=rl->GetHeader();
        AliGenCocktailEventHeader *cock=
            (AliGenCocktailEventHeader*)ah->GenEventHeader();
        TList *headers=cock->GetHeaders();
        const Int_t nvtx=headers->GetEntries();

        hnv->Fill(nvtx);
        Double_t *z=new Double_t[nvtx];
        Double_t *t=new Double_t[nvtx];
        Int_t *idx=new Int_t[nvtx];
        for (Int_t v=0; v<nvtx; v++) {
            AliGenEventHeader *h=(AliGenEventHeader *)headers->At(v);
            t[v]=h->InteractionTime();
            ht->Fill(t[v]);
            Int_t nt=FindContributors(t[v],stack);
            hnt->Fill(nt);

            TArrayF vtx(3); h->PrimaryVertex(vtx);
            z[v]=vtx[2];
            hz->Fill(z[v]);
        }

	TMath::Sort(nvtx,t,idx);
        for (Int_t v=0; v<nvtx-1; v++) {
	  Int_t i1=idx[v];
	  Int_t i2=idx[v+1];
          Float_t dt=t[i1]-t[i2];
          hdt->Fill(dt);
	}        
	TMath::Sort(nvtx,z,idx);
        for (Int_t v=0; v<nvtx-1; v++) {
	  Int_t i1=idx[v];
	  Int_t i2=idx[v+1];
          Float_t dz=z[i1]-z[i2];
          hdz->Fill(dz);
	}        

        delete[] t;
        delete[] z;
        delete[] idx;

    }

    TCanvas *c1=new TCanvas("c1","",0,0,750,1000); c1->Divide(1,2); 
    c1->cd(1); hnv->Draw();
    c1->cd(2); gPad->SetLogy(); hnt->Draw();

    TCanvas *c2=new TCanvas("c2","",0,0,750,1000); c2->Divide(1,2); 
    c2->cd(1); ht->Draw();
    c2->cd(2); hdt->Fit("expo");

    TCanvas *c3=new TCanvas("c3","",0,0,750,1000); c3->Divide(1,2); 
    c3->cd(1); hz->Fit("gaus");
    c3->cd(2); hdz->Draw(); //hdz->Fit("expo");

}
