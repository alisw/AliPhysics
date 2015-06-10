
#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <TMath.h>
  #include <TError.h>
  #include <Riostream.h>
  #include <TH1F.h>
  #include <TH2F.h>
  #include <TTree.h>
  #include <TParticle.h>
  #include <TCanvas.h>
  #include <TLine.h>
  #include <TText.h>
  #include <TBenchmark.h>
  #include <TStyle.h>
  #include <TFile.h>
  #include <TROOT.h>
  #include <TNtuple.h>
  #include <TEllipse.h>
  #include <TGeoManager.h>

  #include "AliStack.h"
  #include "AliHeader.h"
  #include "AliTrackReference.h"
  #include "AliRunLoader.h"
  #include "AliRun.h"
  #include "AliESDEvent.h"
  #include "AliESDtrack.h"

  #include "AliITSUClusterPix.h"
  #include "AliITSULoader.h"
  #include "AliITSUGeomTGeo.h"
#endif


Int_t check_radii(const Char_t *dir=".") {
   TFile *f=TFile::Open("xyz.root","recreate");
   TNtuple *nt=new TNtuple("nt","my ntuple","x:y:z");

   if (gAlice) { 
       delete AliRunLoader::Instance();
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

   AliITSULoader* itsl = (AliITSULoader*)rl->GetLoader("ITSLoader");
   if (itsl == 0x0) {
       ::Error("GoodTracksITS","Can not find the ITSLoader");
       delete rl;
       return 4;
   }
   itsl->LoadRecPoints();
  
   TGeoManager::Import("geometry.root");
   AliITSUGeomTGeo* gm = new AliITSUGeomTGeo(kTRUE,kTRUE);
   AliITSUClusterPix::SetGeom(gm);

   Int_t nev=rl->GetNumberOfEvents();
   ::Info("GoodTracksITS","Number of events : %d\n",nev);  

   //********  Loop over generated events 
   for (Int_t e=0; e<nev; e++) {
     Int_t k;

     rl->GetEvent(e);

     TTree *cTree=itsl->TreeR();
     if (!cTree) {
        ::Error("GoodTracksITS","Can't get the cluster tree !"); 
        delete rl;
        return 8;
     }

     const Int_t nLayers=7;
     TBranch *branch[nLayers];
     TClonesArray clusters[nLayers];
     for (Int_t layer=0; layer<nLayers; layer++) {
       TClonesArray *ptr = 
       new(clusters+layer) TClonesArray("AliITSUClusterPix",1000);
       Char_t bname[33];
       sprintf(bname,"ITSRecPoints%d\0",layer);
       branch[layer]=cTree->GetBranch(bname);
       if (!branch[layer]) {
          ::Error("GoodTracksITS","Can't get the clusters branch !"); 
          delete rl;
          return 9;
       }
       branch[layer]->SetAddress(&ptr);
     }

     Int_t entr=(Int_t)cTree->GetEntries();
     for (k=0; k<entr; k++) {
         cTree->GetEvent(k);
         for (Int_t lay=0; lay<nLayers; lay++) {
             Int_t ncl=clusters[lay].GetEntriesFast(); if (ncl==0) continue;
             while (ncl--) {
                AliITSUClusterPix *pnt=
                (AliITSUClusterPix*)clusters[lay].UncheckedAt(ncl);
                Float_t g[3];
                pnt->GetGlobalXYZ(g);
                //pnt->GetLocalXYZ(g);
                //cout<<g[0]<<' '<<g[1]<<' '<<g[2]<<endl;
                nt->Fill(g);
	     }
             clusters[lay].Clear();
	 }
     }


   } //*** end of the loop over generated events

   nt->Draw("y:x");
   Float_t rmin=19.44;
   TEllipse *ellipse = new TEllipse(0,0,rmin,rmin,0,360,0);
   ellipse->SetFillStyle(0);
   ellipse->SetLineColor(2);
   ellipse->Draw();
   
   Float_t rmax=19.77;
   ellipse = new TEllipse(0,0,rmax,rmax,0,360,0);
   ellipse->SetFillStyle(0);
   ellipse->SetLineColor(4);
   ellipse->Draw();

   delete rl;
   return 0;
}


