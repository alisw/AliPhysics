#if !defined(__CINT__) || defined(__MAKECINT__)
  #include "Riostream.h"
  #include "AliV0vertexer.h"
  #include "TFile.h"
  #include "TStopwatch.h"

  #include "AliRun.h"
  #include "AliRunLoader.h"
  #include "AliITSLoader.h"
#endif

extern AliRun *gAlice;

Int_t AliV0FindVertices(Int_t nev=5) {
   cerr<<"Looking for V0 vertices...\n";

   if (gAlice) {
      delete gAlice->GetRunLoader();
      delete gAlice; 
      gAlice=0;
   } 
   AliRunLoader* rl = AliRunLoader::Open("galice.root");
   if (rl == 0x0) {
      cerr<<"AliV0FindVertices.C : Can not open session RL=NULL"<< endl;
      return 1;
   }
   AliITSLoader* itsl = (AliITSLoader*)rl->GetLoader("ITSLoader");
   if (itsl == 0x0) {
      cerr<<"AliV0FindVertices.C : Can not get ITS loader"<<endl;
      return 2;
   }
   itsl->LoadTracks("read");
   itsl->LoadV0s("recreate");

   Double_t cuts[]={33,  // max. allowed chi2
                    0.16,// min. allowed negative daughter's impact parameter 
                    0.05,// min. allowed positive daughter's impact parameter 
                    0.080,// max. allowed DCA between the daughter tracks
                    0.998,// max. allowed cosine of V0's pointing angle
                    0.9,  // min. radius of the fiducial volume
                    2.9   // max. radius of the fiducial volume
                   };
   TStopwatch timer;
   AliV0vertexer vtxer(cuts);
   Int_t rc=0;
   if (nev>rl->GetNumberOfEvents()) nev=rl->GetNumberOfEvents();
   for (Int_t i=0; i<nev; i++) {
     rl->GetEvent(i);
     //Double_t vtx[3]={0.,0.,0.}; vtxer.SetVertex(vtx); // primary vertex (cm)

     TTree *tTree=itsl->TreeT();
     if (!tTree) {
        cerr<<"AliV0FindVertices.C : Can't get the ITS track tree !"<<endl;
        return 3;
     }
     TTree *vTree=itsl->TreeV0();
     if (!vTree) {
        itsl->MakeTree("V0");
        vTree=itsl->TreeV0();
     }

     rc=vtxer.Tracks2V0vertices(tTree,vTree);

     itsl->WriteV0s("OVERWRITE");
   }
   timer.Stop(); timer.Print();
    
   delete rl;

   return rc;
}
