#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <Riostream.h>
  #include "AliCascadeVertexer.h"
  #include "TFile.h"
  #include "TStopwatch.h"

  #include "AliRun.h"
  #include "AliRunLoader.h"
  #include "AliITSLoader.h"
#endif

Int_t AliCascadeFindVertices(Int_t nev=5) {
   cerr<<"Looking for cascade vertices...\n";

   if (gAlice) {
      delete gAlice->GetRunLoader();
      delete gAlice;
      gAlice=0;
   }
   AliRunLoader* rl = AliRunLoader::Open("galice.root");
   if (rl == 0x0) {
      cerr<<"AliCascadeFindVertices.C : Can not open session RL=NULL"<< endl;
      return 1;
   }
   AliITSLoader* itsl = (AliITSLoader*)rl->GetLoader("ITSLoader");
   if (itsl == 0x0) {
      cerr<<"AliCascadeFindVertices.C : Can not get ITS loader"<<endl;
      return 2;
   }
   itsl->LoadTracks("read");
   itsl->LoadV0s("read");
   itsl->LoadCascades("recreate");

   Double_t cuts[]={33.,    // max. allowed chi2
                    0.05,   // min. allowed V0 impact parameter 
                    0.008,  // window around the Lambda mass 
                    0.035,  // min. allowed bachelor's impact parameter 
                    0.10,   // max. allowed DCA between a V0 and a track
                    0.9985, // max. allowed cosine of the cascade pointing angle
                    0.9,    // min. radius of the fiducial volume
                    2.9     // max. radius of the fiducial volume
                   };
   TStopwatch timer;
   AliCascadeVertexer *vertexer=new AliCascadeVertexer(cuts);
   Int_t rc=0;
   if (nev>rl->GetNumberOfEvents()) nev=rl->GetNumberOfEvents();
   for (Int_t i=0; i<nev; i++) {
     rl->GetEvent(i);

     TTree *tTree=itsl->TreeT();
     if (!tTree) {
       cerr<<"AliCascadeFindVertices.C : Can't get the ITS track tree !"<<endl;
       return 3;
     }
     TTree *vTree=itsl->TreeV0();
     if (!vTree) {
       cerr<<"AliCascadeFindVertices.C : Can't get the V0 tree !"<<endl;
       return 4;
     }
     TTree *xTree=itsl->TreeX();
     if (!xTree) {
        itsl->MakeTree("X");
        xTree=itsl->TreeX();
     }

     rc=vertexer->V0sTracks2CascadeVertices(vTree,tTree,xTree);

     itsl->WriteCascades("OVERWRITE");
   }
   delete vertexer;
   timer.Stop(); timer.Print();
    
   delete rl;

   return rc;
}
