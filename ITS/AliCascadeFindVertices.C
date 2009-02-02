#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <Riostream.h>
  #include "AliCascadeVertexer.h"
  #include "TFile.h"
  #include "TKey.h"
  #include "TStopwatch.h"

  #include "AliRun.h"
  #include "AliESD.h"
  #include "AliTracker.h"
  #include "AliRunLoader.h"
#endif

Int_t AliCascadeFindVertices(Int_t nev=5) {
   cerr<<"Looking for cascade vertices...\n";

   if (gAlice) {
      delete AliRunLoader::Instance();
      delete gAlice;
      gAlice=0;
   }
   AliRunLoader* rl = AliRunLoader::Open("galice.root");
   if (rl == 0x0) {
      cerr<<"AliCascadeFindVertices.C : Can not open session RL=NULL"<< endl;
      return 1;
   }

   if (rl->LoadgAlice()) {
      cerr<<"AliV0FindVertices.C : LoadgAlice returned error"<<endl;
      delete rl;
      return 3;
   }

   // Magnetic field
   AliTracker::SetFieldMap(gAlice->Field(),1); // 1 means uniform magnetic field
   Double_t cuts[]={33.,    // max. allowed chi2
                    0.05,   // min. allowed V0 impact parameter 
                    0.008,  // window around the Lambda mass 
                    0.035,  // min. allowed bachelor's impact parameter 
                    0.10,   // max. allowed DCA between a V0 and a track
                    0.9985,// max. allowed cosine of the cascade pointing angle
                    0.9,    // min. radius of the fiducial volume
                    2.9     // max. radius of the fiducial volume
                   };
   TStopwatch timer;
   AliCascadeVertexer *vertexer=new AliCascadeVertexer(cuts);

   Int_t rc=0;
   if (nev>rl->GetNumberOfEvents()) nev=rl->GetNumberOfEvents();

   TFile *casf=TFile::Open("AliESDcas.root","RECREATE");
   if ((!casf)||(!casf->IsOpen())) {
      cerr<<"Can't AliESDcas.root !\n"; return 1;
   }
   TFile *v0f=TFile::Open("AliESDv0.root");
   if ((!v0f)||(!v0f->IsOpen())) {
      cerr<<"Can't AliESDv0.root !\n"; return 1;
   }

   TKey *key=0;
   TIter next(v0f->GetListOfKeys());
   for (Int_t i=0; i<nev; i++) {
     v0f->cd();
     if ((key=(TKey*)next())==0) break;
     cerr<<"Processing event number: "<<i<<endl;
     AliESD *event=(AliESD*)key->ReadObj();

     rc=vertexer->V0sTracks2CascadeVertices(event);

     if (rc==0) {
        Char_t ename[100]; 
        sprintf(ename,"%d",i);
        casf->cd();
        if (!event->Write(ename)) rc++;
     } 
     if (rc) {
        cerr<<"Something bad happened...\n";
     }
     delete event;
   }
   delete vertexer;
   timer.Stop(); timer.Print();
    
   v0f->Close();
   casf->Close();

   delete rl;

   return rc;
}
