/****************************************************************************
 *           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch                 *
 ****************************************************************************/

#if !defined(__CINT__) || defined(__MAKECINT__)
  #include "Riostream.h"
  #include "AliV0vertexer.h"
  #include "TFile.h"
  #include "TKey.h"
  #include "TStopwatch.h"

  #include "AliRun.h"
  #include "AliTracker.h"
  #include "AliMagF.h"
  #include "AliESD.h"
  #include "AliRunLoader.h"
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

   if (rl->LoadgAlice()) {
      cerr<<"AliV0FindVertices.C : LoadgAlice returned error"<<endl;
      delete rl;
      return 3;
   }

   // Magnetic field
   AliTracker::SetFieldMap(gAlice->Field(),1); // 1 means uniform magnetic field
       
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

   TFile *v0f=TFile::Open("AliESDv0.root","RECREATE");
   if ((!v0f)||(!v0f->IsOpen())) {
      cerr<<"Can't AliESDv0.root !\n"; return 1;
   }
   TFile *itsf=TFile::Open("AliESDits.root");
   if ((!itsf)||(!itsf->IsOpen())) {
      cerr<<"Can't AliESDits.root !\n"; return 1;
   }
   TKey *key=0;
   TIter next(itsf->GetListOfKeys());
   for (Int_t i=0; i<nev; i++) {
     itsf->cd();
     if ((key=(TKey*)next())==0) break;
     cerr<<"Processing event number: "<<i<<endl;
     AliESD *event=(AliESD*)key->ReadObj();

     //Double_t vtx[3]={0.,0.,0.}; vtxer.SetVertex(vtx); // primary vertex (cm)

     rc=vtxer.Tracks2V0vertices(event);

     if (rc==0) {
        Char_t ename[100]; 
        sprintf(ename,"%d",i);
        v0f->cd();
        if (!event->Write(ename)) rc++;
     } 
     if (rc) {
        cerr<<"Something bad happened...\n";
     }
     delete event;
   }
   timer.Stop(); timer.Print();
    
   itsf->Close();
   v0f->Close();

   delete rl;

   return rc;
}
