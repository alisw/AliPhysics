#ifndef __CINT__
  #include "Riostream.h"
  #include "AliV0vertexer.h"
  #include "TFile.h"
  #include "TStopwatch.h"
#endif

Int_t AliV0FindVertices(Int_t nev=1) {
   cerr<<"Looking for V0 vertices...\n";

   TFile *out=TFile::Open("AliV0vertices.root","new");
   if (!out->IsOpen()) {cerr<<"Delete old AliV0vertices.root !\n"; return 1;}

   TFile *in=TFile::Open("AliITStracksV2.root");
   if (!in->IsOpen()) {cerr<<"Can't open AliITStracksV2.root !\n"; return 2;}

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
   for (Int_t i=0; i<nev; i++) {
     //Double_t vtx[3]={0.,0.,0.}; vtxer.SetVertex(vtx); // primary vertex (cm)
     vtxer.SetEvent(i);
     rc=vtxer.Tracks2V0vertices(in,out);
   }
   timer.Stop(); timer.Print();
    
   in->Close();
   out->Close();

   return rc;
}
