#ifndef __CINT__
  #include <iostream.h>
  #include "AliV0vertexer.h"
  #include "TFile.h"
  #include "TStopwatch.h"
#endif

Int_t AliV0FindVertices() {
   cerr<<"Looking for V0 vertices...\n";

   TFile *out=TFile::Open("AliV0vertices.root","new");
   if (!out->IsOpen()) {cerr<<"Delete old AliV0vertices.root !\n"; return 1;}

   TFile *in=TFile::Open("AliITStracksV2.root");
   if (!in->IsOpen()) {cerr<<"Can't open AliITStracksV2.root !\n"; return 2;}

   Double_t cuts[]={33,  // max. allowed chi2
                    0.12,// min. allowed negative daughter's impact parameter 
                    0.06,// min. allowed positive daughter's impact parameter 
                    0.080,// max. allowed DCA between the daughter tracks
                    0.998,// max. allowed cosine of V0's pointing angle
                    0.9,  // min. radius of the fiducial volume
                    2.9   // max. radius of the fiducial volume
                   };
   TStopwatch timer;
   AliV0vertexer *vertexer=new AliV0vertexer(cuts);
   Int_t rc=vertexer->Tracks2V0vertices(in,out);
   delete vertexer;
   timer.Stop(); timer.Print();
    
   in->Close();
   out->Close();

   return rc;
}
