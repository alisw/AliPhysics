#ifndef __CINT__
  #include <iostream.h>
  #include "AliCascadeVertexer.h"
  #include "TFile.h"
  #include "TStopwatch.h"
#endif

Int_t AliCascadeFindVertices() {
   cerr<<"Looking for cascade vertices...\n";

   TFile *out=TFile::Open("AliCascadeVertices.root","new");
   if (!out->IsOpen()) {
      cerr<<"Delete old AliCascadeVertices.root !\n"; return 1;
   }
   TFile *in=TFile::Open("AliITStracksV2.root");
   if (!in->IsOpen()) {cerr<<"Can't open AliITStracksV2.root !\n"; return 2;}

   TFile *file=TFile::Open("AliV0vertices.root");
   if (!file->IsOpen()) {
      cerr<<"Can't open AliV0vertices.root !\n";return 3;
   }
   Double_t cuts[]={33.,  // max. allowed chi2
                    0.015,// min. allowed V0 impact parameter 
                    0.05, // window around the Lambda mass 
                    0.015,// min. allowed track impact parameter 
                    0.060,// max. allowed DCA between a V0 and a track
                    0.997,// max. allowed cosine of the cascade pointing angle
                    0.9,  // min. radius of the fiducial volume
                    2.9   // max. radius of the fiducial volume
                   };
   TStopwatch timer;
   AliCascadeVertexer *vertexer=new AliCascadeVertexer(cuts);
   Int_t rc=vertexer->V0sTracks2CascadeVertices(in,out);
   delete vertexer;
   timer.Stop(); timer.Print();
    
   file->Close();
   in->Close();
   out->Close();

   return rc;
}
