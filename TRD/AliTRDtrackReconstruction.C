#ifndef __CINT__
  #include "alles.h"
  #include "AliMagF.h"
  #include "AliTRDtracker.h"
#endif

Int_t TRDPropagateBack(const Char_t *geoname, const Char_t *clrsname, 
		       const Char_t *inname, const Char_t *outname, Int_t n);

Int_t TRDFindTracks(const Char_t *geoname, const Char_t *clrsname,
		    const Char_t *inname, const Char_t *outname, Int_t n);

Int_t AliTRDtrackReconstruction(Int_t n=1) {
   const Char_t *TRDdigName="galice.root";
   const Char_t *dummyName="dummy.root";     
   const Char_t *TRDclsName="AliTRDclusters.root";
   const Char_t *TRDtrkName="AliTRDtracks.root";
   const Char_t *TPCbkTrkName="AliTPCBackTracks.root";

   AliKalmanTrack::SetConvConst(100/0.299792458/0.2/gAlice->Field()->Factor());

    
// ********** Find TRD tracks from TPC back propagated tracks *********** //
   
   
   if (TRDPropagateBack(TRDclsName, TRDclsName, TPCbkTrkName, TRDtrkName, n)) {
      cerr<<"Failed to propagate back through TRD !\n";
      return 1;
   } 
   
  
// ********** Find TRD tracks and make seeds for TPC *********** //

   /*
   if (TRDFindTracks(TRDclsName,TRDclsName, TRDtrkName, TRDtrkName,n)) {
     cerr<<"Failed to find TRD tracks !\n";
     return 1;
   }
   */
     
   return 0;
}
   

Int_t TRDPropagateBack(const Char_t *geoname, const Char_t *clrsname,
		       const Char_t *inname, const Char_t *outname, Int_t n) 
{
   Int_t rc=0;
   const Char_t *name="TRDPropagateBack";
   cerr<<'\n'<<name<<"...\n";
   gBenchmark->Start(name);

   TFile *geofile =TFile::Open(geoname);
   TFile *out=TFile::Open(outname,"update");
   TFile *in =TFile::Open(inname);
   TFile *clrsfile =TFile::Open(clrsname);

   AliTRDtracker *tracker=new AliTRDtracker(geofile);

   for (Int_t i=0;i<n;i++){
     printf("Processing event %d\n",i);
     tracker->SetEventNumber(i);
     rc=tracker->PropagateBack(in,out);
   }

   delete tracker;
   in->Close();
   out->Close();
   geofile->Close();
   clrsfile->Close();

   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}


Int_t TRDFindTracks(const Char_t *geoname, const Char_t *clrsname,
		    const Char_t *inname, const Char_t *outname, Int_t n)
{
   Int_t rc=0;
   const Char_t *name="TRDFindTracks";
   cerr<<'\n'<<name<<"...\n";
   gBenchmark->Start(name);

   TFile *geofile =TFile::Open(geoname);
   TFile *out=TFile::Open(outname,"update");
   TFile *in =TFile::Open(inname);
   TFile *clrsfile =TFile::Open(clrsname);

   AliTRDtracker *tracker=new AliTRDtracker(geofile);
   tracker->SetAddTRDseeds();

   for (Int_t i=0;i<n;i++){
     printf("Processing event %d\n",i);
     tracker->SetEventNumber(i);
     rc=tracker->Clusters2Tracks(in,out);
   }

   delete tracker;
   in->Close();
   out->Close();
   geofile->Close();
   clrsfile->Close();

   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}

