// The following macro runs the HLT ITS tracker over the HLT
// tracks stored in the ESD and stores the output tracks in a
// separate ESD file AliESDits.root

#if !defined(__CINT__) || defined(__MAKECINT__)
  #include <TMath.h>
  #include <TError.h>
  #include <Riostream.h>
  #include <TH1F.h>
  #include <TH2F.h>
  #include <TNtuple.h>
  #include <TProfile2D.h>
  #include <TF1.h>
  #include <TGraphErrors.h>
  #include <TTree.h>
  #include <TParticle.h>
  #include <TCanvas.h>
  #include <TFile.h>
  #include <TROOT.h>
  #include <TStopwatch.h>
  #include <TSystem.h>

  #include "AliStack.h"
  #include "AliHeader.h"
  #include "AliTrackReference.h"
  #include "AliRunLoader.h"
  #include "AliITS.h"
  #include "AliITSLoader.h"
  #include "AliITSgeom.h"
  #include "AliITStrackerV2.h"
  #include "AliRun.h"
  #include "AliESD.h"
  #include "AliMagF.h"
  #include "AliGenEventHeader.h"

  #include "AliHLTITStrack.h"
  #include "AliHLTITStracker.h"
  #include "AliHLTITSVertexerZ.h"
#endif

//extern TSystem *gSystem;

Int_t RunHLTITS(Int_t nev=1,Int_t run=0) {

  //  gSystem->Load("libAliHLTITS.so");

  TStopwatch timer;
  timer.Start();

   if (gAlice) {
      delete gAlice->GetRunLoader();
      delete gAlice; 
      gAlice=0;
   }

   AliRunLoader *rl = AliRunLoader::Open("galice.root");
   if (rl == 0x0) {
      cerr<<"Can not open session"<<endl;
      return 1;
   }
   Int_t retval = rl->LoadgAlice();
   if (retval) {
      cerr<<"AliESDtest.C : LoadgAlice returned error"<<endl;
      delete rl;
      return 1;
   }
   retval = rl->LoadHeader();
   if (retval) {
      cerr<<"AliESDtest.C : LoadHeader returned error"<<endl;
      delete rl;
      return 2;
   }
   gAlice=rl->GetAliRun();
       
   AliTracker::SetFieldMap(gAlice->Field());

   AliITSLoader* itsl = (AliITSLoader*)rl->GetLoader("ITSLoader");
   if (itsl == 0x0) {
      cerr<<"AliESDtest.C : Can not get the ITS loader"<<endl;
      return 3;
   }
   itsl->LoadRecPoints("read");

   AliITS *dITS = (AliITS*)gAlice->GetDetector("ITS");
   if (!dITS) {
      cerr<<"AliESDtest.C : Can not find the ITS detector !"<<endl;
      return 4;
   }
   //   AliITSgeom *geom = dITS->GetITSgeom();
   AliITSgeom *geom = new AliITSgeom();
   geom->ReadNewFile("$ALICE_ROOT/ITS/ITSgeometry_vPPRasymmFMD.det");

   //An instance of the HLT ITS tracker
   AliHLTITStracker itsTracker(geom);

   TFile *ef=TFile::Open("AliESDs.root");
   if (!ef || !ef->IsOpen()) {cerr<<"Can't AliESDs.root !\n"; return 1;}
   AliESD* event = new AliESD;
   TTree* tree = (TTree*) ef->Get("esdTree");
   if (!tree) {cerr<<"no ESD tree found\n"; return 1;};
   tree->SetBranchAddress("ESD", &event);

   TFile *itsf=TFile::Open("AliESDits.root","RECREATE");
   if ((!itsf)||(!itsf->IsOpen())) {
      cerr<<"Can't AliESDits.root !\n"; return 1;
   }

   Int_t rc=0;
   if (nev>rl->GetNumberOfEvents()) nev=rl->GetNumberOfEvents();
   //The loop over events
   for (Int_t i=0; i<nev; i++) {

     cerr<<"\n\nProcessing event number : "<<i<<endl;
     tree->GetEvent(i);
     rl->GetEvent(i);

     TArrayF v(3);     
     rl->GetHeader()->GenEventHeader()->PrimaryVertex(v);
     Double_t vtx[3]={v[0],v[1],v[2]};
     Double_t cvtx[3]={0.005,0.005,0.010};
     cout<<"MC vertex position: "<<v[2]<<endl;

     AliHLTITSVertexerZ vertexer("null");
     AliESDVertex* vertex = NULL;
     TStopwatch timer2;
     timer2.Start();
     TTree* treeClusters = itsl->TreeR();
     //     vertex = vertexer.FindVertexForCurrentEvent(i);
     //     AliESDVertex *vertex = vertexer.FindVertexForCurrentEvent(geom,treeClusters);
     vertex = new AliESDVertex(vtx,cvtx);
     timer2.Stop();
     timer2.Print();
     if(!vertex){
       cerr<<"Vertex not found"<<endl;
       vertex = new AliESDVertex(vtx,cvtx);
     }
     else {
       vertex->SetTruePos(vtx);  // store also the vertex from MC
     }

     event->SetVertex(vertex);

     Double_t vtxPos[3];
     Double_t vtxErr[3];
     vertex->GetXYZ(vtxPos);
     vertex->GetSigmaXYZ(vtxErr);
     itsTracker.SetVertex(vtxPos,vtxErr);

     TTree *itsTree=itsl->TreeR();
     if (!itsTree) {
       cerr<<"Can't get the ITS cluster tree !\n";
       return 4;
     }     
     itsTracker.LoadClusters(itsTree);
     rc+=itsTracker.Clusters2Tracks(event);
     //     rc+=itsTracker.PropagateBack(event);
     itsTracker.UnloadClusters();

     if (rc==0) {
       TTree* tree = new TTree("esdTree", "Tree with ESD objects");
       tree->Branch("ESD", "AliESD", &event);
       tree->Fill();
       itsf->cd();
       tree->Write();
     } 
     if (rc) {
       cerr<<"Something bad happened...\n";
     }

   }
   delete event;

   itsf->Close();
   ef->Close();

   //   delete rl;

   timer.Stop();
   timer.Print();

   return rc;
}
