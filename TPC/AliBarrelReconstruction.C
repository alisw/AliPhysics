/****************************************************************************
 * This macro is supposed to do reconstruction in the barrel ALICE trackers *
 * (Make sure you have TPC digits and ITS hits before using this macro !!!) *
 * It does the following steps (April 12, 2001):                            *
 *                   1) TPC cluster finding                                 *
 *                   2) TPC track finding                                   *
 *                   3) ITS cluster finding V2 (via fast points !)          *
 *                   4) ITS track finding V2                                *
 *                   5) ITS back track propagation V2                       *
 *                   6) TPC back track propagation                          *
 *                (Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch)          *
 ****************************************************************************/

#ifndef __CINT__
  #include "alles.h"
  #include "AliMagF.h"
  #include "AliTPCtracker.h"

  #include "AliITS.h"
  #include "AliITSgeom.h"
  #include "AliITSRecPoint.h"
  #include "AliITSclusterV2.h"
  #include "AliITSsimulationFastPoints.h"
  #include "AliITStrackerV2.h"
#endif

Int_t TPCFindClusters(const Char_t *inname, const Char_t *outname);
Int_t TPCFindTracks(const Char_t *inname, const Char_t *outname);
Int_t TPCSortTracks(const Char_t *inname, const Char_t *outname);
Int_t TPCPropagateBack(const Char_t *inname, const Char_t *outname);

Int_t ITSFindClusters(const Char_t *inname, const Char_t *outname);
Int_t ITSFindTracks(const Char_t *inname, const Char_t *outname);
Int_t ITSPropagateBack(const Char_t *inname, const Char_t *outname);


Int_t AliBarrelReconstruction() {
   const Char_t *TPCdigName="galice.root";
   const Char_t *TPCclsName="AliTPCclusters.root";
   const Char_t *TPCtrkName="AliTPCtracks.root";

   const Char_t *ITSdigName="galice.root";
   const Char_t *ITSclsName="AliITSclustersV2.root";
   const Char_t *ITStrkName="AliITStracksV2.root";

   AliKalmanTrack::SetConvConst(100/0.299792458/0.2/gAlice->Field()->Factor());

// ********** Find TPC clusters *********** //
   if (TPCFindClusters(TPCdigName,TPCclsName)) {
      cerr<<"Failed to get TPC clusters !\n";
      return 1;
   }      

// ********** Find TPC tracks *********** //
   if (TPCFindTracks(TPCclsName,TPCtrkName)) {
      cerr<<"Failed to get TPC tracks !\n";
      return 1;
   }

// ********** Sort and label TPC tracks *********** //
   if (TPCSortTracks(TPCclsName,TPCtrkName)) {
      cerr<<"Failed to sort TPC tracks !\n";
      return 1;
   }

// ********** Find ITS clusters *********** //
   if (ITSFindClusters(ITSdigName,ITSclsName)) {
      cerr<<"Failed to get ITS clusters !\n";
      return 1;
   }

// ********** Find ITS tracks *********** //
   {TFile *clsFile=TFile::Open(ITSclsName);
   if (ITSFindTracks(TPCtrkName,ITStrkName)) {
      cerr<<"Failed to get ITS tracks !\n";
      return 1;
   }
   clsFile->Close();}

// ********** Back propagation of the ITS tracks *********** //
   {TFile *clsFile=TFile::Open(ITSclsName);
   if (ITSPropagateBack(ITStrkName,TPCtrkName)) {
      cerr<<"Failed to propagate back the ITS tracks !\n";
      return 1;
   }
   clsFile->Close();}


// ********** Back propagation of the TPC tracks *********** //
   {TFile *clsFile=TFile::Open(TPCclsName);
   if (TPCPropagateBack(TPCtrkName,TPCtrkName)) {
      cerr<<"Failed to propagate back the TPC tracks !\n";
      return 1;
   }
   clsFile->Close();}

   return 0;
}


Int_t TPCFindClusters(const Char_t *inname, const Char_t *outname) {
   Int_t rc=0;
   const Char_t *name="TPCFindClusters";
   cerr<<'\n'<<name<<"...\n";
   gBenchmark->Start(name);

   TFile *out=TFile::Open(outname,"recreate");
   TFile *in =TFile::Open(inname);

   AliTPCParam *param=(AliTPCParam *)in->Get("75x40_100x60");
   if (!param) {cerr<<"TPC parameters have not been found !\n"; return 1;}
   AliTPCv2 tpc;
   tpc.SetParam(param);

   tpc.Digits2Clusters(out);

   in->Close();
   out->Close();
   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}

Int_t TPCFindTracks(const Char_t *inname, const Char_t *outname) {
   Int_t rc=0;
   const Char_t *name="TPCFindTracks";
   cerr<<'\n'<<name<<"...\n";
   gBenchmark->Start(name);
   TFile *out=TFile::Open(outname,"recreate");
   TFile *in =TFile::Open(inname);
   AliTPCParam *param=(AliTPCParam *)in->Get("75x40_100x60");
   if (!param) {cerr<<"TPC parameters have not been found !\n"; return 1;}

   AliTPCtracker *tracker=new AliTPCtracker(param);
   rc=tracker->Clusters2Tracks(0,out);
   delete tracker;
 
   in->Close();
   out->Close();
   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}

Int_t TPCSortTracks(const Char_t *inname, const Char_t *outname) {
   Int_t rc=0;
   const Char_t *name="TPCSortTracks";
   cerr<<'\n'<<name<<"...\n";
   gBenchmark->Start(name);

   TFile *out=TFile::Open(outname);
   TFile *in =TFile::Open(inname);
   AliTPCParam *param=(AliTPCParam *)in->Get("75x40_100x60");
   if (!param) {cerr<<"TPC parameters have not been found !\n"; return 1;}

   TObjArray tarray(10000);
   AliTPCtrack *iotrack=0;
   Int_t i;

   out->cd();
   TTree *tracktree=(TTree*)out->Get("TPCf");
   TBranch *tbranch=tracktree->GetBranch("tracks");
   Int_t nentr=(Int_t)tracktree->GetEntries();
   for (i=0; i<nentr; i++) {
       iotrack=new AliTPCtrack;
       tbranch->SetAddress(&iotrack);
       tracktree->GetEvent(i);
       tarray.AddLast(iotrack);
   }   
   tarray.Sort();
   out->Close();
   
   //assign thacks GEANT labels
   in->cd();
   AliTPCtracker *tracker = new AliTPCtracker(param);
   tracker->LoadInnerSectors();
   tracker->LoadOuterSectors();
   for (i=0; i<nentr; i++) {
       iotrack=(AliTPCtrack*)tarray.UncheckedAt(i);
       tracker->CookLabel(iotrack,0.1);
   }   
   delete tracker;
   in->Close();
   //end of GEANT label assignment

   out=TFile::Open(outname,"recreate");
   tracktree=new TTree("TPCf","Tree with TPC tracks");
   tracktree->Branch("tracks","AliTPCtrack",&iotrack,32000,0);
   for (i=0; i<nentr; i++) {
       iotrack=(AliTPCtrack*)tarray.UncheckedAt(i);
       tracktree->Fill();
   }
   tracktree->Write();
   out->Close();

   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}

Int_t ITSFindClusters(const Char_t *inname, const Char_t *outname) {
   Int_t rc=0;
   const Char_t *name="ITSFindClusters";
   cerr<<'\n'<<name<<"...\n";
   gBenchmark->Start(name);
   TFile *out=TFile::Open(outname,"recreate");
   TFile *in =TFile::Open(inname,"update");

   if (!(gAlice=(AliRun*)in->Get("gAlice"))) {
      cerr<<"Can't get gAlice !\n";
      return 1;
   }
   Int_t ev=0;
   gAlice->GetEvent(ev);
   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   if (!ITS) { cerr<<"Can't get the ITS !\n"; return 1;}

   gAlice->MakeTree("R"); ITS->MakeBranch("R",0);
   //////////////// Taken from ITSHitsToFastPoints.C ///////////////////////
   AliITSsimulationFastPoints *sim = new AliITSsimulationFastPoints();
   for (Int_t i=0;i<3;i++) { ITS->SetSimulationModel(i,sim); }
   Int_t nsignal=25;
   Int_t size=-1;
   Int_t bgr_ev=Int_t(ev/nsignal);
   ITS->HitsToFastRecPoints(ev,bgr_ev,size," ","All"," ");
   //////////////////////////////////////////////////////////////////////////

   gAlice->GetEvent(ev);

   AliITSgeom *geom=ITS->GetITSgeom();
   out->cd();
   geom->Write();

   TClonesArray *clusters=new TClonesArray("AliITSclusterV2",10000);
   TTree *cTree=new TTree("cTree","ITS clusters");
   cTree->Branch("Clusters",&clusters);

   TTree *pTree=gAlice->TreeR();
   if (!pTree) { cerr<<"Can't get TreeR !\n"; return 1; }
   TBranch *branch=pTree->GetBranch("ITSRecPoints");
   if (!branch) { cerr<<"Can't get ITSRecPoints branch !\n"; return 1;}
   TClonesArray *points=new TClonesArray("AliITSRecPoint",10000);
   branch->SetAddress(&points);

   TClonesArray &cl=*clusters;
   Int_t nclusters=0;
   Int_t nentr=(Int_t)pTree->GetEntries();
   for (Int_t i=0; i<nentr; i++) {
       if (!pTree->GetEvent(i)) {cTree->Fill(); continue;}
       Int_t lay,lad,det; geom->GetModuleId(i,lay,lad,det);
       Float_t x,y,zshift; geom->GetTrans(lay,lad,det,x,y,zshift); 
       Double_t rot[9];    geom->GetRotMatrix(lay,lad,det,rot);
       Double_t yshift = x*rot[0] + y*rot[1];
       Int_t ndet=(lad-1)*geom->GetNdetectors(lay) + (det-1);
       Int_t ncl=points->GetEntriesFast();
       nclusters+=ncl;
       for (Int_t j=0; j<ncl; j++) {
          AliITSRecPoint *p=(AliITSRecPoint*)points->UncheckedAt(j);
          Float_t lp[5];
          lp[0]=-p->GetX()-yshift; if (lay==1) lp[0]=-lp[0];
          lp[1]=p->GetZ()+zshift;
          lp[2]=p->GetSigmaX2();
          lp[3]=p->GetSigmaZ2();
          lp[4]=p->GetQ();
          Int_t lab[6]; 
          lab[0]=p->GetLabel(0);lab[1]=p->GetLabel(1);lab[2]=p->GetLabel(2);
          lab[3]=ndet;

          Int_t label=lab[0];
          TParticle *part=(TParticle*)gAlice->Particle(label);
          label=-3;
          while (part->P() < 0.005) {
             Int_t m=part->GetFirstMother();
             if (m<0) {cerr<<"Primary momentum: "<<part->P()<<endl; break;}
             label=m;
             part=(TParticle*)gAlice->Particle(label);
          }
          if      (lab[1]<0) lab[1]=label;
          else if (lab[2]<0) lab[2]=label;
          else cerr<<"No empty labels !\n";

          new(cl[j]) AliITSclusterV2(lab,lp);
       }
       cTree->Fill(); clusters->Delete();
       points->Delete();
   }
   cTree->Write();
   cerr<<"Number of clusters: "<<nclusters<<endl;

   delete cTree; delete clusters; delete points;

   delete gAlice; gAlice=0;
   in->Close();
   out->Close();
   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}

Int_t ITSFindTracks(const Char_t *inname, const Char_t *outname) {
   Int_t rc=0;
   const Char_t *name="ITSFindTracks";
   cerr<<'\n'<<name<<"...\n";
   gBenchmark->Start(name);

   AliITSgeom *geom=(AliITSgeom*)gFile->Get("AliITSgeom");
   if (!geom) { cerr<<"can't get ITS geometry !\n"; return 1;}
   AliITStrackerV2 tracker(geom);

   TFile *out=TFile::Open(outname,"recreate");
   TFile *in =TFile::Open(inname);
   rc=tracker.Clusters2Tracks(in,out);
   in->Close();
   out->Close();

   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}

Int_t ITSPropagateBack(const Char_t *inname, const Char_t *outname) {
   Int_t rc=0;
   const Char_t *name="ITSPropagateBack";
   cerr<<'\n'<<name<<"...\n";
   gBenchmark->Start(name);

   AliITSgeom *geom=(AliITSgeom*)gFile->Get("AliITSgeom");
   if (!geom) { cerr<<"can't get ITS geometry !\n"; return 1;}
   AliITStrackerV2 tracker(geom);

   TFile *out=TFile::Open(outname,"update");
   TFile *in =TFile::Open(inname);
   rc=tracker.PropagateBack(in,out);
   in->Close();
   out->Close();

   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}

Int_t TPCPropagateBack(const Char_t *inname, const Char_t *outname) {
   Int_t rc=0;
   const Char_t *name="TPCPropagateBack";
   cerr<<'\n'<<name<<"...\n";
   gBenchmark->Start(name);

   AliTPCParam *param=(AliTPCParam *)gFile->Get("75x40_100x60");
   if (!param) {cerr<<"TPC parameters have not been found !\n"; return 1;}
   AliTPCtracker *tracker=new AliTPCtracker(param);

   TFile *out=TFile::Open(outname,"update");
   TFile *in =TFile::Open(inname);
   rc=tracker->PropagateBack(in,out);
   delete tracker;
   in->Close();
   out->Close();

   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}

