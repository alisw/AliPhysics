/****************************************************************************
 * This macro is supposed to do reconstruction in the barrel ALICE trackers *
 * (Make sure you have TPC digits and ITS hits before using this macro !!!) *
 * It does the following steps (July 9, 2002,Dubna)                         *
 *                   1) TPC cluster finding                                 *
 *                   2) TPC track finding                                   *
 *                   3) ITS cluster finding V2                              *
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
  #include<iostream.h>
#endif

Int_t TPCFindClusters(const Char_t *inname, const Char_t *outname, Int_t n);
Int_t TPCFindTracks(const Char_t *inname, const Char_t *outname, Int_t n);
Int_t TPCSortTracks(const Char_t *inname, const Char_t *inname2, const Char_t *outname,  Int_t n);
Int_t TPCPropagateBack(const Char_t *inname, const Char_t *outname);

Int_t ITSFindClusters(const Char_t *inname,  const Char_t *outname, Int_t n);
Int_t ITSFindTracks(const Char_t *inname, const Char_t *inname2, const Char_t *outname, Int_t n);
Int_t ITSPropagateBack(const Char_t *inname, const Char_t *outname);


Int_t AliBarrelReconstructionV2(Int_t n=1) {
  //const Char_t *TPCdigName="rfio:galice.root";
   const Char_t *TPCdigName="galice.root";
   const Char_t *TPCclsName="AliTPCclusters.root";
   const Char_t *TPCtrkName="AliTPCtracks.root";
   const Char_t *TPCtrkNameS="AliTPCtracksSorted.root";


   //const Char_t *ITSdigName="rfio:galice.root";
   const Char_t *ITSdigName="galice.root";
   const Char_t *ITSclsName="AliITSclustersV2.root";
   const Char_t *ITStrkName="AliITStracksV2.root";

   AliKalmanTrack::SetConvConst(100/0.299792458/0.2/gAlice->Field()->Factor());

    
// ********** Find TPC clusters *********** //
   if (TPCFindClusters(TPCdigName,TPCclsName, n)) {
      cerr<<"Failed to get TPC clusters !\n";
      return 1;
   }      

// ********** Find TPC tracks *********** //
    if (TPCFindTracks(TPCclsName,TPCtrkName,n)) {
      cerr<<"Failed to get TPC tracks !\n";
      return 1;
    }
  
// ********** Sort and label TPC tracks *********** //
   if (TPCSortTracks(TPCclsName,TPCtrkName,TPCtrkNameS,n)) {
      cerr<<"Failed to sort TPC tracks !\n";
      return 1;
    }

// ********** Find ITS clusters *********** //
   if (ITSFindClusters(ITSdigName,ITSclsName,n)) {
      cerr<<"Failed to get ITS clusters !\n";
      return 1;
   }

// ********** Find ITS tracks *********** //
   {
     //TFile *clsFile=TFile::Open(ITSclsName);
   if (ITSFindTracks(TPCtrkNameS,ITSclsName,ITStrkName,n)) {
      cerr<<"Failed to get ITS tracks !\n";
      return 1;
   }
   //clsFile->Close();
   }
   return 1;
// ********** Back propagation of the ITS tracks *********** //
   {TFile *clsFile=TFile::Open(ITSclsName);
   if (ITSPropagateBack(ITStrkName,TPCtrkNameS)) {
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


Int_t TPCFindClusters(const Char_t *inname, const Char_t *outname, Int_t n) {
   Int_t rc=0;
   const Char_t *name="TPCFindClusters";
   cerr<<'\n'<<name<<"...\n";
   gBenchmark->Start(name);

   TFile *out=TFile::Open(outname,"recreate");
   TFile *in =TFile::Open(inname);

   AliTPCParam *param=(AliTPCParam *)in->Get("75x40_100x60_150x60");
   if (!param) {cerr<<"TPC parameters have not been found !\n"; return 1;}
   AliTPCv2 tpc;
   tpc.SetParam(param);

   //tpc.Digits2Clusters(out); //MI change
   cout<<"TPCFindClusters: nev ="<<n<<endl;
   for (Int_t i=0;i<n;i++){
     printf("Processing event %d\n",i);
     tpc.Digits2Clusters(out,i);
     //	 AliTPCclusterer::Digits2Clusters(dig, out, i);
   }
   in->Close();
   out->Close();
   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}

Int_t TPCFindTracks(const Char_t *inname, const Char_t *outname, Int_t n) {
   Int_t rc=0;
   const Char_t *name="TPCFindTracks";
   cerr<<'\n'<<name<<"...\n";
   gBenchmark->Start(name);
   TFile *out=TFile::Open(outname,"recreate");
   TFile *in =TFile::Open(inname);
   AliTPCParam *param=(AliTPCParam *)in->Get("75x40_100x60_150x60");
   if (!param) {cerr<<"TPC parameters have not been found !\n"; return 1;}

   //AliTPCtracker *tracker=new AliTPCtracker(param);
   //rc=tracker->Clusters2Tracks(0,out);
   //delete tracker;
   cout<<"TPCFindTracks: nev ="<<n<<endl;

   for (Int_t i=0;i<n;i++){
     printf("Processing event %d\n",i);
     AliTPCtracker *tracker = new AliTPCtracker(param,i);
     //Int_t rc=
       tracker->Clusters2Tracks(0,out);
     delete tracker;
   }

   in->Close();
   out->Close();
   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}
Int_t TPCSortTracks(const Char_t *inname, const Char_t *inname2, const Char_t * outname,  Int_t eventn){
   Int_t rc=0;
   const Char_t *name="TPCSortTracks";
   cerr<<'\n'<<name<<"...\n";
   gBenchmark->Start(name);

   TFile *out =TFile::Open(outname,"recreate");
   TFile *in1 =TFile::Open(inname);
   TFile *in2 =TFile::Open(inname2);


   AliTPCParam *param=(AliTPCParam *)in1->Get("75x40_100x60_150x60");
   if (!param) {cerr<<"TPC parameters have not been found !\n"; return 1;}

   cout<<"TPCSortedtracks: nev ="<<eventn<<endl;


   // loop over events 
   for (Int_t event=0;event<eventn; event++){
     TObjArray tarray(10000);
     AliTPCtrack *iotrack=0;
     Int_t i;


     in2->cd();
     char   tname[100];
     sprintf(tname,"TreeT_TPC_%d",event);

     TTree *tracktree=(TTree*)in2->Get(tname);
     TBranch *tbranch=tracktree->GetBranch("tracks");
     Int_t nentr=(Int_t)tracktree->GetEntries();
     for (i=0; i<nentr; i++) {
       iotrack=new AliTPCtrack;
       tbranch->SetAddress(&iotrack);
       tracktree->GetEvent(i);
       tarray.AddLast(iotrack);
     }   
     tarray.Sort();
     // in2->Close();
     
     //assign thacks GEANT labels
     in1->cd();
     AliTPCtracker *tracker = new AliTPCtracker(param,event);
     tracker->LoadInnerSectors();
     tracker->LoadOuterSectors();
     for (i=0; i<nentr; i++) {
       iotrack=(AliTPCtrack*)tarray.UncheckedAt(i);
       tracker->CookLabel(iotrack,0.1);
     }   
     delete tracker;
     //in->Close();
     //end of GEANT label assignment
     
     tracktree=new TTree(tname,"Tree with TPC tracks");
     tracktree->Branch("tracks","AliTPCtrack",&iotrack,32000,0);
     for (i=0; i<nentr; i++) {
       iotrack=(AliTPCtrack*)tarray.UncheckedAt(i);
       tracktree->Fill();
     }
     out->cd();
     tracktree->Write();
   }

   out->Close();
   in2->Close();
   in1->Close();

   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}

Int_t ITSFindClusters(const Char_t *inname, const Char_t *outname, Int_t n) {
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

   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   if (!ITS) { cerr<<"Can't get the ITS !\n"; return 1;}
   AliITSgeom *geom=ITS->GetITSgeom();
   out->cd();   
   geom->Write();
     
   cout<<"ITSFindClusters: nev ="<<n<<endl;
   Int_t ev=0;
   for (ev = 0; ev<n; ev++){
     in->cd();   // !!!! MI directory must point to galice. - othervise problem with Tree -connection
     gAlice->GetEvent(ev);
     //gAlice->TreeR()->Reset();   //reset reconstructed tree

     
     TTree *pTree=gAlice->TreeR();
     if (!pTree){
       gAlice->MakeTree("R");
       pTree = gAlice->TreeR();
     }
     TBranch *branch=pTree->GetBranch("ITSRecPoints");
     /*
     if (!branch) {
       //if not reconstructed ITS branch do reconstruction 
       ITS->MakeBranch("R",0);

       //////////////// Taken from ITSHitsToFastPoints.C ///////////////////////
       AliITSsimulationFastPoints *sim = new AliITSsimulationFastPoints();
       for (Int_t i=0;i<3;i++) { ITS->SetSimulationModel(i,sim); }
       Int_t nsignal=25;
       Int_t size=-1;
       Int_t bgr_ev=Int_t(ev/nsignal);
       ITS->HitsToFastRecPoints(ev,bgr_ev,size," ","All"," ");
       //////////////////////////////////////////////////////////////////////////
       gAlice->GetEvent(ev);   //MI comment  - in HitsToFast... they reset treeR to 0 
       //they overwrite full reconstructed event ???? ... so lets connect TreeR one more
       //time
     }
     */

     
     out->cd();
     TClonesArray *clusters=new TClonesArray("AliITSclusterV2",10000);
     char   cname[100];
     sprintf(cname,"TreeC_ITS_%d",ev);
  
     TTree *cTree=new TTree(cname,"ITS clusters");
     cTree->Branch("Clusters",&clusters);
     
     pTree=gAlice->TreeR();
     if (!pTree) { cerr<<"Can't get TreeR !\n"; return 1; }
     branch=pTree->GetBranch("ITSRecPoints");
     if (!branch) { cerr<<"Can't get ITSRecPoints branch !\n"; return 1;}
     TClonesArray *points=new TClonesArray("AliITSRecPoint",10000);
     branch->SetAddress(&points);
     
     TClonesArray &cl=*clusters;
     Int_t nclusters=0;
     Int_t nentr=(Int_t)pTree->GetEntries();
     Int_t lab[6]; 
     Float_t lp[5];
     AliITSgeom *geom=ITS->GetITSgeom();

     cout<<"!!!! nentr ="<<nentr<<endl;
     for (Int_t i=0; i<nentr; i++) {
       if (!pTree->GetEvent(i)) {cTree->Fill(); continue;}
       Int_t lay,lad,det; geom->GetModuleId(i,lay,lad,det);
       Float_t x,y,zshift; geom->GetTrans(lay,lad,det,x,y,zshift); 
       Double_t rot[9];    geom->GetRotMatrix(lay,lad,det,rot);
       Double_t yshift = x*rot[0] + y*rot[1];
       Int_t ndet=(lad-1)*geom->GetNdetectors(lay) + (det-1);
       Int_t ncl=points->GetEntriesFast();
       nclusters+=ncl;

     Float_t kmip=1.;
     if(lay<5&&lay>2){kmip=280.;}; // b.b.
     if(lay<7&&lay>4){kmip=38.;};


       for (Int_t j=0; j<ncl; j++) {
	 AliITSRecPoint *p=(AliITSRecPoint*)points->UncheckedAt(j);
	 lp[0]=-p->GetX()-yshift; if (lay==1) lp[0]=-lp[0];
	 lp[1]=p->GetZ()+zshift;
	 lp[2]=p->GetSigmaX2();
	 lp[3]=p->GetSigmaZ2();
	 lp[4]=p->GetQ(); lp[4]/=kmip;    // b.b.
	 if(ncl==11)cout<<"mod,cl,Q ="<<i<<","<<j<<","<<lp[4]<<endl;
	 lab[0]=p->GetLabel(0);lab[1]=p->GetLabel(1);lab[2]=p->GetLabel(2);
	 lab[3]=ndet;
	 
	 Int_t label=lab[0];
        if(label>=0) { // b.b.
	 TParticle *part=(TParticle*)gAlice->Particle(label);
	 label=-3;
	 while (part->P() < 0.005) {
	   Int_t m=part->GetFirstMother();
	   if (m<0) {cerr<<"Primary momentum: "<<part->P()<<endl; break;}
	   label=m;
	   part=(TParticle*)gAlice->Particle(label);
	 }
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

   }

   delete gAlice; gAlice=0;
   in->Close();
   out->Close();
   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}

Int_t ITSFindTracks(const Char_t * inname, const Char_t *inname2, const Char_t *outname, Int_t n) {
   Int_t rc=0;
   const Char_t *name="ITSFindTracks";
   cerr<<'\n'<<name<<"...\n";
   gBenchmark->Start(name);


   TFile *out=TFile::Open(outname,"recreate");
   TFile *in =TFile::Open(inname);
   TFile *in2 =TFile::Open(inname2);

   AliITSgeom *geom=(AliITSgeom*)gFile->Get("AliITSgeom");
   if (!geom) { cerr<<"can't get ITS geometry !\n"; return 1;}

   cout<<"ITSFindtracks: nev ="<<n<<endl;

   for (Int_t i=0;i<n;i++){
     AliITStrackerV2 tracker(geom,i);
     rc=tracker.Clusters2Tracks(in,out);
   }

   in->Close();
   in2->Close();
   out->Close();

   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}

Int_t ITSPropagateBack(const Char_t *inname, const Char_t *outname) {
   
  Int_t rc=0;
  /*
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
  */
   return rc;
}

Int_t TPCPropagateBack(const Char_t *inname, const Char_t *outname) {
   Int_t rc=0;
   const Char_t *name="TPCPropagateBack";
   cerr<<'\n'<<name<<"...\n";
   gBenchmark->Start(name);

   AliTPCParam *param=(AliTPCParam *)gFile->Get("75x40_100x60_150x60");
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





