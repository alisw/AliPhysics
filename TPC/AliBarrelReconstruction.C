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
  #include "iostream.h"
  #include "AliRun.h"
  #include "AliRunLoader.h"
  #include "AliLoader.h"
  #include "AliTPCLoader.h"
  #include "AliITSLoader.h"
  #include "AliMagF.h"
  #include "AliTPCtracker.h"
  #include "AliITS.h"
  #include "AliITSgeom.h"
  #include "AliITSRecPoint.h"
  #include "AliITSclusterV2.h"
  #include "AliITSsimulationFastPoints.h"
  #include "AliITStrackerV2.h"
#endif



Int_t TPCFindClusters( Int_t n);
Int_t TPCFindTracks(Int_t n);
Int_t TPCSortTracks(const Char_t *outname,  Int_t n);
Int_t TPCPropagateBack();

Int_t ITSFindClusters(Int_t n);
Int_t ITSFindTracks(const Char_t *inname2, Int_t n);
Int_t ITSPropagateBack();

const char* TPCtrkNameS= "TPC.TracksSorted.root";

class AliRunLoader;
class AliTPCLoader;
class AliTPCParam;
AliRunLoader *rl = 0x0;
AliTPCLoader *tpcl = 0x0;   
AliTPCParam  *param = 0x0;
Bool_t debug = kFALSE;
Int_t AliBarrelReconstruction(Int_t n=1) 
 {
 
   if (gAlice)
    {
      delete gAlice->GetRunLoader();
      delete gAlice;//if everything was OK here it is already NULL
      gAlice = 0x0;
    }
    
   rl = AliRunLoader::Open("galice.root");
   if (rl == 0x0)
    {
      cerr<<"Can not open session"<<endl;
      return 1;
    }
   
   if (rl->LoadgAlice())
    {
      cerr<<"Error occured while l"<<endl;
      return 1;
    }
   AliKalmanTrack::SetConvConst(1000/0.299792458/rl->GetAliRun()->Field()->SolenoidField());
   
   tpcl = (AliTPCLoader*)rl->GetLoader("TPCLoader");
   if (tpcl == 0x0)
    {
      cerr<<"Can not get TPC Loader"<<endl;
      return 1;
    }

// ********** Find TPC clusters *********** //
   if ( TPCFindClusters(n) ) 
    {
      cerr<<"Failed to get TPC clusters !\n";
      return 1;
    }
   
// ********** Find TPC tracks *********** //
    if (TPCFindTracks(n)) 
     {
      cerr<<"Failed to get TPC tracks !\n";
      return 1;
     }

//   cout<<"Stopping tracking on TPC\n";
// ********** Sort and label TPC tracks *********** //
   if (TPCSortTracks(TPCtrkNameS,n)) {
      cerr<<"Failed to sort TPC tracks !\n";
      return 1;
    } 

//   cout<<"Stopping tracking on TPC sorted T\n";

   
// ********** Find ITS clusters *********** //
   if (ITSFindClusters(n)) 
    {
      cerr<<"Failed to get ITS clusters !\n";
      return 1;
    }
   
// ********** Find ITS tracks *********** //

   if (ITSFindTracks(TPCtrkNameS,n)) 
    {
      cerr<<"Failed to get ITS tracks !\n";
      return 1;
    }
   //clsFile->Close();
   cout<<"Stopping on ITSFindTracks\n";
   delete rl;
   rl = 0x0;
   return 0;
   return 1;

// ********** Back propagation of the ITS tracks *********** //
   if ( ITSPropagateBack() ) {
      cerr<<"Failed to propagate back the ITS tracks !\n";
      return 1;
   }


// ********** Back propagation of the TPC tracks *********** //
   
   if (TPCPropagateBack()) {
      cerr<<"Failed to propagate back the TPC tracks !\n";
      return 1;
   }

   return 0;
}


Int_t TPCFindClusters(Int_t n) 
 {
   Int_t rc=0;
   const Char_t *name="TPCFindClusters";
   cerr<<'\n'<<name<<"...\n";
   rl->CdGAFile();
   param=(AliTPCParam *)gDirectory->Get("75x40_100x60");
   if (!param) 
    {
     param=(AliTPCParam *)gDirectory->Get("75x40_100x60_150x60");
     if (!param) 
      {
        cerr<<"TPC parameters have not been found !\n";
        return 1;
      }
    }
//   param->Dump();
   param->Update();

   gBenchmark->Start(name);

   AliTPCv2 tpc;
   tpc.SetParam(param);
   tpc.SetLoader(tpcl);
 
   tpcl->LoadDigits("read");
   tpcl->LoadRecPoints("recreate");
  

   //tpc.Digits2Clusters(out); //MI change
   for (Int_t i=0;i<n;i++)
    {
      printf("Processing event %d\n",i);
      tpc.Digits2Clusters(i);
     //     AliTPCclusterer::Digits2Clusters(dig, out, i);
    }
   
   tpcl->UnloadDigits();
   tpcl->UnloadRecPoints();
   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}

Int_t TPCFindTracks(Int_t n) {

   Int_t rc=0;
   AliTPCtracker *tracker = 0x0;
   const Char_t *name="TPCFindTracks";
   rl->GetEvent(0);
   rl->CdGAFile();
   param=(AliTPCParam *)gDirectory->Get("75x40_100x60");
   
   if (!param) 
    {
     param=(AliTPCParam *)gDirectory->Get("75x40_100x60_150x60");
     if (!param) 
      {
        cerr<<"TPC parameters have not been found !\n";
        return 1;
      }
    }
//   param->Dump();
   param->Update();

   gBenchmark->Start(name);

   tpcl->LoadRecPoints("read");
   tpcl->LoadTracks("recreate");

   //AliTPCtracker *tracker=new AliTPCtracker(param);
   //rc=tracker->Clusters2Tracks(0,out);
   //delete tracker;
   
   for (Int_t i=0;i<n;i++){
     rl->GetEvent(i);
     printf("Processing event %d\n",i);
     
     tracker = new AliTPCtracker(param, i, (AliConfig::GetDefaultEventFolderName()).Data());
     //Int_t rc=
     tracker->Clusters2Tracks();
     delete tracker;
   }
   tpcl->UnloadRecPoints();
   tpcl->UnloadTracks();

   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}


Int_t TPCSortTracks(const Char_t * outname,  Int_t eventn){
   Int_t rc=0;
   const Char_t *name="TPCSortTracks";

   cerr<<'\n'<<name<<"...\n";
   rl->CdGAFile();
   param=(AliTPCParam *)gDirectory->Get("75x40_100x60");
   if (!param) 
    {
     param=(AliTPCParam *)gDirectory->Get("75x40_100x60_150x60");
     if (!param) 
      {
        cerr<<"TPC parameters have not been found !\n";
        return 1;
      }
    }
   param->Update();

   gBenchmark->Start(name);
   
 
   AliRunLoader* rl2 = AliRunLoader::Open("galice.root","tmp");
   
   AliLoader* tpcl2 = (AliTPCLoader*)rl2->GetLoader("TPCLoader");
   cout<<"tpcl2->SetTracksFileName("<<outname<<");\n";
   tpcl2->SetTracksFileName(TString(outname));
   tpcl2->LoadTracks("recreate");
   
   // loop over events 
   for (Int_t event=0;event<eventn; event++)
    {
     
     rl->GetEvent(event);
     rl2->GetEvent(event);

     TObjArray tarray(10000);
     AliTPCtrack *iotrack=0;
     Int_t i;
     
     if (tpcl->TreeT() == 0x0)   tpcl->LoadTracks("read");
     TTree *tracktree=tpcl->TreeT();
     if (tracktree == 0x0)
      {
        cerr<<"Can not get TreeT for event "<<event<<endl;
        continue;
      }
     TBranch *tbranch=tracktree->GetBranch("tracks");
     Int_t nentr=(Int_t)tracktree->GetEntries();
     cout<<"IN Tracks nentr = "<<nentr<<endl;
     for (i=0; i<nentr; i++) 
      {
       iotrack=new AliTPCtrack;
       tbranch->SetAddress(&iotrack);
       tracktree->GetEvent(i);
       tarray.AddLast(iotrack);
      }   
     tarray.Sort();
     
     //assign thacks GEANT labels
     cout<<"Running Tracker\n";
     AliTPCtracker *tracker = new AliTPCtracker(param, event);

     cout<<"Load Sectors\n";

     tracker->LoadInnerSectors();
     tracker->LoadOuterSectors();

     cout<<"Cooking Labels\n";

     for (i=0; i<nentr; i++) 
      {
       iotrack=(AliTPCtrack*)tarray.UncheckedAt(i);
       tracker->CookLabel(iotrack,0.1);
      }   
     cout<<"deleting tracker\n";
     delete tracker;
     //in->Close();
     //end of GEANT label assignment
     

     if (tpcl2->TreeT() == 0x0)   tpcl2->MakeTree("T");
     tracktree= tpcl2->TreeT();
     if (tracktree == 0x0)
      {
        cerr<<"Can not get TreeT for Sorted tracks\n";
        return 1;
      }
     tracktree->Branch("tracks","AliTPCtrack",&iotrack,32000,0);
     for (i=0; i<nentr; i++) {
       iotrack=(AliTPCtrack*)tarray.UncheckedAt(i);
       tracktree->Fill();
     }

     cout<<" Tracks File Name is "<<tpcl2->GetTracksFileName()<<endl;
     tpcl2->WriteTracks("OVERWRITE");
   }
  

   delete rl2;

   tpcl->UnloadTracks();
   tpcl->UnloadRecPoints(); 
   
   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}

Int_t ITSFindClusters(Int_t n) 
 {
   Int_t rc=0;
   const Char_t *name="ITSFindClusters";
   cerr<<'\n'<<name<<"...\n";
   gBenchmark->Start(name);

   if (rl->GetEvent(0))
    {
      cerr<<"Problems\n";
      return 1;
    }
    
   AliITSLoader* itsl = (AliITSLoader*)rl->GetLoader("ITSLoader");
   if (itsl == 0x0)
    {
      cerr<<"Can not get ITS Loader from Run Loader\n";
      return 1;
    } 
  
   gAlice=rl->GetAliRun();
   if (!gAlice) {
      cerr<<"Can't get gAlice !\n";
      return 1;
   }
   if (rl->TreeK() == 0x0) rl->LoadKinematics();
   if (rl->TreeE() == 0x0) rl->LoadHeader();
   
   itsl->LoadRawClusters("recreate");
   
   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   if (!ITS) { cerr<<"Can't get the ITS !\n"; return 1;}
   AliITSgeom *geom=ITS->GetITSgeom();

   Int_t ev=0;
   for (ev = 0; ev<n; ev++)
   {

     rl->GetEvent(ev);
     TBranch *branch = 0x0;
     TTree *pTree = 0x0;
     
     pTree=itsl->TreeR();
     if (pTree == 0x0) 
      {
        itsl->LoadRecPoints("read");
        pTree=itsl->TreeR();
        if (pTree == 0x0) 
         {
           ::Error("AliBarrelReonstruction.C::ITSFindClusters",
                   "Can not get TreeR for event %d",ev);
           return 1;
         }
      }
     
     pTree=itsl->TreeR();
     branch=(pTree)?pTree->GetBranch("ITSRecPoints"):0x0;

     if (branch== 0x0) {
       //if not reconstructed ITS branch do reconstruction 
       if (debug) ::Info("AliBarrelReconstruction.C","Did not get ITSRecPoints from TreeR.");
       if (debug) ::Info("AliBarrelReconstruction.C","Making branch and Producing RecPoints");
       itsl->SetRecPointsFileOption("recreate");
       if (itsl->TreeR()==0x0) itsl->MakeTree("R");
       
       ITS->MakeBranch("RF");
       ITS->SetTreeAddress();
       itsl->LoadHits();

       //////////////// Taken from ITSHitsToFastPoints.C ///////////////////////
       for (Int_t i=0;i<3;i++) { 
         ITS->SetSimulationModel(i,new AliITSsimulationFastPoints()); 
       }
       Int_t nsignal=25;
       Int_t size=-1;
       Int_t bgr_ev=Int_t(ev/nsignal);
       
       ITS->HitsToFastRecPoints(ev,bgr_ev,size," ","All"," ");
       //////////////////////////////////////////////////////////////////////////
       //MI comment  - in HitsToFast... they reset treeR to 0 
       //they overwrite full reconstructed event ???? ... so lets connect TreeR one more
       //time
       itsl->UnloadHits();
       ITS->SetTreeAddress();
     }

     TClonesArray *clusters=new TClonesArray("AliITSclusterV2",10000);

     if (itsl->TreeC() == 0x0) itsl->MakeTree("C");
     TTree *cTree = itsl->TreeC();
     TBranch* b = cTree->GetBranch("Clusters");
     if (b == 0x0) 
       cTree->Branch("Clusters",&clusters);
     else b->SetAddress(&clusters);
     
      
     if (itsl->TreeR() == 0x0) itsl->LoadRecPoints();

     if (branch == 0x0)
      {//it means that we produced FastRecPoints above
       pTree=itsl->TreeR();
       if (!pTree) { cerr<<"Can't get TreeR !\n"; return 1; }
       branch=pTree->GetBranch("ITSRecPointsF");
       if (!branch) { cerr<<"Can't get Fast RecPoints branch named ITSRecPointsF  !\n"; return 1;}
      } 
     TClonesArray *points=new TClonesArray("AliITSRecPoint",10000);
     branch->SetAddress(&points);
     
     TClonesArray &cl=*clusters;
     Int_t nclusters=0;
     Int_t nentr=(Int_t)branch->GetEntries();
     AliITSgeom *geom=ITS->GetITSgeom();
     
     cout<<"\n\n Number of entries in TreeR = "<<nentr<<endl;
     for (Int_t i=0; i<nentr; i++) 
      {
       if (!branch->GetEvent(i)) 
        {
          cTree->Fill(); 
          continue;
        }
       
       Int_t lay,lad,det; geom->GetModuleId(i,lay,lad,det);
       Float_t x,y,zshift; geom->GetTrans(lay,lad,det,x,y,zshift); 
       Double_t rot[9];    geom->GetRotMatrix(lay,lad,det,rot);
       Double_t yshift = x*rot[0] + y*rot[1];
       Int_t ndet=(lad-1)*geom->GetNdetectors(lay) + (det-1);
       Int_t ncl=points->GetEntriesFast();
       nclusters+=ncl;

       if (debug) ::Info("AliBarrelReconstruction.C",
              "i=%d lay=%d lad=%d det=%d NRP=%d",
               i,   lay,   lad,   det,   ncl);

       for (Int_t j=0; j<ncl; j++) 
        { 
          Float_t lp[5];
          Int_t lab[4]; 

          AliITSRecPoint *p=(AliITSRecPoint*)points->UncheckedAt(j);
          lp[0]=-p->GetX()-yshift; if (lay==1) lp[0]=-lp[0];
          lp[1]=p->GetZ()+zshift;
          lp[2]=p->GetSigmaX2();
          lp[3]=p->GetSigmaZ2();
          lp[4]=p->GetQ();
          lab[0]=p->GetLabel(0);
          lab[1]=p->GetLabel(1);
          lab[2]=p->GetLabel(2);
          lab[3]=ndet;
     
          Int_t label=lab[0];
          if (label<0) continue;
          TParticle *part=(TParticle*)rl->Stack()->Particle(label);
          label=-3;
          while (part->P() < 0.005) {
            Int_t m=part->GetFirstMother();
            if (m<0) {cerr<<"Primary momentum: "<<part->P()<<endl; break;}
            label=m;
            part=(TParticle*)gAlice->Particle(label);
          }
          if (lab[1]<0) lab[1]=label;
          else if (lab[2]<0) lab[2]=label;
               else cerr<<"No empty labels !\n";
     
          new(cl[j]) AliITSclusterV2(lab,lp);
         }
       cTree->Fill(); clusters->Delete();
       points->Delete();
     }
     itsl->WriteRawClusters("OVERWRITE");
     
     cerr<<"Number of clusters: "<<nclusters<<endl;
     delete clusters; delete points;

   }
   itsl->UnloadRecPoints();
   itsl->UnloadRawClusters();
   rl->UnloadKinematics();

   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}

Int_t ITSFindTracks(const Char_t *inname2, Int_t n) 
{

   Int_t rc=0;
   const Char_t *name="ITSFindTracks";
   cerr<<'\n'<<name<<"...\n";
   gBenchmark->Start(name);

   rl->GetEvent(0);

   AliITSLoader* itsl = (AliITSLoader*)rl->GetLoader("ITSLoader");
   if ( (itsl == 0x0) || (tpcl == 0x0))
    {
      cerr<<"Can not get ITS Loader from Run Loader\n";
      return 1;
    } 
   
   tpcl->SetTracksFileName(TString(inname2));
   if(rl->GetAliRun() == 0x0) rl->LoadgAlice();
   
   AliITS* ITS = rl->GetAliRun()->GetModule("ITS");
   
   rl->CdGAFile();
   AliITSgeom *geom=ITS->GetITSgeom();
   if (!geom) { cerr<<"can't get ITS geometry !\n"; return 1;}

   itsl->LoadTracks("recreate");
   itsl->LoadRawClusters("read");
   for (Int_t i=0;i<n;i++)
    {
      rl->GetEvent(i);
      AliITStrackerV2 tracker(geom,i,AliConfig::GetDefaultEventFolderName());
      rc=tracker.Clusters2Tracks();
    }


   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}

Int_t ITSPropagateBack() 
{
   
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

Int_t TPCPropagateBack() 
{
   Int_t rc=0;
   const Char_t *name="TPCPropagateBack";
   cerr<<'\n'<<name<<"...\n";
   gBenchmark->Start(name);

   param=(AliTPCParam *)gDirectory->Get("75x40_100x60");
   if (!param) 
    {
     param=(AliTPCParam *)gDirectory->Get("75x40_100x60_150x60");
     if (!param) 
      {
        cerr<<"TPC parameters have not been found !\n";
        return 1;
      }
    }
//   param->Dump();
   param->Update();

   AliTPCtracker *tracker = new AliTPCtracker(param, AliConfig::GetDefaultEventFolderName());

//   TFile *out=TFile::Open(outname,"update");
//   TFile *in =TFile::Open(inname);

   rc=tracker->PropagateBack();
   delete tracker;

   gBenchmark->Stop(name);
   gBenchmark->Show(name);

   return rc;
}
