struct GoodTrack {
  Int_t lab;
  Int_t code;
  Float_t px,py,pz;
  Float_t x,y,z;
  Float_t pxg,pyg,pzg,ptg;
  Bool_t flag;
};
Int_t good_tracks(GoodTrack *gt, Int_t max);


Int_t TPCtracks() {
   Int_t i;
   cerr<<"Doing comparison...\n";

  // Connect the Root Galice file containing Geometry, Kine and Hits


   TFile *cf=TFile::Open("AliTPCclusters.root");
   if (!cf->IsOpen()) {cerr<<"Can't open AliTPCclusters.root !\n"; return 1;}
   AliTPCParam *digp= (AliTPCParam*)cf->Get("75x40_100x60");
   if (!digp) { cerr<<"TPC parameters have not been found !\n"; return 2; }

// Load clusters  
   AliTPCClustersArray *ca=new AliTPCClustersArray;
   ca->Setup(digp);
   ca->SetClusterType("AliTPCcluster");
   ca->ConnectTree("Segment Tree");
   Int_t nentr=Int_t(ca->GetTree()->GetEntries());
   for (Int_t i=0; i<nentr; i++) {
       ca->LoadEntry(i);
   }

// Load tracks
   TFile *tf=TFile::Open("AliTPCtracks.root");
   if (!tf->IsOpen()) {cerr<<"Can't open AliTPCtracks.root !\n"; return 3;}
   TObjArray tarray(2000);
   TTree *tracktree=(TTree*)tf->Get("TreeT");
   TBranch *tbranch=tracktree->GetBranch("tracks");
   nentr=(Int_t)tracktree->GetEntries();
   for (i=0; i<nentr; i++) {
       AliTPCtrack *iotrack=new AliTPCtrack;
       tbranch->SetAddress(&iotrack);
       tracktree->GetEvent(i);
       iotrack->CookLabel(ca);
       tarray.AddLast(iotrack);
   }   

   tf->Close();
   delete ca;
   cf->Close();
   //printf("after cf close\n");
   
   cerr<<"Number of found tracks "<<nentr<<endl;
   tarray.Sort();

   TFile *pfile = new TFile("tpctracks.root","RECREATE"); 

   TTree tracktree1("TreeT","Tree with TPC tracks");
   AliTPCtrack *iotrack=0;   
   tracktree1.Branch("tracks","AliTPCtrack",&iotrack,32000,0);

   for (i=0; i<nentr; i++) {
       iotrack=(AliTPCtrack*)tarray.UncheckedAt(i);
       tracktree1.Fill();
   }   

   tracktree1.Write();
   pfile->Close();
   tarray.Clear();


/////////////////////////////////////////////////////////////////////////
   GoodTrack gt[7000];
   Int_t ngood=0;

   cerr<<"Marking good tracks (this will take a while)...\n";
   ngood=good_tracks(gt,7000);
   ofstream out("good_tracks");
   if (out) {
         for (Int_t ngd=0; ngd<ngood; ngd++)            
	    out<<gt[ngd].lab<<' '<<gt[ngd].code<<' '
	       <<gt[ngd].px <<' '<<gt[ngd].py<<' '<<gt[ngd].pz<<' '
	       <<gt[ngd].x  <<' '<<gt[ngd].y <<' '<<gt[ngd].z <<' '
	       <<gt[ngd].pxg  <<' '<<gt[ngd].pyg <<' '<<gt[ngd].pzg <<' '
	       <<gt[ngd].ptg <<gt[ngd].flag<<endl;
   } else cerr<<"Can not open file (good_tracks) !\n";
   out.close();

   cerr<<"Number of good tracks : "<<ngood<<endl;

   printf("before return in tpctracks\n");
   return 0;
}

//---------------------------------

Int_t good_tracks(GoodTrack *gt, Int_t max) {
   Int_t nt=0;

   //TFile *file=TFile::Open("rfio:galice.root");  // via server
   TFile *file=TFile::Open("galice.root");
   if (!file->IsOpen()) {cerr<<"Can't open galice.root !\n"; exit(4);}

   if (!(gAlice=(AliRun*)file->Get("gAlice"))) {
     cerr<<"gAlice have not been found on galice.root !\n";
     exit(5);
   }

   gAlice->GetEvent(0);   

   AliTPC *TPC=(AliTPC*)gAlice->GetDetector("TPC");
   Int_t ver = TPC->IsVersion(); 
   cerr<<"TPC version "<<ver<<" has been found !\n";

   AliTPCParam *digp=(AliTPCParam*)file->Get("75x40_100x60");
   if (!digp) { cerr<<"TPC parameters have not been found !\n"; exit(6); }
   TPC->SetParam(digp);

   Int_t nrow_up=digp->GetNRowUp();
   Int_t nrows=digp->GetNRowLow()+nrow_up;
   Int_t zero=digp->GetZeroSup();
   Int_t gap=Int_t(0.125*nrows);
   Int_t good_number=Int_t(0.4*nrows);

   TClonesArray *particles=gAlice->Particles(); 
   Int_t np=particles->GetEntriesFast();
   Int_t *good=new Int_t[np];
   for (Int_t ii=0; ii<np; ii++) good[ii]=0;
   
   //MI change to be possible compile macro
   //definition out of the switch statement
   Int_t sectors_by_rows=0;
   TTree *TD=0;
   AliSimDigits da, *digits=&da;
   Int_t *count=0;

   switch (ver) {
   case 1:
     {
      TFile *cf=TFile::Open("AliTPCclusters.root");
      if (!cf->IsOpen()){cerr<<"Can't open AliTPCclusters.root !\n";exit(5);}
      AliTPCClustersArray *ca=new AliTPCClustersArray;
      ca->Setup(digp);
      ca->SetClusterType("AliTPCcluster");
      ca->ConnectTree("Segment Tree");
      Int_t nrows=Int_t(ca->GetTree()->GetEntries());
      for (Int_t n=0; n<nrows; n++) {
          AliSegmentID *s=ca->LoadEntry(n);
          Int_t sec,row;
          digp->AdjustSectorRow(s->GetID(),sec,row);
          AliTPCClustersRow &clrow = *ca->GetRow(sec,row);
          Int_t ncl=clrow.GetArray()->GetEntriesFast();
          while (ncl--) {
              AliTPCcluster *c=(AliTPCcluster*)clrow[ncl];
              Int_t lab=c->GetLabel(0);
              if (lab<0) continue; //noise cluster
              lab=TMath::Abs(lab);
              if (sec>=digp->GetNInnerSector())
              if (row==nrow_up-1    ) good[lab]|=0x1000;
              if (sec>=digp->GetNInnerSector())
              if (row==nrow_up-1-gap) good[lab]|=0x800;
              good[lab]++;
          }
          ca->ClearRow(sec,row);
      }
      cf->Close();
     }
      break;
   case 2:
      TD=(TTree*)gDirectory->Get("TreeD_75x40_100x60");
      TD->GetBranch("Segment")->SetAddress(&digits);
      count = new Int_t[np];
      Int_t i;
      for (i=0; i<np; i++) count[i]=0;
      Int_t sectors_by_rows=(Int_t)TD->GetEntries();
      for (i=0; i<sectors_by_rows; i++) {
          if (!TD->GetEvent(i)) continue;
          Int_t sec,row;
          digp->AdjustSectorRow(digits->GetID(),sec,row);
          cerr<<sec<<' '<<row<<"                                     \r";
          digits->First();
          while (digits->Next()) {
              Int_t it=digits->CurrentRow(), ip=digits->CurrentColumn();
              Short_t dig = digits->GetDigit(it,ip);
              Int_t idx0=digits->GetTrackID(it,ip,0); 
              Int_t idx1=digits->GetTrackID(it,ip,1);
              Int_t idx2=digits->GetTrackID(it,ip,2);
              if (idx0>=0 && dig>=zero) count[idx0]+=1;
              if (idx1>=0 && dig>=zero) count[idx1]+=1;
              if (idx2>=0 && dig>=zero) count[idx2]+=1;
          }
          for (Int_t j=0; j<np; j++) {
              if (count[j]>1) {
                 if (sec>=digp->GetNInnerSector())
                 if (row==nrow_up-1    ) good[j]|=0x1000;
                 if (sec>=digp->GetNInnerSector())
                 if (row==nrow_up-1-gap) good[j]|=0x800;
                 good[j]++;
              }
              count[j]=0;
          }
      }
      delete[] count;
      break;
   default:
      cerr<<"Invalid TPC version !\n";
      file->Close();
      exit(7);
   }

   TTree *TH=gAlice->TreeH();
   //TClonesArray *hits=TPC->Hits();
   Int_t npart=TH->GetEntries();

   while (npart--) {
      AliTPChit *hit0=0;
      TPC->ResetHits();
      TH->GetEvent(npart);
      AliTPChit *hit = (AliTPChit*) TPC->FirstHit(-1);
      
      while(hit) {
         if(hit->fQ==0.) break;
	 hit = (AliTPChit*) TPC->NextHit();
      }
      if(hit) {
         hit0 = new AliTPChit(*hit); // Make copy of hit
	 hit=hit0;
      }
      else continue;
      AliTPChit *hit1=(AliTPChit*) TPC->NextHit();
      if(hit1==0) continue;
      
      if (hit1->fQ != 0.) continue;
     // Int_t i=hit->fTrack;
	   Int_t i=hit->Track();  //modificata in accordo a nuovo AliTPCComparison
      TParticle *p = (TParticle*)particles->UncheckedAt(i);
      if (p->GetFirstMother()>=0) continue;  //secondary particle
      if (good[i] < 0x1000+0x800+2+good_number) continue;
      if (p->Pt()<0.100) continue;
      if (TMath::Abs(p->Pz()/p->Pt())>0.999) continue;

      gt[nt].lab=i;
      gt[nt].code=p->GetPdgCode();
//**** px py pz - in global coordinate system, x y z - in local !
     // gt[nt].px=hit->fX; gt[nt].py=hit->fY; gt[nt].pz=hit->fZ;  //modificato tenendo conto di AliTPCComparison
      gt[nt].px=hit->X(); gt[nt].py=hit->Y(); gt[nt].pz=hit->Z();	  
      Float_t cs,sn; digp->AdjustCosSin(hit1->fSector,cs,sn);
      //gt[nt].x = hit1->fX*cs + hit1->fY*sn;
     // gt[nt].y =-hit1->fX*sn + hit1->fY*cs;  //modificato tenedo conto di AliTPCComaprison
      //gt[nt].z = hit1->fZ;
      gt[nt].x = hit1->X()*cs + hit1->Y()*sn;
      gt[nt].y =-hit1->X()*sn + hit1->Y()*cs;
      gt[nt].z = hit1->Z();		
      gt[nt].pxg = p->Px();
      gt[nt].pyg = p->Py();
      gt[nt].pzg = p->Pz();				
      gt[nt].ptg = p->Pt();		
      gt[nt].flag = 0;		
      nt++;
        
	if(hit0) delete hit0;
      cerr<<i<<"                \r";
      if (nt==max) {cerr<<"Too many good tracks !\n"; break;}
   }
   delete[] good;

   printf("before delete gAlice\n");

   delete gAlice; gAlice=0;

   printf("after delete gAlice\n");
   file->Close();
   printf("after file close\n");
   return nt;
}

//--------------------------------------

