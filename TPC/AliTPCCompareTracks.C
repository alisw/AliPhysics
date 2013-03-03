#ifndef __CINT__
#include "alles.h"
#include "AliComplexCluster.h"
//#include "AliTPCclusterM.h"
#include "AliTPCclusterMI.h"
#endif


Int_t AliTPCCompareTracks(Int_t eventn, Bool_t all = kFALSE) {
   cerr<<"Comparing tracks...\n";
   //CONNECT FILES
   TFile *file=TFile::Open("galice.root");
   if (!file->IsOpen()) {cerr<<"Can't open galice.root !\n"; return 1;}
   //
   TFile *ftracks=TFile::Open("AliTPCtracks.root","update");
   if (!ftracks->IsOpen()){cerr<<"Can't open AliTPCtracks.root !\n"; return 3;}
   //
   AliTPCParam *param=(AliTPCParam *)file->Get("75x40_100x60_150x60");
   if (!param) {cerr<<"TPC parameters have not been found !\n"; return 2;}   
   //
   TFile *fout=TFile::Open("AliTPCTracksDif.root","new"); 
   if (!fout->IsOpen()) fout=TFile::Open("AliTPCClustersDif.root","recreate");  
   //
   //connect exact clusters 
   file->cd();
   AliTPCClustersArray *ce=new AliTPCClustersArray;
   ce->Setup(param);
   ce->SetClusterType("AliComplexCluster");
   char  cname[100];
   sprintf(cname,"TreeCExact_%d",eventn);
   ce->ConnectTree(cname);
   //
   //connect reconstructed tracks
   ftracks->cd();
   sprintf(cname,"Seeds");
   TTree * treetracks = (TTree*)ftracks->Get(cname);
   TBranch * branchtracks = treetracks->GetBranch("seeds");   
   //
   //load seeds to the memory
   Int_t trackmap[500000][4];  // map of tracks corresponding to given track number
   memset(trackmap,0,sizeof(Int_t)*4*500000);
   Int_t ntracks = treetracks->GetEntries();
   TObjArray * trackarr= new TObjArray(ntracks);
   Int_t nproces = TMath::Min(ntracks,4000);

   //   for (Int_t itrack =0; itrack<ntracks; itrack++){
   for (Int_t itrack =0; itrack<nproces; itrack++){
     AliTPCseed * seed = new AliTPCseed;  
     //
     seed->fPoints = new TClonesArray("AliTPCTrackPoint",200);
     branchtracks->SetAddress(&seed); 
     branchtracks->GetEntry(itrack);

     //if (seed->fRemoval>0 && (itrack%4) ) continue;
     trackarr->AddLast(seed);    

     //crete array with exact position information
     seed->fEPoints = new TClonesArray("AliTPCTrackPointRef",1);
     seed->fEPoints->ExpandCreateFast(200);

     Int_t label = TMath::Abs(seed->GetLabel());
     Int_t i;
     if (label>500000) {
       //too big track label
     }else{
       for (i=0;i<4;i++) {
	 if ( trackmap[label][i]==0) 
	   break;
       }
       if(i<4) trackmap[label][i]=itrack;
     }
   }

   //add information about exact positions   
   Int_t nrows=Int_t(ce->GetTree()->GetEntries());

   for (Int_t n=0; n<nrows; n++) {
       AliSegmentID *se=ce->LoadEntry(n);
       Int_t sec,row;
       param->AdjustSectorRow(se->GetID(),sec,row);
       //
       AliTPCClustersRow* clrowe = ce->GetRow(sec,row);
       //
       Int_t ncl=clrowe->GetArray()->GetEntriesFast();
       const Int_t kNIS=param->GetNInnerSector(), kNOS=param->GetNOuterSector();
       Int_t npads, sign;
       if (sec < kNIS) {
          npads = param->GetNPadsLow(row);
          sign = (sec < kNIS/2) ? 1 : -1;
       } else {
          npads = param->GetNPadsUp(row);
          sign = ((sec-kNIS) < kNOS/2) ? 1 : -1;
       }
       
       Int_t trrow = row;
       if (sec>=param->GetNInnerSector()) trrow += param->GetNRowLow(); 
       
       while (ncl--) {
	   AliComplexCluster *cl=(AliComplexCluster*)(clrowe->GetArray()->At(ncl));
           Double_t x=param->GetPadRowRadii(sec,row);
	   Double_t y=cl->fY;
	   y = ( y + 0.5 -  0.5*npads) *param->GetPadPitchWidth(sec);
	   Double_t z=cl->fX*param->GetZWidth();
	   z = sign*(param->GetZLength() - z);
           Float_t cs, sn, tmp;
           param->AdjustCosSin(sec,cs,sn);
	   for (Int_t i=0;i<3;i++){
	     if (cl->fTracks[0]<500000) if (trackmap[cl->fTracks[0]][i]) {
	       AliTPCseed * seed = (AliTPCseed*)trackarr->At(trackmap[cl->fTracks[0]][i]);
	       TClonesArray * clarr = seed->fPoints;
	       if (!clarr){
		 //printf("Seed %d without info",trackmap[cl->fTracks[0]][i]);
		 continue;
	       }
	       AliTPCTrackPoint    * trcluster = (AliTPCTrackPoint*)(seed->fPoints->At(trrow));
	       AliTPCTrackPointRef * ecluster  = (AliTPCTrackPointRef*)(seed->fEPoints->At(trrow));
	       //
	       ecluster->GetCPoint() = trcluster->GetCPoint();
	       ecluster->GetTPoint() = trcluster->GetTPoint();
	       //
	       AliTPCExactPoint & epoint =  ecluster->GetExactPoint();
	       /*
		 trcluster->fEZ = z;
		 trcluster->fEY = y;
		 trcluster->fAngleY = cl->fSigmaY2*param->GetPadPitchLength(sec);
		 trcluster->fAngleZ = cl->fSigmaX2*param->GetPadPitchLength(sec);
		 trcluster->fEAmp = cl->fQ;
		 trcluster->fEPrim = cl->fMax;
	       */
	       epoint.fEZ = z;
	       epoint.fEY = y;
	       epoint.fEAngleY = cl->fSigmaY2*param->GetPadPitchLength(sec);
	       epoint.fEAngleZ = cl->fSigmaX2*param->GetPadPitchLength(sec);
	       epoint.fEAmp = cl->fQ;
	       epoint.fEPrim = cl->fMax;
	     }
	   } 
       }
       //       cr->ClearRow(sec,row);
       ce->ClearRow(sec,row);
   }

   // make new tree - with tracks - exact position
   fout->cd();
   TTree * treenew = new TTree("Seedref","Seedref");
   AliTPCseed *  ioseed  =   (AliTPCseed*)trackarr->At(0);   
   TBranch * br = treenew->Branch("seeds","AliTPCseed",&ioseed,32000,99);

   for (Int_t itrack =0; itrack<ntracks; itrack++){
     ioseed  =  (AliTPCseed*) trackarr->At(itrack); 
     br->SetAddress(&ioseed);
     treenew->Fill();
   }

   treenew->Write();
   printf("1\n");
   delete ce;
   printf("2\n");
   //delete cr;
   //   carray->GetTree()->Write();
   printf("3\n");
   //   delete carray;
   printf("4\n");
   //   cf->Close(); 
   printf("5\n");
   fout->Close(); 
   printf("6\n");
   ftracks->Close();
   printf("7\n");
   return 0;
}
