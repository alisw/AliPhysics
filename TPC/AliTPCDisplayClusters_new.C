#ifndef __CINT__
  #include "alles.h"
  #include "AliTPCtracker.h"
  #include "TView.h"
  #include "TPolyMarker3D.h"
  #include "AliSimDigits.h"

#endif
Int_t AliTPCDisplayClusters(Int_t eventn=0, Int_t noiseth=15) {
   cerr<<"Displaying clusters...\n";

   TFile *file=TFile::Open("galice.root");
   if (!file->IsOpen()) {cerr<<"Can't open galice.root !\n"; return 1;}

   TFile *cf=TFile::Open("AliTPCclusters.root");
   if (!cf->IsOpen()){cerr<<"Can't open AliTPCclusters.root !\n"; return 3;}

   AliTPCParam *dig=(AliTPCParam *)cf->Get("75x40_100x60");
   if (!dig) {cerr<<"TPC parameters have not been found !\n"; return 2;}


   // some "constants"
   Int_t markerColorSignal = 5;
   Int_t markerColorBgr = 2;
   Int_t MASK = 10000000;


   TCanvas *c1=new TCanvas("cdisplay", "Cluster display",0,0,700,730);
   TView *v=new TView(1);
   v->SetRange(-430,-560,-430,430,560,1710);
   c1->Clear();
   c1->SetFillColor(1);
   c1->SetTheta(90.);
   c1->SetPhi(0.);

   AliTPCClustersArray *ca=new AliTPCClustersArray;
   ca->Setup(dig);
   ca->SetClusterType("AliTPCcluster");
   char  cname[100];
   sprintf(cname,"TreeC_TPC_%d",eventn);

   ca->ConnectTree(cname);
   Int_t nrows=Int_t(ca->GetTree()->GetEntries());
   for (Int_t n=0; n<nrows; n++) {
       AliSegmentID *s=ca->LoadEntry(n);
       Int_t sec,row;
       dig->AdjustSectorRow(s->GetID(),sec,row);
       AliTPCClustersRow &clrow = *ca->GetRow(sec,row);
       Int_t ncl=clrow.GetArray()->GetEntriesFast();
       TPolyMarker3D *pm=new TPolyMarker3D(ncl);
       TPolyMarker3D *pmSignal=new TPolyMarker3D(ncl); // polymarker for signal
       Int_t imarBgr=0;
       Int_t imarSignal=0;
       while (ncl--) {
           AliTPCcluster *cl=(AliTPCcluster*)clrow[ncl];
           Double_t x=dig->GetPadRowRadii(sec,row), y=cl->GetY(), z=cl->GetZ();
	   if (cl->GetQ()<noiseth) continue;
           Float_t cs, sn, tmp;
           dig->AdjustCosSin(sec,cs,sn);
           tmp = x*cs-y*sn; y= x*sn+y*cs; x=tmp; 
	   if (TMath::Abs(z)>50) continue;
	   Int_t trackId = cl->GetLabel(0);
	   if (trackId<0) continue;  //bug in tracks ID
	   if (trackId < MASK-1) {
	     pmSignal->SetPoint(imarSignal,x,y,z);
	     imarSignal++;
	   } else {
	     pm->SetPoint(imarBgr,x,y,z);
	     imarBgr++;
	   }          
       }
       ca->ClearRow(sec,row);
       
      // change color for signal
       pm->SetMarkerSize(1); 
       pm->SetMarkerColor(markerColorBgr);
       pm->SetMarkerStyle(1);
       pm->Draw();

       pmSignal->SetMarkerSize(1); 
       pmSignal->SetMarkerColor(markerColorSignal);
       pmSignal->SetMarkerStyle(1);
       pmSignal->Draw();      
   }
   delete ca;
   cf->Close();

   TGeometry *geom=(TGeometry*)file->Get("AliceGeom");
   //   TList *list = geom->GetListOfNodes();
   TNode * main = (TNode*)((geom->GetListOfNodes())->First());
   TIter next(main->GetListOfNodes());
   TNode  *module=0;
   while((module = (TNode*)next())) {
     char ch[100];
     sprintf(ch,"%s\n",module->GetTitle());
     //printf("%s\n",module->GetTitle());
     if (ch[0]=='T'&&ch[1]=='P' && ch[2]=='C')  //if TPC draw
       module->SetVisibility(3);
     else
       module->SetVisibility(-1);
   }
     
   
   geom->Draw("same");
   //v->Draw();
   c1->Modified(); c1->Update(); 

   file->Close();
   return 0;
}

