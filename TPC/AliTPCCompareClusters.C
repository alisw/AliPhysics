/****************************************************************************
 *           Origin: M.Ivanov, CERN, Marian.Ivanov@cern.ch                 *
 ****************************************************************************/



#ifndef __CINT__
#include "alles.h"
#include "AliComplexCluster.h"
#endif

Int_t AliTPCCompareClusters(Int_t eventn) {
   cerr<<"Comparing clusters...\n";

   TFile *file=TFile::Open("galice.root");
   if (!file->IsOpen()) {cerr<<"Can't open galice.root !\n"; return 1;}

   //   TFile *cf=TFile::Open("AliTPCclusters.root");
   TFile *cf=TFile::Open("AliTPCclusters.root");
   if (!cf->IsOpen()){cerr<<"Can't open AliTPCclusters.root !\n"; return 3;}

   AliTPCParam *dig=(AliTPCParam *)cf->Get("75x40_100x60");
   if (!dig) {cerr<<"TPC parameters have not been found !\n"; return 2;}
   
   TFile *fout=TFile::Open("AliTPCClustersDif.root","new"); 
   if (!fout->IsOpen()) fout=TFile::Open("AliTPCClustersDif.root","recreate");  

   //connect exact clusters
   file->cd();
   AliTPCClustersArray *ce=new AliTPCClustersArray;
   ce->Setup(dig);
   ce->SetClusterType("AliComplexCluster");
   char  cname[100];
   sprintf(cname,"TreeCExact_%d",eventn);
   ce->ConnectTree(cname);
   //connect reconstructed clusters
   cf->cd();
   AliTPCClustersArray *cr=new AliTPCClustersArray;
   cr->Setup(dig);
   cr->SetClusterType("AliTPCcluster");  
   sprintf(cname,"TreeC_TPC_%d",eventn);
   cr->ConnectTree(cname);
   
   //create differantial clusters
   fout->cd();
   AliTPCClustersArray carray;
   carray.Setup(dig);
   carray.SetClusterType("AliDifCluster");
   carray.MakeTree();
   
   Int_t nrows=Int_t(cr->GetTree()->GetEntries());
   for (Int_t n=0; n<nrows; n++) {
       AliSegmentID *sr=cr->LoadEntry(n);
       Int_t sec,row;
       dig->AdjustSectorRow(sr->GetID(),sec,row);
       ce->LoadRow(sec,row);
       cr->LoadRow(sec,row);
       AliTPCClustersRow* clrowe = ce->GetRow(sec,row);
       AliTPCClustersRow* clrowr = cr->GetRow(sec,row);
       AliTPCClustersRow* clrowdiff = carray.CreateRow(sec,row);

       Int_t ncl=clrowe->GetArray()->GetEntriesFast();
       const Int_t kNIS=dig->GetNInnerSector(), kNOS=dig->GetNOuterSector();
       Int_t npads, sign;
       if (sec < kNIS) {
          npads = dig->GetNPadsLow(row);
          sign = (sec < kNIS/2) ? 1 : -1;
       } else {
          npads = dig->GetNPadsUp(row);
          sign = ((sec-kNIS) < kNOS/2) ? 1 : -1;
       }

       while (ncl--) {
	 //           AliTPCcluster *cl=(AliTPCcluster*)clrow[ncl];
           AliComplexCluster *cl=(AliComplexCluster*)(clrowe->GetArray()->At(ncl));
	   AliDifCluster * cldif = new(clrowdiff->Append()) AliDifCluster;

           Double_t x=dig->GetPadRowRadii(sec,row);
	   Double_t y=cl->fY;
	   y = ( y + 0.5 -  0.5*npads) *dig->GetPadPitchWidth(sec);
	   Double_t z=cl->fX*dig->GetZWidth();
	   z = sign*(dig->GetZLength() - z);
           Float_t cs, sn, tmp;
           dig->AdjustCosSin(sec,cs,sn);
	   Int_t ncr = clrowr->GetArray()->GetEntriesFast();
	   AliTPCcluster *clrec =0;
	   while (ncr--) {
	     clrec = (AliTPCcluster*)(clrowr->GetArray()->At(ncr));
	     if (clrec->GetLabel(0) == cl->fTracks[0]){
	       //cout<<cl->fTracks[0]<<"\n";
	       cldif->fTracks[0]= cl->fTracks[0];
	       cldif->fTracks[1]= clrec->GetLabel(1);
	       cldif->fTracks[2]= clrec->GetLabel(2);

	       cldif->fX = clrec->GetZ();
	       cldif->fY = clrec->GetY();
	       cldif->fQ = clrec->GetQ();
	       cldif->fOrigQ = cl->fQ;
	       cldif->fDx = cldif->fX - z;
	       cldif->fAngleX= cl->fSigmaX2;
	       cldif->fAngleY= cl->fSigmaY2;	       
	       cldif->fDy = cldif->fY - y;	       	       
	     }
	     
	   }
           //tmp = x*cs-y*sn; y= x*sn+y*cs; x=tmp;

	   // pm->SetPoint(ncl,x,y,z);
       }
       cr->ClearRow(sec,row);
       ce->ClearRow(sec,row);
       carray.StoreRow(sec,row);
       carray.ClearRow(sec,row);
   }
   delete ce;
   delete cr;
   carray.GetTree()->Write("Diff");
   cf->Close();
   fout->Close();
   file->Close();
   return 0;
}
