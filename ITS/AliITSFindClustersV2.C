#ifndef __CINT__
  #include <iostream.h>

  #include "AliRun.h"
  #include "AliITS.h"
  #include "AliITSgeom.h"
  #include "AliITSRecPoint.h"
  #include "AliITSclusterV2.h"

  #include "TFile.h"
  #include "TTree.h"
  #include "TParticle.h"
#endif

Int_t AliITSFindClustersV2(Char_t SlowOrFast='f') {
/****************************************************************
 *  This macro converts AliITSRecPoint(s) to AliITSclusterV2(s) *
 ****************************************************************/
   cerr<<"AliITSRecPoint(s) -> AliITSclusterV2(s)...\n";

   if (gAlice) {delete gAlice; gAlice=0;}

   TFile *in=TFile::Open("galice.root");
   if (!in->IsOpen()) { cerr<<"Can't open galice.root !\n"; return 1; }

   if (!(gAlice=(AliRun*)in->Get("gAlice"))) {
      cerr<<"Can't find gAlice !\n";
      return 2;
   }
   gAlice->GetEvent(0);

   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   if (!ITS) { cerr<<"Can't find the ITS !\n"; return 3; }
   AliITSgeom *geom=ITS->GetITSgeom();
 
   TFile *out=TFile::Open("AliITSclustersV2.root","new");
   if (!out->IsOpen()) {
      cerr<<"Delete old AliITSclustersV2.root !\n"; 
      return 4;
   }
   geom->Write();

   TClonesArray *clusters=new TClonesArray("AliITSclusterV2",10000);
   TTree *cTree=new TTree("TreeC_ITS_0","ITS clusters");
   cTree->Branch("Clusters",&clusters);

   TTree *pTree=gAlice->TreeR();
   if (!pTree) { cerr<<"Can't get TreeR !\n"; return 5; }
   TBranch *branch = 0;
   if (SlowOrFast=='f') {
     branch = pTree->GetBranch("ITSRecPointsF");
   }
   else {
     branch = pTree->GetBranch("ITSRecPoints");
   }
   if (!branch) { cerr<<"Can't get ITSRecPoints branch !\n"; return 6; }
   TClonesArray *points=new TClonesArray("AliITSRecPoint",10000);
   branch->SetAddress(&points);

   TClonesArray &cl=*clusters;
   Int_t nclusters=0;
   Int_t nentr=(Int_t)branch->GetEntries();

   cerr<<"Number of entries: "<<nentr<<endl;

   //   Float_t lp[5]; Int_t lab[6]; //Why can't it be inside a loop ?
   Float_t * lp = new Float_t[5]; 
   Int_t * lab = new Int_t[6];

   for (Int_t i=0; i<nentr; i++) {
       points->Clear();
       branch->GetEvent(i);
       Int_t ncl=points->GetEntriesFast(); if (ncl==0){cTree->Fill();continue;}
       Int_t lay,lad,det; geom->GetModuleId(i,lay,lad,det);
       Float_t x,y,zshift; geom->GetTrans(lay,lad,det,x,y,zshift); 
       Double_t rot[9];    geom->GetRotMatrix(lay,lad,det,rot);
       Double_t yshift = x*rot[0] + y*rot[1];
       Int_t ndet=(lad-1)*geom->GetNdetectors(lay) + (det-1);
       nclusters+=ncl;

       Float_t kmip=1; // ADC->mip normalization factor for the SDD and SSD 
       if(lay==4 || lay==3){kmip=280.;};
       if(lay==6 || lay==5){kmip=38.;};

       for (Int_t j=0; j<ncl; j++) {
          AliITSRecPoint *p=(AliITSRecPoint*)points->UncheckedAt(j);
          //Float_t lp[5];
          lp[0]=-p->GetX()-yshift; if (lay==1) lp[0]=-lp[0];
          lp[1]=p->GetZ()+zshift;
          lp[2]=p->GetSigmaX2();
          lp[3]=p->GetSigmaZ2();
          lp[4]=p->GetQ(); lp[4]/=kmip;
          //Int_t lab[6]; 
          lab[0]=p->GetLabel(0);lab[1]=p->GetLabel(1);lab[2]=p->GetLabel(2);
          lab[3]=ndet;

          Int_t label=lab[0];
          if (label>=0) {
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
	  }

          new(cl[j]) AliITSclusterV2(lab,lp);
       }
       cTree->Fill(); clusters->Delete();
       points->Delete();
   }
   cTree->Write();

   cerr<<"Number of clusters: "<<nclusters<<endl;

   delete [] lp;
   delete [] lab;

   delete cTree; delete clusters; delete points;

   delete gAlice; gAlice=0;

   in->Close();
   out->Close();

   return 0;

}














