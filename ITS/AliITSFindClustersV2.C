/****************************************************************
*  This macro converts AliITSRecPoint(s) to AliITSclusterV2(s)  *
*           Origin: I.Belikov, CERN, Jouri.Belikov@cern.ch      *
*****************************************************************/

#ifndef __CINT__
  #include <Riostream.h>

  #include "AliRun.h"
  #include "AliRunLoader.h"
  #include "AliLoader.h"
  #include "AliITS.h"
  #include "AliITSgeom.h"
  #include "AliITSclustererV2.h"

  #include "TStopwatch.h"
#endif

Int_t AliITSFindClustersV2(Char_t SlowOrFast='f')
{

    cerr<<"AliITSRecPoint(s) -> AliITSclusterV2(s)...\n";
    
    if (gAlice) 
     {
      delete gAlice->GetRunLoader();
      delete gAlice; 
      gAlice=0;
     }
 
    AliRunLoader* rl = AliRunLoader::Open("galice.root");
    if (rl == 0x0)
     {
      cerr<<"AliITSHits2DigitsDefault.C : Can not open session RL=NULL"
           << endl;
       return 3;
     }
     
    Int_t retval = rl->LoadgAlice();
    if (retval)
     {
      cerr<<"AliITSHits2DigitsDefault.C : LoadgAlice returned error"
           << endl;
       delete rl;
       return 3;
     }
    gAlice=rl->GetAliRun();
    rl->LoadHeader();
    retval = rl->LoadKinematics();
    if (retval)
     {
      cerr<<"AliITSHits2DigitsDefault.C : LoadKinematics returned error"
           << endl;
       delete rl;
       return 3;
     }
    
    AliITSLoader* gime = (AliITSLoader*)rl->GetLoader("ITSLoader");
    if (gime == 0x0)
     {
      cerr<<"AliITSHits2DigitsDefault.C : can not get ITS loader"
           << endl;
     }

   rl->GetEvent(0);

   AliITS *ITS  = (AliITS*)gAlice->GetModule("ITS");
   if (!ITS) { cerr<<"Can't find the ITS !\n"; delete rl; return 3; }
   AliITSgeom *geom=ITS->GetITSgeom();
 
   TClonesArray *clusters=new TClonesArray("AliITSclusterV2",10000);
  
   gime->LoadRawClusters("recreate");

   if (SlowOrFast=='f') 
    {
       gime->SetRecPointsFileName("ITS.FastRecPoints.root");
    }
   if (gime->LoadRecPoints())
    {
      cerr<<"Load Rec Pints returned error !\n"; 
      delete rl;
      return 4;
    }

   TClonesArray *points = new TClonesArray("AliITSRecPoint",10000);

   Float_t lp[5]; 
   Int_t lab[6];
   
   Int_t iEvent;
   for (iEvent = 0; iEvent< rl->GetNumberOfEvents() ; iEvent++) 
    {

     rl->GetEvent(iEvent);

     TTree *cTree = gime->TreeC();
     if (cTree == 0x0)  
      {
       gime->MakeTree("C");
       cTree = gime->TreeC();
      }
   
     cTree->Branch("Clusters",&clusters);
     
     TTree *pTree=gime->TreeR();
     if (pTree == 0x0) 
      {
        cerr<<"Can not get TreeR !\n"; delete rl;return 5;
      }
     TBranch *branch = 0;
     if (SlowOrFast=='f') {
       branch = pTree->GetBranch("ITSRecPointsF");
     }
     else {
       branch = pTree->GetBranch("ITSRecPoints");
     }
     if (!branch) 
      { 
       cerr<<"Can't get ITSRecPoints branch !\n"; 
       delete rl;
       return 6; 
      }
    
     branch->SetAddress(&points);
  
     AliStack* stack = rl->Stack();
     if (stack == 0x0)
      {
       cerr<<"AliITSFindClustersV2.C : Can not get stack"
           << endl;
       delete rl;
       return 3;
      }

     TClonesArray &cl=*clusters;
     Int_t nclusters=0;
     Int_t nentr=(Int_t)branch->GetEntries();

     cerr<<"Number of entries: "<<nentr<<endl;

     for (Int_t i=0; i<nentr; i++) 
      {
       points->Clear();
       clusters->Clear();
       branch->GetEvent(i);
       Int_t ncl=points->GetEntriesFast(); if (ncl==0){cTree->Fill();continue;}
       Int_t lay,lad,det; geom->GetModuleId(i,lay,lad,det);
       if ( (lay<0) || (lad<0) || (det<0))
        {
          ::Error("AliITSFindClustersV2.C","No such a module %d",i);
          continue;
        }
       Float_t x,y,zshift; geom->GetTrans(lay,lad,det,x,y,zshift); 
       Double_t rot[9];    geom->GetRotMatrix(lay,lad,det,rot);
       Double_t yshift = x*rot[0] + y*rot[1];
       Int_t ndet=(lad-1)*geom->GetNdetectors(lay) + (det-1);
       nclusters+=ncl;

       Float_t kmip=1; // ADC->mip normalization factor for the SDD and SSD 
       if(lay==4 || lay==3){kmip=280.;};
       if(lay==6 || lay==5){kmip=38.;};

       for (Int_t j=0; j<ncl; j++) 
        {
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
             TParticle *part=(TParticle*)stack->Particle(label);
             if (part == 0x0)
              cerr<<"Can not get particle with label "<<label<<endl;
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
       cTree->Fill(); 
//       clusters->Delete(); points->Delete();
     }
    gime->WriteRawClusters("OVERWRITE");
    cerr<<"Number of clusters: "<<nclusters<<endl;
    
   }
   
   delete clusters;
   delete points;
   delete rl;
   return 0;

}














