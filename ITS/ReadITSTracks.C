#include "iostream.h"


  void ReadITSTracks(Int_t evNumber1=0,Int_t evNumber2=0) {

  const char *filename="itstracks.root";
  
  ///////////////// Dynamically link some shared libs ////////////////////////////////
  
  if (gClassTable->GetID("AliRun") < 0) {
    gROOT->LoadMacro("loadlibs.C");
    loadlibs();
  } else {
    delete gAlice;
    gAlice=0;
  }

// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *file = (TFile*)gROOT->GetListOfFiles()->FindObject(filename);
   if (!file) file = new TFile(filename);

//
//   Loop over events 
//

   char tname[30];

   for (int nev=0; nev<= evNumber2; nev++) {

   sprintf(tname,"TreeT%d",nev);
   TTree *tracktree=(TTree*)file->Get(tname);
   TBranch *tbranch=tracktree->GetBranch("ITStracks");
	//cout<<" open the file \n"; 
	
   Int_t nentr=tracktree->GetEntries();

   TObjArray tarray(nentr);
   AliITSiotrack *iotrack=0;
   printf("nentr %d\n",nentr);
	
   for (Int_t i=0; i<nentr; i++) {
      AliITSiotrack *iotrack=new AliITSiotrack;
      // tarray.AddAt(new AliITSiotrack,i);
      // iotrack=(AliITSiotrack*)tarray.UncheckedAt(i);
       tbranch->SetAddress(&iotrack);
       tracktree->GetEvent(i);
		 tarray.AddLast(iotrack);
   }
	file->Close();		 
	
	  AliITSiotrack *iotrack;
   for (Int_t i=0; i<nentr; i++) {	
	 iotrack=(AliITSiotrack*)tarray.UncheckedAt(i);
	 if(!iotrack) continue;
     Int_t label=iotrack->GetLabel();
	  Double_t phistate=iotrack->GetStatePhi();
	  Double_t tgl=iotrack->GetStateTgl();	
	  Double_t Zstate=iotrack->GetStateZ();
	  Double_t Dr=iotrack->GetStateD();		  
	  Double_t C=iotrack->GetStateC();
	  
      cout<<" track label = "<<label<<"\n";
      cout<<" phi z D tanl C = "<<phistate<<" "<<Zstate<<" "<<Dr<<" "<<tgl<<" "<<C<<"\n"; 	  
	  
	  		  		    
	    

 /*
       Int_t label=iotrack->GetLabel();
       Float_t x=iotrack->GetX();
       Float_t y=iotrack->GetY();
       Float_t z=iotrack->GetZ();
       printf(" i %d label x y z %d %f %f %f\n",i,label,x,y,z);
       // delete iotrack;
		 */
		 
   }  

   }   // event loop    
}

