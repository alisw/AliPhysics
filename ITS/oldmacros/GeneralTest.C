Int_t GeneralTest(Int_t verpoint=2) {
   Int_t rc=0;

//Test TPC simulation

   gROOT->LoadMacro("$(ALICE_ROOT)/macros/grun.C");
   grun();

   Int_t ver=gAlice->GetDetector("TPC")->IsVersion();
   delete gAlice; gAlice=0;

   if ((ver!=1)&&(ver!=2)) {
      cerr<<"Invalid TPC version: "<<ver<<" ! (must be 1 or 2)\n";
      return 12345;
   }

   if (ver==2) {
     gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCHits2Digits.C");
     if (rc=AliTPCHits2Digits()) return rc;

   }

//Test TPC reconstruction
   gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCFindClusters.C");
   if (rc=AliTPCFindClusters()) return rc;

   gROOT->LoadMacro("$(ALICE_ROOT)/TPC/AliTPCFindTracks.C");
   if (rc=AliTPCFindTracks()) return rc;
   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/TPCtracks.C");
   if (rc=TPCtracks()) return rc; 


   TFile *file=TFile::Open("galice.root");
   if (!file->IsOpen()) {cerr<<"Can't open galice.root !\n"; exit(4);}

   if (!(gAlice=(AliRun*)file->Get("gAlice"))) {
     cerr<<"gAlice have not been found on galice.root !\n";
     exit(5);
   }
           
   AliITSgeom *gm = ((AliITS*)gAlice->GetDetector("ITS"))->GetITSgeom();
   delete gAlice; gAlice=0;

   if (!gm) {
       cerr  << "This version of the ITS geometry does not have a AliITSgeom"
	     << " defined" << endl;
      return 12345;
   }
cout<<" verpoint = "<<verpoint<<"\n";
   switch(verpoint){
   case 1:
   {
   printf(" Start Fast Points calculation \n");
   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/ITSHitsToFastPoints.C");
   if (rc=ITSHitsToFastPoints()) return rc;
    }
    break;  
  case 2:
  printf("Start digitization \n");

  gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSHits2DigitsDefault.C");
  if (rc=AliITSHits2Digits()) return rc;

   printf("start reconstruction\n");

//Test ITS reconstruction
   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/AliITSFindClusters.C");
   if (rc=AliITSFindClusters()) return rc;   
  
  break;
  default:
     cerr<<"Invalid Recpoint version !\n";
      file->Close();
      exit(7);  
  } 
  
   file->Close();
   
   
   printf(" Start ITS tracking \n");     
   gROOT->LoadMacro("$(ALICE_ROOT)/ITS/ITStrackingGeneral.C");       
   if (rc=ITStrackingGeneral()) return rc;
   
   return rc;   
}
