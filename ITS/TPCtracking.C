void TPCtracking( Int_t N =-1 )
{
  // Dynamically link some shared libs                    
/*    gROOT->LoadMacro("loadlibs.C");
    loadlibs();

// Connect the Root Galice file containing Geometry, Kine and Hits
   TFile *in=TFile::Open("galice.root");
   if (!in->IsOpen()) {cerr<<"Can't open galice.root !\n"; return;}
 
   if (!(gAlice=(AliRun*)in->Get("gAlice"))) {
     cerr<<"gAlice have not been found on galice.root !\n";
     return;
   }              
  
  if (N < 0 ) N = (Int_t)gAlice->TreeE()->GetEntries();
  cout<<endl<<N<<" EVENTS FOUND"<<endl;
  */

  cout<<"Running AliBarrelReconstruction.C  ...\n";
  gROOT->LoadMacro("AliBarrelReconstruction.C");
  AliBarrelReconstruction(N);
  
  cout<<"Running AliTPCComparison.C  ...\n";
  gROOT->LoadMacro("AliTPCComparison.C");
  AliTPCComparison(N);
  
  
 }
