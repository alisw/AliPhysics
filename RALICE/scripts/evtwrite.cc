//////////////////////////////////////////////////////////////////////////
// Test of the AliEvent, AliTrack, AliVertex and AliJet functionality.
// Various tracks are created and some momenta and masses are entered
// 'by hand'. These tracks are then grouped into vertices and jets
// using the RALICE facilities to build events.
// In the program several events of different size are created to enable
// testing of the multiple event structure on the output file.
// The data structures are written onto the output file events.root
// 
//--- NvE 28-may-2001 UU-SAP Utrecht ***
//////////////////////////////////////////////////////////////////////////
{
 gSystem->Load("ralice");

 // Create an output file and a Tree 
 TFile* f=new TFile("events.root","RECREATE","AliEvent test output");
 TTree* tree=new TTree("T","Test tree for AliEvent output");

 // Define an event (this is also the main vertex)
 AliEvent* vmain=new AliEvent();

 // Branch in the tree for trk output
 Int_t split=0;
 Int_t bsize=32000;
 tree->Branch("Events","AliEvent",&vmain,bsize,split); 

 // Define some particle data
 Float_t p1[3]={0,0,1};
 Float_t p2[3]={0,0,-1};
 Float_t p3[3]={1,0,0};
 Float_t p4[3]={-1,0,0};

 Ali3Vector p;

 AliTrack t1,t2,t3,t4;

 AliJet j1,j2;

 AliVertex vt,vj,v1,v2;

 Float_t pos[3]={0,0,0};

 Double_t m1,m2,thcms,phicms;
 Double_t pi=acos(-1.);

 Int_t nev=5; // The number of events to be created for the output file
 for (Int_t iev=1; iev<=nev; iev++)
 {
  cout << endl;
  cout << " === Going for event number : " << iev << endl;
  cout << endl;

  vmain->SetRunNumber(12345);
  vmain->SetEventNumber(iev);

  if ((iev != 1) && (iev != 4)) // Create first and 4th event as empty
  {

   cout << " Event : " << iev << endl;

   p.SetVector(p1,"car");
   t1.Set3Momentum(p);
   t1.SetMass(0.135);
   t1.SetCharge(1);

   p.SetVector(p2,"car");
   t2.Set3Momentum(p);
   t2.SetMass(0.135);
   t2.SetCharge(2);

   p.SetVector(p3,"car");
   t3.Set3Momentum(p);
   t3.SetMass(0.135);
   t3.SetCharge(3);

   p.SetVector(p4,"car");
   t4.Set3Momentum(p);
   t4.SetMass(0.135);
   t4.SetCharge(4);

   if (iev != 3) // No decays for event number 3
   { 
    // Testing decay of t3
    m1=0.01;
    m2=0.02;
    thcms=0;
    phicms=0;
    t3.Decay(m1,m2,thcms,phicms);

    // Test second level decay
    m1=0.001;
    m2=0.002;
    thcms=0;
    phicms=0;

    AliTrack* tx;
    for (Int_t i=1; i<=t3.GetNdecay(); i++)
    {
     tx=t3.GetDecayTrack(i);
     tx->Decay(m1,m2,thcms,phicms);
    }
   }

   j1.AddTrack(t1);
   j1.AddTrack(t2);

   j2.AddTrack(t3);
   j2.AddTrack(t4);

   // Vertex testing with Tracks
   vt.AddTrack(t1);
   vt.AddTrack(t2);
   vt.AddTrack(t3);
   vt.AddTrack(t4);
   pos[0]=1;
   pos[1]=2;
   pos[2]=3;
   vt.SetPosition(pos,"car");

   // Vertex testing with Jets
   vj.AddJet(j1);
   vj.AddJet(j2);
   pos[0]=4;
   pos[1]=5;
   pos[2]=6;
   vj.SetPosition(pos,"car");

   // Secondary Vertex testing
   v1.AddTrack(t1);
   v1.AddTrack(t2);
   v1.AddTrack(t4);
   pos[0]=11;
   pos[1]=12;
   pos[2]=13;
   v1.SetPosition(pos,"car");

   v2.AddTrack(t2);
   v2.AddTrack(t3);
   pos[0]=21;
   pos[1]=22;
   pos[2]=23;
   v2.SetPosition(pos,"car");

   v1.AddVertex(v2); // Make v2 a secondary vertex of v1

   // Define the main vertex vmain and connect all tracks and vertices to it
   pos[0]=0.01;
   pos[1]=0.02;
   pos[2]=0.03;
   vmain->SetPosition(pos,"car");
   vmain->AddTrack(t1);
   vmain->AddTrack(t2);
   vmain->AddTrack(t3);
   vmain->AddTrack(t4);
   if (iev != 2) // These vertices are absent in event number 2
   {
    vmain->AddVertex(vt);
    vmain->AddVertex(vj);
   }
   vmain->AddVertex(v1);
  }
   
  // Print the full event data structure
  vmain.ListAll();

  // Write the complete structure to the output Tree
  tree->Fill();

  // Reset the complete main Vertex structure
  vmain->Reset();
  vt.Reset();
  vj.Reset();
  v1.Reset();
  v2.Reset();
  t1.Reset();
  t2.Reset();
  t3.Reset();
  t4.Reset();
  j1.Reset();
  j2.Reset();
 }

 // Provide overview of the Tree contents
 cout << endl;
 tree->Print();

 // Close output file
 f->Write();
 f->Close();
}
