void testsim (Int_t nevent=100, const char *config="testconfig.C")
{
  //
  // Simple macro to run aliroot in a batch mode
  //
  cerr<<" ___________________________________________________________________ "<<endl;
  cerr<<" "<<endl;
  cerr<<"             MESS ==> Beginning of the PHOS simulation."<<endl;
  cerr<<" ___________________________________________________________________ "<<endl;
  TStopwatch timer;
  timer.Start();

  gAlice->Init(config); 
  cerr<<" ___________________________________________________________________ "<<endl;
  cerr<<" "<<endl;
  cerr<<"             MESS ==> " << nevent << " : The simulation will last about 4 minutes."<<endl;
   cerr<<" ___________________________________________________________________ "<<endl;
  gAlice->Run(nevent);
  timer.Stop();
  timer.Print();
   cerr<<" ___________________________________________________________________ "<<endl;
  cerr<<" "<<endl;
  cerr << "             MESS ==> Simulation ended successfully. " << endl ; 
  cerr<<" ___________________________________________________________________ "<<endl;
}
