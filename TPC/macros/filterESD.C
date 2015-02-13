/// \file filterESD.C

void filterESD(){
  //
  //
  //
  //gSystem->Setenv("alien_CLOSE_SE","ALICE::GSI::SE");
  //TGrid * alien =     TGrid::Connect("alien://",0,0,"t");
  printf("Filtering esd.C\n");
  gSystem->AddIncludePath("-I$ALICE_ROOT/TPC/macros");
  gROOT->LoadMacro("$ALICE_ROOT/TPC/macros/AliXRDPROOFtoolkit.cxx+");
  //  AliXRDPROOFtoolkit::FilterList("esd.txt","AliESDs.root esdTree AliESDfriends.root *",0);
  AliXRDPROOFtoolkit::FilterList("esd.txt","AliESDs.root esdTree",0);
}
