void display(char *file="")
{
  AliL3Logger l;
//  l.UnSet(AliL3Logger::kDebug);
//  l.UnSet(AliL3Logger::kAll);
//  l.Set(AliL3Logger::kInformational);
   l.UseStdout();
   l.UseStream();
   
  gStyle->SetOptStat(0);

  int slice[2] = {0,35};

  a = new AliL3Display(slice);
  a->Setup("tracks.raw","./");
  a->DisplayAll();
  //a->DisplayTracks();
  //a->DisplayClusters();
  //a->DisplayClusterRow(1,151,file,"colz");

}

