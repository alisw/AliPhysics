void display(int minslice,int maxslice,char *file="")
{
  AliL3Logger l;
//  l.UnSet(AliL3Logger::kDebug);
//  l.UnSet(AliL3Logger::kAll);
//  l.Set(AliL3Logger::kInformational);
  //l.UseStdout();
   l.UseStream();
   
  gStyle->SetOptStat(0);

  int slice[2] = {minslice,maxslice};

  a = new AliL3Display(slice);
  a->Setup("tracks.raw","./");
  //a->Setup("/prog/alice/data/Rawdata/1_patch/pp/recon_6/tracks.raw","/prog/alice/data/Rawdata/1_patch/pp/recon_6/");
  //a->DisplayAll();
  a->DisplayTracks();
  //a->DisplayClusters();
  //a->DisplayClusterRow(1,151,file,"colz");

}


