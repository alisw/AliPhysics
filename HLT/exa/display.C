// $Id$

/**
   Macro for displaying L3 cluster/track or raw data. 
   Uses the AliL3Display class. 
*/

void display(int minslice,int maxslice,char *file="tracks.raw",char *path="./",char *gfile="$(LEVEL3)/GEO/alice.geom")
{
  AliL3Logger l;
  //l.UnSet(AliL3Logger::kDebug);
  //l.UnSet(AliL3Logger::kAll);
  //l.Set(AliL3Logger::kInformational);
  //l.UseStdout();
  l.UseStream();

  gStyle->SetOptStat(0);
  int slice[2] = {minslice,maxslice};

  a = new AliL3Display(slice,gfile);
  a->Setup(file,path);
  //a->Setup("/prog/alice/data/Rawdata/1_patch/pp/recon_6/tracks.raw","/prog/alice/data/Rawdata/1_patch/pp/recon_6/");

  //a->DisplayAll();
  a->DisplayTracks();
  //a->DisplayClusters();

  //a->DisplayClusterRow(1,151,file,"colz");
}


