// $Id$

/**
   Macro for displaying L3 cluster/track or raw data. 
   Uses the AliHLTDisplay class. 
*/

void display(Int_t minslice,Int_t maxslice,Char_t *file="tracks.raw",Char_t *path="./",Char_t *gfile="$(ALIGEOPATH)/alice.geom")
{
  gStyle->SetOptStat(0);
  Int_t slice[2] = {minslice,maxslice};

  AliHLTDisplay *a = new AliHLTDisplay(slice,gfile);
  //a->Setup(file,path);
  a->Setup(file,path,0,kTRUE);
  
  
  /* Choose one of the following */
  a->DisplayAll();
  //a->DisplayTracks();
  //a->DisplayClusters();

  //a->DisplayClusterRow(1,151,file,"colz");
}


void display_cl(Int_t ev=0,Char_t *path="./",Char_t *gfile="$(LEVEL3)/GEO/alice.geom")
{
  gStyle->SetOptStat(0);
  Int_t slice[2] = {0,35};

  Char_t file[1024];
  sprintf(file,"%s/tracks_%d.raw",path,ev);

  a = new AliHLTDisplay(slice,gfile);
  //a->Setup(file,path);
  a->Setup(file,path,ev,kTRUE);
  
  
  /* Choose one of the following */
  a->DisplayAll();
  //a->DisplayTracks();
  //a->DisplayClusters();

  //a->DisplayClusterRow(1,151,file,"colz");
}
