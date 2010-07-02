const Char_t *name[]={
  "clResY", "clPullY", 
  "trkltResY", "trkltPullY", "trkltResZ", "trkltPullZ", "trkltResPhi"
};
const Char_t *title[]={
  "cluster resolution r-phi"
 ,"cluster pulls r-phi"
 ,"tracklet resolution r-phi"
 ,"tracklet pulls r-phi"
 ,"tracklet resolution z"
 ,"tracklet pulls z"
 ,"tracklet resolution phi"
};
void AliTRDmakeTrendDB()
{
  gSystem->Load("libANALYSIS.so");
  gSystem->Load("libANALYSISalice.so");
  gSystem->Load("libTENDER.so");
  gSystem->Load("libPWG1.so");


/*  AliCDBManager *cdb(AliCDBManager::Instance());
  cdb->SetDefaultStorage("local://$ALICE_ROOT/OCDB");
  cdb->SetRun(0);*/
  AliTRDtrendingManager *tm(AliTRDtrendingManager::Instance());
  AliTRDtrendValue *val(NULL);

  for(Int_t iv(0); iv<7; iv++){
    Double_t limits[]={-1., 1., -2., 2., -3., 3., -4., 4., -5., 5.};
    char *messages[]={
      "acceptable"
     ,"worrisome"
     ,"far from normal"
     ,"off the limits"
     ,"totally off"
    };
    tm->AddValue("resolution", name[iv], title[iv], limits, messages, "Alex Bercuci/a.bercuci@gsi.de", "Markus Fasel/m.fasel@gsi.de,Ionut Arsene/I.C.Arsene@gsi.de,Christoph Blume/blume@ikf.uni-frankfurt.de");
    val=tm->GetValue("resolution", name[iv]);
    val->Set(0.5+iv);
  }
  tm->Print();
  //tm->Save();
}