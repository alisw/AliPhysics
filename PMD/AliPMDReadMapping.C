void AliPMDReadMapping()
{
  TFile * f = TFile::Open("/Users/basanta/Misc/Public/PMD/OCDBMacro/CDB_IDEAL/PMD/Calib/Mapping/Run0_999999999_v0_s0.root");

  f->ls();
  if (!AliCDBEntry)
    {
      printf("Something is wrong ************ \n");
    }
  else if(AliCDBEntry)
    {
      AliCDBEntry->PrintId(); 
      AliCDBEntry->PrintMetaData();
      AliPMDMappingData * mapdata = 0;
      mapdata = (AliPMDMappingData*)AliCDBEntry->GetObject();
      mapdata->Print(); 
    }
}
